#!/bin/bash
# Main orchestrator for HPC aligner benchmark
# Usage: ./run_benchmark.sh [--targets targets.conf] [--large] [--dry-run]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config/global.env"

# Parse arguments
TARGETS_FILE="$SCRIPT_DIR/config/targets.conf"
INCLUDE_LARGE=0
DRY_RUN=0

while [[ $# -gt 0 ]]; do
    case $1 in
        --targets)
            TARGETS_FILE="$2"
            shift 2
            ;;
        --large)
            INCLUDE_LARGE=1
            shift
            ;;
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        -h|--help)
            cat <<EOF
Usage: $0 [OPTIONS]

Options:
  --targets FILE    Use custom targets file (default: config/targets.conf)
  --large           Include large targets (T8b 17Gb, T8c 4.3Gb)
  --dry-run         Show what would be run without submitting jobs
  -h, --help        Show this help message

Examples:
  $0                    # Run on default targets (small + medium)
  $0 --large            # Include large wheat genomes
  $0 --dry-run          # Preview array size and configuration
EOF
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            echo "Use --help for usage information" >&2
            exit 1
            ;;
    esac
done

# Verify targets file exists
if [[ ! -f "$TARGETS_FILE" ]]; then
    echo "ERROR: Targets file not found: $TARGETS_FILE" >&2
    exit 1
fi

# Generate run timestamp
RUN_ID="run_$(date +%Y%m%d_%H%M%S)"
RUN_DIR="$RESULTS_ROOT/$RUN_ID"
mkdir -p "$RUN_DIR"

echo "=== Phraya Aligner Benchmark ==="
echo "Run ID:     $RUN_ID"
echo "Output:     $RUN_DIR"
echo "Targets:    $TARGETS_FILE"
echo "Threads:    $THREADS"
echo "Replicates: $REPLICATES"
echo

# Discover available nodes for rotation
echo "Discovering available nodes..."
NODELIST=$("$SCRIPT_DIR/utils/nodelist.sh" "$REPLICATES" 2>&1) || {
    echo "ERROR: Node discovery failed" >&2
    echo "$NODELIST" >&2
    echo
    echo "This may indicate SLURM is not available or no nodes are idle/mixed." >&2
    echo "On non-HPC systems, you may need to modify nodelist.sh or skip node rotation." >&2
    exit 1
}
echo "Using nodes: $NODELIST"
echo

# Count targets (excluding comments and large if not requested)
if [[ $INCLUDE_LARGE -eq 0 ]]; then
    TARGET_COUNT=$(grep -v '^#' "$TARGETS_FILE" | grep -v '|large|' | wc -l)
else
    TARGET_COUNT=$(grep -v '^#' "$TARGETS_FILE" | wc -l)
fi

if [[ $TARGET_COUNT -eq 0 ]]; then
    echo "ERROR: No targets found in $TARGETS_FILE" >&2
    exit 1
fi

# Array dimensions: targets × aligners × replicates
# bwa-mem2 and minimap2 are alignment-only baselines; bwa-pipeline is the fair
# phraya comparison (alignment + sort + mpileup + bcftools call).
ALIGNERS=("bwa-mem2" "minimap2" "bwa-pipeline" "phraya")
NUM_ALIGNERS=${#ALIGNERS[@]}
ARRAY_SIZE=$((TARGET_COUNT * NUM_ALIGNERS * REPLICATES))

echo "Benchmark configuration:"
echo "  Targets:    $TARGET_COUNT"
echo "  Aligners:   $NUM_ALIGNERS (${ALIGNERS[*]})"
echo "  Replicates: $REPLICATES"
echo "  Array size: $ARRAY_SIZE tasks"
echo

if [[ $DRY_RUN -eq 1 ]]; then
    echo "DRY RUN: Would submit $ARRAY_SIZE array job tasks"
    echo
    echo "Tasks would be:"
    grep -v '^#' "$TARGETS_FILE" | grep -v '|large|' | while IFS='|' read -r tid tpath tclass tsize; do
        for aligner in "${ALIGNERS[@]}"; do
            for rep in $(seq 0 $((REPLICATES - 1))); do
                echo "  - $tid / $aligner / rep_$rep (genome: ${tsize}GB)"
            done
        done
    done
    exit 0
fi

# Step 1: Run STREAM Triad characterization (if not already cached)
echo "=== Step 1: STREAM Triad Platform Characterization ==="
if [[ ! -f "$RUN_DIR/stream_triad.txt" ]]; then
    echo "Submitting STREAM Triad job..."
    STREAM_JOB=$(sbatch --parsable \
        --job-name="benchmark_stream_$RUN_ID" \
        --output="$RUN_DIR/stream_%N.log" \
        --nodelist="$NODELIST" \
        "$SCRIPT_DIR/stream.slurm" "$RUN_DIR")

    echo "  STREAM job ID: $STREAM_JOB"
    echo "  Waiting for completion..."

    # Wait for STREAM to complete
    while squeue -j "$STREAM_JOB" -h &>/dev/null; do
        sleep 5
    done

    if [[ ! -f "$RUN_DIR/stream_triad.txt" ]]; then
        echo "ERROR: STREAM job completed but stream_triad.txt not found" >&2
        exit 1
    fi

    echo "  STREAM Triad measurement complete"
else
    echo "  Using cached STREAM Triad results"
fi
echo

# Step 2: Submit main benchmark array job
echo "=== Step 2: Submit Benchmark Array Job ==="
echo "Submitting $ARRAY_SIZE array tasks..."
BENCHMARK_JOB=$(sbatch --parsable \
    --job-name="phraya_benchmark_$RUN_ID" \
    --array="0-$((ARRAY_SIZE - 1))" \
    --output="$RUN_DIR/slurm-%A_%a.log" \
    --export=ALL,RUN_DIR="$RUN_DIR",TARGETS_FILE="$TARGETS_FILE",NODELIST="$NODELIST",INCLUDE_LARGE="$INCLUDE_LARGE" \
    "$SCRIPT_DIR/benchmark.slurm")

echo "  Benchmark job ID: $BENCHMARK_JOB"
echo "  Monitor: squeue -j $BENCHMARK_JOB"
echo "  Logs:    $RUN_DIR/slurm-*.log"
echo

# Step 3: Submit aggregation job (depends on benchmark completion)
echo "=== Step 3: Submit Aggregation Job ==="
AGGREGATE_JOB=$(sbatch --parsable \
    --job-name="benchmark_aggregate_$RUN_ID" \
    --dependency=afterok:$BENCHMARK_JOB \
    --output="$RUN_DIR/aggregate.log" \
    --export=ALL,RUN_DIR="$RUN_DIR" \
    "$SCRIPT_DIR/utils/aggregate.slurm")

echo "  Aggregation job ID: $AGGREGATE_JOB (runs after $BENCHMARK_JOB)"
echo

echo "=== Benchmark Submitted ==="
echo "Run ID: $RUN_ID"
echo "Status: squeue -j $BENCHMARK_JOB,$AGGREGATE_JOB"
echo "Cancel: scancel $BENCHMARK_JOB $AGGREGATE_JOB"
echo
echo "Results will appear in: $RUN_DIR/results.json"
echo
echo "After completion, score results:"
echo "  python ~/data-commons/test/benchmarking/alignment/score.py $RUN_DIR/results.json --sensitivity"
