#!/bin/bash
# Main orchestrator for HPC aligner benchmark
# Usage: ./run_benchmark.sh [--targets targets.conf] [--large] [--dry-run] [--alphabet {dna|protein|both}]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config/global.env"

# Parse arguments
TARGETS_FILE_OVERRIDE=""
INCLUDE_LARGE=0
DRY_RUN=0
ALPHABET="dna"

while [[ $# -gt 0 ]]; do
    case $1 in
        --targets)
            TARGETS_FILE_OVERRIDE="$2"
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
        --alphabet)
            ALPHABET="$2"
            shift 2
            ;;
        -h|--help)
            cat <<EOF
Usage: $0 [OPTIONS]

Options:
  --targets FILE       Use custom targets file (default: config/targets.conf for
                        dna, config/targets_protein.conf for protein). Not valid
                        with --alphabet both (ambiguous which alphabet it targets)
                        — run dna and protein separately with their own --targets.
  --alphabet MODE      dna (default) | protein | both. See ADR-0013;
                        --alphabet protein exercises DIAMOND/BLASTP + phraya
                        --alphabet protein instead of the DNA aligner set, on
                        a 4-arg <proteome> <query> <out_dir> <threads> wrapper
                        contract (no paired reads) — see BENCHMARK_EXPANSION.md.
  --large              Include large targets (T8b 17Gb, T8c 4.3Gb) (dna only)
  --dry-run            Show what would be run without submitting jobs
  -h, --help           Show this help message

Examples:
  $0                        # Run DNA aligners on default targets (small + medium)
  $0 --large                # Include large wheat genomes
  $0 --alphabet protein     # Run DIAMOND/BLASTP/phraya-protein on protein targets
  $0 --alphabet both        # Submit both DNA and protein runs (separate array jobs)
  $0 --dry-run              # Preview array size and configuration
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

case "$ALPHABET" in
    dna|protein|both) ;;
    *)
        echo "ERROR: --alphabet must be dna, protein, or both (got: $ALPHABET)" >&2
        exit 1
        ;;
esac

if [[ -n "$TARGETS_FILE_OVERRIDE" && "$ALPHABET" == "both" ]]; then
    echo "ERROR: --targets is ambiguous with --alphabet both (which alphabet does it target?)." >&2
    echo "  Run --alphabet dna and --alphabet protein separately, each with its own --targets." >&2
    exit 1
fi

# DNA and protein each get an independent aligner set and wrapper-arg-count
# contract (5-arg paired-reads vs 4-arg single-query-FASTA — see
# BENCHMARK_EXPANSION.md). Kept as separate arrays rather than one combined
# list so a protein-only run never asks a DNA-only tool (or vice versa) to run
# on data its wrapper contract can't accept.
ALIGNERS_DNA=("bwa-mem2" "minimap2" "bwa-pipeline" "phraya" "phraya-sensitive" "phraya-sensitive-linear" "phraya-fast" "bowtie2" "minibwa" "rammap")
ALIGNERS_PROTEIN=("diamond-blastp" "blastp" "phraya-protein" "phraya-protein-sensitive")

# Generate run timestamp (shared across both alphabets when --alphabet both)
RUN_ID="run_$(date +%Y%m%d_%H%M%S)"
RUN_DIR_ROOT="$RESULTS_ROOT/$RUN_ID"
mkdir -p "$RUN_DIR_ROOT"

echo "=== Phraya Aligner Benchmark ==="
echo "Run ID:     $RUN_ID"
echo "Output:     $RUN_DIR_ROOT"
echo "Alphabet:   $ALPHABET"
echo "Threads:    $THREADS"
echo "Replicates: $REPLICATES"
echo

# Discover available nodes for rotation (shared across both alphabets)
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

# Step 1: STREAM Triad characterization (hardware property, alphabet-independent;
# shared cache and shared per-run copy regardless of how many alphabets run).
echo "=== Step 1: STREAM Triad Platform Characterization ==="
STREAM_CACHE="$SCRIPT_DIR/cache/stream_triad_c064.txt"
if [[ -f "$STREAM_CACHE" ]]; then
    cp "$STREAM_CACHE" "$RUN_DIR_ROOT/stream_triad.txt"
    echo "  Using cached STREAM Triad: $(cat "$STREAM_CACHE") MB/s"
elif [[ ! -f "$RUN_DIR_ROOT/stream_triad.txt" ]]; then
    echo "  Submitting STREAM Triad job..."
    STREAM_JOB=$(sbatch --parsable \
        --job-name="benchmark_stream_$RUN_ID" \
        --output="$RUN_DIR_ROOT/stream_%N.log" \
        --nodelist="$NODELIST" \
        "$SCRIPT_DIR/stream.slurm" "$RUN_DIR_ROOT")

    echo "  STREAM job ID: $STREAM_JOB"
    echo "  Waiting for completion..."

    while squeue -j "$STREAM_JOB" -h 2>/dev/null | grep -q "$STREAM_JOB"; do
        sleep 5
    done

    if [[ ! -f "$RUN_DIR_ROOT/stream_triad.txt" ]]; then
        echo "ERROR: STREAM job completed but stream_triad.txt not found" >&2
        exit 1
    fi
    echo "  STREAM Triad measurement complete"
    cp "$RUN_DIR_ROOT/stream_triad.txt" "$STREAM_CACHE"
    echo "  Cached result for future runs"
fi
echo

# submit_alphabet_run <alphabet> <targets_file> <aligners_array_name>
#
# Submits one benchmark array job + its dependent aggregation job for one
# alphabet. When --alphabet is a single value, RUN_DIR == RUN_DIR_ROOT (today's
# exact layout, unchanged). When --alphabet both, each alphabet gets its own
# RUN_DIR_ROOT/<alphabet> subdirectory so the two runs' outputs never collide.
submit_alphabet_run() {
    local alphabet="$1" targets_file="$2" aligners_ref="$3"
    local -n aligners="$aligners_ref"
    local run_dir="$RUN_DIR_ROOT"
    if [[ "$ALPHABET" == "both" ]]; then
        run_dir="$RUN_DIR_ROOT/$alphabet"
        mkdir -p "$run_dir"
        cp "$RUN_DIR_ROOT/stream_triad.txt" "$run_dir/stream_triad.txt"
    fi

    if [[ ! -f "$targets_file" ]]; then
        echo "ERROR: Targets file not found: $targets_file" >&2
        exit 1
    fi

    local target_count
    if [[ "$alphabet" == "dna" && $INCLUDE_LARGE -eq 0 ]]; then
        target_count=$(grep -v '^#' "$targets_file" | grep -v '^[[:space:]]*$' | grep -v '|large|' | wc -l)
    else
        target_count=$(grep -v '^#' "$targets_file" | grep -v '^[[:space:]]*$' | wc -l)
    fi

    if [[ $target_count -eq 0 ]]; then
        echo "ERROR: No targets found in $targets_file" >&2
        exit 1
    fi

    local num_aligners=${#aligners[@]}
    local array_size=$((target_count * num_aligners * REPLICATES))

    echo "=== [$alphabet] Benchmark configuration ==="
    echo "  Targets file: $targets_file"
    echo "  Targets:      $target_count"
    echo "  Aligners:     $num_aligners (${aligners[*]})"
    echo "  Replicates:   $REPLICATES"
    echo "  Array size:   $array_size tasks"
    echo

    if [[ $DRY_RUN -eq 1 ]]; then
        echo "DRY RUN [$alphabet]: Would submit $array_size array job tasks"
        local large_filter="grep -v '|large|'"
        if [[ "$alphabet" != "dna" || $INCLUDE_LARGE -eq 1 ]]; then
            large_filter="cat"
        fi
        grep -v '^#' "$targets_file" | grep -v '^[[:space:]]*$' | eval "$large_filter" | while IFS='|' read -r tid tpath tclass tsize; do
            for aligner in "${aligners[@]}"; do
                for rep in $(seq 0 $((REPLICATES - 1))); do
                    echo "  - [$alphabet] $tid / $aligner / rep_$rep"
                done
            done
        done
        echo
        return 0
    fi

    echo "=== [$alphabet] Submit Benchmark Array Job ==="
    local benchmark_job
    benchmark_job=$(sbatch --parsable \
        --job-name="phraya_benchmark_${RUN_ID}_${alphabet}" \
        --array="0-$((array_size - 1))" \
        --output="$run_dir/slurm-%A_%a.log" \
        --export=ALL,SCRIPT_DIR="$SCRIPT_DIR",RUN_DIR="$run_dir",TARGETS_FILE="$targets_file",NODELIST="$NODELIST",INCLUDE_LARGE="$INCLUDE_LARGE",ALPHABET="$alphabet" \
        "$SCRIPT_DIR/benchmark.slurm")

    echo "  Benchmark job ID: $benchmark_job"
    echo "  Monitor: squeue -j $benchmark_job"
    echo "  Logs:    $run_dir/slurm-*.log"
    echo

    echo "=== [$alphabet] Submit Aggregation Job ==="
    local aggregate_job
    aggregate_job=$(sbatch --parsable \
        --job-name="benchmark_aggregate_${RUN_ID}_${alphabet}" \
        --dependency=afterok:$benchmark_job \
        --output="$run_dir/aggregate.log" \
        --export=ALL,SCRIPT_DIR="$SCRIPT_DIR",RUN_DIR="$run_dir" \
        "$SCRIPT_DIR/utils/aggregate.slurm")

    echo "  Aggregation job ID: $aggregate_job (runs after $benchmark_job)"
    echo
    echo "  Results will appear in: $run_dir/results.json"
    echo
}

TARGETS_DNA="${TARGETS_FILE_OVERRIDE:-$SCRIPT_DIR/config/targets.conf}"
TARGETS_PROTEIN="${TARGETS_FILE_OVERRIDE:-$SCRIPT_DIR/config/targets_protein.conf}"

case "$ALPHABET" in
    dna)
        submit_alphabet_run dna "$TARGETS_DNA" ALIGNERS_DNA
        ;;
    protein)
        submit_alphabet_run protein "$TARGETS_PROTEIN" ALIGNERS_PROTEIN
        ;;
    both)
        submit_alphabet_run dna "$TARGETS_DNA" ALIGNERS_DNA
        submit_alphabet_run protein "$TARGETS_PROTEIN" ALIGNERS_PROTEIN
        ;;
esac

if [[ $DRY_RUN -eq 1 ]]; then
    exit 0
fi

echo "=== Benchmark Submitted ==="
echo "Run ID: $RUN_ID"
echo "After completion, score results:"
if [[ "$ALPHABET" == "both" ]]; then
    echo "  python ~/data-commons/test/benchmarking/alignment/score.py $RUN_DIR_ROOT/dna/results.json --sensitivity"
    echo "  python ~/data-commons/test/benchmarking/alignment/score.py $RUN_DIR_ROOT/protein/results.json --sensitivity"
else
    echo "  python ~/data-commons/test/benchmarking/alignment/score.py $RUN_DIR_ROOT/results.json --sensitivity"
fi
