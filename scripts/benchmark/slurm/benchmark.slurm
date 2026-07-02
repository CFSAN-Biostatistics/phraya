#!/bin/bash
#SBATCH --job-name=phraya_bench
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=4:00:00

# Main benchmark array job
# Decodes task ID → (target, aligner, replicate), runs alignment+calling with timing
#
# Aligners (4 total):
#   bwa-mem2      — alignment only (throughput baseline, SAM output)
#   minimap2      — alignment only (throughput baseline, SAM output)
#   bwa-pipeline  — bwa mem + sort + mpileup + bcftools call (fair phraya comparison)
#   phraya        — full alignment + variant observation pipeline

set -euo pipefail

# Load environment (uses absolute tool paths — module system unreliable in SLURM)
SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "$SCRIPT_DIR/config/global.env"

# Verify required variables from environment
if [[ -z "${RUN_DIR:-}" ]] || [[ -z "${TARGETS_FILE:-}" ]]; then
    echo "ERROR: RUN_DIR and TARGETS_FILE must be exported by orchestrator" >&2
    exit 1
fi

# Decode array task ID → (target_idx, aligner_idx, replicate_idx)
TASK_ID=${SLURM_ARRAY_TASK_ID}
NUM_ALIGNERS=4  # bwa-mem2, minimap2, bwa-pipeline, phraya
NUM_REPLICATES=$REPLICATES

TARGET_IDX=$((TASK_ID / (NUM_ALIGNERS * NUM_REPLICATES)))
ALIGNER_IDX=$(((TASK_ID / NUM_REPLICATES) % NUM_ALIGNERS))
REPLICATE_IDX=$((TASK_ID % NUM_REPLICATES))

# Map aligner index to name
ALIGNERS=("bwa-mem2" "minimap2" "bwa-pipeline" "phraya")
ALIGNER="${ALIGNERS[$ALIGNER_IDX]}"

# Read target from targets.conf (skip comments, filter large if needed)
if [[ "${INCLUDE_LARGE:-0}" -eq 0 ]]; then
    TARGET_LINE=$(grep -v '^#' "$TARGETS_FILE" | grep -v '|large|' | sed -n "$((TARGET_IDX + 1))p")
else
    TARGET_LINE=$(grep -v '^#' "$TARGETS_FILE" | sed -n "$((TARGET_IDX + 1))p")
fi

if [[ -z "$TARGET_LINE" ]]; then
    echo "ERROR: Invalid target index $TARGET_IDX" >&2
    exit 1
fi

IFS='|' read -r TARGET_ID TARGET_PATH SIZE_CLASS GENOME_SIZE_GB <<< "$TARGET_LINE"

# Construct paths
DATA_DIR="$DATA_ROOT/$TARGET_PATH/data"
REF_FASTA="$DATA_DIR/reference/reference.fasta"
READS_1="$DATA_DIR/reads/reads_1_30k.fastq.gz"
READS_2="$DATA_DIR/reads/reads_2_30k.fastq.gz"

# Output directory
OUT_DIR="$RUN_DIR/$TARGET_ID/$ALIGNER/rep_$REPLICATE_IDX"
mkdir -p "$OUT_DIR"

echo "=== Benchmark Task ==="
echo "Task ID:    $TASK_ID"
echo "Target:     $TARGET_ID ($TARGET_PATH)"
echo "Aligner:    $ALIGNER"
echo "Replicate:  $REPLICATE_IDX"
echo "Node:       $SLURM_JOB_NODELIST"
echo "Output:     $OUT_DIR"
echo

# Verify input files exist
for f in "$REF_FASTA" "$READS_1" "$READS_2"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Input file not found: $f" >&2
        echo "Did you run acquire.py for target $TARGET_ID?" >&2
        exit 1
    fi
done

# Pre-flight check for phraya binary
if [[ "$ALIGNER" == "phraya" ]]; then
    PHRAYA="${PHRAYA_ROOT}/target/release/phraya"
    if [[ ! -f "$PHRAYA" ]]; then
        echo "ERROR: phraya binary not found at $PHRAYA" >&2
        exit 1
    fi
fi

# Run aligner via wrapper script (wrappers handle index build internally)
WRAPPER="$SCRIPT_DIR/wrappers/${ALIGNER}.sh"
if [[ ! -f "$WRAPPER" ]]; then
    echo "ERROR: Wrapper not found: $WRAPPER" >&2
    exit 1
fi

echo "Running $ALIGNER..."
"$WRAPPER" "$REF_FASTA" "$READS_1" "$READS_2" "$OUT_DIR" "$THREADS"

# Write timing.json (read wall_seconds from timing.txt, add metadata)
if [[ -f "$OUT_DIR/timing.txt" ]]; then
    WALL=$(grep 'wall_seconds=' "$OUT_DIR/timing.txt" | head -1 | cut -d= -f2 | tr -d ' ')
    cat > "$OUT_DIR/timing.json" <<JSON
{"aligner":"$ALIGNER","target":"$TARGET_ID","replicate":$REPLICATE_IDX,"wall_seconds":${WALL:-null},"threads":$THREADS,"node":"$SLURM_JOB_NODELIST"}
JSON
else
    echo "ERROR: timing.txt not found after run" >&2
    exit 1
fi

echo "=== Task Complete: ${WALL:-?}s ==="
