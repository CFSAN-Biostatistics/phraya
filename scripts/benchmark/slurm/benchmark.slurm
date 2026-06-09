#!/bin/bash
#SBATCH --job-name=phraya_bench
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00

# Main benchmark array job
# Decodes task ID → (target, aligner, replicate), runs alignment with timing

set -euo pipefail

# Load environment
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config/global.env"

# Load modules
module load "$MODULE_BWA" || echo "WARNING: Failed to load $MODULE_BWA"
module load "$MODULE_MINIMAP" || echo "WARNING: Failed to load $MODULE_MINIMAP"
module load "$MODULE_SAMTOOLS" || echo "WARNING: Failed to load $MODULE_SAMTOOLS"

# Verify required variables from environment
if [[ -z "${RUN_DIR:-}" ]] || [[ -z "${TARGETS_FILE:-}" ]]; then
    echo "ERROR: RUN_DIR and TARGETS_FILE must be exported by orchestrator" >&2
    exit 1
fi

# Decode array task ID → (target_idx, aligner_idx, replicate_idx)
TASK_ID=${SLURM_ARRAY_TASK_ID}
NUM_ALIGNERS=2  # bwa-mem2, minimap2
NUM_REPLICATES=$REPLICATES

TARGET_IDX=$((TASK_ID / (NUM_ALIGNERS * NUM_REPLICATES)))
ALIGNER_IDX=$(((TASK_ID / NUM_REPLICATES) % NUM_ALIGNERS))
REPLICATE_IDX=$((TASK_ID % NUM_REPLICATES))

# Map aligner index to name
ALIGNERS=("bwa-mem2" "minimap2")
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
READS_1="$DATA_DIR/reads/reads_1.fastq.gz"
READS_2="$DATA_DIR/reads/reads_2.fastq.gz"

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

# Build index if not present (one-time per target, flock-protected)
INDEX_DIR="$DATA_DIR/reference"
case $ALIGNER in
    bwa-mem2)
        if [[ ! -f "$INDEX_DIR/reference.fasta.bwt.2bit.64" ]]; then
            (
                flock -x 200
                if [[ ! -f "$INDEX_DIR/reference.fasta.bwt.2bit.64" ]]; then
                    echo "Building BWA-MEM2 index for $TARGET_ID..."
                    bwa-mem2 index "$REF_FASTA"
                fi
            ) 200>"$INDEX_DIR/.bwa_index.lock"
        fi
        ;;
    minimap2)
        if [[ ! -f "$INDEX_DIR/reference.mmi" ]]; then
            (
                flock -x 200
                if [[ ! -f "$INDEX_DIR/reference.mmi" ]]; then
                    echo "Building minimap2 index for $TARGET_ID..."
                    minimap2 -d "$INDEX_DIR/reference.mmi" "$REF_FASTA"
                fi
            ) 200>"$INDEX_DIR/.minimap_index.lock"
        fi
        ;;
esac

# Run aligner via wrapper script
WRAPPER="$SCRIPT_DIR/wrappers/${ALIGNER}.sh"
if [[ ! -f "$WRAPPER" ]]; then
    echo "ERROR: Wrapper not found: $WRAPPER" >&2
    exit 1
fi

echo "Running $ALIGNER alignment..."
"$WRAPPER" "$REF_FASTA" "$READS_1" "$READS_2" "$OUT_DIR" "$THREADS"

# Parse timing output
if [[ -f "$OUT_DIR/timing.txt" ]]; then
    python3 "$SCRIPT_DIR/utils/parse_time.py" "$OUT_DIR/timing.txt" > "$OUT_DIR/timing.json"
else
    echo "ERROR: timing.txt not found after alignment" >&2
    exit 1
fi

echo "=== Task Complete ==="
