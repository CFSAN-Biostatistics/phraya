#!/bin/bash
#SBATCH --job-name=phraya_bench
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=4:00:00

# Main benchmark array job
# Decodes task ID → (target, aligner, replicate), runs alignment+calling with timing
#
# Aligners (9 total):
#   bwa-mem2      — alignment only (throughput baseline, SAM output)
#   minimap2      — alignment only (throughput baseline, SAM output)
#   bwa-pipeline  — bwa mem + sort + mpileup + bcftools call (fair phraya comparison)
#   phraya          — phraya --strategy balanced (Myers ≤500bp, WFA fallback, K=2, ±50bp)
#   phraya-sensitive — phraya --strategy sensitive (WFA, K=50 chains, ±25bp)
#   phraya-fast     — phraya --strategy fast (K=1 chain + divergence cutoff, ±150bp)
#   bowtie2         — industry-standard gapped aligner (--very-sensitive preset)
#   minibwa         — hybrid bwa-mem/minimap2 (3× faster than bwa-mem)
#   rammap          — pure Rust minimap2 clone (SIMD DP, -x sr preset)

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
NUM_ALIGNERS=9  # bwa-mem2, minimap2, bwa-pipeline, phraya×3, bowtie2, minibwa, rammap
NUM_REPLICATES=$REPLICATES

TARGET_IDX=$((TASK_ID / (NUM_ALIGNERS * NUM_REPLICATES)))
ALIGNER_IDX=$(((TASK_ID / NUM_REPLICATES) % NUM_ALIGNERS))
REPLICATE_IDX=$((TASK_ID % NUM_REPLICATES))

# Map aligner index to name
ALIGNERS=("bwa-mem2" "minimap2" "bwa-pipeline" "phraya" "phraya-sensitive" "phraya-fast" "bowtie2" "minibwa" "rammap")
ALIGNER="${ALIGNERS[$ALIGNER_IDX]}"

# Read target from targets.conf (skip comments, filter large if needed)
if [[ "${INCLUDE_LARGE:-0}" -eq 0 ]]; then
    TARGET_LINE=$(grep -v '^#' "$TARGETS_FILE" | grep -v '^[[:space:]]*$' | grep -v '|large|' | sed -n "$((TARGET_IDX + 1))p")
else
    TARGET_LINE=$(grep -v '^#' "$TARGETS_FILE" | grep -v '^[[:space:]]*$' | sed -n "$((TARGET_IDX + 1))p")
fi

if [[ -z "$TARGET_LINE" ]]; then
    echo "ERROR: Invalid target index $TARGET_IDX" >&2
    exit 1
fi

IFS='|' read -r TARGET_ID TARGET_PATH SIZE_CLASS GENOME_SIZE_GB <<< "$TARGET_LINE"

# Construct paths (prefer 30k-subsampled reads for speed; fall back to full reads)
DATA_DIR="$DATA_ROOT/$TARGET_PATH/data"
REF_FASTA="$DATA_DIR/reference/reference.fasta"
if [[ -f "$DATA_DIR/reads/reads_1_30k.fastq.gz" ]]; then
    READS_1="$DATA_DIR/reads/reads_1_30k.fastq.gz"
    READS_2="$DATA_DIR/reads/reads_2_30k.fastq.gz"
else
    READS_1="$DATA_DIR/reads/reads_1.fastq.gz"
    READS_2="$DATA_DIR/reads/reads_2.fastq.gz"
fi

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

# Pre-flight check for phraya binary (all phraya-* variants share the same binary)
if [[ "$ALIGNER" == phraya* ]]; then
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

# Write timing.json: parse all fields from timing.txt and add SLURM metadata
if [[ -f "$OUT_DIR/timing.txt" ]]; then
    _get() { grep "^$1=" "$OUT_DIR/timing.txt" | head -1 | cut -d= -f2 | tr -d ' '; }
    WALL=$(_get wall_seconds)
    PEAK_RSS=$(_get peak_rss_gb)
    TOTAL_READS=$(_get total_reads)
    N_ALIGNED=$(_get n_aligned)
    UNALIGNED_FRAC=$(_get unaligned_frac)
    cat > "$OUT_DIR/timing.json" <<JSON
{
  "aligner":         "$ALIGNER",
  "target":          "$TARGET_ID",
  "replicate":       $REPLICATE_IDX,
  "wall_seconds":    ${WALL:-null},
  "wall_time_s":     ${WALL:-null},
  "peak_rss_gb":     ${PEAK_RSS:-null},
  "total_reads":     ${TOTAL_READS:-0},
  "n_aligned":       ${N_ALIGNED:-0},
  "unaligned_frac":  ${UNALIGNED_FRAC:-null},
  "threads":         $THREADS,
  "node":            "$SLURM_JOB_NODELIST"
}
JSON
else
    echo "ERROR: timing.txt not found after run" >&2
    exit 1
fi

echo "=== Task Complete: ${WALL:-?}s ==="
