#!/bin/bash
#SBATCH --job-name=phraya_bench
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=4:00:00

# Main benchmark array job
# Decodes task ID → (target, aligner, replicate), runs alignment+calling with timing
#
# DNA aligners (10 total, ALPHABET=dna):
#   bwa-mem2               — alignment only (throughput baseline, SAM output)
#   minimap2                — alignment only (throughput baseline, SAM output)
#   bwa-pipeline             — bwa mem + sort + mpileup + bcftools call (fair phraya comparison)
#   phraya                   — phraya --strategy balanced (Myers ≤500bp, WFA fallback, K=2, ±50bp)
#   phraya-sensitive          — phraya --strategy sensitive (WFA, K=50 chains, ±25bp; gap-affine
#                               by default once ADR-0014 ships)
#   phraya-sensitive-linear   — phraya --strategy sensitive --gap-model linear (ADR-0014: the old
#                               linear-cost WFA path, kept for direct affine-vs-linear comparison)
#   phraya-fast               — phraya --strategy fast (K=1 chain + divergence cutoff, ±150bp)
#   bowtie2                   — industry-standard gapped aligner (--very-sensitive preset)
#   minibwa                   — hybrid bwa-mem/minimap2 (3× faster than bwa-mem)
#   rammap                    — pure Rust minimap2 clone (SIMD DP, -x sr preset)
#
# Protein aligners (4 total, ALPHABET=protein — ADR-0013):
#   diamond-blastp            — fast seed-and-extend protein search (DNA harness's minimap2 analog)
#   blastp                    — maximally-sensitive protein search reference (DNA harness's
#                               sensitive/BWA analog)
#   phraya-protein             — phraya --alphabet protein --strategy balanced
#   phraya-protein-sensitive   — phraya --alphabet protein --strategy sensitive
#
# Wrapper contract differs by alphabet (see BENCHMARK_EXPANSION.md): DNA wrappers
# take 5 positional args (<ref> <reads_1> <reads_2> <out_dir> <threads>, paired
# FASTQ); protein wrappers take 4 (<proteome> <query> <out_dir> <threads>, no
# read pairing).

set -euo pipefail

# Load environment (uses absolute tool paths — module system unreliable in SLURM)
SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "$SCRIPT_DIR/config/global.env"

# Verify required variables from environment
if [[ -z "${RUN_DIR:-}" ]] || [[ -z "${TARGETS_FILE:-}" ]]; then
    echo "ERROR: RUN_DIR and TARGETS_FILE must be exported by orchestrator" >&2
    exit 1
fi
ALPHABET="${ALPHABET:-dna}"

# Map alphabet to its aligner array (must match run_benchmark.sh's ALIGNERS_DNA /
# ALIGNERS_PROTEIN exactly — kept in sync by hand, same as the pre-alphabet code).
case "$ALPHABET" in
    dna)
        ALIGNERS=("bwa-mem2" "minimap2" "bwa-pipeline" "phraya" "phraya-sensitive" "phraya-sensitive-linear" "phraya-fast" "bowtie2" "minibwa" "rammap")
        ;;
    protein)
        ALIGNERS=("diamond-blastp" "blastp" "phraya-protein" "phraya-protein-sensitive")
        ;;
    *)
        echo "ERROR: Unknown ALPHABET: $ALPHABET (expected dna or protein)" >&2
        exit 1
        ;;
esac
NUM_ALIGNERS=${#ALIGNERS[@]}
NUM_REPLICATES=$REPLICATES

# Decode array task ID → (target_idx, aligner_idx, replicate_idx)
TASK_ID=${SLURM_ARRAY_TASK_ID}
TARGET_IDX=$((TASK_ID / (NUM_ALIGNERS * NUM_REPLICATES)))
ALIGNER_IDX=$(((TASK_ID / NUM_REPLICATES) % NUM_ALIGNERS))
REPLICATE_IDX=$((TASK_ID % NUM_REPLICATES))
ALIGNER="${ALIGNERS[$ALIGNER_IDX]}"

# Read target from targets.conf (skip comments; DNA additionally filters |large|
# unless requested — protein targets.conf has no size tiers to filter).
if [[ "$ALPHABET" == "dna" && "${INCLUDE_LARGE:-0}" -eq 0 ]]; then
    TARGET_LINE=$(grep -v '^#' "$TARGETS_FILE" | grep -v '^[[:space:]]*$' | grep -v '|large|' | sed -n "$((TARGET_IDX + 1))p")
else
    TARGET_LINE=$(grep -v '^#' "$TARGETS_FILE" | grep -v '^[[:space:]]*$' | sed -n "$((TARGET_IDX + 1))p")
fi

if [[ -z "$TARGET_LINE" ]]; then
    echo "ERROR: Invalid target index $TARGET_IDX" >&2
    exit 1
fi

# Both targets.conf and targets_protein.conf share the same 4-field pipe-delimited
# shape (TARGET_ID|ORGANISM_PATH|SIZE_CLASS|SIZE_METRIC); only the meaning of the
# 4th field differs (genome GB vs. proteome residues), which this script doesn't
# need to interpret.
IFS='|' read -r TARGET_ID TARGET_PATH SIZE_CLASS SIZE_METRIC <<< "$TARGET_LINE"

DATA_DIR="$DATA_ROOT/$TARGET_PATH/data"
OUT_DIR="$RUN_DIR/$TARGET_ID/$ALIGNER/rep_$REPLICATE_IDX"
mkdir -p "$OUT_DIR"

echo "=== Benchmark Task ==="
echo "Task ID:    $TASK_ID"
echo "Alphabet:   $ALPHABET"
echo "Target:     $TARGET_ID ($TARGET_PATH)"
echo "Aligner:    $ALIGNER"
echo "Replicate:  $REPLICATE_IDX"
echo "Node:       $SLURM_JOB_NODELIST"
echo "Output:     $OUT_DIR"
echo

WRAPPER="$SCRIPT_DIR/wrappers/${ALIGNER}.sh"
if [[ ! -f "$WRAPPER" ]]; then
    echo "ERROR: Wrapper not found: $WRAPPER" >&2
    exit 1
fi

if [[ "$ALPHABET" == "dna" ]]; then
    # 5-arg contract: <ref.fasta> <reads_1.fq.gz> <reads_2.fq.gz> <out_dir> <threads>
    REF_FASTA="$DATA_DIR/reference/reference.fasta"
    if [[ -f "$DATA_DIR/reads/reads_1_30k.fastq.gz" ]]; then
        READS_1="$DATA_DIR/reads/reads_1_30k.fastq.gz"
        READS_2="$DATA_DIR/reads/reads_2_30k.fastq.gz"
    else
        READS_1="$DATA_DIR/reads/reads_1.fastq.gz"
        READS_2="$DATA_DIR/reads/reads_2.fastq.gz"
    fi

    for f in "$REF_FASTA" "$READS_1" "$READS_2"; do
        if [[ ! -f "$f" ]]; then
            echo "ERROR: Input file not found: $f" >&2
            echo "Did you run acquire.py for target $TARGET_ID?" >&2
            exit 1
        fi
    done

    if [[ "$ALIGNER" == phraya* ]]; then
        PHRAYA="${PHRAYA_ROOT}/target/release/phraya"
        [[ -f "$PHRAYA" ]] || { echo "ERROR: phraya binary not found at $PHRAYA" >&2; exit 1; }
    fi

    echo "Running $ALIGNER..."
    "$WRAPPER" "$REF_FASTA" "$READS_1" "$READS_2" "$OUT_DIR" "$THREADS"
else
    # 4-arg contract: <proteome.fasta> <query.fasta> <out_dir> <threads> (ADR-0013:
    # no read pairing for protein search). Staged by gen_synthetic_protein.py per
    # BENCHMARK_EXPANSION.md, not acquire.py — see targets_protein.conf.
    PROTEOME="$DATA_DIR/protein/proteome.fasta"
    QUERY="$DATA_DIR/protein/queries.fasta"

    for f in "$PROTEOME" "$QUERY"; do
        if [[ ! -f "$f" ]]; then
            echo "ERROR: Input file not found: $f" >&2
            echo "Did you stage protein data for target $TARGET_ID with gen_synthetic_protein.py?" >&2
            exit 1
        fi
    done

    if [[ "$ALIGNER" == phraya* ]]; then
        PHRAYA="${PHRAYA_ROOT}/target/release/phraya"
        [[ -f "$PHRAYA" ]] || { echo "ERROR: phraya binary not found at $PHRAYA" >&2; exit 1; }
    fi

    echo "Running $ALIGNER..."
    "$WRAPPER" "$PROTEOME" "$QUERY" "$OUT_DIR" "$THREADS"
fi

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
  "alphabet":        "$ALPHABET",
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
