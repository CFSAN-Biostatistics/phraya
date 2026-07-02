#!/bin/bash
# BWA full variant-calling pipeline (fair phraya comparison)
# Steps: bwa mem → samtools sort → samtools mpileup → bcftools call → VCF
# Usage: bwa-pipeline.sh <ref.fasta> <reads_1.fq.gz> <reads_2.fq.gz> <out_dir> <threads>
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$SCRIPT_DIR/config/global.env"

REF=$1; READS_1=$2; READS_2=$3; OUT_DIR=$4; THREADS=$5

for f in "$REF" "$READS_1" "$READS_2"; do
    [[ -f "$f" ]] || { echo "ERROR: not found: $f" >&2; exit 1; }
done

# Build BWA index (flock-protected)
if [[ ! -f "${REF}.bwt" ]]; then
    (flock -x 200; [[ -f "${REF}.bwt" ]] || $BWA_BIN index "$REF") 200>"${REF}.bwa_index.lock"
fi

START=$SECONDS

# Step 1: align + sort to BAM
$BWA_BIN mem -t "$THREADS" "$REF" "$READS_1" "$READS_2" \
    | $SAMTOOLS_BIN sort -@ "$THREADS" -o "$OUT_DIR/alignment.bam" -
T_ALIGN=$((SECONDS - START))

# Step 2: index BAM
$SAMTOOLS_BIN index "$OUT_DIR/alignment.bam"

# Step 3: pileup → VCF (bcftools mpileup | call pipeline; -Ou passes BCF through pipe)
$BCFTOOLS_BIN mpileup -f "$REF" -d 10000 -q 20 -Q 20 -Ou \
    "$OUT_DIR/alignment.bam" \
    | $BCFTOOLS_BIN call -mv -Oz -o "$OUT_DIR/variants.vcf.gz"
$BCFTOOLS_BIN index "$OUT_DIR/variants.vcf.gz"

ELAPSED=$((SECONDS - START))

# Variant count for sanity check
N_VARIANTS=$($BCFTOOLS_BIN stats "$OUT_DIR/variants.vcf.gz" | grep "^SN.*number of SNPs" | cut -f4)

cat > "$OUT_DIR/timing.txt" <<TIMING
wall_seconds=$ELAPSED threads=$THREADS aligner=bwa-pipeline
t_align_sort=${T_ALIGN}s
n_variants=${N_VARIANTS:-unknown}
TIMING

echo "Elapsed: ${ELAPSED}s (align+sort: ${T_ALIGN}s), variants: ${N_VARIANTS:-?}" >&2
