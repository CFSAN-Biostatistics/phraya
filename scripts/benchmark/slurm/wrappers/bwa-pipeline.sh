#!/bin/bash
# BWA full variant-calling pipeline (fair phraya comparison)
# Steps: bwa mem → samtools sort → bcftools mpileup | bcftools call → VCF
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

# Step 1: align + sort to BAM; capture peak RSS via measure_rss.py
PYTHON="${PYTHON3_BIN:-python3}"
MEASURE="$SCRIPT_DIR/utils/measure_rss.py"
"$PYTHON" "$MEASURE" "$OUT_DIR/time_verbose.txt" -- \
    bash -c "$BWA_BIN mem -t $THREADS $REF $READS_1 $READS_2 2>$OUT_DIR/bwa.log \
        | $SAMTOOLS_BIN sort -@ $THREADS -o $OUT_DIR/alignment.bam -"

T_ALIGN=$((SECONDS - START))

# Step 2: index BAM
$SAMTOOLS_BIN index "$OUT_DIR/alignment.bam"

# Step 3: pileup → VCF (-Ou passes BCF through pipe)
$BCFTOOLS_BIN mpileup -f "$REF" -d 10000 -q 20 -Q 20 -Ou \
    "$OUT_DIR/alignment.bam" \
    | $BCFTOOLS_BIN call -mv -Oz -o "$OUT_DIR/variants.vcf.gz"
$BCFTOOLS_BIN index "$OUT_DIR/variants.vcf.gz"

ELAPSED=$((SECONDS - START))

PEAK_RSS_KB=$(grep 'Maximum resident' "$OUT_DIR/time_verbose.txt" | grep -oP '\d+' | tail -1)
PEAK_RSS_GB=$(awk "BEGIN{printf \"%.3f\", ${PEAK_RSS_KB:-0}/1048576}")

# Variant count for sanity check
N_VARIANTS=$($BCFTOOLS_BIN stats "$OUT_DIR/variants.vcf.gz" | grep "^SN.*number of SNPs" | cut -f4)

# Unaligned fraction from BAM
N_TOTAL=$($SAMTOOLS_BIN view -c "$OUT_DIR/alignment.bam" 2>/dev/null || echo 0)
N_MAPPED=$($SAMTOOLS_BIN view -c -F4 "$OUT_DIR/alignment.bam" 2>/dev/null || echo 0)
N_UNMAPPED=$(( N_TOTAL - N_MAPPED ))
UNALIGNED_FRAC=$(awk "BEGIN{if($N_TOTAL>0) printf \"%.4f\", $N_UNMAPPED/$N_TOTAL; else print \"0.0000\"}")

cat > "$OUT_DIR/timing.txt" <<TIMING
wall_seconds=$ELAPSED
threads=$THREADS
aligner=bwa-pipeline
peak_rss_gb=${PEAK_RSS_GB}
t_align_sort=${T_ALIGN}
n_variants=${N_VARIANTS:-unknown}
total_reads=${N_TOTAL}
n_aligned=${N_MAPPED}
n_unaligned=${N_UNMAPPED}
unaligned_frac=${UNALIGNED_FRAC}
TIMING

echo "Elapsed: ${ELAPSED}s (align+sort: ${T_ALIGN}s), variants: ${N_VARIANTS:-?}, RSS=${PEAK_RSS_GB}GB, unaligned=${UNALIGNED_FRAC}" >&2
