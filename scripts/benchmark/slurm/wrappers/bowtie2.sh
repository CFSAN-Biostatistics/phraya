#!/bin/bash
# Bowtie2 alignment-only wrapper (throughput baseline — SAM output, no variant calling)
# Usage: bowtie2.sh <ref.fasta> <reads_1.fq.gz> <reads_2.fq.gz> <out_dir> <threads>
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$SCRIPT_DIR/config/global.env"

REF=$1; READS_1=$2; READS_2=$3; OUT_DIR=$4; THREADS=$5

for f in "$REF" "$READS_1" "$READS_2"; do
    [[ -f "$f" ]] || { echo "ERROR: not found: $f" >&2; exit 1; }
done

# Use bowtie2-align-s binary directly (avoids Perl wrapper Sys::Hostname issue)
BT2_ALIGN="/nfs/software/apps/micromamba/1.5.8/envs/bowtie2-v2.5.4/bin/bowtie2-align-s"
BT2_BUILD="${BOWTIE2_BUILD_BIN:-/nfs/software/apps/micromamba/1.5.8/envs/bowtie2-v2.5.4/bin/bowtie2-build}"

# Build Bowtie2 index (flock-protected, one-time per reference)
INDEX_BASE="${REF%.fasta}"
INDEX_BASE="${INDEX_BASE%.fa}"
if [[ ! -f "${INDEX_BASE}.1.bt2" ]]; then
    (flock -x 200; [[ -f "${INDEX_BASE}.1.bt2" ]] || $BT2_BUILD --threads $THREADS "$REF" "$INDEX_BASE" 2>"$OUT_DIR/bowtie2-build.log") 200>"${INDEX_BASE}.bt2_index.lock"
fi

PYTHON="${PYTHON3_BIN:-python3}"
MEASURE="$SCRIPT_DIR/utils/measure_rss.py"
START=$SECONDS
"$PYTHON" "$MEASURE" "$OUT_DIR/time_verbose.txt" -- \
    bash -c "$BT2_ALIGN --very-sensitive -x $INDEX_BASE -1 $READS_1 -2 $READS_2 -p $THREADS 2>$OUT_DIR/bowtie2.log" | $SAMTOOLS_BIN view -bS - > "$OUT_DIR/alignment.bam"
ELAPSED=$((SECONDS - START))

PEAK_RSS_KB=$(grep 'Maximum resident' "$OUT_DIR/time_verbose.txt" | grep -oP '\d+' | tail -1)
PEAK_RSS_GB=$(awk "BEGIN{printf \"%.3f\", ${PEAK_RSS_KB:-0}/1048576}")

# Unaligned fraction from BAM flagstat
N_TOTAL=$($SAMTOOLS_BIN view -c "$OUT_DIR/alignment.bam" 2>/dev/null || echo 0)
N_MAPPED=$($SAMTOOLS_BIN view -c -F4 "$OUT_DIR/alignment.bam" 2>/dev/null || echo 0)
N_UNMAPPED=$(( N_TOTAL - N_MAPPED ))
UNALIGNED_FRAC=$(awk "BEGIN{if($N_TOTAL>0) printf \"%.4f\", $N_UNMAPPED/$N_TOTAL; else print \"0.0000\"}")

cat > "$OUT_DIR/timing.txt" <<EOF
wall_seconds=$ELAPSED
threads=$THREADS
aligner=bowtie2
peak_rss_gb=${PEAK_RSS_GB}
total_reads=${N_TOTAL}
n_aligned=${N_MAPPED}
n_unaligned=${N_UNMAPPED}
unaligned_frac=${UNALIGNED_FRAC}
EOF

echo "Elapsed: ${ELAPSED}s, RSS=${PEAK_RSS_GB}GB, aligned=${N_MAPPED}/${N_TOTAL}" >&2
