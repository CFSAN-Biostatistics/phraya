#!/bin/bash
# minimap2 alignment-only wrapper (throughput baseline — SAM output)
# Usage: minimap2.sh <ref.fasta> <reads_1.fq.gz> <reads_2.fq.gz> <out_dir> <threads>
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$SCRIPT_DIR/config/global.env"

REF=$1; READS_1=$2; READS_2=$3; OUT_DIR=$4; THREADS=$5

for f in "$REF" "$READS_1" "$READS_2"; do
    [[ -f "$f" ]] || { echo "ERROR: not found: $f" >&2; exit 1; }
done

# Build .mmi index (flock-protected, one-time per reference)
INDEX="${REF%.fasta}.mmi"
if [[ ! -f "$INDEX" ]]; then
    (flock -x 200; [[ -f "$INDEX" ]] || $MINIMAP2_BIN -d "$INDEX" "$REF") 200>"${REF}.mmi.lock"
fi

START=$SECONDS
{ /usr/bin/time -v \
    bash -c "$MINIMAP2_BIN -ax sr -t $THREADS $INDEX $READS_1 $READS_2 > $OUT_DIR/alignment.sam 2>$OUT_DIR/minimap2.log" \
  ; } 2> "$OUT_DIR/time_verbose.txt"
ELAPSED=$((SECONDS - START))

PEAK_RSS_KB=$(grep 'Maximum resident' "$OUT_DIR/time_verbose.txt" | grep -oP '\d+' | tail -1)
PEAK_RSS_GB=$(awk "BEGIN{printf \"%.3f\", ${PEAK_RSS_KB:-0}/1048576}")

N_TOTAL=$($SAMTOOLS_BIN view -c "$OUT_DIR/alignment.sam" 2>/dev/null || echo 0)
N_MAPPED=$($SAMTOOLS_BIN view -c -F4 "$OUT_DIR/alignment.sam" 2>/dev/null || echo 0)
N_UNMAPPED=$(( N_TOTAL - N_MAPPED ))
UNALIGNED_FRAC=$(awk "BEGIN{if($N_TOTAL>0) printf \"%.4f\", $N_UNMAPPED/$N_TOTAL; else print \"0.0000\"}")

cat > "$OUT_DIR/timing.txt" <<EOF
wall_seconds=$ELAPSED
threads=$THREADS
aligner=minimap2
peak_rss_gb=${PEAK_RSS_GB}
total_reads=${N_TOTAL}
n_aligned=${N_MAPPED}
n_unaligned=${N_UNMAPPED}
unaligned_frac=${UNALIGNED_FRAC}
EOF

echo "Elapsed: ${ELAPSED}s, RSS=${PEAK_RSS_GB}GB, aligned=${N_MAPPED}/${N_TOTAL}" >&2
