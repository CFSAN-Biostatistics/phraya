#!/bin/bash
# minimap2 PacBio CLR long-read alignment wrapper (throughput baseline — SAM output)
# Usage: minimap2-pb.sh <ref.fasta> <reads.fq> <out_dir> <threads>
#
# Uses minimap2's -ax map-pb preset for PacBio CLR reads (~10-15% error).
# For ONT R10.4 (~15% error), use -x map-ont instead.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config/global.env"

REF=$1; READS=$2; OUT_DIR=$3; THREADS=$4

for f in "$REF" "$READS"; do
    [[ -f "$f" ]] || { echo "ERROR: not found: $f" >&2; exit 1; }
done

# Build .mmi index (flock-protected, one-time per reference)
INDEX="${REF%.fasta}.mmi"
if [[ ! -f "$INDEX" ]]; then
    (flock -x 200; [[ -f "$INDEX" ]] || $MINIMAP2_BIN -d "$INDEX" "$REF") 200>"${REF}.mmi.lock"
fi

PYTHON="${PYTHON3_BIN:-python3}"
MEASURE="$SCRIPT_DIR/utils/measure_rss.py"
START=$(now_s)
"$PYTHON" "$MEASURE" "$OUT_DIR/time_verbose.txt" -- \
    bash -c "$MINIMAP2_BIN -ax map-pb -t $THREADS $INDEX $READS > $OUT_DIR/alignment.sam 2>$OUT_DIR/minimap2-pb.log"
ELAPSED=$(elapsed_s "$START")

PEAK_RSS_KB=$(grep 'Maximum resident' "$OUT_DIR/time_verbose.txt" | grep -oP '\d+' | tail -1)
PEAK_RSS_GB=$(awk "BEGIN{printf \"%.3f\", ${PEAK_RSS_KB:-0}/1048576}")

N_TOTAL=$($SAMTOOLS_BIN view -c "$OUT_DIR/alignment.sam" 2>/dev/null || echo 0)
N_MAPPED=$($SAMTOOLS_BIN view -c -F4 "$OUT_DIR/alignment.sam" 2>/dev/null || echo 0)
N_UNMAPPED=$(( N_TOTAL - N_MAPPED ))
UNALIGNED_FRAC=$(awk "BEGIN{if($N_TOTAL>0) printf \"%.4f\", $N_UNMAPPED/$N_TOTAL; else print \"0.0000\"}")

cat > "$OUT_DIR/timing.txt" <<EOF
wall_seconds=$ELAPSED
threads=$THREADS
aligner=minimap2-pb
peak_rss_gb=${PEAK_RSS_GB}
total_reads=${N_TOTAL}
n_aligned=${N_MAPPED}
n_unmapped=${N_UNMAPPED}
unaligned_frac=${UNALIGNED_FRAC}
