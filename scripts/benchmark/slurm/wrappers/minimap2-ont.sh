#!/bin/bash
# minimap2 ONT long-read alignment wrapper (throughput baseline — SAM output)
# Usage: minimap2-ont.sh <ref.fasta> <reads.fq.gz> <out_dir> <threads>
#
# Uses minimap2's -ax map-ont preset for ONT R10.4 reads (~15% error).
# For PacBio CLR (~10% error), use -ax map-pb; for HiFi (~1% error), use -ax map-pb-fast.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$SCRIPT_DIR/config/global.env"

REF=$1; READS=$2; OUT_DIR=$3; THREADS=$4

for f in "$REF" "$READS"; do
    [[ -f "$f" ]] || { echo "ERROR: not found: $f" >&2; exit 1; }
done

# Build a preset-specific .mmi index (flock-protected, one-time per reference). The preset
# MUST be passed at index time: minimap2 overrides -k/-w with the prebuilt index's values
# ("Indexing parameters (-k, -w or -H) overridden by parameters used in the prebuilt index"),
# so an index built without -x silently runs this wrapper with default k/w.
INDEX="${REF%.fasta}.map-ont.mmi"
# Build to a temp name and rename: a concurrent task that sees $INDEX must never read a
# half-written file (observed: a partial index loaded as "#seq: 0" gave PA 0.0).
if [[ ! -f "$INDEX" ]]; then
    (flock -x 200; [[ -f "$INDEX" ]] || { $MINIMAP2_BIN -x map-ont -d "$INDEX.tmp.$$" "$REF" && mv "$INDEX.tmp.$$" "$INDEX"; }) 200>"${INDEX}.lock"
fi

PYTHON="${PYTHON3_BIN:-python3}"
MEASURE="$SCRIPT_DIR/utils/measure_rss.py"
START=$(now_s)
"$PYTHON" "$MEASURE" "$OUT_DIR/time_verbose.txt" -- \
    bash -c "$MINIMAP2_BIN -ax map-ont -t $THREADS $INDEX $READS > $OUT_DIR/alignment.sam 2>$OUT_DIR/minimap2-ont.log"
ELAPSED=$(elapsed_s "$START")

PEAK_RSS_KB=$(grep 'Maximum resident' "$OUT_DIR/time_verbose.txt" | grep -oP '\d+' | tail -1)
PEAK_RSS_GB=$(awk "BEGIN{printf \"%.3f\", ${PEAK_RSS_KB:-0}/1048576}")

# Primary records only (-F 0x900): long reads routinely produce secondary/supplementary
# records, which would inflate total_reads and skew unaligned_frac.
N_TOTAL=$($SAMTOOLS_BIN view -c -F 0x900 "$OUT_DIR/alignment.sam" 2>/dev/null || echo 0)
N_MAPPED=$($SAMTOOLS_BIN view -c -F 0x904 "$OUT_DIR/alignment.sam" 2>/dev/null || echo 0)
N_UNMAPPED=$(( N_TOTAL - N_MAPPED ))
UNALIGNED_FRAC=$(awk "BEGIN{if($N_TOTAL>0) printf \"%.4f\", $N_UNMAPPED/$N_TOTAL; else print \"0.0000\"}")
TOTAL_BASES=$(zcat -f "$READS" | awk 'NR%4==2{n+=length($0)} END{print n+0}')

cat > "$OUT_DIR/timing.txt" <<EOF
wall_seconds=$ELAPSED
threads=$THREADS
aligner=minimap2-ont
peak_rss_gb=${PEAK_RSS_GB}
total_reads=${N_TOTAL}
total_bases=${TOTAL_BASES}
n_aligned=${N_MAPPED}
n_unmapped=${N_UNMAPPED}
unaligned_frac=${UNALIGNED_FRAC}
EOF

echo "Elapsed: ${ELAPSED}s, RSS=${PEAK_RSS_GB}GB, aligned=${N_MAPPED}/${N_TOTAL}" >&2