#!/bin/bash
# rammap alignment-only wrapper (pure Rust minimap2 clone with SIMD DP)
# Usage: rammap.sh <ref.fasta> <reads_1.fq.gz> <reads_2.fq.gz> <out_dir> <threads>
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$SCRIPT_DIR/config/global.env"

REF=$1; READS_1=$2; READS_2=$3; OUT_DIR=$4; THREADS=$5

for f in "$REF" "$READS_1" "$READS_2"; do
    [[ -f "$f" ]] || { echo "ERROR: not found: $f" >&2; exit 1; }
done

# rammap tool path (check explicit path, then PATH)
if [[ -n "${RAMMAP_BIN:-}" ]] && [[ -x "$RAMMAP_BIN" ]]; then
    RAMMAP="$RAMMAP_BIN"
elif command -v rammap &>/dev/null; then
    RAMMAP="rammap"
else
    echo "ERROR: rammap not found. Set RAMMAP_BIN or ensure rammap is in PATH." >&2
    echo "Install: cargo install --git https://github.com/jwanglab/rammap" >&2
    exit 1
fi

# rammap does not require separate index build (inline with alignment, like minimap2)
# No index pre-build step needed.

PYTHON="${PYTHON3_BIN:-python3}"
MEASURE="$SCRIPT_DIR/utils/measure_rss.py"

# Concatenate paired reads for rammap (it expects single input for paired-end in sr preset)
# Actually, rammap -x sr can take two read files like minimap2. Let's use that.
START=$SECONDS
"$PYTHON" "$MEASURE" "$OUT_DIR/time_verbose.txt" -- \
    bash -c "$RAMMAP -x sr -a -t $THREADS $REF $READS_1 $READS_2 > $OUT_DIR/alignment.sam 2>$OUT_DIR/rammap.log"
ELAPSED=$((SECONDS - START))

PEAK_RSS_KB=$(grep 'Maximum resident' "$OUT_DIR/time_verbose.txt" | grep -oP '\d+' | tail -1)
PEAK_RSS_GB=$(awk "BEGIN{printf \"%.3f\", ${PEAK_RSS_KB:-0}/1048576}")

# Unaligned fraction from SAM flagstat
N_TOTAL=$($SAMTOOLS_BIN view -c "$OUT_DIR/alignment.sam" 2>/dev/null || echo 0)
N_MAPPED=$($SAMTOOLS_BIN view -c -F4 "$OUT_DIR/alignment.sam" 2>/dev/null || echo 0)
N_UNMAPPED=$(( N_TOTAL - N_MAPPED ))
UNALIGNED_FRAC=$(awk "BEGIN{if($N_TOTAL>0) printf \"%.4f\", $N_UNMAPPED/$N_TOTAL; else print \"0.0000\"}")

cat > "$OUT_DIR/timing.txt" <<EOF
wall_seconds=$ELAPSED
threads=$THREADS
aligner=rammap
peak_rss_gb=${PEAK_RSS_GB}
total_reads=${N_TOTAL}
n_aligned=${N_MAPPED}
n_unaligned=${N_UNMAPPED}
unaligned_frac=${UNALIGNED_FRAC}
EOF

echo "Elapsed: ${ELAPSED}s, RSS=${PEAK_RSS_GB}GB, aligned=${N_MAPPED}/${N_TOTAL}" >&2
