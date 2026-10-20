#!/bin/bash
# Phraya --strategy long-read wrapper (ONT/PacBio CLR)
# Usage: phraya-longread.sh <ref.fasta> <reads.fq> <out_dir> <threads>
#
# NOTE: This uses unbanded WFA via the standard wfa_extend() path.
# The 512 MB wavefront memory cap (default_max_s_cap) causes alignment
# abandonment for reads with edit distance >~1% divergence (~700 edits on 50kb reads).
# The LongReadAligner module in long_read.rs provides correct chunked logic but
# is NOT wired into the alignment pipeline.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$SCRIPT_DIR/config/global.env"

REF=$1; READS=$2; OUT_DIR=$3; THREADS=$4

for f in "$REF" "$READS"; do
    [[ -f "$f" ]] || { echo "ERROR: not found: $f" >&2; exit 1; }
done

PHRAYA="${PHRAYA_ROOT}/target/release/phraya"
[[ -f "$PHRAYA" ]] || PHRAYA="${PHRAYA_ROOT}/target/release/phraya.exe"
[[ -f "$PHRAYA" ]] || { echo "ERROR: phraya binary not found at ${PHRAYA_ROOT}/target/release/phraya" >&2; exit 1; }

echo "=== Phraya Alignment (strategy=long-read) ===" >&2
echo "Reference: $REF" >&2
echo "Reads: $READS" >&2
echo "Threads: $THREADS" >&2

PYTHON="${PYTHON3_BIN:-python3}"
MEASURE="$SCRIPT_DIR/utils/measure_rss.py"

# Ensure output directory exists
mkdir -p "$OUT_DIR"

# Build plan first (required for align --reference)
PLAN_FILE="$OUT_DIR/plan.phrayaplan"
"$PHRAYA" plan --inputs "$READS" --reference "$REF" --output "$PLAN_FILE" 2>&1 | tee "$OUT_DIR/plan.log"
[[ -f "$PLAN_FILE" ]] || { echo "ERROR: Plan file not created" >&2; exit 1; }

# Run alignment with reference-palette mode
START=$(date +%s.%N)
"$PYTHON" "$MEASURE" "$OUT_DIR/time_verbose.txt" -- \
    bash -c "RAYON_NUM_THREADS=$THREADS '$PHRAYA' align --strategy long-read --reference '$REF' --output '$OUT_DIR' '$PLAN_FILE' >'$OUT_DIR/align.log' 2>&1"
ELAPSED=$(awk -v s="$START" -v e="$(date +%s.%N)" 'BEGIN{printf "%.3f", e - s}')

# Parse timing output
PEAK_RSS_KB=$(grep 'Maximum resident' "$OUT_DIR/time_verbose.txt" 2>/dev/null | grep -oP '\d+' | tail -1 || echo "0")
PEAK_RSS_GB=$(awk "BEGIN{printf \"%.3f\", ${PEAK_RSS_KB:-0}/1048576}")

# Count reads with >=1 placement (and total reads) from the queries sidecar.
# `align --output $OUT_DIR` writes the sidecar directly under $OUT_DIR.
QUERIES="$OUT_DIR/cross_space.phraya.queries"
[[ -f "$QUERIES" ]] || { echo "ERROR: queries sidecar not found: $QUERIES" >&2; exit 1; }
N_ALIGNED=$("$PYTHON" "$SCRIPT_DIR/utils/count_phraya_aligned.py" "$QUERIES")
N_TOTAL=$(awk 'NR%4==1' "$READS" | wc -l)
UNALIGNED_FRAC=$(awk -v t="$N_TOTAL" -v a="$N_ALIGNED" 'BEGIN{if(t>0) printf "%.4f", (t-a)/t; else print "0.0000"}')

cat > "$OUT_DIR/timing.txt" <<EOF
wall_seconds=$ELAPSED
threads=$THREADS
aligner=phraya-longread
peak_rss_gb=${PEAK_RSS_GB}
total_reads=${N_TOTAL}
n_aligned=${N_ALIGNED}
unaligned_frac=${UNALIGNED_FRAC}
EOF

echo "Done: ${ELAPSED}s, RSS=${PEAK_RSS_GB}GB, aligned=${N_ALIGNED}/${N_TOTAL} (unaligned=${UNALIGNED_FRAC})" >&2