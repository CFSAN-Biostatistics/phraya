#!/bin/bash
# Shared phraya runner — called by phraya.sh, phraya-sensitive.sh, phraya-fast.sh
# Usage: _phraya_run.sh <strategy> <ref.fasta> <reads_1.fq.gz> <reads_2.fq.gz> <out_dir> <threads>
set -euo pipefail

STRATEGY=$1; REF=$2; READS_1=$3; READS_2=$4; OUT_DIR=$5; THREADS=$6

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$SCRIPT_DIR/config/global.env"

for f in "$REF" "$READS_1" "$READS_2"; do
    [[ -f "$f" ]] || { echo "ERROR: not found: $f" >&2; exit 1; }
done

PHRAYA="${PHRAYA_ROOT}/target/release/phraya"
[[ -f "$PHRAYA" ]] || { echo "ERROR: phraya binary not found at $PHRAYA" >&2; exit 1; }

echo "=== Phraya Alignment (strategy=$STRATEGY) ==="
echo "Reference: $REF"
echo "Reads: $READS_1, $READS_2"
echo "Threads: $THREADS"

PLAN_FILE="$OUT_DIR/plan.phrayaplan"
"$PHRAYA" plan \
    --inputs "$READS_1" \
    --inputs "$READS_2" \
    --reference "$REF" \
    --output "$PLAN_FILE" \
    --batch-to 1 \
    --batch-output-pattern "$OUT_DIR/alignment.phraya"

[[ -f "$PLAN_FILE" ]] || { echo "ERROR: Plan file not created" >&2; exit 1; }

# Count total reads in both FASTQ files for unaligned fraction
TOTAL_READS=$(( $(zcat "$READS_1" | wc -l) / 4 + $(zcat "$READS_2" | wc -l) / 4 )) || TOTAL_READS=0

START_SECS=$SECONDS

# measure_rss.py polls /proc/PID/status for peak RSS; phraya stdout+stderr → align.log
PYTHON="${PYTHON3_BIN:-python3}"
MEASURE="$SCRIPT_DIR/utils/measure_rss.py"
"$PYTHON" "$MEASURE" "$OUT_DIR/time_verbose.txt" -- \
    bash -c "RAYON_NUM_THREADS=$THREADS \"$PHRAYA\" align --strategy \"$STRATEGY\" --worker 0 \"$PLAN_FILE\" >\"$OUT_DIR/align.log\" 2>&1"

ALIGN_EXIT=$?
ELAPSED=$((SECONDS - START_SECS))

if [[ $ALIGN_EXIT -ne 0 ]]; then
    echo "ERROR: Alignment failed (exit $ALIGN_EXIT)" >&2
    tail -20 "$OUT_DIR/align.log" >&2
    exit $ALIGN_EXIT
fi

# Parse peak RSS from /usr/bin/time -v output (KB → GB)
PEAK_RSS_KB=$(grep 'Maximum resident' "$OUT_DIR/time_verbose.txt" | grep -oP '\d+' | tail -1)
PEAK_RSS_GB=$(awk "BEGIN{printf \"%.3f\", ${PEAK_RSS_KB:-0}/1048576}")

# Count aligned reads: observations in .phraya are per-position not per-read;
# use .phraya.queries which has one key per query/read that placed ≥1 alignment
QUERIES_FILE="$OUT_DIR/alignment.phraya.queries"
if [[ -f "$QUERIES_FILE" ]]; then
    PYTHON="${PYTHON3_BIN:-python3}"
    N_ALIGNED=$("$PYTHON" "$SCRIPT_DIR/utils/count_phraya_aligned.py" "$QUERIES_FILE" 2>/dev/null || echo 0)
else
    N_ALIGNED=0
fi

N_UNALIGNED=$(( TOTAL_READS > N_ALIGNED ? TOTAL_READS - N_ALIGNED : 0 ))
UNALIGNED_FRAC=$(awk "BEGIN{if($TOTAL_READS>0) printf \"%.4f\", $N_UNALIGNED/$TOTAL_READS; else print \"0.0000\"}")

cat > "$OUT_DIR/timing.txt" <<TIMING_EOF
wall_seconds=${ELAPSED}
threads=${THREADS}
aligner=phraya-${STRATEGY}
peak_rss_gb=${PEAK_RSS_GB}
total_reads=${TOTAL_READS}
n_aligned=${N_ALIGNED}
n_unaligned=${N_UNALIGNED}
unaligned_frac=${UNALIGNED_FRAC}
TIMING_EOF

echo "Done: ${ELAPSED}s, RSS=${PEAK_RSS_GB}GB, aligned=${N_ALIGNED}/${TOTAL_READS} (unaligned=${UNALIGNED_FRAC})" >&2
