#!/bin/bash
# Shared phraya protein runner — called by phraya-protein.sh, phraya-protein-sensitive.sh
# Usage: _phraya_run_protein.sh <strategy> <proteome.fasta> <query.fasta> <out_dir> <threads>
#
# 4-arg contract (no paired reads) — see _phraya_run.sh for the DNA/5-arg analog.
# NOT YET RUNNABLE: `phraya plan --alphabet protein` doesn't exist until ADR-0013
# ships. Written to the target interface now; errors clearly at the `plan` step
# until then (see BENCHMARK_EXPANSION.md, "What ships now vs. what waits").
set -euo pipefail

STRATEGY=$1; PROTEOME=$2; QUERY=$3; OUT_DIR=$4; THREADS=$5

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$SCRIPT_DIR/config/global.env"

for f in "$PROTEOME" "$QUERY"; do
    [[ -f "$f" ]] || { echo "ERROR: not found: $f" >&2; exit 1; }
done

PHRAYA="${PHRAYA_ROOT}/target/release/phraya"
[[ -f "$PHRAYA" ]] || { echo "ERROR: phraya binary not found at $PHRAYA" >&2; exit 1; }

echo "=== Phraya Protein Alignment (strategy=$STRATEGY) ===" >&2
echo "Proteome: $PROTEOME" >&2
echo "Query: $QUERY" >&2
echo "Threads: $THREADS" >&2

PLAN_FILE="$OUT_DIR/plan.phrayaplan"
"$PHRAYA" plan \
    --inputs "$QUERY" \
    --reference "$PROTEOME" \
    --alphabet protein \
    --output "$PLAN_FILE" \
    --batch-to 1 \
    --batch-output-pattern "$OUT_DIR/alignment.phraya"

[[ -f "$PLAN_FILE" ]] || { echo "ERROR: Plan file not created" >&2; exit 1; }

TOTAL_QUERIES=$(grep -c '^>' "$QUERY" || echo 0)

START_SECS=$SECONDS
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

PEAK_RSS_KB=$(grep 'Maximum resident' "$OUT_DIR/time_verbose.txt" | grep -oP '\d+' | tail -1)
PEAK_RSS_GB=$(awk "BEGIN{printf \"%.3f\", ${PEAK_RSS_KB:-0}/1048576}")

QUERIES_FILE="$OUT_DIR/alignment.phraya.queries"
if [[ -f "$QUERIES_FILE" ]]; then
    PYTHON="${PYTHON3_BIN:-python3}"
    N_ALIGNED=$("$PYTHON" "$SCRIPT_DIR/utils/count_phraya_aligned.py" "$QUERIES_FILE" 2>/dev/null || echo 0)
else
    N_ALIGNED=0
fi

N_UNALIGNED=$(( TOTAL_QUERIES > N_ALIGNED ? TOTAL_QUERIES - N_ALIGNED : 0 ))
UNALIGNED_FRAC=$(awk "BEGIN{if($TOTAL_QUERIES>0) printf \"%.4f\", $N_UNALIGNED/$TOTAL_QUERIES; else print \"0.0000\"}")

cat > "$OUT_DIR/timing.txt" <<TIMING_EOF
wall_seconds=${ELAPSED}
threads=${THREADS}
aligner=phraya-protein-${STRATEGY}
peak_rss_gb=${PEAK_RSS_GB}
total_reads=${TOTAL_QUERIES}
n_aligned=${N_ALIGNED}
n_unaligned=${N_UNALIGNED}
unaligned_frac=${UNALIGNED_FRAC}
TIMING_EOF

echo "Done: ${ELAPSED}s, RSS=${PEAK_RSS_GB}GB, aligned=${N_ALIGNED}/${TOTAL_QUERIES} (unaligned=${UNALIGNED_FRAC})" >&2
