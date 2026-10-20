#!/usr/bin/env bash
# Local benchmark for long-read alignment (ONT/PacBio style)
# Usage: run_local_bench_long.sh <label> [genome_size=5M] [num_reads=50] [min_len=10000] [max_len=30000]
#
# NOTE: Requires minimap2. If it is not already on PATH, the script loads it via
# the environment-modules system (Reedling2: MODULE_MINIMAP2=minimap2/2.28; list
# others with `module avail`). Override MINIMAP2 to use an explicit binary.
# Both aligners are scored against the generator's ground truth with the same
# script (slurm/utils/long_read_accuracy.py), tolerance TOLERANCE bp (default 100).
# RESCORE_ONLY=1 reuses existing alignment outputs and re-runs scoring only.
set -euo pipefail

LABEL="${1:?usage: run_local_bench_long.sh <label> [genome_size=5M] [num_reads=50] [min_len=10000] [max_len=30000]}"
GENOME_SIZE="${2:-5000000}"
NUM_READS="${3:-50}"
MIN_LEN="${4:-10000}"
MAX_LEN="${5:-30000}"

BENCH_DIR="${BENCH_DIR:-${TMPDIR:-/tmp}/phraya-bench-long}"
if [[ -z "${PHRAYA:-}" ]]; then
    PHRAYA="$(pwd)/target/release/phraya"
    [[ -x "$PHRAYA" ]] || PHRAYA="$PHRAYA.exe"
fi
MODULE_MINIMAP2="${MODULE_MINIMAP2:-minimap2/2.28}"
MINIMAP2="${MINIMAP2:-minimap2}"
TOLERANCE="${TOLERANCE:-100}"

# Tools on Reedling2 are provided via `module`; non-interactive shells (ssh
# command, batch) don't define it until /etc/profile is sourced.
load_tool() {  # load_tool <binary> <module>
    command -v "$1" >/dev/null 2>&1 && return 0
    set +u  # module/profile scripts reference unset variables
    type module >/dev/null 2>&1 || source /etc/profile >/dev/null 2>&1 || true
    type module >/dev/null 2>&1 && module load "$2" >/dev/null 2>&1 || true
    set -u
    command -v "$1" >/dev/null 2>&1 || { echo "ERROR: $1 not found and 'module load $2' failed (see 'module avail')" >&2; exit 1; }
}
load_tool "$MINIMAP2" "$MODULE_MINIMAP2"
PYTHON="${PYTHON3:-python3}"
GEN_LONG="$(dirname "$0")/gen_synthetic_long.py"
ACCURACY_PY="$(dirname "$0")/../slurm/utils/long_read_accuracy.py"

SEED="${SEED:-42}"
REF="$BENCH_DIR/$LABEL/ref.fa"
READS="$BENCH_DIR/$LABEL/reads.fq"
TRUTH="$BENCH_DIR/$LABEL/reads.fq.truth.tsv"
OUT_DIR="$BENCH_DIR/$LABEL/output"
PLAN="$OUT_DIR/plan.phrayaplan"

mkdir -p "$OUT_DIR"

echo "=== Long-Read Benchmark: $LABEL ===" >&2

PHRAYA_FAILED=0
if [[ "${RESCORE_ONLY:-0}" == "1" ]]; then
    # Reuse existing data + alignment outputs; only re-run scoring. Walls come from the prior results.txt.
    echo "=== RESCORE_ONLY: skipping data generation, plan, and alignment ===" >&2
    for f in "$TRUTH" "$OUT_DIR/results.txt" "$OUT_DIR/phraya/cross_space.phraya.queries" "$OUT_DIR/minimap2.sam"; do
        [[ -f "$f" ]] || { echo "ERROR: RESCORE_ONLY needs existing $f (run once without RESCORE_ONLY)" >&2; exit 1; }
    done
    ELAPSED_PHRAYA=$(sed -n 's/^phraya_wall_seconds=//p' "$OUT_DIR/results.txt")
    ELAPSED_MM2=$(sed -n 's/^minimap2_wall_seconds=//p' "$OUT_DIR/results.txt")
else
    echo "Generating long-read data..."
    python3 "$GEN_LONG" \
        --genome-size "$GENOME_SIZE" \
        --num-reads "$NUM_READS" \
        --min-len "$MIN_LEN" \
        --max-len "$MAX_LEN" \
        --seed "$SEED" \
        --ref-out "$REF" \
        --reads-out "$READS"

    # Always rebuild (~0.2s for 5 Mbp): the reference is regenerated each run, so a cached index could be stale.
    echo "Building minimap2 index..."
    MMI="$REF.mmi"
    "$MINIMAP2" -d "$MMI" "$REF" 2> "$OUT_DIR/minimap2_index.log"

    echo "=== Phase 1: Plan Creation ===" >&2
    "$PHRAYA" plan --inputs "$READS" --reference "$REF" --output "$PLAN" 2>&1 | tee "$OUT_DIR/plan.log"
    [[ -f "$PLAN" ]] || { echo "ERROR: Plan file not created" >&2; exit 1; }

    echo "=== Phase 2: Align with phraya --strategy long-read ===" >&2
    START_PHRAYA=$(date +%s.%N)
    # Capture phraya's own exit status, not tee's; don't let set -e abort before scoring minimap2.
    set +e
    "$PHRAYA" align --strategy long-read --reference "$REF" --output "$OUT_DIR/phraya" "$PLAN" 2>&1 | tee "$OUT_DIR/phraya_align.log"
    [[ "${PIPESTATUS[0]}" == "0" ]] || PHRAYA_FAILED=1
    set -e
    ELAPSED_PHRAYA=$(awk -v s="$START_PHRAYA" -v e="$(date +%s.%N)" 'BEGIN{printf "%.3f", e - s}')

    echo "=== Phase 3: Align with minimap2 -ax map-ont ===" >&2
    START_MM2=$(date +%s.%N)
    "$MINIMAP2" -ax map-ont -t 1 "$MMI" "$READS" > "$OUT_DIR/minimap2.sam" 2> "$OUT_DIR/minimap2.log"
    ELAPSED_MM2=$(awk -v s="$START_MM2" -v e="$(date +%s.%N)" 'BEGIN{printf "%.3f", e - s}')
fi

echo "=== Phase 4: Placement Accuracy (vs ground truth, tolerance ${TOLERANCE}bp) ===" >&2
# Both aligners are scored by the same script against the same truth. No stderr
# suppression or fallbacks: a scoring failure aborts the run (set -e).
if [[ "${PHRAYA_FAILED:-0}" == "1" ]]; then
    SCORE_PHRAYA="FAILED"
else
    SCORE_PHRAYA=$("$PYTHON" "$ACCURACY_PY" "$OUT_DIR/phraya/cross_space.phraya.queries" "$TRUTH" --tolerance "$TOLERANCE")
fi
SCORE_MM2=$("$PYTHON" "$ACCURACY_PY" "$OUT_DIR/minimap2.sam" "$TRUTH" --tolerance "$TOLERANCE")

# score fields: pa precision n_truth n_mapped n_correct n_correct_any
field() { echo "$1" | cut -f"$2"; }

echo "=== Results ===" >&2
printf '%-16s %10s %8s %10s %9s %9s %12s\n' aligner wall_s PA precision mapped correct correct_any
if [[ "$SCORE_PHRAYA" == "FAILED" ]]; then
    echo "phraya-longread  FAILED"
else
    printf '%-16s %10s %8s %10s %6s/%-3s %9s %12s\n' phraya-longread "$ELAPSED_PHRAYA" \
        "$(field "$SCORE_PHRAYA" 1)" "$(field "$SCORE_PHRAYA" 2)" "$(field "$SCORE_PHRAYA" 4)" \
        "$(field "$SCORE_PHRAYA" 3)" "$(field "$SCORE_PHRAYA" 5)" "$(field "$SCORE_PHRAYA" 6)"
fi
printf '%-16s %10s %8s %10s %6s/%-3s %9s %12s\n' minimap2-ont "$ELAPSED_MM2" \
    "$(field "$SCORE_MM2" 1)" "$(field "$SCORE_MM2" 2)" "$(field "$SCORE_MM2" 4)" \
    "$(field "$SCORE_MM2" 3)" "$(field "$SCORE_MM2" 5)" "$(field "$SCORE_MM2" 6)"
echo ""
echo "Results file: $OUT_DIR/results.txt"
{
    echo "label=$LABEL"
    echo "genome_size=$GENOME_SIZE"
    echo "num_reads=$NUM_READS"
    echo "read_len_range=$MIN_LEN-$MAX_LEN"
    echo "tolerance=$TOLERANCE"
    echo "phraya_wall_seconds=$ELAPSED_PHRAYA"
    echo "phraya_score=$SCORE_PHRAYA"
    echo "minimap2_wall_seconds=$ELAPSED_MM2"
    echo "minimap2_score=$SCORE_MM2"
    echo "score_fields=pa precision n_truth n_mapped n_correct n_correct_any"
} > "$OUT_DIR/results.txt"
