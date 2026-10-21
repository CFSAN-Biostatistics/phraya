#!/bin/bash
# Phraya --strategy long-read wrapper (ONT/PacBio CLR)
# Usage: phraya-longread.sh <ref.fasta> <reads.fq> <out_dir> <threads>
#
# Uses the minimap2-style seed→chain→align long-read path in phraya-align/src/executor.rs:
# minimizer seeding → chain-to-locus dedup → exact-segment gap-fill (global Myers fills
# between seeds, fitting WFA end-extensions) with chain cap and end caps.
#
# Timing contract (same as _phraya_run.sh / minimap2-ont.sh):
#   wall_seconds = `phraya align` only. `phraya plan` (read + reference sketching) is the
#   analogue of minimap2's index build, which minimap2-ont.sh also excludes; it is recorded
#   separately as plan_seconds so the write-up can report both.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$SCRIPT_DIR/config/global.env"

REF=$1; READS=$2; OUT_DIR=$3; THREADS=$4

for f in "$REF" "$READS"; do
    [[ -f "$f" ]] || { echo "ERROR: not found: $f" >&2; exit 1; }
done

PHRAYA="${PHRAYA_ROOT}/target/release/phraya"
[[ -f "$PHRAYA" ]] || { echo "ERROR: phraya binary not found at $PHRAYA" >&2; exit 1; }

mkdir -p "$OUT_DIR"
PYTHON="${PYTHON3_BIN:-python3}"
MEASURE="$SCRIPT_DIR/utils/measure_rss.py"

echo ">>> plan" >&2
PLAN_FILE="$OUT_DIR/plan.phrayaplan"
PLAN_START=$(now_s)
RAYON_NUM_THREADS="$THREADS" "$PHRAYA" plan --inputs "$READS" --reference "$REF" --output "$PLAN_FILE" >"$OUT_DIR/plan.log" 2>&1 \
    || { echo "ERROR: phraya plan failed" >&2; tail -20 "$OUT_DIR/plan.log" >&2; exit 1; }
PLAN_SECONDS=$(elapsed_s "$PLAN_START")
[[ -f "$PLAN_FILE" ]] || { echo "ERROR: Plan file not created" >&2; exit 1; }

# Input size, from the (uncompressed or gzipped) single-end FASTQ.
TOTAL_READS=$(( $(zcat -f "$READS" | wc -l) / 4 ))
TOTAL_BASES=$(zcat -f "$READS" | awk 'NR%4==2{n+=length($0)} END{print n+0}')

echo ">>> align" >&2
START_SECS=$(now_s)
# `set -e` is active; capture the status explicitly so a failure still reports align.log.
ALIGN_EXIT=0
"$PYTHON" "$MEASURE" "$OUT_DIR/time_verbose.txt" -- \
    bash -c "RAYON_NUM_THREADS=$THREADS \"$PHRAYA\" align --strategy long-read --reference \"$REF\" --output \"$OUT_DIR\" \"$PLAN_FILE\" >\"$OUT_DIR/align.log\" 2>&1" \
    || ALIGN_EXIT=$?
ELAPSED=$(elapsed_s "$START_SECS")

if [[ $ALIGN_EXIT -ne 0 ]]; then
    echo "ERROR: phraya align failed (exit $ALIGN_EXIT)" >&2
    tail -20 "$OUT_DIR/align.log" >&2
    exit "$ALIGN_EXIT"
fi

PEAK_RSS_KB=$(grep 'Maximum resident' "$OUT_DIR/time_verbose.txt" | grep -oP '\d+' | tail -1)
PEAK_RSS_GB=$(awk "BEGIN{printf \"%.3f\", ${PEAK_RSS_KB:-0}/1048576}")

# n_aligned = reads with >=1 retained placement (see count_phraya_aligned.py, issue #193).
QUERIES_FILE="$OUT_DIR/cross_space.phraya.queries"
[[ -f "$QUERIES_FILE" ]] || { echo "ERROR: $QUERIES_FILE not written by phraya align" >&2; exit 1; }
N_ALIGNED=$("$PYTHON" "$SCRIPT_DIR/utils/count_phraya_aligned.py" "$QUERIES_FILE")

N_UNALIGNED=$(( TOTAL_READS > N_ALIGNED ? TOTAL_READS - N_ALIGNED : 0 ))
UNALIGNED_FRAC=$(awk "BEGIN{if($TOTAL_READS>0) printf \"%.4f\", $N_UNALIGNED/$TOTAL_READS; else print \"0.0000\"}")

cat > "$OUT_DIR/timing.txt" <<TIMING_EOF
wall_seconds=${ELAPSED}
plan_seconds=${PLAN_SECONDS}
threads=${THREADS}
aligner=phraya-longread
peak_rss_gb=${PEAK_RSS_GB}
total_reads=${TOTAL_READS}
total_bases=${TOTAL_BASES}
n_aligned=${N_ALIGNED}
n_unaligned=${N_UNALIGNED}
unaligned_frac=${UNALIGNED_FRAC}
TIMING_EOF

echo "Done: align ${ELAPSED}s (plan ${PLAN_SECONDS}s), RSS=${PEAK_RSS_GB}GB, aligned=${N_ALIGNED}/${TOTAL_READS}" >&2
