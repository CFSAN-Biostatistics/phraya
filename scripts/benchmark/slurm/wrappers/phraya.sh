#!/bin/bash
# Phraya aligner wrapper with timing capture
# Usage: phraya.sh <ref.fasta> <reads_1.fq.gz> <reads_2.fq.gz> <out_dir> <threads>

set -euo pipefail

# Load environment
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$SCRIPT_DIR/config/global.env"

REF=$1
READS_1=$2
READS_2=$3
OUT_DIR=$4
THREADS=$5

# Verify inputs exist
for f in "$REF" "$READS_1" "$READS_2"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Input file not found: $f" >&2
        exit 1
    fi
done

# Ensure phraya binary is available
PHRAYA="${PHRAYA_ROOT}/target/release/phraya"
if [[ ! -f "$PHRAYA" ]]; then
    echo "ERROR: phraya binary not found at $PHRAYA" >&2
    exit 1
fi

echo "=== Phraya Alignment ==="
echo "Reference: $REF"
echo "Reads: $READS_1, $READS_2"
echo "Output: $OUT_DIR"
echo "Threads: $THREADS"
echo

# Step 1: Create plan
PLAN_FILE="$OUT_DIR/plan.phrayaplan"
echo "Step 1/2: Creating alignment plan..."
"$PHRAYA" plan \
    --inputs "$READS_1" \
    --inputs "$READS_2" \
    --reference "$REF" \
    --output "$PLAN_FILE" \
    --batch-to 1 \
    --batch-output-pattern "$OUT_DIR/alignment.phraya"

if [[ ! -f "$PLAN_FILE" ]]; then
    echo "ERROR: Plan file not created" >&2
    exit 1
fi
echo "Plan created: $PLAN_FILE"
echo

# Step 2: Align with timing
echo "Step 2/2: Running alignment (worker 0)..."
START_SECS=$SECONDS

RAYON_NUM_THREADS="$THREADS" "$PHRAYA" align \
    --strategy balanced \
    --worker 0 \
    "$PLAN_FILE" \
    > "$OUT_DIR/align.log" 2>&1

ALIGN_EXIT=$?
ELAPSED=$((SECONDS - START_SECS))

# Write timing information
cat > "$OUT_DIR/timing.txt" <<TIMING_EOF
Phraya Alignment Timing
Command: phraya align --strategy balanced --worker 0 $PLAN_FILE
Threads: $THREADS (RAYON_NUM_THREADS)
Exit code: $ALIGN_EXIT
wall_seconds=${ELAPSED}
TIMING_EOF

if [[ $ALIGN_EXIT -ne 0 ]]; then
    echo "ERROR: Alignment failed with exit code $ALIGN_EXIT" >&2
    echo "Log: $OUT_DIR/align.log" >&2
    tail -20 "$OUT_DIR/align.log" >&2
    exit $ALIGN_EXIT
fi

echo "Alignment complete!"
echo "Output: $OUT_DIR/alignment.phraya"
ls -lh "$OUT_DIR"/alignment.phraya* 2>/dev/null || echo "(No output files found)"

# Note: phraya outputs .phraya format, not SAM
# For this benchmark we're measuring alignment time, not PA computation
