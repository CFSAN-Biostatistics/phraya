#!/bin/bash
# phraya aligner wrapper with timing capture
# Usage: phraya.sh <ref.fasta> <reads_1.fq.gz> <reads_2.fq.gz> <out_dir> <threads>

set -euo pipefail

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
PHRAYA="${PHRAYA_ROOT:-$HOME/phraya}/target/release/phraya"
if [[ ! -f "$PHRAYA" ]]; then
    echo "ERROR: phraya binary not found at $PHRAYA" >&2
    echo "Run: cd $PHRAYA_ROOT && cargo build --release" >&2
    exit 1
fi

# Create plan (batch mode: single worker to process all reads)
PLAN_FILE="$OUT_DIR/plan.phrayaplan"
"$PHRAYA" plan \
    --inputs "$READS_1" "$READS_2" \
    --reference "$REF" \
    --output "$PLAN_FILE" \
    --batch-to 1 \
    --batch-output-pattern "$OUT_DIR/alignment.phraya"

# Align with timing (batch worker 0 processes all reads)
# Note: phraya doesn't have a -t/--threads flag yet, uses rayon default (all cores)
# For fair comparison, we should limit via RAYON_NUM_THREADS env var
RAYON_NUM_THREADS="$THREADS" /usr/bin/time -v "$PHRAYA" align "$PLAN_FILE" --worker 0 \
    2> "$OUT_DIR/timing.txt"

# phraya outputs .phraya format, not SAM
# For PA computation in aggregation, we'd need SAM output
# TODO: Add SAM export to phraya or convert .phraya → SAM
echo "WARNING: phraya outputs .phraya format, not SAM. PA computation may fail." >&2
