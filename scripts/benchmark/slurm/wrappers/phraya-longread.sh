#!/bin/bash
# Phraya --strategy long-read wrapper (ONT/PacBio CLR)
# Usage: phraya-longread.sh <ref.fasta> <reads.fq> <out_dir> <threads>
#
# Uses the minimap2-style seed→chain→align long-read path in phraya-align/src/executor.rs:
# minimizer seeding → chain-to-locus dedup → exact-segment gap-fill (global Myers fills
# between seeds, fitting WFA end-extensions) with chain cap and end caps.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$SCRIPT_DIR/config/global.env"

REF=$1; READS=$2; OUT_DIR=$3; THREADS=$4

for f in "$REF" "$READS"; do
    [[ -f "$f" ]] || { echo "ERROR: not found: $f" >&2; exit 1; }
done

PHRAYA="${PHRAYA_ROOT}/target/release/phraya"
[[ -f "$PHRAYA" ]] || PHRAYA="${PHRAYA_ROOT}/target/release/phraya.exe"

mkdir -p "$OUT_DIR"

echo ">>> plan" >&2
"$PHRAYA" plan --inputs "$READS" --reference "$REF" --output "$OUT_DIR/plan.phrayaplan" 2>&1 | tee "$OUT_DIR/plan.log"

echo ">>> align" >&2
"$PHRAYA" align --strategy long-read --reference "$REF" --output "$OUT_DIR" "$OUT_DIR/plan.phrayaplan" 2>&1 | tee "$OUT_DIR/align.log"

echo ">>> done. outputs in $OUT_DIR" >&2
