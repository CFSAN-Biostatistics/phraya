#!/bin/bash
# minimap2 alignment-only wrapper (throughput baseline — SAM output)
# Usage: minimap2.sh <ref.fasta> <reads_1.fq.gz> <reads_2.fq.gz> <out_dir> <threads>
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$SCRIPT_DIR/config/global.env"

REF=$1; READS_1=$2; READS_2=$3; OUT_DIR=$4; THREADS=$5

for f in "$REF" "$READS_1" "$READS_2"; do
    [[ -f "$f" ]] || { echo "ERROR: not found: $f" >&2; exit 1; }
done

# Build .mmi index (flock-protected, one-time per reference)
INDEX="${REF%.fasta}.mmi"
if [[ ! -f "$INDEX" ]]; then
    (flock -x 200; [[ -f "$INDEX" ]] || $MINIMAP2_BIN -d "$INDEX" "$REF") 200>"${REF}.mmi.lock"
fi

START=$SECONDS
$MINIMAP2_BIN -ax sr -t "$THREADS" "$INDEX" "$READS_1" "$READS_2" > "$OUT_DIR/alignment.sam"
ELAPSED=$((SECONDS - START))

echo "wall_seconds=$ELAPSED threads=$THREADS aligner=minimap2" > "$OUT_DIR/timing.txt"
echo "Elapsed: ${ELAPSED}s" >&2
