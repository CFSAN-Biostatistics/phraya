#!/bin/bash
# minimap2 aligner wrapper with timing capture
# Usage: minimap2.sh <ref.fasta> <reads_1.fq.gz> <reads_2.fq.gz> <out_dir> <threads>

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

# Use pre-built .mmi index if available, else use FASTA directly
INDEX="${REF%.fasta}.mmi"
if [[ -f "$INDEX" ]]; then
    REF_ARG="$INDEX"
else
    REF_ARG="$REF"
fi

# Run minimap2 with timing (-x sr = short read preset)
/usr/bin/time -v minimap2 -ax sr -t "$THREADS" "$REF_ARG" "$READS_1" "$READS_2" \
    > "$OUT_DIR/alignment.sam" 2> "$OUT_DIR/timing.txt"
