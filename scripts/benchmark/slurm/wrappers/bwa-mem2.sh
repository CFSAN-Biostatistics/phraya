#!/bin/bash
# BWA-MEM2 aligner wrapper with timing capture
# Usage: bwa-mem2.sh <ref.fasta> <reads_1.fq.gz> <reads_2.fq.gz> <out_dir> <threads>

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

# Run BWA-MEM2 with timing
/usr/bin/time -v bwa-mem2 mem -t "$THREADS" "$REF" "$READS_1" "$READS_2" \
    > "$OUT_DIR/alignment.sam" 2> "$OUT_DIR/timing.txt"
