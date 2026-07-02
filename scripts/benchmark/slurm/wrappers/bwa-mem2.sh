#!/bin/bash
# BWA alignment-only wrapper (throughput baseline — SAM output, no variant calling)
# Usage: bwa-mem2.sh <ref.fasta> <reads_1.fq.gz> <reads_2.fq.gz> <out_dir> <threads>
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$SCRIPT_DIR/config/global.env"

REF=$1; READS_1=$2; READS_2=$3; OUT_DIR=$4; THREADS=$5

for f in "$REF" "$READS_1" "$READS_2"; do
    [[ -f "$f" ]] || { echo "ERROR: not found: $f" >&2; exit 1; }
done

# Build BWA index (flock-protected, one-time per reference)
if [[ ! -f "${REF}.bwt" ]]; then
    (flock -x 200; [[ -f "${REF}.bwt" ]] || $BWA_BIN index "$REF") 200>"${REF}.bwa_index.lock"
fi

START=$SECONDS
$BWA_BIN mem -t "$THREADS" "$REF" "$READS_1" "$READS_2" > "$OUT_DIR/alignment.sam"
ELAPSED=$((SECONDS - START))

echo "wall_seconds=$ELAPSED threads=$THREADS aligner=bwa-mem" > "$OUT_DIR/timing.txt"
echo "Elapsed: ${ELAPSED}s" >&2
