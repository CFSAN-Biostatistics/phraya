#!/bin/bash
# NCBI BLASTP wrapper (protein accuracy reference — slow, maximally sensitive; the
# protein-mode analog of the DNA harness's `sensitive`/BWA role, ADR-0013)
# Usage: blastp.sh <proteome.fasta> <query.fasta> <out_dir> <threads>
#
# 4-arg contract — see diamond-blastp.sh for why protein wrappers drop the
# paired-reads slot the DNA 5-arg contract uses.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$SCRIPT_DIR/config/global.env"

PROTEOME=$1; QUERY=$2; OUT_DIR=$3; THREADS=$4

for f in "$PROTEOME" "$QUERY"; do
    [[ -f "$f" ]] || { echo "ERROR: not found: $f" >&2; exit 1; }
done

# BLAST+ tool paths (check explicit path, then PATH)
if [[ -n "${BLASTP_BIN:-}" ]] && [[ -x "$BLASTP_BIN" ]]; then
    BLASTP="$BLASTP_BIN"
elif command -v blastp &>/dev/null; then
    BLASTP="blastp"
else
    echo "ERROR: blastp not found. Set BLASTP_BIN or ensure blastp is in PATH." >&2
    echo "Install: conda install -c bioconda blast" >&2
    exit 1
fi
if [[ -n "${MAKEBLASTDB_BIN:-}" ]] && [[ -x "$MAKEBLASTDB_BIN" ]]; then
    MAKEBLASTDB="$MAKEBLASTDB_BIN"
elif command -v makeblastdb &>/dev/null; then
    MAKEBLASTDB="makeblastdb"
else
    echo "ERROR: makeblastdb not found. Set MAKEBLASTDB_BIN or ensure it is in PATH." >&2
    exit 1
fi

# Build BLAST protein DB (flock-protected, one-time per proteome)
if [[ ! -f "${PROTEOME}.pin" ]]; then
    (flock -x 200
     [[ -f "${PROTEOME}.pin" ]] || $MAKEBLASTDB -in "$PROTEOME" -dbtype prot -out "$PROTEOME" 2>"$OUT_DIR/makeblastdb.log"
    ) 200>"${PROTEOME}.blastdb_index.lock"
fi

PYTHON="${PYTHON3_BIN:-python3}"
MEASURE="$SCRIPT_DIR/utils/measure_rss.py"
OUT_TSV="$OUT_DIR/alignment.tsv"
START=$SECONDS
"$PYTHON" "$MEASURE" "$OUT_DIR/time_verbose.txt" -- \
    bash -c "$BLASTP -db '$PROTEOME' -query '$QUERY' -out '$OUT_TSV' -num_threads $THREADS \
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
        2>$OUT_DIR/blastp.log"
ELAPSED=$((SECONDS - START))

PEAK_RSS_KB=$(grep 'Maximum resident' "$OUT_DIR/time_verbose.txt" | grep -oP '\d+' | tail -1)
PEAK_RSS_GB=$(awk "BEGIN{printf \"%.3f\", ${PEAK_RSS_KB:-0}/1048576}")

N_TOTAL=$(grep -c '^>' "$QUERY" || echo 0)
N_MAPPED=$(cut -f1 "$OUT_TSV" 2>/dev/null | sort -u | wc -l || echo 0)
N_UNMAPPED=$(( N_TOTAL > N_MAPPED ? N_TOTAL - N_MAPPED : 0 ))
UNALIGNED_FRAC=$(awk "BEGIN{if($N_TOTAL>0) printf \"%.4f\", $N_UNMAPPED/$N_TOTAL; else print \"0.0000\"}")

cat > "$OUT_DIR/timing.txt" <<EOF
wall_seconds=$ELAPSED
threads=$THREADS
aligner=blastp
peak_rss_gb=${PEAK_RSS_GB}
total_reads=${N_TOTAL}
n_aligned=${N_MAPPED}
n_unaligned=${N_UNMAPPED}
unaligned_frac=${UNALIGNED_FRAC}
EOF

echo "Elapsed: ${ELAPSED}s, RSS=${PEAK_RSS_GB}GB, aligned=${N_MAPPED}/${N_TOTAL}" >&2
