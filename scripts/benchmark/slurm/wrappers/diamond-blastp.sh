#!/bin/bash
# DIAMOND blastp wrapper (protein throughput baseline — tabular output, ADR-0013)
# Usage: diamond-blastp.sh <proteome.fasta> <query.fasta> <out_dir> <threads>
#
# 4-arg contract, not the DNA harness's 5-arg <ref> <reads_1> <reads_2> <out_dir>
# <threads> — protein search is single-FASTA-vs-single-FASTA, no read pairing.
# See BENCHMARK_EXPANSION.md for why the wrapper contract differs by alphabet.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$SCRIPT_DIR/config/global.env"

PROTEOME=$1; QUERY=$2; OUT_DIR=$3; THREADS=$4

for f in "$PROTEOME" "$QUERY"; do
    [[ -f "$f" ]] || { echo "ERROR: not found: $f" >&2; exit 1; }
done

# DIAMOND tool path (check explicit path, then PATH)
if [[ -n "${DIAMOND_BIN:-}" ]] && [[ -x "$DIAMOND_BIN" ]]; then
    DIAMOND="$DIAMOND_BIN"
elif command -v diamond &>/dev/null; then
    DIAMOND="diamond"
else
    echo "ERROR: diamond not found. Set DIAMOND_BIN or ensure diamond is in PATH." >&2
    echo "Install: conda install -c bioconda diamond" >&2
    exit 1
fi

# Build DIAMOND protein DB (flock-protected, one-time per proteome)
DB="${PROTEOME}.dmnd"
if [[ ! -f "$DB" ]]; then
    (flock -x 200
     [[ -f "$DB" ]] || $DIAMOND makedb --in "$PROTEOME" -d "${PROTEOME}" --threads "$THREADS" 2>"$OUT_DIR/diamond-makedb.log"
    ) 200>"${PROTEOME}.dmnd_index.lock"
fi

PYTHON="${PYTHON3_BIN:-python3}"
MEASURE="$SCRIPT_DIR/utils/measure_rss.py"
OUT_TSV="$OUT_DIR/alignment.tsv"
START=$SECONDS
"$PYTHON" "$MEASURE" "$OUT_DIR/time_verbose.txt" -- \
    bash -c "$DIAMOND blastp -d '$DB' -q '$QUERY' -o '$OUT_TSV' --threads $THREADS \
        -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
        2>$OUT_DIR/diamond.log"
ELAPSED=$((SECONDS - START))

PEAK_RSS_KB=$(grep 'Maximum resident' "$OUT_DIR/time_verbose.txt" | grep -oP '\d+' | tail -1)
PEAK_RSS_GB=$(awk "BEGIN{printf \"%.3f\", ${PEAK_RSS_KB:-0}/1048576}")

# Query count from FASTA (count of '>' header lines); aligned = distinct qseqid
# in the tabular output (DIAMOND reports 0+ hit rows per query, no row at all
# for queries with no hit above its default e-value threshold).
N_TOTAL=$(grep -c '^>' "$QUERY" || echo 0)
N_MAPPED=$(cut -f1 "$OUT_TSV" 2>/dev/null | sort -u | wc -l || echo 0)
N_UNMAPPED=$(( N_TOTAL > N_MAPPED ? N_TOTAL - N_MAPPED : 0 ))
UNALIGNED_FRAC=$(awk "BEGIN{if($N_TOTAL>0) printf \"%.4f\", $N_UNMAPPED/$N_TOTAL; else print \"0.0000\"}")

cat > "$OUT_DIR/timing.txt" <<EOF
wall_seconds=$ELAPSED
threads=$THREADS
aligner=diamond-blastp
peak_rss_gb=${PEAK_RSS_GB}
total_reads=${N_TOTAL}
n_aligned=${N_MAPPED}
n_unaligned=${N_UNMAPPED}
unaligned_frac=${UNALIGNED_FRAC}
EOF

echo "Elapsed: ${ELAPSED}s, RSS=${PEAK_RSS_GB}GB, aligned=${N_MAPPED}/${N_TOTAL}" >&2
