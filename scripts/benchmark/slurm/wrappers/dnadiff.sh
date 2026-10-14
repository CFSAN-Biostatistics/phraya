#!/bin/bash
# dnadiff (MUMmer 4) wrapper: the CSP2 baseline for assembly-pair benchmarking.
# Usage: dnadiff.sh <ref.fasta> <qry.fasta> <out_dir> <threads>
#
# Runs exactly what CSP2's runMUMmer runs (`dnadiff -p PREFIX REF QRY`; dnadiff has no
# thread option, so <threads> is recorded but unused). Output: PREFIX.{report,snps,1coords,...}
# in <out_dir>, consumed by utils/snpdiffs_concordance.py --dnadiff <out_dir>/PREFIX.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$SCRIPT_DIR/config/global.env"

REF=$1; QRY=$2; OUT_DIR=$3; THREADS=$4

for f in "$REF" "$QRY"; do
    [[ -f "$f" ]] || { echo "ERROR: not found: $f" >&2; exit 1; }
done

# Ensure output directory exists
mkdir -p "$OUT_DIR"

# Find dnadiff binary - check PATH, MUMMER environment, or micromamba
DNADIFF_BIN="${DNADIFF_BIN:-dnadiff}"
if ! command -v "$DNADIFF_BIN" &> /dev/null; then
    DNADIFF_BIN="/nfs/software/apps/micromamba/1.5.8/envs/mummer-v4.0.0/bin/dnadiff"
    if [[ ! -f "$DNADIFF_BIN" ]]; then
        DNADIFF_BIN="/usr/local/bin/dnadiff"
    fi
fi

echo "=== dnadiff (MUMmer) Assembly Comparison ===" >&2
echo "Reference: $REF" >&2
echo "Query: $QRY" >&2
echo "Threads: $THREADS" >&2
echo "Output: $OUT_DIR" >&2

PYTHON="${PYTHON3_BIN:-python3}"
MEASURE="$SCRIPT_DIR/utils/measure_rss.py"

# Build prefix for dnadiff output (basename of reference)
REF_BASENAME=$(basename "$REF")
REF_BASENAME="${REF_BASENAME%.*}"  # Remove extension
PREFIX="$REF_BASENAME"

START=$(now_s)
"$PYTHON" "$MEASURE" "$OUT_DIR/time_verbose.txt" -- \
    bash -c "$DNADIFF_BIN -p '$OUT_DIR/$PREFIX' '$REF' '$QRY' > '$OUT_DIR/dnadiff.log' 2>&1"
ELAPSED=$(elapsed_s "$START")

# Check if dnadiff produced output
if [[ ! -f "$OUT_DIR/$PREFIX.report" ]]; then
    echo "ERROR: dnadiff did not produce expected output files" >&2
    exit 1
fi

# Parse timing output
PEAK_RSS_KB=$(grep 'Maximum resident' "$OUT_DIR/time_verbose.txt" 2>/dev/null | grep -oP '\d+' | tail -1 || echo "0")
PEAK_RSS_GB=$(awk "BEGIN{printf \"%.3f\", ${PEAK_RSS_KB:-0}/1048576}")

# Extract summary metrics from .report (columns: Measure REF QRY; first occurrence = 1-to-1
# section, as compileMUMmer.py reads it). Values look like 2687317(90.9964%).
REF_TOTAL_BASES=$(grep '^TotalBases' "$OUT_DIR/$PREFIX.report" | head -1 | awk '{print $2}' | grep -oP '^\d+')
REF_ALIGNED_BASES=$(grep '^AlignedBases' "$OUT_DIR/$PREFIX.report" | head -1 | awk '{print $2}' | grep -oP '^\d+')
REF_TOTAL_SNPS=$(grep '^TotalSNPs' "$OUT_DIR/$PREFIX.report" | head -1 | awk '{print $2}')
REF_TOTAL_INDELS=$(grep '^TotalIndels' "$OUT_DIR/$PREFIX.report" | head -1 | awk '{print $2}')

# Compute percentages
if [[ -n "$REF_TOTAL_BASES" && "$REF_TOTAL_BASES" -gt 0 ]]; then
    REF_BREADTH_PCT=$(awk "BEGIN{printf \"%.2f\", ${REF_ALIGNED_BASES:-0}/${REF_TOTAL_BASES}*100}")
else
    REF_BREADTH_PCT="0.00"
fi

cat > "$OUT_DIR/timing.txt" <<EOF
wall_seconds=$ELAPSED
threads=$THREADS
aligner=dnadiff
peak_rss_gb=${PEAK_RSS_GB}
total_bases=${REF_TOTAL_BASES:-0}
aligned_bases=${REF_ALIGNED_BASES:-0}
snps=${REF_TOTAL_SNPS:-0}
indels=${REF_TOTAL_INDELS:-0}
ref_breadth_pct=${REF_BREADTH_PCT}
EOF

echo "Done: ${ELAPSED}s, RSS=${PEAK_RSS_GB}GB, ${REF_TOTAL_SNPS:-0} SNPs, ${REF_TOTAL_INDELS:-0} indels" >&2