#!/usr/bin/env bash
# Local benchmark for long-read alignment (ONT/PacBio style)
# Usage: run_local_bench_long.sh <label> [genome_size=5M] [num_reads=50] [min_len=10000] [max_len=30000]
#
# NOTE: Requires minimap2 and samtools to be installed and accessible in PATH.
# Run on Reedling2 HPC for best results.
set -euo pipefail

LABEL="${1:?usage: run_local_bench_long.sh <label> [genome_size=5M] [num_reads=50] [min_len=10000] [max_len=30000]}"
GENOME_SIZE="${2:-5000000}"
NUM_READS="${3:-50}"
MIN_LEN="${4:-10000}"
MAX_LEN="${5:-30000}"

BENCH_DIR="${BENCH_DIR:-${TMPDIR:-/tmp}/phraya-bench-long}"
PHRAYA="${PHRAYA:-$(pwd)/target/release/phraya.exe}"
MINIMAP2="${MINIMAP2:-minimap2}"
SAMTOOLS="${SAMTOOLS:-samtools}"
PYTHON="${PYTHON3:-python3}"
GEN_LONG="$(dirname "$0")/gen_synthetic_long.py"
ACCURACY_PY="$(dirname "$0")/../slurm/utils/long_read_accuracy.py"

SEED="${SEED:-42}"
REF="$BENCH_DIR/$LABEL/ref.fa"
READS="$BENCH_DIR/$LABEL/reads.fq"
TRUTH="$BENCH_DIR/$LABEL/reads.fq.truth.tsv"
OUT_DIR="$BENCH_DIR/$LABEL/output"
PLAN="$OUT_DIR/plan.phrayaplan"

mkdir -p "$OUT_DIR"

echo "=== Long-Read Benchmark: $LABEL ===" >&2

echo "Generating long-read data..."
python3 "$GEN_LONG" \
    --genome-size "$GENOME_SIZE" \
    --num-reads "$NUM_READS" \
    --min-len "$MIN_LEN" \
    --max-len "$MAX_LEN" \
    --seed "$SEED" \
    --ref-out "$REF" \
    --reads-out "$READS"

echo "Building minimap2 index..."
MMI="${REF%.fasta}.mmi"
if [[ ! -f "$MMI" ]]; then
    "$MINIMAP2" -d "$MMI" "$REF"
fi

echo "=== Phase 1: Plan Creation ===" >&2
"$PHRAYA" plan --inputs "$READS" --reference "$REF" --output "$PLAN" 2>&1 | tee "$OUT_DIR/plan.log"
[[ -f "$PLAN" ]] || { echo "ERROR: Plan file not created" >&2; exit 1; }

echo "=== Phase 2: Align with phraya --strategy long-read ===" >&2
START_PHRAYA=$(date +%s.%N)
"$PHRAYA" align --strategy long-read --reference "$REF" --output "$OUT_DIR/phraya" "$PLAN" 2>&1 | tee "$OUT_DIR/phraya_align.log" || PHRAYA_FAILED=1
ELAPSED_PHRAYA=$(awk "BEGIN{printf \"%.3f\", \$(date +%s.%N) - $START_PHRAYA}")

echo "=== Phase 3: Align with minimap2 -ax map-ont ===" >&2
START_MM2=$(date +%s.%N)
"$MINIMAP2" -ax map-ont -t 1 "$MMI" "$READS" > "$OUT_DIR/minimap2.sam" 2>&1
ELAPSED_MM2=$(awk "BEGIN{printf \"%.3f\", \$(date +%s.%N) - $START_MM2}")

echo "=== Phase 4: Placement Accuracy ===" >&2

# Count total reads from truth file (subtract header line)
N_TOTAL=$(tail -n +2 "$TRUTH" | wc -l)

# Get mapped counts from minimap2
N_MAPPED_MM2=$("$SAMTOOLS" view -c -F4 "$OUT_DIR/minimap2.sam" 2>/dev/null || echo 0)

# Get phraya mapped count from queries file
N_MAPPED_PHRAYA=0
if [[ -f "$OUT_DIR/phraya/cross_space.phraya.queries" ]]; then
    N_MAPPED_PHRAYA=$("$PYTHON" -c "
import msgpack, zstandard
dctx = zstandard.ZstdDecompressor()
data = dctx.decompress(open('$OUT_DIR/phraya/cross_space.phraya.queries','rb').read())
result = msgpack.unpackb(data, raw=False)
print(sum(len(v) for v in result.values()))
" 2>/dev/null || echo 0)
fi

# Compute PA for each aligner independently
if [[ "${PHRAYA_FAILED:-0}" == "1" ]]; then
    PA_PHRAYA="FAILED"
else
    PA_PHRAYA=$("$PYTHON" "$ACCURACY_PY" \
        "$OUT_DIR/phraya/cross_space.phraya.queries" "$TRUTH" --tolerance 100 2>/dev/null | cut -f1 || echo "0.0000")
fi
PA_MM2=$(awk "BEGIN{if($N_TOTAL>0) printf \"%.4f\", $N_MAPPED_MM2/$N_TOTAL; else print \"0.0000\"}")

echo "=== Results ===" >&2
echo "phraya-longread: wall=${ELAPSED_PHRAYA}s, mapped=${N_MAPPED_PHRAYA}/${N_TOTAL}, PA=${PA_PHRAYA}"
echo "minimap2-ont:    wall=${ELAPSED_MM2}s, mapped=${N_MAPPED_MM2}/${N_TOTAL}, PA=${PA_MM2}"
echo ""
echo "Results file: $BENCH_DIR/$LABEL/results.txt"
{
    echo "label=$LABEL"
    echo "genome_size=$GENOME_SIZE"
    echo "num_reads=$NUM_READS"
    echo "read_len_range=$MIN_LEN-$MAX_LEN"
    echo "phraya_wall_seconds=$ELAPSED_PHRAYA"
    echo "phraya_mapped=$N_MAPPED_PHRAYA"
    echo "phraya_pa=$PA_PHRAYA"
    echo "minimap2_wall_seconds=$ELAPSED_MM2"
    echo "minimap2_mapped=$N_MAPPED_MM2"
    echo "minimap2_pa=$PA_MM2"
} > "$OUT_DIR/results.txt"
