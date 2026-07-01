#!/usr/bin/env bash
# Local before/after benchmark for the `phraya align` hot path.
#
# Generates a seeded synthetic reference + reads (once, reused across runs), then
# times `phraya align` under /usr/bin/time -v, single-threaded for a clean
# algorithmic signal. Appends wall-clock + peak RSS to a results TSV tagged with a
# label, and leaves the output .phraya in place so runs can be diffed for
# correctness (a pure perf change must not alter output).
#
# Usage:
#   scripts/benchmark/local/run_local_bench.sh <label> [genome_size] [num_reads] [read_len] [divergence]
#
# Env:
#   BENCH_DIR   where data/outputs live (default: ${TMPDIR:-/tmp}/phraya-bench)
#   PHRAYA      phraya binary (default: target/release/phraya)
set -euo pipefail

LABEL="${1:?usage: run_local_bench.sh <label> [genome_size] [num_reads] [read_len] [divergence]}"
GENOME_SIZE="${2:-2000000}"
NUM_READS="${3:-20000}"
READ_LEN="${4:-150}"
DIVERGENCE="${5:-0.01}"
SEED="${SEED:-1}"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
BENCH_DIR="${BENCH_DIR:-${TMPDIR:-/tmp}/phraya-bench}"
PHRAYA="${PHRAYA:-$REPO_ROOT/target/release/phraya}"
GEN="$REPO_ROOT/scripts/benchmark/local/gen_synthetic.py"

# Data is keyed by scale params so identical params reuse identical inputs.
KEY="g${GENOME_SIZE}_n${NUM_READS}_l${READ_LEN}_d${DIVERGENCE}_s${SEED}"
DATA_DIR="$BENCH_DIR/$KEY"
REF="$DATA_DIR/ref.fa"
READS="$DATA_DIR/reads.fq"
PLAN="$DATA_DIR/plan.phrayaplan"
OUT="$DATA_DIR/${LABEL}.phraya"
RESULTS="$BENCH_DIR/results.tsv"

mkdir -p "$DATA_DIR"

if [[ ! -f "$PHRAYA" ]]; then
    echo "ERROR: phraya binary not found at $PHRAYA" >&2
    echo "  build with: RUSTFLAGS=\"-C target-cpu=native\" cargo build --release" >&2
    exit 1
fi

if [[ ! -f "$REF" || ! -f "$READS" ]]; then
    echo ">> generating synthetic data ($KEY)"
    python3 "$GEN" --genome-size "$GENOME_SIZE" --num-reads "$NUM_READS" \
        --read-len "$READ_LEN" --divergence "$DIVERGENCE" --seed "$SEED" \
        --ref-out "$REF" --reads-out "$READS"
fi

echo ">> planning"
"$PHRAYA" plan --inputs "$READS" --reference "$REF" --output "$PLAN" \
    --batch-to 1 --batch-output-pattern "$OUT" >/dev/null

echo ">> aligning (label=$LABEL, single-threaded)"
TIMER="$REPO_ROOT/scripts/benchmark/local/time_run.py"
MEASURE="$(RAYON_NUM_THREADS=1 python3 "$TIMER" "$PHRAYA" align "$PLAN" --worker 0)"
WALL="$(printf '%s' "$MEASURE" | cut -f1)"
RSS_MB="$(printf '%s' "$MEASURE" | cut -f2)"

if [[ ! -f "$RESULTS" ]]; then
    printf 'label\tkey\twall_s\tpeak_rss_mb\toutput\n' > "$RESULTS"
fi
printf '%s\t%s\t%s\t%s\t%s\n' "$LABEL" "$KEY" "$WALL" "$RSS_MB" "$OUT" >> "$RESULTS"

echo ">> $LABEL: wall=${WALL}s peak_rss=${RSS_MB}MB"
echo ">> results appended to $RESULTS"

# Correctness oracle: the .phraya embeds a wall-clock timestamp, so hashing it
# directly is useless across runs. Instead hash the sorted variant TSV (no
# timestamp; sorting neutralises any HashMap iteration order) plus the query
# sidecar. These digests MUST match across code revisions for a pure perf change.
TSV_DIGEST="$("$PHRAYA" filter "$OUT" --format tsv 2>/dev/null | sort | sha256sum | cut -d' ' -f1)"
QUERIES_DIGEST="$(sha256sum "${OUT}.queries" 2>/dev/null | cut -d' ' -f1)"
echo ">> correctness: variants_tsv=$TSV_DIGEST queries=$QUERIES_DIGEST"
