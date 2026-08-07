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
#   BENCH_DIR      where data/outputs live (default: ${TMPDIR:-/tmp}/phraya-bench)
#   PHRAYA         phraya binary (default: target/release/phraya)
#   INDEL_RATE     per-site indel probability for gen_synthetic.py (default: 0, no indels)
#   INDEL_MAX_LEN  max indel length in bp (default: 5)
#   STRATEGY       phraya --strategy override (default: unset -> phraya's own default)
#   GAP_MODEL      phraya --gap-model override; requires STRATEGY=sensitive (ADR-0014).
#                  NOT YET SUPPORTED by `phraya align` — set only once ADR-0014 ships;
#                  until then this errors at the align step with "unrecognized argument",
#                  which is expected, not a bug in this script.
#
# When INDEL_RATE > 0, also runs compute_indel_recovery.py against the generated
# truth sidecar and reports Indel Event Concordance (IEC) alongside timing.
set -euo pipefail

LABEL="${1:?usage: run_local_bench.sh <label> [genome_size] [num_reads] [read_len] [divergence]}"
GENOME_SIZE="${2:-2000000}"
NUM_READS="${3:-20000}"
READ_LEN="${4:-150}"
DIVERGENCE="${5:-0.01}"
SEED="${SEED:-1}"
INDEL_RATE="${INDEL_RATE:-0}"
INDEL_MAX_LEN="${INDEL_MAX_LEN:-5}"
STRATEGY="${STRATEGY:-}"
GAP_MODEL="${GAP_MODEL:-}"

if [[ -n "$GAP_MODEL" && "$STRATEGY" != "sensitive" ]]; then
    echo "ERROR: GAP_MODEL requires STRATEGY=sensitive (ADR-0014: gap-affine is sensitive-only)" >&2
    exit 1
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
BENCH_DIR="${BENCH_DIR:-${TMPDIR:-/tmp}/phraya-bench}"
PHRAYA="${PHRAYA:-$REPO_ROOT/target/release/phraya}"
GEN="$REPO_ROOT/scripts/benchmark/local/gen_synthetic.py"

# Data is keyed by scale params (including the indel axis) so identical params
# reuse identical inputs; INDEL_RATE=0 keeps the original key shape unchanged.
KEY="g${GENOME_SIZE}_n${NUM_READS}_l${READ_LEN}_d${DIVERGENCE}_s${SEED}"
if [[ "$INDEL_RATE" != "0" ]]; then
    KEY="${KEY}_indel${INDEL_RATE}_max${INDEL_MAX_LEN}"
fi
DATA_DIR="$BENCH_DIR/$KEY"
REF="$DATA_DIR/ref.fa"
READS="$DATA_DIR/reads.fq"
TRUTH="$DATA_DIR/reads.fq.truth.tsv"
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
        --indel-rate "$INDEL_RATE" --indel-max-len "$INDEL_MAX_LEN" \
        --ref-out "$REF" --reads-out "$READS"
fi

echo ">> planning"
"$PHRAYA" plan --inputs "$READS" --reference "$REF" --output "$PLAN" \
    --batch-to 1 --batch-output-pattern "$OUT" >/dev/null

ALIGN_ARGS=(align "$PLAN" --worker 0)
[[ -n "$STRATEGY" ]] && ALIGN_ARGS+=(--strategy "$STRATEGY")
[[ -n "$GAP_MODEL" ]] && ALIGN_ARGS+=(--gap-model "$GAP_MODEL")

echo ">> aligning (label=$LABEL, strategy=${STRATEGY:-default}, gap_model=${GAP_MODEL:-default}, single-threaded)"
TIMER="$REPO_ROOT/scripts/benchmark/local/time_run.py"
MEASURE="$(RAYON_NUM_THREADS=1 python3 "$TIMER" "$PHRAYA" "${ALIGN_ARGS[@]}")"
WALL="$(printf '%s' "$MEASURE" | cut -f1)"
RSS_MB="$(printf '%s' "$MEASURE" | cut -f2)"

if [[ ! -f "$RESULTS" ]]; then
    printf 'label\tkey\twall_s\tpeak_rss_mb\tstrategy\tgap_model\toutput\n' > "$RESULTS"
fi
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$LABEL" "$KEY" "$WALL" "$RSS_MB" "${STRATEGY:-default}" "${GAP_MODEL:-default}" "$OUT" >> "$RESULTS"

echo ">> $LABEL: wall=${WALL}s peak_rss=${RSS_MB}MB"
echo ">> results appended to $RESULTS"

# Correctness oracle: the .phraya embeds a wall-clock timestamp, so hashing it
# directly is useless across runs. Instead hash the variant TSV with allele-token
# order normalized (normalize_tsv.py) and lines sorted — this is invariant to the
# per-process HashMap iteration order, so it changes ONLY if the actual variant
# data changes. This digest MUST match across code revisions for a pure perf change.
NORM="$REPO_ROOT/scripts/benchmark/local/normalize_tsv.py"
FILTER_TSV="$("$PHRAYA" filter "$OUT" --format tsv 2>/dev/null)"
TSV_DIGEST="$(printf '%s\n' "$FILTER_TSV" | python3 "$NORM" | sort | sha256sum | cut -d' ' -f1)"
echo ">> correctness: variants_tsv=$TSV_DIGEST"

# Indel Event Concordance (see BENCHMARK_EXPANSION.md): only meaningful when the
# generated data actually has indel events, and only informative once a real
# gap-affine/linear comparison is possible (ADR-0014) — computed regardless so
# `--strategy sensitive` linear-mode runs get a baseline IEC too.
if [[ "$INDEL_RATE" != "0" && -f "$TRUTH" ]]; then
    IEC_SCRIPT="$REPO_ROOT/scripts/benchmark/local/compute_indel_recovery.py"
    echo ">> indel event concordance:"
    printf '%s\n' "$FILTER_TSV" | python3 "$IEC_SCRIPT" --truth "$TRUTH" | sed 's/^/   /'
fi
