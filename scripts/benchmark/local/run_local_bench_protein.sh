#!/usr/bin/env bash
# Local before/after benchmark for `phraya align` in protein mode (ADR-0013).
#
# Mirrors run_local_bench.sh's design (seeded synthetic data reused across runs,
# timed single-threaded alignment, correctness-oracle TSV digest) but for a
# proteome + query-protein workload instead of a reference + FASTQ-reads workload.
#
# Uses reference-palette mode (`align --reference <proteome> --output <dir>`), not
# batch mode: a generated proteome is always multi-record (a database of many
# proteins is the realistic shape of a protein-search reference), and batch mode
# hard-errors on any multi-record reference/centroid file. This is the same
# restriction the SLURM DNA harness hit on multi-chromosome targets and fixed in
# 58c7c2d by switching to --reference; this script had the identical bug from the
# day it was written (BENCHMARK_EXPANSION.md), just never exercised past the
# ADR-0013 plan-time panic until that was fixed.
#
# Usage:
#   scripts/benchmark/local/run_local_bench_protein.sh <label> [num_proteins] [protein_len] [num_queries] [divergence]
#
# Env:
#   BENCH_DIR      where data/outputs live (default: ${TMPDIR:-/tmp}/phraya-bench)
#   PHRAYA         phraya binary (default: target/release/phraya)
#   QUERY_LEN      query length in aa (default: same as protein_len -> whole-protein queries)
#   INDEL_RATE     per-site indel probability (default: 0, no indels)
#   INDEL_MAX_LEN  max indel length in residues (default: 3)
#   STRATEGY       phraya --strategy override (default: unset -> phraya's own default)
set -euo pipefail

LABEL="${1:?usage: run_local_bench_protein.sh <label> [num_proteins] [protein_len] [num_queries] [divergence]}"
NUM_PROTEINS="${2:-2000}"
PROTEIN_LEN="${3:-300}"
NUM_QUERIES="${4:-500}"
DIVERGENCE="${5:-0.02}"
SEED="${SEED:-1}"
QUERY_LEN="${QUERY_LEN:-$PROTEIN_LEN}"
INDEL_RATE="${INDEL_RATE:-0}"
INDEL_MAX_LEN="${INDEL_MAX_LEN:-3}"
STRATEGY="${STRATEGY:-}"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
BENCH_DIR="${BENCH_DIR:-${TMPDIR:-/tmp}/phraya-bench}"
PHRAYA="${PHRAYA:-$REPO_ROOT/target/release/phraya}"
GEN="$REPO_ROOT/scripts/benchmark/local/gen_synthetic_protein.py"

KEY="protein_p${NUM_PROTEINS}_l${PROTEIN_LEN}_q${NUM_QUERIES}_ql${QUERY_LEN}_d${DIVERGENCE}_s${SEED}"
if [[ "$INDEL_RATE" != "0" ]]; then
    KEY="${KEY}_indel${INDEL_RATE}_max${INDEL_MAX_LEN}"
fi
DATA_DIR="$BENCH_DIR/$KEY"
PROTEOME="$DATA_DIR/proteome.fa"
QUERIES="$DATA_DIR/queries.fa"
TRUTH="$DATA_DIR/queries.fa.truth.tsv"
PLAN="$DATA_DIR/plan.phrayaplan"
OUT_DIR="$DATA_DIR/${LABEL}_out"
RESULTS="$BENCH_DIR/results_protein.tsv"

mkdir -p "$DATA_DIR"

if [[ ! -f "$PHRAYA" ]]; then
    echo "ERROR: phraya binary not found at $PHRAYA" >&2
    echo "  build with: RUSTFLAGS=\"-C target-cpu=native\" cargo build --release" >&2
    exit 1
fi

if [[ ! -f "$PROTEOME" || ! -f "$QUERIES" ]]; then
    echo ">> generating synthetic protein data ($KEY)"
    python3 "$GEN" --num-proteins "$NUM_PROTEINS" --protein-len "$PROTEIN_LEN" \
        --num-queries "$NUM_QUERIES" --query-len "$QUERY_LEN" --divergence "$DIVERGENCE" \
        --indel-rate "$INDEL_RATE" --indel-max-len "$INDEL_MAX_LEN" --seed "$SEED" \
        --proteome-out "$PROTEOME" --queries-out "$QUERIES"
fi

echo ">> planning (--alphabet protein)"
"$PHRAYA" plan --inputs "$QUERIES" --reference "$PROTEOME" --output "$PLAN" \
    --alphabet protein >/dev/null

# Reference-palette mode, not batch mode: PROTEOME is always multi-record (a
# proteome is a database of many proteins), and batch mode hard-errors on any
# multi-record reference. --output is a directory: one <protein_id>.phraya per
# proteome record plus a cross_space.phraya.queries union sidecar.
rm -rf "$OUT_DIR"
ALIGN_ARGS=(align --reference "$PROTEOME" --output "$OUT_DIR" "$PLAN")
[[ -n "$STRATEGY" ]] && ALIGN_ARGS+=(--strategy "$STRATEGY")

echo ">> aligning (label=$LABEL, strategy=${STRATEGY:-default}, single-threaded)"
TIMER="$REPO_ROOT/scripts/benchmark/local/time_run.py"
MEASURE="$(RAYON_NUM_THREADS=1 python3 "$TIMER" "$PHRAYA" "${ALIGN_ARGS[@]}")"
WALL="$(printf '%s' "$MEASURE" | cut -f1)"
RSS_MB="$(printf '%s' "$MEASURE" | cut -f2)"

if [[ ! -f "$RESULTS" ]]; then
    printf 'label\tkey\twall_s\tpeak_rss_mb\tstrategy\toutput\n' > "$RESULTS"
fi
printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$LABEL" "$KEY" "$WALL" "$RSS_MB" "${STRATEGY:-default}" "$OUT_DIR" >> "$RESULTS"

echo ">> $LABEL: wall=${WALL}s peak_rss=${RSS_MB}MB"
echo ">> results appended to $RESULTS"

# Correctness oracle: palette mode writes one .phraya per proteome record, not one
# named file, so concatenate every record's filtered TSV (skipping the
# cross_space.phraya.queries sidecar, a different format `filter` doesn't read)
# before normalizing/hashing — same convention as run_local_bench.sh, generalized
# to N reference spaces.
NORM="$REPO_ROOT/scripts/benchmark/local/normalize_tsv.py"
FILTER_TSV=""
for space_file in "$OUT_DIR"/*.phraya; do
    FILTER_TSV+="$("$PHRAYA" filter "$space_file" --format tsv 2>/dev/null)"$'\n'
done
TSV_DIGEST="$(printf '%s' "$FILTER_TSV" | python3 "$NORM" | sort | sha256sum | cut -d' ' -f1)"
echo ">> correctness: variants_tsv=$TSV_DIGEST"

if [[ "$INDEL_RATE" != "0" && -f "$TRUTH" ]]; then
    IEC_SCRIPT="$REPO_ROOT/scripts/benchmark/local/compute_indel_recovery.py"
    echo ">> indel event concordance:"
    printf '%s' "$FILTER_TSV" | python3 "$IEC_SCRIPT" --truth "$TRUTH" | sed 's/^/   /'
fi
