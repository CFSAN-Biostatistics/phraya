# Local Benchmark Scripts

Fast, dependency-light before/after benchmarking for the `phraya align` hot
path — no SLURM, no staged HPC data, no external aligners. Everything here
runs on a laptop against seeded synthetic data. See `../BENCHMARK_EXPANSION.md`
for the design rationale behind the indel/protein/gap-model additions.

## DNA

```bash
scripts/benchmark/local/run_local_bench.sh <label> [genome_size] [num_reads] [read_len] [divergence]
```

Generates (or reuses) a seeded synthetic reference + FASTQ reads via
`gen_synthetic.py`, plans + aligns, times it, and hashes the resulting variant
TSV as a correctness oracle (must match across revisions for a pure perf
change — see `normalize_tsv.py`).

Env overrides:
- `INDEL_RATE` / `INDEL_MAX_LEN` — generate indel-bearing reads (default: no
  indels, byte-identical to pre-indel-support output). When set, also runs
  `compute_indel_recovery.py` and reports Indel Event Concordance (IEC) — see
  "Indel Event Concordance" below.
- `STRATEGY` / `GAP_MODEL` — pass `--strategy`/`--gap-model` through to
  `phraya align` (ADR-0014; `GAP_MODEL` requires `STRATEGY=sensitive`,
  enforced by the script). **`GAP_MODEL` isn't runnable yet** — `phraya align`
  doesn't accept `--gap-model` until ADR-0014 ships; setting it fails clearly
  at the align step until then, which is expected, not a bug.

## Protein (ADR-0013)

```bash
scripts/benchmark/local/run_local_bench_protein.sh <label> [num_proteins] [protein_len] [num_queries] [divergence]
```

Same shape as the DNA script, over `gen_synthetic_protein.py`'s reference
proteome + query proteins instead of a reference genome + reads. No
reverse-complement step (protein has no reverse strand). Supports the same
`INDEL_RATE`/`INDEL_MAX_LEN`/`STRATEGY` env overrides.

**Not runnable yet** — `phraya plan --alphabet protein` doesn't exist until
ADR-0013 ships; the script fails clearly at the `plan` step until then. Data
generation (`gen_synthetic_protein.py`) works today and is independently
useful for validating the generator.

## Scripts

| Script | Role |
|---|---|
| `gen_synthetic.py` | Seeded DNA reference + FASTQ read generator. Substitutions + indels; writes a `<reads>.truth.tsv` sidecar (read_id, start, strand, n_subs, indel_events). |
| `gen_synthetic_protein.py` | Seeded protein proteome + query generator. Imports `mutate_with_truth`/`format_indel_events` from `gen_synthetic.py` (one mutation implementation, two alphabets) — no reimplementation, no drift. |
| `compute_indel_recovery.py` | Computes Indel Event Concordance: joins a `phraya filter --format tsv` dump against a `.truth.tsv` on read/query ID, checks whether each read's true indel-event count matches its reported CIGAR's I/D op count. |
| `run_local_bench.sh` / `run_local_bench_protein.sh` | End-to-end plan+align+time+hash harness, DNA and protein respectively. |
| `time_run.py` | Portable `wall_seconds\tpeak_rss_mb` measurement (no `/usr/bin/time` dependency). |
| `normalize_tsv.py` | Sorts the `all_alleles` HashMap-rendered column so the variant TSV digest is stable across runs regardless of HashMap iteration order — the correctness-oracle convention both bench scripts rely on. |

## Indel Event Concordance (IEC)

Placement accuracy answers "did the read land near the right place?" — nothing
about whether the *edit script* is right. Gap-affine's (ADR-0014) whole value
proposition is scoring indels better, so IEC checks CIGARs instead: for each
simulated read with `N` true indel events, is the reported CIGAR's indel
(`I`/`D`) op count also `N`? IEC = concordant reads / reads with `N > 0`.

This already works against today's linear-cost default and gives a real
baseline: on a `div=0.01, indel_rate=0.01, indel_max_len=4` run (300 reads,
20 kb reference), `phraya align`'s default strategy scored **IEC = 0.147** —
most true indels get fragmented into several CIGAR ops instead of one (one
3-event read was reported as 20 separate ops). That's the exact failure mode
ADR-0014's gap-affine WFA exists to fix; once it ships, rerun the same command
with `STRATEGY=sensitive GAP_MODEL=affine` for a direct before/after.

```bash
INDEL_RATE=0.01 INDEL_MAX_LEN=4 scripts/benchmark/local/run_local_bench.sh baseline 20000 300 100 0.01
```
