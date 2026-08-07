# Benchmark harness expansion: protein alphabets + gap-affine scoring

Design spec for extending `scripts/benchmark/` ahead of implementing ADR-0013
(protein-space alphabet support) and ADR-0014 (gap-affine scoring in `sensitive`).
Written first, per team direction, so both features land with a benchmark story
already in place rather than bolted on after.

## Why the existing harness doesn't cover either feature

Read end-to-end before writing anything:

- **Wrapper contract is fixed at 5 positional args**: `<ref.fasta> <reads_1.fq.gz>
  <reads_2.fq.gz> <out_dir> <threads>` (`benchmark.slurm`, all `wrappers/*.sh`).
  Paired FASTQ is baked in. Protein tools (DIAMOND, BLASTP) take a proteome FASTA
  and a query-protein FASTA — no pairing, no FASTQ.
- **`NUM_ALIGNERS=9` and the `ALIGNERS=(...)` array are hard-coded** in both
  `run_benchmark.sh` and `benchmark.slurm`, kept in sync by hand. A protein-only
  tool set (DIAMOND/BLASTP) and a DNA-only tool set (bwa-mem2/minimap2/...) can't
  share one array without every DNA tool also being asked to run on protein input.
- **`targets.conf`'s 4th column is `GENOME_SIZE_GB`**, feeding BNT
  (bandwidth-normalized throughput) and MEI (memory efficiency index). Both
  formulas normalize by genome gigabytes — meaningless for a proteome, which is
  orders of magnitude smaller and measured in residues, not bases.
- **Placement accuracy (`compute_pa_sam`/`compute_pa_phraya`) decodes wgsim/dwgsim
  read-name position encoding** — a DNA short-read simulator convention. It also
  only tells you *where* a read landed, never whether an indel was scored
  correctly, so gap-affine's actual value proposition (correctly consolidating a
  multi-base indel into one gap-affine-cheap event instead of scattering it as
  ties-with-mismatches) is currently unmeasurable even in the existing DNA suite.
- **`gen_synthetic.py`'s `mutate()` only does substitutions.** There is no indel
  axis in the synthetic data generator at all — so today's harness could not
  benchmark gap-affine's advertised benefit (better modeling of multi-base
  indels) even if the flag existed. This is the sharpest concrete gap:
  fix it regardless of the protein work, because ADR-0014 is unbenchmarkable
  without it.

## Two independent extensions, not one

Protein-alphabet and gap-affine are orthogonal features (ADR-0013/0014 both say
so); the harness work has the same shape and can proceed in either order, but
they touch different files:

| Concern | Protein (ADR-0013) | Gap-affine (ADR-0014) |
|---|---|---|
| New synthetic data need | protein-alphabet generator (new script) | indel events in the *existing* DNA generator |
| New comparison baseline | DIAMOND blastp, BLASTP (real external mapper-class tools) | **none externally** — see below |
| New accuracy metric | protein placement accuracy (position-based, reusable logic) | Indel Event Concordance (op-count based, new) |
| Wrapper contract | 4 args (no paired reads) | same 5-arg contract, new `--gap-model` dimension |
| Targets | new proteome targets (our own generator, staged once) | reuse existing DNA targets unchanged |

### Why no third-party gap-affine baseline in the *timing* harness

WFA2-lib (the algorithm's reference implementation, cited in ADR-0014) is a
**pairwise DP library**, not a mapper — it has no seeding/indexing stage and
expects pre-paired sequences, not a FASTQ-vs-reference workload. Wedging it into
a 5-arg `<ref> <reads_1> <reads_2> <out_dir> <threads>` contract would require
building a seeding/mapping stage around it that doesn't exist, which is out of
scope for a benchmark wrapper. WFA2-lib's proper role is as the **correctness
oracle for the differential test suite** ADR-0014 already calls for (issue #230
acceptance criteria) — a `cargo test` concern, not a `scripts/benchmark/`
concern. The benchmark harness's job here is narrower and better-defined:
compare `phraya --strategy sensitive --gap-model linear` against
`--gap-model affine` on the **same tool**, on indel-enriched data, on timing,
RSS, and indel-recovery accuracy. That is a same-tool, cross-gap-model
comparison — exactly what a perf/behavior harness should measure, leaving
algorithmic correctness to differential tests. Both are proposed as parallel
`phraya-*` wrapper variants (`phraya-sensitive-linear.sh`,
`phraya-sensitive-affine.sh`), the same pattern the harness already uses for
`phraya`/`phraya-sensitive`/`phraya-fast`.

## New metric: Indel Event Concordance (IEC)

Placement accuracy answers "did the read land near the right place?" — it says
nothing about whether the *edit script* matches what was actually simulated.
Gap-affine's whole value proposition is scoring indels better, so the harness
needs a metric that looks at CIGARs, not just positions.

**Definition**: for each simulated read with `N` true indel events (from the
generator's truth sidecar, see below), let `K` be the number of indel
(`I`/`D`) operations in the read's reported CIGAR. The read is **concordant**
if `K == N`. IEC for a run is `concordant_reads / reads_with_indel_events`.

This is deliberately simpler than trying to match indel *positions* between
truth and CIGAR: position-matching requires reconstructing each read's
reference-coordinate alignment start from the `.phraya` TSV, which is a
variant-position sidecar (one row per (position, supporting read), no explicit
per-read alignment-start column) — reconstructing it reliably is more
machinery than the metric is worth. Op-count concordance needs only two columns
that already exist per read: `Provenance` (read ID, joins to truth) and `Cigar`
(same string on every TSV row for that read, since CIGAR is a whole-read
property) from `phraya filter --format tsv`. It is not scored, but it isn't
useless either — with the low indel rates a realistic bacterial-comparison
benchmark uses, "wrong count" is dominated by "gap not consolidated" (linear's
biggest weakness) far more than by "gap consolidated at the wrong position."
It is intentionally *not* folded into a composite score (CAS is an accepted
formula elsewhere and this ADR doesn't touch it) — IEC is a new, independently
reported field, consistent with treating it as a diagnostic for the specific
question gap-affine exists to answer.

## Concrete file plan

### New / changed — synthetic data (implemented in this pass, verifiable locally)

- `scripts/benchmark/local/gen_synthetic.py` — add `--indel-rate` (probability
  per surviving site of an indel event instead of a substitution — the two are
  mutually exclusive per site so total divergence stays interpretable) and
  `--indel-max-len` (uniform `1..=max` length, insertion/deletion 50/50). Write
  a `<reads>.truth.tsv` sidecar: `read_id, start, strand, n_subs, indel_events`
  where `indel_events` is `;`-separated `offset:I|D:len`. This is additive —
  `--indel-rate 0` (default) reproduces today's substitution-only output
  byte-for-byte, so the existing DNA differential/perf harness is unaffected.
- `scripts/benchmark/local/gen_synthetic_protein.py` (new) — same seeded,
  dependency-free design, over the 20-letter amino-acid alphabet. Generates a
  reference "proteome" (concatenated random ORFs) and query proteins at
  controlled substitution+indel divergence, no reverse-complement step (protein
  has no reverse strand — matches ADR-0013's decision to skip RC entirely).
  Emits the same truth-sidecar shape as the DNA generator for the IEC metric to
  reuse without a second implementation.
- `scripts/benchmark/local/compute_indel_recovery.py` (new) — joins a
  `phraya filter --format tsv` dump against a `.truth.tsv` on `Provenance` /
  `read_id`, computes IEC per the definition above.

### New / changed — local quick-bench (implemented in this pass; the phraya
invocation steps reference `--gap-model`/`--alphabet`, which don't exist until
ADR-0014/0013 ship — **these scripts will error at the `phraya align` step
until then**, by design; they exist now so the harness needs zero further
authoring work once the flags land)

- `scripts/benchmark/local/run_local_bench.sh` — add optional
  `INDEL_RATE`/`GAP_MODEL` env overrides, pass `--gap-model` through to
  `phraya align` when set, run `compute_indel_recovery.py` after alignment when
  the generated data has indel events.
- `scripts/benchmark/local/run_local_bench_protein.sh` (new) — mirrors
  `run_local_bench.sh` using `gen_synthetic_protein.py` and
  `phraya plan --alphabet protein`.

### New / changed — SLURM harness (written to existing conventions in this
pass; **not executable here** — no SLURM cluster, no DIAMOND/BLASTP binaries,
no staged proteome data in this sandbox. Structurally reviewed against the
existing wrapper/orchestrator contracts, not run.)

- `scripts/benchmark/slurm/config/targets_protein.conf` (new) — schema
  `TARGET_ID|ORGANISM_PATH|SIZE_CLASS|PROTEOME_RESIDUES`. Unlike the DNA
  targets (real staged genomes, reads simulated externally by dwgsim), there
  is no standard "real proteome + query set with known ground truth" corpus
  to point at — so protein targets are generated once by
  `gen_synthetic_protein.py` at a chosen scale and staged like STREAM Triad's
  one-time platform characterization, not sourced externally. `PROTEOME_RESIDUES`
  is a placeholder (`0`) until populated from the actual generated proteome,
  matching `targets.conf`'s existing placeholder convention for
  `GENOME_SIZE_GB`.
- **`targets.conf` is unchanged.** Gap-affine reuses the existing DNA targets
  as-is: dwgsim's default read simulation already includes realistic indels
  (it's not a substitution-only simulator), so `phraya-sensitive-linear` vs
  `phraya-sensitive-affine` on today's staged data already exercises real
  indel content, no new target axis required. What's *not* in scope for this
  pass: an IEC-equivalent op-count-accuracy metric at SLURM scale, because
  that needs dwgsim's true-indel encoding parsed out of its read names (a
  different, unverified format from the `wgsim` encoding
  `compute_pa_sam`/`compute_pa_phraya` already parse) — flagged as follow-up
  work requiring dwgsim format research, not implemented speculatively here.
  IEC itself is fully implemented and immediately usable at local-bench scale
  (`compute_indel_recovery.py`, verified above) — SLURM-scale IEC is a
  data-format problem, not a metric-design problem.
- `scripts/benchmark/slurm/wrappers/diamond-blastp.sh`,
  `wrappers/blastp.sh` (new) — 4-arg contract:
  `<proteome.fasta> <query.fasta> <out_dir> <threads>`. DIAMOND is the
  DNA-benchmark's `minimap2`/`bwa-mem2` analog (fast, seed-and-extend,
  industry-standard for large-scale protein search); BLASTP is the slow,
  maximally-sensitive accuracy reference (same role BWA/`sensitive` plays on
  the DNA side) — both installed via bioconda, same pattern as the existing
  `bowtie2`/`minibwa` wrapper installation docs.
- `scripts/benchmark/slurm/wrappers/phraya-protein.sh`,
  `phraya-protein-sensitive.sh` (new) — `phraya --alphabet protein` variants,
  4-arg contract, mirroring `phraya.sh`/`phraya-sensitive.sh`.
- `scripts/benchmark/slurm/wrappers/phraya-sensitive-linear.sh`,
  `phraya-sensitive-affine.sh` (new) — `phraya --strategy sensitive
  --gap-model {linear,affine}`, 5-arg (DNA) contract, replacing
  `phraya-sensitive.sh` once `sensitive` defaults to affine (kept as a
  redirect to `-affine` for continuity, since that becomes the true default).
- `run_benchmark.sh` / `benchmark.slurm` — add `--alphabet {dna|protein|both}`
  (default `dna`, so existing invocations are unchanged). Split the hard-coded
  `ALIGNERS` array into `ALIGNERS_DNA` / `ALIGNERS_PROTEIN`, select
  `targets.conf` vs `targets_protein.conf` accordingly, and branch the wrapper
  invocation on arg count (4 vs 5) by alphabet rather than by aligner name.
- `utils/aggregate_results.py` — add `alphabet` and `gap_model` fields to each
  target entry (parsed from the aligner-variant name, e.g.
  `phraya-sensitive-affine` → `gap_model=affine`); branch PA computation for
  protein targets onto a position-only mapeval reusing
  `gen_synthetic_protein.py`'s own truth sidecar (no wgsim-format parsing
  needed, since we control the data); replace `genome_size_gb`/BNT with
  `proteome_residues`/RNT (`queries / (wall_time_s × threads ×
  stream_triad_gbps)`, i.e. BNT's formula with `queries` swapped for `reads`
  — same units, same interpretation, applied to the axis that's actually
  meaningful for protein workloads) for protein target entries.

### Docs

- `scripts/benchmark/slurm/README.md` — new "Protein targets" and "Gap-affine
  targets" sections documenting the two new config files, wrapper contracts,
  and metrics (IEC, RNT), following the existing structure.
- `scripts/benchmark/local/README.md` (new — none existed) — the local
  scripts had no top-level doc; add one now that there are twice as many of
  them, covering the DNA/protein/indel entry points and the correctness-oracle
  convention (`normalize_tsv.py`, now joined by `compute_indel_recovery.py`).

## What ships now vs. what waits

Everything above ships in this pass **except** the two SLURM invocation flags
the wrappers assume (`phraya --alphabet`, `phraya --gap-model`) — those are
ADR-0013/0014's own deliverables and are explicitly out of scope until the ADRs
are accepted and implemented (per the hold the team put on both issues). The
harness is built to the target interface now so no further harness-authoring
work blocks either feature's implementation; the local quick-bench scripts
will visibly fail at the `phraya align` step with a clear "unrecognized
argument" error until then, which is the correct and expected state, not a
bug to silently work around.
