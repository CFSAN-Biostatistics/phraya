# Changelog

All notable changes to Phraya are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Phraya uses [semantic versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **Protein-space alignment** (ADR-0013): `phraya plan --alphabet {auto|dna|protein}` — auto-detects amino-acid input from content (presence of E/F/I/L/P/Q, letters with no IUPAC nucleotide meaning) and switches minimizer seeding to non-canonical (protein has no reverse strand, so no canonicalization and no reverse-complement search) at protein-scale k=6/w=5 defaults. WFA/Myers extension is unchanged — both were already alphabet-blind byte comparators. DNA behavior and performance are unaffected (same code path; protein does strictly less work per query, since it skips the dual-strand search DNA always performs).
- **Gap-affine scoring, defaulted in `sensitive`** (ADR-0014): a gap-affine WFA extension (separate M/I/D wavefronts, `gap_open + L*gap_extend` per L-base gap) that consolidates a real multi-base indel into one CIGAR gap op instead of leaving the search indifferent between one gap and several mismatch-cost-tied substitutions. `--strategy sensitive` uses it by default; `--gap-model {linear|affine}` overrides (sensitive-only — `balanced`/`fast` are Myers-primary with no affine mode and reject the flag). Reported `edit_distance` stays on the traditional mismatches+indel-bases definition regardless of gap model, so `score_alignments`'s 0.95 threshold and all downstream consumers are unaffected. Measured on a synthetic indel-enriched dataset: Indel Event Concorda…
- **QC-layer filter fields**, closing the gap between Phraya's filterable field space and what a downstream SNP pipeline (e.g. CSP2) filters on: `identity` (`1 − edit_distance/query_aligned_len`, the same quantity ADR-0011 stores in `.phraya.queries`) and `match_fraction` (`M/(M+X+I+D)`, BLAST/MUMmer convention) and `aligned_length`, all derived from one CIGAR walk (`phraya_core::cigar::CigarStats`, replacing the private `parse_cigar` in `phraya-align`); `edge_distance` — distance to the nearest end of the read/contig that produced a variant, via a new `VariantObservation::query_position` field stamped at all three CIGAR-op emission sites; `snp_density_15`/`snp_density_125`/`snp_density_1000` — count of other distinct variant positions within a 15/125/1000bp…
- **`phraya plan --use-case {auto|reads-with-ref|contigs-with-reads|contigs-only}`**: escape hatch for use-case auto-detection (CSP2 spike finding B1). Bypasses `detect_use_case` entirely when set.
- **`phraya plan --min-homology <f64>` / `--no-homology-gate`**: gates Case 4's (contigs-only, no reference) all-pairs task generation by minimizer-sketch Jaccard similarity (CSP2 spike finding B3), instead of a dense `i<j` all-pairs list over every contig regardless of relatedness. Default threshold 0.1; `--no-homology-gate` restores the full dense list.
- **`PhrayaPlan::sequence_ids`**: ordered index→sequence-ID table, positionally matching `task_list`'s index space (CSP2 spike finding B2, see Fixed below).

### Changed
- **Breaking**: Filter presets renamed: `conservative` → `strict`, `sensitive` → `tolerant` (ADR-0010). Threshold values unchanged; this is a pure rename to avoid overloading "sensitive" with the alignment strategy layer.

### Fixed
- **`--reference` no longer truncates a multi-record FASTA to its first sequence** (#233). `phraya plan`, `phraya align --reference` (ADR-0011 palette mode), and `phraya plan`'s Case-2/4 task generation now treat every record in a `--reference` file as its own content-hashed reference space, aligning against all of them (N×M for N reads/contigs × M reference records) instead of silently dropping every contig after the first. A multi-contig or multi-chromosome reference (draft assemblies, genomes with plasmids) is the normal case, not the edge case. Batch mode (`--worker`/`--ensure`), whose per-worker output is single-target, now hard-errors on a multi-record reference/centroid file instead of silently aligning against only the first record — use `--reference` mode for multi-contig references.
- **`merge` silently discarded most per-observation QC fields**: rebuilding each merged observation via `VariantObservation::new(...)` reset `in_tandem_repeat` to `false`, `variant_type` to `Snp`, `kmer_uniqueness` to `1.0`, `strand` to its default, and `mate_info` to `None` — every time, on every merge. Consequences on the documented `plan → align → merge → filter` path: `exclude_tandem_repeats` (including in the `strict` preset) was a silent no-op, `--min-kmer-uniqueness` was a silent no-op, indels in a merged file were emitted as SNPs with the wrong REF/ALT in VCF output, and `--require-both-mates-mapped` rejected every variant in any merged file. Fixed by cloning and adjusting only the fields merge legitimately changes (local coverage, pair/insert-size a…
- **CSP2 spike finding B1 (use-case misclassification)**: `detect_use_case`'s with-reference branch required *every* non-reference record to be contig-length (`.all(|seq| seq.len() >= 5000)`) — one short contig (a plasmid, an assembly tail) flipped the whole comparison to Case 2. Replaced with a base-weighted majority rule (≥90% of total bases from contig-length records). The no-reference branch treated "more than one input file" as sufficient evidence for Case 3 (centroid selection) — two draft assemblies passed as plain positional inputs, neither flagged `--reference` (CSP2's actual invocation shape), was always misclassified as Case 3 instead of the direct pairwise Case 4 comparison it should be. Now classifies per input file (contig-like vs read-like by the same base-weighted rule) and only falls to Case 3 when the files are genuinely mixed.
- **CSP2 spike finding B2 (broken `plan-tasks`↔`align` id contract)**: `phraya plan-tasks` printed raw `u32` plan-internal array indices (e.g. `0\t1`); `phraya align`'s traditional mode resolves `QUERY_ID`/`TARGET_ID` by FASTA/FASTQ header string. The documented `phraya plan-tasks | parallel phraya align` pipeline in the README could never work, for any use case. Fixed by adding an ordered `sequence_ids` table to `.phrayaplan` (positionally matching `task_list`'s index space) and having `plan-tasks` resolve through it. Surfaced a second, independent `.phrayaplan` (de)serialization bug while wiring this in: rmp_serde's struct-as-array positional encoding treats a run of `skip_serializing_if` fields as safe only when every member of that run skips or writes *as a unit*; `sequence_ids` (non-empty on essentially every real plan) placed after `sparse_mode`/`dense_kmer_index`/`w11_membership` (routinely skipped) desynced the array — `sequence_ids`'s bytes were read into `insert_size_distribution`'s slot, surfacing as a msgpack marker error rather than a wrong value. Fixed by moving `sequence_ids` into the always-serialized head block alongside `alphabet`, with a regression test covering exactly this "non-empty head field, fully-skipped tail" shape.
- **CSP2 spike finding B3 (no homology gate on Case 4 all-pairs)**: the contigs-only, no-reference task list was a dense `i<j` all-pairs list over every input contig with no pre-filter — two draft assemblies with a few hundred contigs each produced tens of thousands of tasks, almost all between contigs sharing no real homology. Gated by minimizer-sketch Jaccard similarity (already computed at plan time), configurable via `--min-homology`/`--no-homology-gate`.
- **CSP2 spike finding B4 (OOM on ~100kb contigs)**: Case 4 always compares two *comparable-length* sequences, which both the linear (`fill_wfa_fitting_impl`) and affine (`fill_wfa_affine_generic`) WFA engines route through their "global" alignment mode. The linear global path (`fill_wfa`) had no cap parameter at all — it explored up to `s = qn + tn` wavefront phases, each allocating an `O(qn + tn)` array retained for traceback, before any cap was checked (existing `max_s_cap` support only rejected the already-fully-computed result afterward). Production entry points (`wfa_extend`, `wfa_extend_affine`) also never supplied a cap. Combined: a divergent or unrelated ~100kb pair could allocate tens of GB. Fixed by making `fill_wfa` enforce a cap live inside its loop, and giving every production entry point a derived default cap (`default_max_s_cap`). Two further bugs surfaced only once these regression tests were actually executed (not just compiled), both now fixed: (1) `fill_wfa_affine_generic` swept the *full* `-tn..=qn` diagonal range on every wavefront phase regardless of `s`, unlike the linear engine's `[-(s.min(tn)), s.min(qn)]` band — turning O(cap · s) work into O(cap · (qn + tn)), catastrophic once `cap << qn + tn` at contig scale; banded to match the linear engine, which is a sound bound since `gap_open ≥ 0` only ever makes the affine model's true reachable diagonal range a *subset* of the linear model's. (2) every `cap as i32` cast of a caller-supplied `usize` cap silently truncated `usize::MAX` (a common "uncapped" sentinel) to `-1`, making the capped loop's range empty and the call abandon immediately; fixed by saturating through `cap.min(i32::MAX as usize)` before the cast.

## [v0.1.0] - 2026-06-06

### Added

- Phase 1 MVP: Cases 2 (reads + reference), 3 (contigs + reads, auto-centroid), and 4 (contigs only) working end-to-end
- WFA O(s·n) alignment — wavefront-based, not diagonal DP
- SIMD-accelerated diagonal fill via SSE4.2/NEON (`wide` crate)
- K-mer sketching via `simd-minimizers` (AVX2/NEON, k=21, w=11)
- `.phrayaplan` v2 format: MessagePack + zstd, sketch reuse, task list
- `.phraya` position index: VariantObservation with CIGAR, mapq, coverage track, multi-mapping
- `.phraya.queries` query index: multi-mapping alternatives per read
- BAM/CRAM input via `noodles` (pure Rust, no htslib)
- `phraya filter`: threshold-based filtering + named presets (strict/tolerant)
- VCF, TSV, and `.phraya` output formats
- Tandem repeat detection and annotation on variants
- Local coverage computed from alignment (±50bp window)
- Real mapq and avg_base_quality derived from input data
- `phraya-filter` crate: public library API for custom pipelines
- Parallel execution via `rayon`; plan tasks exported for GNU Parallel/SLURM/WDL/Nextflow
- Paired-end filtering with mate info and insert size distribution

### Architecture

- Zero binary dependencies: all alignment in Rust
- Library-first: `phraya-filter` exposes API; CLI is a thin wrapper
- Deferred filtering: alignment produces rich `.phraya`; filter parameters applied post-hoc

---

*Release notes template for future versions:*

```markdown
## [vX.Y.Z] - YYYY-MM-DD

### Added
- ...

### Changed
- ...

### Fixed
- ...

### Removed
- ...
```

[Unreleased]: https://github.com/CFSAN-Biostatistics/phraya/compare/v0.1.0...HEAD
[v0.1.0]: https://github.com/CFSAN-Biostatistics/phraya/releases/tag/v0.1.0
