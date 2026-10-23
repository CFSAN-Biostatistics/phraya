# Changelog

All notable changes to Phraya are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Phraya uses [semantic versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **`--strategy assembly`** (contig-vs-contig, ADR-0015): every chain in either orientation becomes a local alignment block (its exact seed segments joined by global WFA gap fills) instead of one whole-query fitting extension per (query, reference space). Long contigs are no longer abandoned at the WFA cap; a query contig that overlaps a reference contig only partially, or spans several, aligns to each overlap; reverse-complemented contigs align on the reverse strand; large indels split a block, as in dnadiff. Blocks shorter than 100 bp or below 80% identity are dropped, variants come from every kept block, and `AlignmentResult.blocks` carries the blocks (empty for other strategies). On synthetic assembly pairs (`scripts/benchmark/local/gen_assembly_pair.py`, 0.5% SNPs) callable SNP recall goes from 0.44 (identical breakpoints) and 0.00 (shifted breakpoints) under `balanced` to 1.00, with zero false positives.
- **Phase 3 Tier 1 — AVX2/AVX-512 SIMD runtime dispatch**: `count_matching_prefix` now dispatches at runtime from AVX-512F (64 bytes/step) → AVX2 (32 bytes/step) → SSE2 (16 bytes/step) → u64 → scalar. Previously SSE2/NEON was compile-time only. Runtime detection via `is_x86_feature_detected!()` with graceful fallback. Measured 1.5% wall-time improvement on balanced strategy with low-divergence reads on Intel Xeon Gold 5416S (AVX-512F available).
- **Phase 3 Tier 1 — `Strategy::Exact`**: New alignment strategy offering highest precision (K=250 chain cap, ±5bp coverage window, full WFA with no Myers fallback). `--strategy exact` on CLI. ~2× faster than balanced on 150bp reads at 1% divergence (0.614s vs 1.223s).
- **Phase 3 Tier 2 — ShortReadAligner API**: New `phraya_align::ShortReadAligner` module with paired-end configuration (`PairedEndConfig`), scoring parameters (`ShortReadScoring`), mate position prediction (`MateRescueResult`), and adaptive k-mer size selection. 7 unit tests covering single-end/paired-end modes, insert size validation, orientation checks, and k-mer tuning.
- **Core SNP FASTA output** (`--format core-snp`): Multi-sequence FASTA at positions present in ALL samples (indels excluded, missing positions excluded). Includes `--sample-ids` override for multi-sample sample name customization.
- **Distance matrix TSV output** (`--format distance-matrix`): Square TSV Hamming distance matrix on core SNP bases (excludes N/missing).
- **Enhanced QC summary**: `phraya qc` appends 8 columns (strategy, elapsed_ms, min_coverage, min_mapq, min_identity, kmer_uniqueness_min/mean/max). k-mer uniqueness computed from observations; strategy/timing/filter params emit NA (not stored in .phraya format).
- **Gap-affine scoring, defaulted in `sensitive`** (ADR-0014): a gap-affine WFA extension (separate M/I/D wavefronts, `gap_open + L*gap_extend` per L-base gap) that consolidates a real multi-base indel into one CIGAR gap op instead of leaving the search indifferent between one gap and several mismatch-cost-tied substitutions. `--strategy sensitive` uses it by default; `--gap-model {linear|affine}` overrides (sensitive-only — `balanced`/`fast` are Myers-primary with no affine mode and reject the flag). Reported `edit_distance` stays on the traditional mismatches+indel-bases definition regardless of gap model, so `score_alignments`'s 0.95 threshold and all downstream consumers are unaffected. Measured on a synthetic indel-enriched dataset: Indel Event Concorda…
- **QC-layer filter fields**, closing the gap between Phraya's filterable field space and what a downstream SNP pipeline (e.g. CSP2) filters on: `identity` (`1 − edit_distance/query_aligned_len`, the same quantity ADR-0011 stores in `.phraya.queries`) and `match_fraction` (`M/(M+X+I+D)`, BLAST/MUMmer convention) and `aligned_length`, all derived from one CIGAR walk (`phraya_core::cigar::CigarStats`, replacing the private `parse_cigar` in `phraya-align`); `edge_distance` — distance to the nearest end of the read/contig that produced a variant, via a new `VariantObservation::query_position` field stamped at all three CIGAR-op emission sites; `snp_density_15`/`snp_density_125`/`snp_density_1000` — count of other distinct variant positions within a 15/125/1000bp…
- **`phraya plan --min-homology <f64>` / `--no-homology-gate`**: gates Case 4's (contigs-only, no reference) all-pairs task generation by minimizer-sketch Jaccard similarity (CSP2 spike finding B3), instead of a dense `i<j` all-pairs list over every contig regardless of relatedness. Default threshold 0.1; `--no-homology-gate` restores the full dense list.
- **`PhrayaPlan::sequence_ids`**: ordered index→sequence-ID table, positionally matching `task_list`'s index space (CSP2 spike finding B2, see Fixed below).
- **`phraya filter --format snpdiffs` accepts multiple `.phraya` files**: one per reference contig from `phraya align --reference ... --strategy assembly`, merged into a single snpdiffs block (`phraya filter out/*.phraya --format snpdiffs`). Previously only `inputs[0]` was used; all additional files were silently dropped.
### Removed
- **`UseCase` enum removed from `PhrayaPlan`** (breaking for any consumer reading
  `plan.use_case`): input-shape detection (`Case 2`/`3`/`4`) is teaching language per
  AGENTS.md, not a stored-format type. `--use-case {auto|reads-with-ref|contigs-with-reads|contigs-only}`
  CLI flag removed; `detect_use_case` and `use_case_number` in `phraya-cli/src/main.rs`
  consolidated into a single `generate_task_list` that derives tasks directly from input
  characteristics (read vs contig classification, file count, reference presence). The
  `--use-case` escape hatch for CSP2 spike finding B1 is no longer needed — the B1
  heuristic fix (majority-by-bases, per-file classification) makes auto-detection
  reliable; override was a workaround, not a feature.
- **`phraya-io` `use_case` module deleted** (`classify_input`, `InputType`, `UseCaseError`):
  a second, divergent file-level Read/Contig classifier (first-record <5kb vs base-weighted
  majority) that had no callers outside its own tests — the CLI uses `is_contig_like` in
  `main.rs`. Removing it eliminates the risk of two classifiers disagreeing.

### Changed
- **`.phrayaplan` format v8** (breaking): `SharedFrame` no longer serializes `use_case`,
  dropping the second positional msgpack array slot. Since rmp_serde encodes structs as
  positional arrays (no field-name keys), old v7 plans are rejected with
  `VersionMismatch{expected:8, got:7}` rather than silently mis-parsing (`PHRAYAPLAN_VERSION`
  bumped 7→8). Existing plan files must be regenerated with `phraya plan`.
- **Reference-palette mode (`align --reference`) was ~10.6s/record regardless of target size, and single-threaded**: `run_align_reference` aligned reads with a serial loop while every other alignment path uses rayon, and `align_read` recomputed each read's content hash, forward sketch, reverse complement, and reverse-complement sketch on every call — all of it space-independent, but repeated once per (space, read) rather than once per read. On a 214-record reference this made a 227x-slower-than-minibwa benchmark result. Split the space-independent work into `prepare_query`/`align_read_prepared` (computed once per read) and parallelized the per-space read loop with the same batched-rayon pattern the batch worker already used.
- **`phraya plan` never stored a reference space's FASTA record ID as its name**, so every `align --reference` output and cross-space sidecar entry fell back to a 16-character content-hash-prefix label (`reference_space_label`'s hash-prefix branch). Reference spaces now take their name from the record ID at plan time; collision-safe (`<label>__<hash16>`) only when two records in the plan's palette sanitize to the same label, so composability (`align({A,B})` writes byte-identical `A.phraya` to `align({A})`) is unaffected.
- **`phraya_accuracy.py` (benchmark scoring) crashed on the cross-space sidecar's 3-element `(space, pos, identity)` placements** (`ValueError: invalid literal for int()`), assuming the legacy 2-element `(pos, score)` schema — leaving placement accuracy `null` for every multi-record reference in the benchmark. Now detects the schema per file and, for cross-space placements, requires both position and (sanitized) space label to match the read's true chromosome; falls back to position-only scoring with a warning for a sidecar written before reference spaces carried names.
- **`scripts/benchmark/slurm/utils/aggregate.slurm` self-clobbered `SCRIPT_DIR`** by always recomputing it from `${BASH_SOURCE[0]}`, which resolves to the SLURM spool copy under `sbatch` — the chained aggregation job could never find `aggregate_results.py`. Now respects an orchestrator-exported `SCRIPT_DIR`, matching `benchmark.slurm`.
- **Every read was force-aligned against every reference it did not come from, fabricating variant calls**: `anchors_from_chains` fell back to a `(0,0)` anchor whenever chaining found nothing, unconditionally. For a read that shares no seed with a target far longer than itself there is no candidate locus, so this aligned the read against `target[0..2*query_len)`, scored ~50-60% identity, and emitted ~`query_len/2` `VariantObservation`s — in a reference the read never came from. `write_cross_space_queries` filters placements at identity ≥ 0.95, so the cross-space sidecar hid this, but the per-space `.phraya` has no such filter and kept every fabricated call. Measured on a 20-record palette where all reads came from one record: 1.87 MB of fabricated calls in each of the 19 non-homologous spaces against 42 KB of real signal in the one true space, at 256µs per (read, space) pair — 2.5x the cost of a genuine placement. Anchor selection is now evidence-ordered (`select_anchors`): chains when they exist; otherwise the seeds' own diagonals, which is what makes a repetitive query that chains to nothing recoverable at its real locus rather than by luck at offset 0; otherwise `(0,0)`, but only where query and target are comparable in length (the contig-vs-contig case issue #146 added it for). Real output is unchanged — the space reads do come from, and the cross-space sidecar, are both byte-identical to before.
- **Reference-palette mode paid full per-space setup for spaces no read could ever place in**: a palette space whose sketch shares no minimizer value with any read cannot produce a seed, hence no chain, hence no placement — but `align --reference` still built its minimizer index and scanned it for tandem repeats. On a 214-record vertebrate palette only ~23 records ever receive a placement. Such spaces are now detected up front (exact, not heuristic: a shared minimizer is a *necessary* condition for a seed) and written directly as an empty track via the new O(1) `CoverageTrack::uncovered`. Combined with the anchor-selection fix above, a 20-record / 2000-read palette alignment went from ~24.8s to ~1.4s (~18x) with byte-identical real output.
- **`detect_tandem_repeats` copied the entire sequence twice per call** — `String::from_utf8_lossy` at the call site, then `str::to_uppercase` inside — ~400 MB of transient allocation per call on a 197 Mb chromosome, for a scan that only reads bytes. Now takes `&[u8]` and folds case per comparison. `TargetContext` also defers detection until a variant actually needs it (`OnceLock`), so a target that places no read never runs the O(target) scan at all.
- **`CoverageTrack::new` forced an 8-byte-per-base copy of the whole reference**: it took `Vec<usize>`, so all three `phraya align` call sites materialized a `Vec<usize>` (1.6 GB on a 197 Mb chromosome) purely to be quantized and dropped, plus an intermediate `Vec<u8>` built inside. Now takes `&[u32]` — the coverage type the rest of the pipeline already uses — and run-length encodes in a single pass. The `merge` path's accumulator narrowed from `usize` to `u32` for the same reason.
- **`phraya align QUERY TARGET` hard-errored when a pair did not align.** Case 3/4 plans enumerate every query×target pair precisely because homology is unknown up front, so most pairs in a multi-contig plan are expected to produce nothing; this only ever "succeeded" because the `(0,0)` fallback fabricated an alignment for every pair. A pair that does not align now writes an empty `.phraya` and empty sidecar entry, matching what palette mode writes for an unreachable space.
- **Extension scored a mismatch for every soft-masked base, collapsing identity on repeat-annotated eukaryotic genomes**: NCBI RefSeq marks repeat-annotated regions with lowercase (gallus_gallus, plasmodium_falciparum, and other repeat-rich assemblies ship this way by default — a normal convention, not an edge case), but every WFA/Myers extension engine compared raw bytes with no case normalization, while minimizer sketching already treated case as non-semantic. A read seeded and chained to its exact correct position and strand scored one mismatch per lowercase target base against its own near-perfect match, collapsing identity from ~99% to ~50% and dropping it below the 0.95 report threshold — identically across `fast`/`balanced`/`sensitive`, since all three route through one of the affected engines. This is the actual root cause behind issue #194's T2/T5 placement collapse; #195's occurrence-cap masking and ADR-0009's dense seeding were real fixes for real (different) problems, but neither touched this. Fixed by normalizing case once, at the `Sequence::new` chokepoint every sequence is built through. Measured on real *P. falciparum* chr1 (61.8% soft-masked): single-chromosome placement 4.3% → 99.96%; full 15-chromosome palette 21.7× more reads placed.
- **(0,0)` fallback anchor emitted fabricated `VariantObservation`s for comparable-length non-homologous pairs**: the existing scope fix (line 38) only suppressed the fallback when the target exceeded 10× the query (read-vs-long-contig). In reference-palette mode with comparable-length contigs (e.g. two bacterial assemblies of similar contig sizes), `fallback_anchor_applies` returns true and the `(0,0)` anchor fires for every non-homologous (query, reference-space) pair, producing a full-length ~50-60% identity alignment whose ~query_len/2 variant calls are pure noise — deposited into the per-space `.phraya` (the cross-space sidecar already discards them at its 0.95 identity filter). Added `FALLBACK_ONLY_MIN_IDENTITY` (0.55): when an alignment relied solely on the `(0,0)` fallback (no shared minimizer seed) and scored below this floor, variant extraction is suppressed while coverage depth and the placement record are still emitted (the latter is filtered by the sidecar anyway). Seeded alignments are unaffected. Resolves the 5.2× SNP inflation in the CSP2 snpdiffs escalation (173 840 → ~33 000 noise-free).

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
