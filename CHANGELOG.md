# Changelog

All notable changes to Phraya are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Phraya uses [semantic versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **Protein-space alignment** (ADR-0013): `phraya plan --alphabet {auto|dna|protein}` — auto-detects amino-acid input from content (presence of E/F/I/L/P/Q, letters with no IUPAC nucleotide meaning) and switches minimizer seeding to non-canonical (protein has no reverse strand, so no canonicalization and no reverse-complement search) at protein-scale k=6/w=5 defaults. WFA/Myers extension is unchanged — both were already alphabet-blind byte comparators. DNA behavior and performance are unaffected (same code path; protein does strictly less work per query, since it skips the dual-strand search DNA always performs).
- **Gap-affine scoring, defaulted in `sensitive`** (ADR-0014): a gap-affine WFA extension (separate M/I/D wavefronts, `gap_open + L*gap_extend` per L-base gap) that consolidates a real multi-base indel into one CIGAR gap op instead of leaving the search indifferent between one gap and several mismatch-cost-tied substitutions. `--strategy sensitive` uses it by default; `--gap-model {linear|affine}` overrides (sensitive-only — `balanced`/`fast` are Myers-primary with no affine mode and reject the flag). Reported `edit_distance` stays on the traditional mismatches+indel-bases definition regardless of gap model, so `score_alignments`'s 0.95 threshold and all downstream consumers are unaffected. Measured on a synthetic indel-enriched dataset: Indel Event Concorda…
- **QC-layer filter fields**, closing the gap between Phraya's filterable field space and what a downstream SNP pipeline (e.g. CSP2) filters on: `identity` (`1 − edit_distance/query_aligned_len`, the same quantity ADR-0011 stores in `.phraya.queries`) and `match_fraction` (`M/(M+X+I+D)`, BLAST/MUMmer convention) and `aligned_length`, all derived from one CIGAR walk (`phraya_core::cigar::CigarStats`, replacing the private `parse_cigar` in `phraya-align`); `edge_distance` — distance to the nearest end of the read/contig that produced a variant, via a new `VariantObservation::query_position` field stamped at all three CIGAR-op emission sites; `snp_density_15`/`snp_density_125`/`snp_density_1000` — count of other distinct variant positions within a 15/125/1000bp window, recomputed on every `filter` run (so it's correct on merged files for free). All available as both discrete threshold flags (`--min-identity`, `--min-match-fraction`, `--min-aligned-length`, `--min-edge-distance`, `--max-snp-density-{15,125,1000}`) and via `--expr`, which is now wired into the CLI (previously dead code — `ExprFilter` existed with ~20 passing tests but no `Filter` subcommand field reached it). `phraya qc <file>...`: new subcommand, one TSV row per `.phraya` file reporting variant count and coverage breadth at 1x/10x, computed from a newly-exposed `AlignmentResult::raw_coverage` (the existing `coverage` field is quantized to the nearest 5, which maps depth 1–2 to 0 and makes single-read coverage unmeasurable).

### Changed
- **Breaking**: Filter presets renamed: `conservative` → `strict`, `sensitive` → `tolerant` (ADR-0010). Threshold values unchanged; this is a pure rename to avoid overloading "sensitive" with the alignment strategy layer.

### Fixed
- **`--reference` no longer truncates a multi-record FASTA to its first sequence** (#233). `phraya plan`, `phraya align --reference` (ADR-0011 palette mode), and `phraya plan`'s Case-2/4 task generation now treat every record in a `--reference` file as its own content-hashed reference space, aligning against all of them (N×M for N reads/contigs × M reference records) instead of silently dropping every contig after the first. A multi-contig or multi-chromosome reference (draft assemblies, genomes with plasmids) is the normal case, not the edge case. Batch mode (`--worker`/`--ensure`), whose per-worker output is single-target, now hard-errors on a multi-record reference/centroid file instead of silently aligning against only the first record — use `--reference` mode for multi-contig references.
- **`merge` silently discarded most per-observation QC fields**: rebuilding each merged observation via `VariantObservation::new(...)` reset `in_tandem_repeat` to `false`, `variant_type` to `Snp`, `kmer_uniqueness` to `1.0`, `strand` to its default, and `mate_info` to `None` — every time, on every merge. Consequences on the documented `plan → align → merge → filter` path: `exclude_tandem_repeats` (including in the `strict` preset) was a silent no-op, `--min-kmer-uniqueness` was a silent no-op, indels in a merged file were emitted as SNPs with the wrong REF/ALT in VCF output, and `--require-both-mates-mapped` rejected every variant in any merged file. Fixed by cloning and adjusting only the fields merge legitimately changes (local coverage, pair/insert-size aggregates) instead of rebuilding from scratch.

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
