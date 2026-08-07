# Changelog

All notable changes to Phraya are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Phraya uses [semantic versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **Protein-space alignment** (ADR-0013): `phraya plan --alphabet {auto|dna|protein}` — auto-detects amino-acid input from content (presence of E/F/I/L/P/Q, letters with no IUPAC nucleotide meaning) and switches minimizer seeding to non-canonical (protein has no reverse strand, so no canonicalization and no reverse-complement search) at protein-scale k=6/w=5 defaults. WFA/Myers extension is unchanged — both were already alphabet-blind byte comparators. DNA behavior and performance are unaffected (same code path; protein does strictly less work per query, since it skips the dual-strand search DNA always performs).
- **Gap-affine scoring, defaulted in `sensitive`** (ADR-0014): a gap-affine WFA extension (separate M/I/D wavefronts, `gap_open + L*gap_extend` per L-base gap) that consolidates a real multi-base indel into one CIGAR gap op instead of leaving the search indifferent between one gap and several mismatch-cost-tied substitutions. `--strategy sensitive` uses it by default; `--gap-model {linear|affine}` overrides (sensitive-only — `balanced`/`fast` are Myers-primary with no affine mode and reject the flag). Reported `edit_distance` stays on the traditional mismatches+indel-bases definition regardless of gap model, so `score_alignments`'s 0.95 threshold and all downstream consumers are unaffected. Measured on a synthetic indel-enriched dataset: Indel Event Concordance improved from 0.34 (linear) to 0.54 (affine) under `sensitive`.

### Changed
- **Breaking**: Filter presets renamed: `conservative` → `strict`, `sensitive` → `tolerant` (ADR-0010). Threshold values unchanged; this is a pure rename to avoid overloading "sensitive" with the alignment strategy layer.

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
