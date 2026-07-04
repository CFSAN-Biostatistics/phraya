# Changelog

All notable changes to Phraya are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Phraya uses [semantic versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
