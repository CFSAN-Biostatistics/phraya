# Phraya - Claude Context

## Project

Unified SNP caller for bacterial genomics. Reads/assemblies/hybrid input. Native Rust alignment algorithms (no BWA/MUMmer deps). Platform-optimized (SIMD, GPU).

## Status

Pre-dev. PRD done. No implementation yet. Workspace set up, crates empty.

## Workspace Structure

```
phraya-core/     # Evidence types, scoring, errors
phraya-align/    # Aligners: ContigAligner, ShortReadAligner, HybridAligner, LongReadAligner
phraya-index/    # Reference indexes: FM-index, FMD-index, k-mer sketch
phraya-io/       # I/O: FASTQ/FASTA/VCF/BAM, .phraya binary format
phraya-filter/   # Post-alignment filtering, presets
phraya-cli/      # Binary CLI
```

## Key Design Points

- **Zero binary deps**: All alignment in Rust. No runtime dep on external aligners.
- **Strategy continuum**: `exact|balanced|fast|sketch` selectable at runtime, same binary.
- **Platform-native**: x86 (SSE4.2/AVX2/AVX-512 runtime dispatch), ARM64 (NEON), GPU (CUDA opt-in).
- **Evidence transparency**: Structured `VariantEvidence` objects with full source-aware confidence metadata. Users filter with explicit thresholds, not opaque defaults.

## Alignment Algorithms

- **ContigAligner**: FM-index or minimizer seed → WFA extension
- **ShortReadAligner**: FMD-index seed → banded Smith-Waterman → paired-end rescue
- **LongReadAligner**: Minimizer seed → adaptive-band WFA
- **HybridAligner**: Contig align first, targeted read deployment at uncertain sites

## SIMD Strategy

- x86: `multiversion` crate for runtime dispatch (SSE4.2/AVX2/AVX-512 compiled into same binary)
- ARM64: NEON unconditional (mandatory on aarch64)
- Kernels: WFA diagonal extension, Smith-Waterman cell fill, minimizer hash

## Phases

1. **Phase 1** (8wk): `ContigAligner`, assembly-only, SSE/NEON baseline, `balanced` strategy
2. **Phase 2** (10wk): `ShortReadAligner`, AVX2/AVX-512 dispatch, `exact`/`fast`, `phraya filter`
3. **Phase 3** (12wk): `HybridAligner`, CUDA batch, mixed manifest, `sketch`, long read scaffold
4. **Phase 4** (14wk): `LongReadAligner` complete, docs, benchmarks, v1.0

## Coding Conventions

- Stable Rust (not nightly for library crates)
- `thiserror` for errors
- `rayon` for parallelism
- `serde` for serialization
- No unsafe unless SIMD intrinsics or FFI (CUDA), then document safety invariants
- Tests alongside code (`#[cfg(test)] mod tests`)
- Benchmarks in `benches/` (criterion)

## Don't Implement Yet

This is project setup only. No algorithm implementation until Phase 1 explicitly starts.
