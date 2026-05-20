# Phraya

**Probabilistic Heuristic for Running Awesome Yield-Agnostic Alignments**

Unified SNP calling for bacterial genomics. Accepts reads, assemblies, or any mixture. Zero binary dependencies. Native Rust implementation of alignment algorithms optimized for 2026 hardware (x86 AVX2/AVX-512, ARM64/Apple Silicon, optional GPU acceleration).

## Status

Pre-development. PRD draft complete. Implementation not started.

## Quick Start

```bash
# Not yet implemented
phraya snp --manifest samples.tsv --reference ref.fasta --strategy balanced
```

## Features (Planned)

- **Input agnostic**: Illumina reads (PE/SE), assemblies, long reads (ONT/PacBio), or hybrid
- **Zero dependencies**: No BWA, MUMmer, samtools, or external aligners required
- **Speed-sensitivity continuum**: `exact`, `balanced`, `fast`, `sketch` strategies in one binary
- **Platform-native performance**: SIMD (AVX2/AVX-512/NEON), multi-core (Rayon), GPU (CUDA opt-in)
- **Structured evidence**: Full transparency for downstream filtering

## Architecture

Workspace with 6 crates:
- `phraya-core`: Evidence types, scoring model, error types
- `phraya-align`: Alignment algorithms (FM-index, WFA, Smith-Waterman, SIMD kernels)
- `phraya-index`: Reference indexing (BWT, FMD-index, k-mer sketches)
- `phraya-io`: I/O (FASTQ/FASTA/VCF/BAM, .phraya native format)
- `phraya-filter`: Post-alignment filtering, presets
- `phraya-cli`: Binary CLI

## Building

```bash
cargo build --release
```

Optional GPU acceleration:
```bash
cargo build --release --features cuda
```

## License

Unlicense. As a work product of the US Government (17 USC 105), Phraya is in the public domain.

## Contributing

PRD and architecture review feedback welcome. Implementation contributions once Phase 1 starts.
