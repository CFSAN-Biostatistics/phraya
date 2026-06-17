# Phraya

**Pairwise Haplotype-Resolved Alignment, Yield Afterward**

*The confluence of the two great streams.*

General-purpose pairwise sequence aligner for bacterial genomics. Short reads, long reads, contigs - Phraya aligns all your data. Produces rich alignment superpositions with deferred filtering for SNP calling, in-silico typing, classification, and other downstream analyses. Zero binary dependencies. Native Rust implementation with SIMD optimization (AVX2/NEON).

## Status

**Phase 1 MVP in development.** Architecture revision completed 2026-05-27. See [issue #58](https://github.com/CFSAN-Biostatistics/phraya/issues/58) for PRD.

## Installation

```bash
cargo install --git https://github.com/CFSAN-Biostatistics/phraya --locked phraya-cli
```

This installs the `phraya` binary using Rust's portable SIMD path. On ARM64 (Graviton, Apple Silicon), NEON is always active. On x86-64, a scalar fallback is used — portable builds run at approximately 40–60% the speed of a native SIMD build.

For full AVX2 acceleration on x86-64:

```bash
RUSTFLAGS="-C target-cpu=native" cargo install --git https://github.com/CFSAN-Biostatistics/phraya --locked phraya-cli
```

Requires Rust 1.75+. No external binary dependencies (BWA, minimap2, samtools, htslib).

## Philosophy

Most aligners force you to choose filtering parameters (mapping quality, coverage thresholds, multi-mapping behavior) before seeing results. Wrong assumptions mean expensive re-alignment.

Phraya separates **alignment computation** from **filtering decisions**:

1. **Plan**: Analyze inputs, detect use case, build k-mer evidence index
2. **Align**: Execute alignments, store rich metadata (multi-mapping, CIGAR, coverage)
3. **Filter**: Experiment with different parameters without re-aligning

Cache alignment results. Try different filters. Reuse for multiple downstream analyses.

## Pipeline

```bash
# 1. Create alignment plan (detects use case, builds k-mer evidence)
phraya plan --inputs reads/*.fastq --reference ref.fasta --output cohort.phrayaplan

# 2. Extract task list for parallel execution
phraya plan-tasks cohort.phrayaplan > tasks.tsv
cat tasks.tsv | parallel --colsep '\t' phraya align cohort.phrayaplan {1} {2}

# 3. Merge results from multiple samples
phraya merge sample_*.phraya --output cohort_merged.phraya

# 4. Filter with different parameters (no re-alignment needed)
phraya filter cohort_merged.phraya --min-coverage 10 --min-mapq 30 --format vcf > variants.vcf
phraya filter cohort_merged.phraya --min-coverage 5 --max-multi-map-fraction 0.3 --format tsv > variants.tsv
```

## Use Cases

Phraya automatically detects your workflow:

- **Case 2** (main use case): N reads + reference → traditional BWA-like alignment
- **Case 3** (key innovation): M contigs + N reads, no reference → selects centroid, aligns all to it
- **Case 4**: M contigs ± reference → minimap2-like contig alignment

## Key Features

- **Multi-mapping storage**: Tracks alternative alignment positions (score ratio ≥ 0.95). Filter ambiguous variants post-hoc.
- **Evidence-informed**: K-mer uniqueness and variation hotspots computed before alignment.
- **Rich metadata**: Every variant observation includes CIGAR, mapping quality, edit distance, local coverage (±50bp), all alleles, provenance.
- **Coverage tracks**: Quantized to nearest 5, RLE-compressed, full reference length.
- **Mergeable format**: Combine samples with order-independent merge preserving provenance.
- **Library-first filtering**: `phraya-filter` crate exposes public API for custom tools.
- **Parallel-ready**: Plan files emit task lists for GNU Parallel, SLURM, WDL, Nextflow.

## Architecture

Workspace with 5 crates:

- **phraya-core**: Core types (Sequence, VariantObservation, EvidenceLayer, CoverageTrack, MinimizerSketch), errors, k-mer sketching (via simd-minimizers), centroid selection, k-mer uniqueness
- **phraya-io**: FASTA/FASTQ/BAM/CRAM parsing, `.phrayaplan`/`.phraya`/`.phraya.queries` formats (MessagePack + zstd)
- **phraya-align**: WFA extension, SIMD diagonal fill (SSE4.2/NEON), seeding from minimizer sketches
- **phraya-filter**: Filtering library (threshold/expression/preset), output formatters (VCF/TSV/phraya)
- **phraya-cli**: Binary CLI (plan/plan-tasks/align/filter subcommands)

## File Formats

- **`.phrayaplan`** (v2): Plan file (k-mer evidence + task list). Read-only during alignment. Binary MessagePack + zstd.
- **`.phraya`**: Position index (variant observations + coverage track). Mergeable. Binary MessagePack + zstd.
- **`.phraya.queries`**: Query index (multi-mapping alternatives per read). Sidecar file. Binary MessagePack + zstd.

## Installation

### Prebuilt binaries (recommended)

Download the tarball for your platform from [GitHub Releases](https://github.com/CFSAN-Biostatistics/phraya/releases):

| Tarball | OS | Arch | SIMD | Use when |
|---------|----|------|------|----------|
| `phraya-*-x86_64-linux-gnu-native.tar.gz` | Linux | x86_64 | AVX2 | Modern x86_64 Linux (≥2013 CPUs — Haswell/Excavator or newer) |
| `phraya-*-x86_64-linux-gnu-portable.tar.gz` | Linux | x86_64 | SSE4.2 | Any x86_64 Linux; broadest compatibility |
| `phraya-*-aarch64-linux-gnu.tar.gz` | Linux | ARM64 | NEON | AWS Graviton, Ampere Altra, ARM servers |
| `phraya-*-x86_64-darwin.tar.gz` | macOS | Intel | AVX2 | Intel Mac |
| `phraya-*-aarch64-darwin.tar.gz` | macOS | Apple Silicon | NEON | M1/M2/M3/M4 Mac |

```bash
tar xzf phraya-*-x86_64-linux-gnu-native.tar.gz
./phraya --version
```

**Portable vs native (x86_64 Linux):** The native build uses AVX2 via
`-C target-cpu=x86-64-v3` and is **~2× faster for k-mer sketching** thanks to the
simd-minimizers AVX2 path. Use it on any CPU from ~2013 onward. If it exits with
`Illegal instruction`, fall back to the portable build (SSE4.2 baseline, runs on
every x86_64 CPU since ~2008).

**ARM builds** (Linux ARM64 and Apple Silicon) always use NEON — there is no
portable/native split because NEON is mandatory on AArch64 and always available.

## Docker Quick Start

```bash
# Pull the latest image (amd64 and arm64 supported)
docker pull ghcr.io/cfsan-biostatistics/phraya:latest

# Verify installation
docker run --rm ghcr.io/cfsan-biostatistics/phraya:latest --version

# Run with your data (mount current directory as /data)
docker run --rm -v $(pwd):/data ghcr.io/cfsan-biostatistics/phraya:latest \
    plan --inputs /data/reads/*.fastq --reference /data/ref.fasta --output /data/cohort.phrayaplan

docker run --rm -v $(pwd):/data ghcr.io/cfsan-biostatistics/phraya:latest \
    align /data/cohort.phrayaplan query_id target_id

docker run --rm -v $(pwd):/data ghcr.io/cfsan-biostatistics/phraya:latest \
    filter /data/cohort.phraya --min-coverage 10 --min-mapq 30 --format vcf > variants.vcf
```

### Available tags

| Tag | Description |
|-----|-------------|
| `latest` | Most recent release |
| `v1.2.3` | Exact version |
| `v1.2` | Latest patch for minor version |

### SIMD in Docker

The Docker image is built with the **SSE4.2 baseline** (`-C target-feature=+sse4.2`) rather than `-C target-cpu=native`. This ensures the image runs on any modern x86-64 CPU but does not use AVX2 acceleration for k-mer sketching.

**For HPC workloads** where you control the hardware, building from source with `-C target-cpu=native` will enable AVX2 (x86-64) or NEON (ARM64) and improve sketching throughput:

```bash
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

### Building from source

```bash
cargo build --release
```

Requires Rust 1.75+. No external binary dependencies (BWA, minimap2, samtools).

For best k-mer sketching performance, enable native SIMD:

```bash
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

Without `-C target-cpu=native`, simd-minimizers falls back to a scalar path and is slower. On ARM64 (Graviton, Apple Silicon), NEON is always enabled.

## Dependencies

Phraya depends on [**simd-minimizers**](https://github.com/ragnargrootkoerkamp/simd-minimizers) for k-mer sketching and seeding. This library implements SIMD-accelerated canonical minimizers using AVX2 (x86-64) or NEON (ARM64), and is described in:

> Ragnar Groot Koerkamp, Igor Martayan. **SimdMinimizers: Computing random minimizers, fast.** *SEA 2025.* doi:[10.4230/LIPIcs.SEA.2025.20](https://doi.org/10.4230/LIPIcs.SEA.2025.20)

We use canonical minimizers with default parameters k=21, w=11 (appropriate for bacterial genomics) and ntHash rolling hashes. Sketches are computed once during `phraya plan` and reused during `phraya align`, eliminating redundant computation.

## Design Decisions

- **Score ratio threshold**: 0.95 (hard-coded). Stores alternatives within 95% of best identity. Opinionated choice for storage efficiency.
- **K-mer parameters**: k=21, w=11 (canonical minimizers, standard for bacterial genomes). l = w+k-1 = 31 satisfies the odd-l canonicality requirement of simd-minimizers.
- **Coverage quantization**: Nearest 5. Enables RLE compression, negligible precision loss for variant calling decisions.
- **Sketch reuse**: Plan-time sketches stored in `.phrayaplan` (v2) keyed by sequence ID; alignment reuses them rather than recomputing.

## Inspired By

- **Deacon** (https://github.com/bede/deacon): General-purpose aligner with flexible post-processing.

Phraya differentiates on:
- Richer intermediate format (more cacheable/reusable)
- More deferred parameters (multi-mapping, coverage computed during alignment, filtered post-hoc)
- Case 3 (contigs + reads without reference via centroid selection)

## Phase 1 MVP Scope

**In scope:**
- Cases 2 (reads + ref), 3 (contigs + reads), 4 (contigs only)
- K-mer evidence (uniqueness only)
- Threshold-based filtering
- VCF/TSV/phraya output formats
- Library API (phraya-filter)

**Phase 2+:**
- Case 1 (read MSA without reference)
- Expression-based filters (`--expr`)
- Named presets (`--preset conservative`)
- Variation hotspot estimation
- Python/R bindings
- GPU acceleration

## For Maintainers

### Release Process

Phraya uses a tag-triggered release workflow. Push a `v*` tag and all automated channels update within ~15 minutes.

```bash
git tag v0.2.0
git push origin v0.2.0
```

#### Automated channels (triggered by tag push)

| Channel | What happens | Secret required |
|---------|-------------|-----------------|
| **GitHub Releases** | 5 prebuilt binaries uploaded (portable + native Linux x86_64, Linux ARM64, macOS Intel, macOS M1) | `GITHUB_TOKEN` (automatic) |
| **Docker** | Multi-arch image pushed to `ghcr.io/cfsan-biostatistics/phraya` with `:latest` + versioned tags | `GITHUB_TOKEN` (automatic) |
| **crates.io** | All 5 crates published in dependency order | `CARGO_REGISTRY_TOKEN` |

Pre-releases (tags containing `-rc`, `-alpha`, `-beta`) skip crates.io publish and do not update the `:latest` Docker tag.

#### Manual channels (require external PRs)

**Bioconda** (`bioconda-recipes` repo):
1. Fork [bioconda/bioconda-recipes](https://github.com/bioconda/bioconda-recipes)
2. Update `recipes/phraya/meta.yaml` — bump `version`, update `sha256` from the GitHub Release SHA256SUMS.txt
3. Open PR to `bioconda/bioconda-recipes`

**Homebrew** (if using a tap rather than homebrew-core):
1. Update `Formula/phraya.rb` in the tap repo — bump `version` and `sha256`
2. Test locally: `brew install --build-from-source Formula/phraya.rb`
3. Commit and push; Homebrew users get the update on next `brew update`

#### Required secrets

- `CARGO_REGISTRY_TOKEN`: crates.io API token for publishing. Set in repo Settings → Secrets → Actions.
- `GITHUB_TOKEN`: Automatically provided by GitHub Actions. No setup needed.

#### Verifying a release

After the workflow completes, verify all channels:

```bash
# GitHub Releases: check all 5 binaries exist
gh release view v0.2.0 --json assets --jq '[.assets[].name]'

# Docker
docker pull ghcr.io/cfsan-biostatistics/phraya:v0.2.0
docker run --rm ghcr.io/cfsan-biostatistics/phraya:v0.2.0 --version

# crates.io: package page should show new version
# https://crates.io/crates/phraya-cli
```

#### Platform binary selection guide

| Platform | Recommended binary | Notes |
|----------|--------------------|-------|
| Linux x86_64, HPC cluster | `phraya-linux-x86_64-native` | AVX2, ~2× faster k-mer sketching |
| Linux x86_64, older hardware | `phraya-linux-x86_64-portable` | SSE4.2 baseline, runs everywhere |
| Linux ARM64 (Graviton, Raspberry Pi) | `phraya-linux-aarch64` | NEON, always enabled |
| macOS Intel | `phraya-macos-x86_64` | AVX2 |
| macOS Apple Silicon | `phraya-macos-aarch64` | NEON |
| Container / unknown CPU | Docker image | Portable SSE4.2 build |

## License

Unlicense. As a work product of the US Government (17 USC 105), Phraya is in the public domain.

## Contributing

See [issue #58](https://github.com/CFSAN-Biostatistics/phraya/issues/58) for Phase 1 PRD. Implementation contributions welcome.
