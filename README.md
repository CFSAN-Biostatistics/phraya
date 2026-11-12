# Phraya

**Pairwise Haplotype-Resolved Alignment, Yield Afterward**

*The confluence of the two great streams.*

General-purpose pairwise sequence aligner for bacterial genomics. Short reads, long reads, contigs - Phraya aligns all your data. Produces rich alignment superpositions with deferred filtering for SNP calling, in-silico typing, classification, and other downstream analyses. Zero binary dependencies. Native Rust implementation with SIMD optimization (AVX2/NEON).

## Status

**Phase 1 MVP shipped (2026-06-06).** Architecture revision completed 2026-05-27. See [issue #58](https://github.com/CFSAN-Biostatistics/phraya/issues/58) for PRD.

## Installation

```bash
cargo install --git https://github.com/CFSAN-Biostatistics/phraya --locked phraya
```

This installs the `phraya` binary using Rust's portable SIMD path. On ARM64 (Graviton, Apple Silicon), NEON is always active. On x86-64, a scalar fallback is used — portable builds run at approximately 40–60% the speed of a native SIMD build.

For full AVX2 acceleration on x86-64:

```bash
RUSTFLAGS="-C target-cpu=native" cargo install --git https://github.com/CFSAN-Biostatistics/phraya --locked phraya
```

Requires Rust 1.75+. No external binary dependencies (BWA, minimap2, samtools, htslib).

## Philosophy

Most aligners force you to choose filtering parameters (mapping quality, coverage thresholds, multi-mapping behavior) before seeing results. Wrong assumptions mean expensive re-alignment.

Phraya separates **alignment computation** from **filtering decisions**:

1. **Plan**: Analyze inputs, detect input shape, build k-mer evidence index against content-hashed reference spaces
2. **Align**: Execute alignments, store rich metadata (multi-mapping, CIGAR, coverage)
3. **Filter**: Experiment with different parameters without re-aligning

Cache alignment results. Try different filters. Reuse for multiple downstream analyses.

## Pipeline

```bash
# 1. Create alignment plan (detects input shape, builds k-mer evidence)
phraya plan --inputs reads/*.fastq --reference ref.fasta --output cohort.phrayaplan

# 2. Extract task list for parallel execution
phraya plan-tasks cohort.phrayaplan > tasks.tsv
cat tasks.tsv | parallel --colsep '\t' phraya align cohort.phrayaplan {1} {2}

# 3. Merge results from multiple samples
phraya merge sample_*.phraya --output cohort_merged.phraya

# 4. Filter with different parameters (no re-alignment needed)
phraya filter cohort_merged.phraya --min-coverage 10 --min-mapq 30 --format vcf > variants.vcf
phraya filter cohort_merged.phraya --min-coverage 5 --min-identity 0.98 --format tsv > variants.tsv
```

## Inputs

Phraya aligns reads or contigs against one or more content-hashed reference spaces (a single
reference is the degenerate case). The input shape — reads vs contigs, presence of a reference —
is auto-detected at plan time.

- **Reads + reference** (`reads-with-ref`): traditional BWA-like alignment
- **Contigs + reads** (`contigs-with-reads`, no reference): selects a centroid contig, aligns all inputs to it
- **Contigs** (`contigs-only`, with or without reference): minimap2-like contig-to-contig alignment

All shapes work over DNA or protein sequences — `phraya plan --alphabet {auto|dna|protein}`
auto-detects amino-acid input from content by default (ADR-0013); no separate command or
workflow for protein alignment.

### Assembly vs assembly (dnadiff-style SNP comparison)

```bash
phraya plan --reference ref_assembly.fasta --inputs query_assembly.fasta --output pair.phrayaplan
phraya align --reference ref_assembly.fasta pair.phrayaplan --output aln/ --strategy assembly
phraya filter aln/*.phraya --format snpdiffs \
    --reference-fasta ref_assembly.fasta --query-fasta query_assembly.fasta \
    --reference-id REF --query-id QRY --output pair.snpdiffs
```

`--strategy assembly` aligns every homologous region of each query contig as its own local
block (one `.phraya` per reference contig). `--format snpdiffs` selects blocks one-to-one
across all reference contigs (like `delta-filter -1`) and writes CSP2's `.snpdiffs` format:
one BED row per kept block, `.` rows for uncovered regions, and one row per variant inside a
kept block, with header coverage computed from those blocks (ADR-0015).

### Depletion and typing against several reference spaces

```bash
phraya plan --inputs reads.fastq --reference host.fasta --reference target.fasta --output dep.phrayaplan
phraya align dep.phrayaplan --reference host.fasta --reference target.fasta --output aln/
# aln/<space>.phraya per space + aln/cross_space.phraya.queries
phraya filter aln/target.phraya --queries-sidecar aln/cross_space.phraya.queries --deplete-space host
phraya filter --queries-sidecar aln/cross_space.phraya.queries --classify-by-space > typing.tsv
```

Align keeps every read's placements in every space (ADR-0011); `filter` decides. Reads whose
best two spaces are within `--min-margin` identity are kept as ambiguous rather than depleted.

## CLI Reference

`phraya <COMMAND> --help` is authoritative; this section mirrors it.

### `phraya plan`

| Flag | Description |
|---|---|
| `--inputs <FILE>` | Input files (FASTA/FASTQ/BAM/CRAM). Required |
| `--reference <FILE>` | Reference file, repeatable. Every record is its own reference space; every query is tasked against every reference record |
| `--output <FILE>` | Output `.phrayaplan`. Required |
| `--alphabet {auto\|dna\|protein}` | Sequence alphabet (ADR-0013). Default `auto` (content-sniffed) |
| `--min-homology <F>` | Minimum sketch Jaccard for a contigs-only all-pairs task. Default `0.1`; ignored with a reference |
| `--no-homology-gate` | Generate the full dense all-pairs contig task list |
| `--chunks <N>` | Pre-split read sketches into N seekable chunk frames. Default `1` |
| `--batch-to <N>` / `--batch-by <N>` | Batch mode: split reads into N chunks / N reads per chunk |
| `--batch-output-pattern <PATTERN>` | Per-worker output pattern with `{worker}`; required with `--batch-to`/`--batch-by` |

### `phraya plan-tasks <PLAN_FILE>`

Prints the plan's task list as TSV (query ID, target ID) for GNU Parallel, SLURM, WDL or Nextflow.

### `phraya align <PLAN_FILE> [QUERY_ID] [TARGET_ID]`

Three mutually exclusive modes: a single task (`QUERY_ID TARGET_ID`), batch (`--worker`/`--ensure`), or reference palette (`--reference`).

| Flag | Description |
|---|---|
| `--output <FILE>` | Output `.phraya`; a directory in `--reference` mode; unused with `--worker` |
| `--strategy <STRATEGY>` | Strategy preset, default `balanced` (see below) |
| `--coverage-window <N>` | Override the per-variant coverage-window radius (bp) |
| `--gap-model {linear\|affine}` | Gap cost model (ADR-0014). Only with `--strategy sensitive`, where it defaults to `affine` |
| `--worker <N>` | Batch mode: process chunk N (0-indexed) |
| `--ensure` | Batch mode: process all missing chunks sequentially |
| `--threads <N>` | Threads for `--ensure`. Default: auto-detect |
| `--reference <FILE>` | Palette mode (ADR-0011), repeatable: align the plan's reads against each presented reference space, resolved by content hash. Writes `<output>/<space>.phraya` plus `<output>/cross_space.phraya.queries` |
| `--sealed` | With `--reference`: a reference whose content hash is not in the plan is a hard error instead of a warning |

| Strategy | Chain cap K | Engine | Coverage window | Use |
|---|---|---|---|---|
| `fast` | 1 | Myers, divergence cutoff | ±150 bp | Confident low-divergence reads; low sensitivity |
| `balanced` | 2 | Myers fitting + WFA fallback | ±50 bp | Default |
| `sensitive` | 50 | Seeded gap-affine WFA | ±25 bp | Divergent reads, indels |
| `exact` | 250 | Full WFA | ±5 bp | Highest precision; small genomes, validation |
| `long-read` | 100 | WFA + chunking | ±25 bp | ONT/PacBio reads, 10–50 kb |
| `assembly` | every chain | Local blocks, WFA gap fills | ±5 bp | Contig vs contig (ADR-0015) |

### `phraya merge --output <FILE> <FILE>...`

Merges `.phraya` files (order-independent, provenance and assembly blocks preserved).

### `phraya filter [INPUT]...`

| Flag | Description |
|---|---|
| `--format <FORMAT>` | `vcf` (default), `tsv`, `phraya`, `snpdiffs`, `core-snp`, `distance-matrix` |
| `--output <FILE>` | Output file; required for `--format phraya` |
| `--preset {strict\|tolerant}` | Named threshold preset; individual flags override it |
| `--expr <EXPR>` | Extra conjunctive predicate, e.g. `"identity >= 0.98 && coverage >= 10"` |
| `--min-coverage` / `--max-coverage <N>` | Coverage bounds |
| `--min-mapq` / `--max-mapq <N>` | MAPQ bounds |
| `--min-kmer-uniqueness <F>` | Minimum k-mer uniqueness (0.0–1.0) |
| `--min-identity <F>` | Minimum `1 − edit_distance / query_aligned_len` |
| `--min-match-fraction <F>` | Minimum `M / (M+X+I+D)` |
| `--min-aligned-length <N>` | Minimum alignment columns |
| `--max-snp-density-15` / `-125` / `-1000 <N>` | Maximum other variant positions in a 15/125/1000 bp window |
| `--min-edge-distance <N>` | Minimum distance to the query's nearest end |
| `--exclude-discordant-pairs` | Drop pairs with insert size beyond mean ± `--discordant-sigma` (default `3.0`) σ |
| `--require-proper-pairs <N>` | Minimum proper-pair fraction of covering reads (`0.9` or `90`) |
| `--min-insert-size` / `--max-insert-size <N>` | Absolute insert-size bounds |
| `--require-both-mates-mapped` | Require both mates mapped |
| `--reference-fasta` / `--query-fasta <FILE>` | `snpdiffs`: assemblies for header metadata, coverage denominators, uncovered rows |
| `--reference-id` / `--query-id <NAME>` | `snpdiffs`: header labels; default the FASTA file names |
| `--sample-ids <ID1,ID2,...>` | `core-snp`/`distance-matrix`: one sample ID per input file |
| `--queries-sidecar <FILE>` | Cross-space `.phraya.queries`; required by the two flags below |
| `--deplete-space <SPACE>` | Drop reads that confidently align to this space; near-ties are kept |
| `--min-margin <F>` | Identity gap below which two placements are near-ties. Default `0.02` |
| `--classify-by-space` | Typing mode: TSV of `read_id, best_space, identity, margin` |

### `phraya qc <FILE>...`

One TSV row per `.phraya` file: variant count, coverage breadth, k-mer uniqueness summary.

## Key Features

- **Multi-mapping storage**: Tracks alternative alignment positions (score ratio ≥ 0.95). Filter ambiguous variants post-hoc.
- **Evidence-informed**: K-mer uniqueness and variation hotspots computed before alignment.
- **Rich metadata**: Every variant observation includes CIGAR, mapping quality, edit distance, local coverage (±50bp), all alleles, provenance.
- **Coverage tracks**: Depths below 5 stored exact; 5 and above quantized to nearest 5. RLE-compressed, full reference length.
- **Mergeable format**: Combine samples with order-independent merge preserving provenance.
- **Library-first filtering**: `phraya-filter` crate exposes public API for custom tools.
- **Parallel-ready**: Plan files emit task lists for GNU Parallel, SLURM, WDL, Nextflow.

## Architecture

Workspace with 5 crates:

- **phraya-core**: Core types (Sequence, VariantObservation, EvidenceLayer, CoverageTrack, MinimizerSketch), errors, k-mer sketching (via simd-minimizers), centroid selection, k-mer uniqueness
- **phraya-io**: FASTA/FASTQ/BAM/CRAM parsing, `.phrayaplan`/`.phraya`/`.phraya.queries` formats (MessagePack + zstd)
- **phraya-align**: WFA (linear and gap-affine) and Myers extension, SIMD match extension (SSE2/NEON baseline, AVX2/AVX-512 runtime dispatch), seeding from minimizer sketches, assembly blocks
- **phraya-filter**: Filtering library (threshold/expression/preset), one-to-one block selection, output formatters (VCF/TSV/phraya/snpdiffs/core-SNP/distance matrix)
- **phraya** (`phraya-cli/`): Binary CLI (plan/plan-tasks/align/merge/filter/qc subcommands)

## File Formats

- **`.phrayaplan`** (v8): Plan file (k-mer evidence + task list + read/contig sketches + reference palette). Read-only during alignment. Chunk-addressable, seekable zstd frames; workers load only the chunk they need.
- **`.phraya`** (v2): Position index (variant observations + coverage track, plus the local alignment blocks of `--strategy assembly`). Mergeable. Binary MessagePack + zstd. v1 files still read.
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

We use canonical minimizers with default parameters k=21, w=11 (appropriate for bacterial genomics) and ntHash rolling hashes. Sketches are computed once during `phraya plan` and reused during `phraya align`, eliminating redundant computation. Protein-alphabet input (ADR-0013) uses non-canonical minimizers instead — there is no reverse complement for amino acids, so canonicalization is meaningless — at protein-scale k=6, w=5 defaults.

## Design Decisions

- **Score ratio threshold**: 0.95 (hard-coded). Stores alternatives within 95% of best identity. Opinionated choice for storage efficiency.
- **K-mer parameters**: k=21, w=11 (canonical minimizers, standard for bacterial genomes). l = w+k-1 = 31 satisfies the odd-l canonicality requirement of simd-minimizers. Protein alphabet: k=6, w=5, non-canonical.
- **Alphabet auto-detection**: content is classified `Protein` iff it contains any of E/F/I/L/P/Q — one-letter amino-acid codes with no IUPAC nucleotide meaning (every other amino-acid letter doubles as a valid DNA base or ambiguity code). `--alphabet {dna|protein}` overrides detection for the one genuine ambiguity: a short peptide spelled entirely in residues that are also valid nucleotide letters (e.g. only A/C/G/T-coding residues) is indistinguishable from DNA by content alone.
- **Gap-affine scoring in `sensitive`** (ADR-0014): `gap_open=4, gap_extend=1, mismatch=2` (hard-coded). `gap_open+gap_extend=5 > mismatch` keeps an isolated SNP cheaper than opening any gap; `gap_extend=1` makes a run of ≥2 gap bases cheaper than the equivalent run of mismatches once a gap is open, so real multi-base indels consolidate into one CIGAR op instead of scattering across cost-tied mismatches. The affine cost only steers *which* alignment the search picks — reported `edit_distance` stays the traditional mismatches+indel-bases count either way. `balanced`/`fast` are unaffected (Myers-primary, no affine mode); `--gap-model linear` opts `sensitive` back into the old uniform-cost path.
- **Coverage quantization**: Depths below 5 stored exact (the boundary that matters most for coverage-breadth thresholds); 5 and above round to nearest 5. Enables RLE compression with negligible precision loss for variant calling decisions.
- **Sketch reuse**: Plan-time sketches stored in `.phrayaplan` (v8) keyed by sequence ID; alignment reuses them rather than recomputing.

## Inspired By

- **Deacon** (https://github.com/bede/deacon): General-purpose aligner with flexible post-processing.

Phraya differentiates on:
- Richer intermediate format (more cacheable/reusable)
- More deferred parameters (multi-mapping, coverage computed during alignment, filtered post-hoc)
- Contigs + reads without an explicit reference: auto-selects a centroid to establish the coordinate space

## Implementation Status

**Phase 1 MVP** (shipped 2026-06-06):

- Reads + reference, contigs + reads (auto-centroid), contigs-only — all working end-to-end
- Automatic centroid selection for contigs + reads without a reference; `--reference` overrides it
- `.phrayaplan` v8 with sketch reuse (chunk-addressable, seekable zstd frames)
- Real WFA O(s·n) alignment (wavefront-based)
- SIMD-accelerated match extension (`count_matching_prefix`: SSE2/NEON/AVX2/AVX-512 with runtime dispatch)
- Local coverage (±50bp window), mapq, avg base quality, confidence from input data
- Tandem repeat detection: annotation on variants, `exclude_tandem_repeats` filter option
- Variation hotspot detection (low k-mer-uniqueness regions) computed at plan time
- `phraya filter`: threshold-based + named presets (strict / tolerant) + expression filters (`--expr`)
- QC-layer filter fields: `identity`, `match_fraction`, `aligned_length`, `edge_distance`, `snp_density_15/125/1000`; `phraya qc` coverage-breadth report
- `phraya merge` preserves per-observation QC fields (tandem repeat, variant type, k-mer uniqueness, strand, mate info)

**Since Phase 1** (unreleased):

- Strategy ladder: `fast`, `balanced`, `sensitive` (gap-affine), `exact`, `long-read`, `assembly`
- Multiple content-hashed reference spaces per plan, `align --reference` palette mode, `--sealed`
- `filter`: depletion (`--deplete-space`) and typing (`--classify-by-space`) over the cross-space sidecar
- `filter` formats: `snpdiffs` (CSP2/dnadiff replacement), `core-snp`, `distance-matrix`
- Packaging: GitHub Releases binaries, multi-arch Docker on ghcr.io, crates.io; Bioconda and Homebrew recipes in-repo, channel submissions pending (#167, #168)

**Planned:**

- Reads-only MSA without a reference
- Two-tier evidence (k-mer → alignment refinement)
- Python / R bindings
- Advanced coverage: per-base breakdown (not just ±50bp window)
- CRAM parsing with external reference (currently only reference-free/unmapped CRAM records)
- GPU acceleration — deferred (#242): profiling showed the alignment DP kernels are too small a share of runtime for a GPU to pay off

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
# https://crates.io/crates/phraya
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

### Commit conventions

Commits that implement an issue (fix a bug, add a feature, or complete Phase work) **must
close the corresponding GitHub issue**. GitHub auto-closes issues when a commit on the
default branch references one using a closing keyword:

```
feat: gap-affine scoring for sensitive strategy (closes #230)
fix: batch-mode worker panics with index out of bounds (fixes #234)
Phase 1: enhanced QC summary — add strategy/filter/k-mer stats columns (closes #246)
```

**Supported keywords**: `close`, `closes`, `closed`, `fix`, `fixes`, `fixed`, `resolve`,
`resolves`, `resolved`.

If you commit without the keyword, **you must manually close the issue** on GitHub before
considering the work complete. Stale open issues that were actually implemented waste
reviewer and CI time — checking for unclosed-but-implemented issues is a required step of
every task's verification checklist.
