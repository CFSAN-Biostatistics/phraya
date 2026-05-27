# Phraya - Claude Context

## Project

General-purpose pairwise sequence aligner for bacterial genomics. Reads/assemblies/hybrid input. Native Rust alignment algorithms (no BWA/MUMmer deps). Platform-optimized SIMD (SSE4.2/AVX2/NEON via SimdMinimizers). Produces rich alignment superpositions with deferred filtering for SNP calling, in-silico typing, classification, and other downstream analyses.

## Status

Architecture revision complete (2026-05-27). Breaking redesign supersedes original issues #1-20. Ready for Phase 1 MVP implementation.

## Workspace Structure

```
phraya-core/     # Core types (Sequence, VariantObservation, EvidenceLayer), errors
phraya-align/    # WFA extension, SIMD diagonal fill (SSE4.2/NEON)
phraya-index/    # SimdMinimizers sketches, seeding, k-mer uniqueness
phraya-io/       # I/O: FASTQ/FASTA, .phraya/.phrayaplan binary formats
phraya-filter/   # Post-alignment filtering library (threshold/expression/preset)
phraya-cli/      # Binary CLI: plan/align/filter subcommands
```

## Key Design Points

- **Zero binary deps**: All alignment in Rust using SimdMinimizers + WFA. No runtime dep on BWA/minimap2/MUMmer.
- **Multi-mapping storage**: Query index tracks alternative alignment positions (score ratio ≥ 0.95 normalized edit distance). Enables filtering on mapping uniqueness post-alignment.
- **Deferred filtering**: Alignment produces rich `.phraya` files (multi-mapping, CIGAR, coverage tracks, k-mer uniqueness). Filter parameters applied post-hoc, not during alignment.
- **Library-first**: `phraya-filter` crate exposes filtering API. CLI is thin wrapper. Enables Python bindings, R integration, custom pipelines.
- **Evidence-informed alignment**: `.phrayaplan` files contain k-mer landscape + variation hotspots estimated from input sequences before alignment begins.
- **Platform-native SIMD**: SimdMinimizers handles AVX2/NEON dispatch. WFA diagonal fill uses SSE4.2/NEON intrinsics.

## Pipeline

```
phraya plan   → .phrayaplan (k-mer evidence + task list)
              ↓
phraya align  → .phraya (position index) + .phraya.queries (query index with multi-mapping)
              ↓
phraya filter → VCF | TSV | filtered .phraya
```

### Use Cases (detected by `phraya plan`)

1. **Case 1**: N reads, no reference → MSA (N×(N-1)/2 alignments) [Phase 2+]
2. **Case 2**: N reads + reference → N alignments (BWA-like, main use case)
3. **Case 3**: M contigs + N reads, no ref → centroid selection + M+N alignments (key innovation)
4. **Case 4**: M contigs ± reference → M or M×(M-1)/2 alignments (minimap2-like)

### Alignment Algorithm

- **Seeding**: SimdMinimizers sketches (k=21, w=11 default) → find shared minimizers
- **Extension**: WFA with SIMD-accelerated diagonal fill (SSE4.2/NEON)
- **Scoring**: Multi-mapping score ratio = (1 - edit_dist/query_len) for primary vs alternatives
- **Output**: Store all alignments with score_ratio ≥ 0.95 (hard-coded opinion)

## File Formats

### `.phrayaplan` (binary MessagePack + zstd)
- **Read-only** during alignment, transmitted to all workers
- Contains: metadata (use case, input files), k-mer index (SimdMinimizers), variation hotspots, task list
- Evidence estimated from k-mers only (no expensive alignment first)
- CLI tool: `phraya plan-tasks` dumps task list for GNU Parallel/xargs/WDL/Nextflow

### `.phraya` (position index, binary MessagePack + zstd)
- **VariantObservation** fields: position, ref_base, all_alleles (with counts), confidence, CIGAR, mapq, edit_distance, local_coverage (±50bp window), avg_base_quality, provenance
- **Coverage track**: quantized to nearest 5, RLE-compressed, full reference length
- **Mergeable**: combine multiple samples via position-centric merge (order-independent)
- Requires `libphraya` to parse (data-dense for caching/long-term storage)

### `.phraya.queries` (query index, sidecar)
- For each query: list of alignment positions + scores above threshold (score_ratio ≥ 0.95)
- Enables multi-mapping analysis: "exclude variants where >50% supporting reads multi-map"
- Separate file to keep merge fast (not needed for typical variant calling)

## Evidence Layer

### Plan-time (k-mer estimated)
- K-mer uniqueness: which k-mers appear in 1 vs many sequences
- Jaccard similarity: pairwise contig similarity from shared minimizers
- Divergence hotspots: regions with few shared k-mers between similar sequences
- Coverage estimate: k-mer depth across reference positions

### Alignment-time (ground truth)
- Per-position: multi_map_fraction, avg_score_ratio_gap to next-best alignment
- Polymorphic sites: all alleles with counts, reference base
- Invariant positions: where all samples match reference

## Filter Operations

### Supported Styles (all chainable)
1. **Threshold-based**: `--min-coverage 10 --min-mapq 30 --max-multi-map-fraction 0.3`
2. **Expression-based**: `--expr "coverage >= 10 && mapq > 30"`
3. **Named presets**: `--preset conservative|sensitive` with optional `--override`
4. **Pipeline composition**: `phraya filter --min-coverage 10 | phraya filter --snps-only`

### Output Formats
- Filtered `.phraya` (subset of records, same format for chaining)
- VCF (standard variant calling output)
- TSV/CSV (arbitrary column selection)

### Feature Space (available for filtering)
- Alignment: coverage, mapq, CIGAR complexity, edit_distance, multi_map_fraction, score_ratio_gap
- Context: edge_distance, local_gc, k-mer_uniqueness, in_homopolymer, in_tandem_repeat, snp_density (15bp/125bp/1000bp windows)
- Alleles: allele_frequency, ref_base, alt_bases
- Quality: avg_base_quality, confidence

## Phase 1 MVP (Current Target)

**Goal**: Validate architecture end-to-end for all 4 use cases with core primitives.

**In scope:**
- `.phrayaplan` format (binary MessagePack + zstd)
- `phraya plan` CLI (detects 4 cases, emits plan with k-mer evidence + tasks)
- `phraya plan-tasks` CLI (dumps task list for parallel execution)
- Query index with multi-mapping (score ratio ≥ 0.95)
- Richer `VariantObservation` (CIGAR, mapq, local coverage, all alleles)
- Coverage tracks (quantized to 5, RLE-compressed)
- `phraya align` reads plan, writes `.phraya` + `.phraya.queries`
- `phraya filter` threshold-based only, outputs VCF/TSV/phraya
- `phraya-filter` crate (library-first API)
- Case 2 (reads + ref) fully working
- Case 3 (contigs + reads) with centroid selection
- Case 4 (contigs only) basic support

**Simplified for MVP:**
- Evidence = k-mer uniqueness only (no hotspot estimation yet)
- Filter = thresholds only (no expressions/presets yet)
- No MSA case 1 (N reads, no ref) - requires assembly or centroid from reads

**Deferred to Phase 2+:**
- Expression-based filters (`--expr`)
- Named presets (`--preset conservative`)
- Variation hotspot estimation in plan
- Case 1 (MSA without reference)
- Two-tier evidence (k-mer → alignment refinement)
- Python bindings / R integration
- Advanced coverage: per-base breakdown (not just ±50bp window)

## Dependencies

- **SimdMinimizers** (`simd-minimizers` crate): AVX2/NEON minimizer sketching via two-stacks sliding window
- **MessagePack** (`rmp-serde`): binary serialization for `.phraya`/`.phrayaplan` formats
- **Compression** (`zstd`): for plan/alignment file storage
- **Parallelism** (`rayon`): for embarrassingly parallel tasks within single-node execution
- Standard: `thiserror` (errors), `serde` (serialization), `clap` (CLI)

## Coding Conventions

- Stable Rust (not nightly for library crates)
- No unsafe unless SIMD intrinsics (WFA diagonal fill), then document SAFETY invariants
- Tests alongside code (`#[cfg(test)] mod tests`)
- Benchmarks in `benches/` (criterion)
- Library-first: `phraya-filter` crate exposes public API, CLI is thin wrapper

## Implementation Notes

- **Breaking redesign**: Original issues #1-20 superseded by this architecture
- **Minimizer library swap**: Replace custom implementation with `simd-minimizers` crate (tactical improvement, doesn't affect architecture)
- **Score ratio hard-coded**: 0.95 threshold for multi-mapping is Phraya's opinion, not user-configurable (aligns with "defer parameters" philosophy for alignment, but this one we enforce)
- **Centroid selection**: For case 3 (contigs + reads, no ref), select contig closest to k-mer space center (median Jaccard) as reference coordinate space
