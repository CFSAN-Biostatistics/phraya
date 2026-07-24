# Phraya - Claude Context

## Project

General-purpose pairwise sequence aligner for bacterial genomics. Reads/assemblies/hybrid input. Native Rust alignment algorithms (no BWA/MUMmer deps). Platform-optimized SIMD (AVX2/NEON sketching via simd-minimizers; portable SIMD alignment via `wide` crate). Produces rich alignment superpositions with deferred filtering for SNP calling, in-silico typing, classification, and other downstream analyses.

## Status

Phase 1 MVP shipped (2026-06-06). Cases 2, 3, 4 working end-to-end. Real WFA O(s·n) alignment, named filter presets, tandem repeat wiring, real mapq/base-quality/confidence all complete.

## Workspace Structure

```
phraya-core/     # Core types (Sequence, VariantObservation, EvidenceLayer, MinimizerSketch),
                 # errors, k-mer sketching (via simd-minimizers), centroid selection, uniqueness
phraya-align/    # WFA extension, SIMD diagonal fill (SSE4.2/NEON), seeding (Seed, find_seeds)
phraya-io/       # I/O: FASTA/FASTQ/BAM/CRAM parsing, .phraya/.phrayaplan/.phraya.queries formats
phraya-filter/   # Post-alignment filtering library (threshold/expression/preset)
phraya-cli/      # Binary CLI: plan/plan-tasks/align/merge/filter subcommands
```

**No `phraya-index` crate.** K-mer sketching was moved from a custom implementation into `phraya-core` (backed by `simd-minimizers`). The crate was deleted in the simd-minimizers refactor.

## Key Design Points

- **Zero binary deps**: All alignment in Rust using simd-minimizers + WFA. No runtime dep on BWA/minimap2/MUMmer.
- **Multi-mapping storage**: Query index tracks alternative alignment positions (score ratio ≥ 0.95 normalized edit distance). Enables filtering on mapping uniqueness post-alignment.
- **Deferred filtering**: Alignment produces rich `.phraya` files (multi-mapping, CIGAR, coverage tracks, k-mer uniqueness). Filter parameters applied post-hoc, not during alignment.
- **Library-first**: `phraya-filter` crate exposes filtering API. CLI is thin wrapper. Enables Python bindings, R integration, custom pipelines.
- **Evidence-informed alignment**: `.phrayaplan` files contain k-mer landscape + variation hotspots estimated from input sequences before alignment begins.
- **Platform-native SIMD**: simd-minimizers handles AVX2/NEON dispatch for sketching. The WFA/Myers inner loop (match extension) is SIMD-accelerated via `count_matching_prefix` (SSE2 on x86_64 / NEON on aarch64, both arch baselines — selected at compile time, no runtime dispatch; portable u64-XOR fallback elsewhere). ~5–9× faster than scalar on long match runs.
- **Sketch reuse**: `phraya plan` computes `MinimizerSketch` per sequence and stores them in `.phrayaplan` (v2) keyed by sequence ID. `phraya align` reuses them instead of recomputing; falls back to recomputing if sketch not in plan.

## Pipeline

```
phraya plan   → .phrayaplan v2 (k-mer evidence + task list + sketches)
              ↓
phraya align  → .phraya (position index) + .phraya.queries (query index with multi-mapping)
              ↓
phraya filter → VCF | TSV | filtered .phraya
```

### Use Cases (detected by `phraya plan`)

1. **Case 2**: N reads + reference → N alignments (BWA-like, main use case)
2. **Case 3**: M contigs + N reads, no ref → centroid selection + M+N-1 alignments (key innovation)
3. **Case 4**: M contigs ± reference → M or M×(M-1)/2 alignments (minimap2-like)

**Not supported**: reads-only with no reference and no contigs. All-vs-all short-read pairwise is de novo assembly (PRD §Won't Have). Phraya requires a coordinate space — provide a reference or mix in contigs and Case 3 auto-selects a centroid.

### Alignment Algorithm

- **Seeding**: `sketch_sequence_default()` (phraya-core, via simd-minimizers, k=21 w=11) → `find_seeds()` (phraya-align/seeding.rs)
- **Strand**: `align_read` aligns *both* the read's forward bytes and its reverse complement (`phraya_core::types::reverse_complement`) and keeps the better-scoring orientation (ties → forward). Canonical minimizers seed either strand, but extension is strand-naïve, so a reverse-strand read only places when its RC is tried. Extension always runs against the forward target, so CIGAR and alleles come out reference-forward regardless of strand; the chosen orientation is recorded on each `VariantObservation` (`strand()` → `Strand::Forward`/`Reverse`). The reverse orientation re-sketches the RC bytes (a forward sketch's seed query positions are wrong for RC).
- **Extension**: WFA O(s·n) (`wfa_extend`, fitting mode) or Myers bit-parallel O(nm/w) (`myers_extend`, fitting mode) — both share the SIMD `count_matching_prefix` match-extension primitive and emit the identical CIGAR convention (M/X consume both, `I` = target-only, `D` = query-only). Differential-tested to agree on edit distance.
- **Scoring**: Multi-mapping score ratio = (1 - edit_dist/query_len) for primary vs alternatives
- **Output**: Store all alignments with score_ratio ≥ 0.95 (hard-coded opinion)

#### Strategy ladder (`--strategy`, `executor::Strategy`)

Each preset selects an algorithm **and** a default coverage-window radius; `--coverage-window N` (or `AlignConfig::with_coverage_window_radius`) overrides the radius orthogonally.

| Strategy | K | Algorithm | Anchors | Coverage radius | Role |
|----------|---|-----------|---------|-----------------|------|
| `sensitive` | ∞ | seeded WFA, all anchors | all distinct seed target-starts + (0,0) | ±25bp | canonical reference path |
| `balanced` (default) | 2 | Myers fitting ≤500bp, WFA fallback | top 2 chain anchors | ±50bp | recall-balanced, faster engine |
| `fast` | 1 | Myers/WFA + divergence cutoff | single best-voted anchor | ±150bp | low-recall survey |

- The anchor cap K controls recall and ambiguity preservation: `sensitive` (K=∞) reports all seeds, `balanced` (K=2) reports top 2 by chain score, `fast` (K=1) reports only the best. Myers and WFA give identical edit distances (proven by differential sweep), so `balanced` is safe as the default; `sensitive` remains the trusted reference.
- `fast` trades recall for speed: K=1 anchor capping collapses per-read anchor count to O(1) in repeats (under-reports multi-mapping), and reads whose best alignment exceeds `FAST_MAX_DIVERGENCE` (0.20, hard-coded opinion) are dropped.
- Myers/WFA both run in **fitting** mode (query fully consumed, target end free) above the `n ≤ m + m/2 + 10` threshold; below it both fall back to global, so they stay consistent.
- **Score-bounded branch-and-bound (ADR-0007, `sensitive` only, `align_read`)**: `sensitive` extends its top-voted anchor first to seed an incumbent `d_best`, then extends each remaining alternate with `wfa_extend_capped` at `max_s = score_bound_max_s(L, d_best) = floor(0.05·L + 0.95·d_best)` — the exact edit distance at which the 0.95 multi-mapping score ratio is hit. Alternates whose wavefront exceeds the cap abandon early instead of running to completion; the incumbent (and cap) tightens as better alignments appear. This is lossless because WFA is boundable and the cap is always ≥ `d_best`, so only alternates `score_alignments` would have filtered are skipped. **Not applied to `balanced`/`fast`**: Myers is fixed-cost (`O(n·⌈m/64⌉)`, not boundable), and `balanced` is already anchor-bounded by K=2, so capping its alternates onto WFA would be *slower* and would break variant-level agreement on degenerate reads (Myers and WFA can pick different co-optimal CIGARs). Note: `.phraya`/`.phraya.queries` output is not byte-deterministic run-to-run (HashMap iteration order leaks into serialization); the regression gate is variant-level equivalence on well-aligned reads, not byte-identity.

## File Formats

### `.phrayaplan` v2 (binary MessagePack + zstd)
- **Read-only** during alignment, transmitted to all workers
- Contains: metadata (use case, input files), `kmer_index: HashMap<String, MinimizerSketch>` (sketches keyed by sequence ID), k-mer uniqueness scores, task list
- Version field checked on read; v1 files are rejected (incompatible kmer_index type)
- CLI tool: `phraya plan-tasks` dumps task list for GNU Parallel/xargs/WDL/Nextflow

### `.phraya` (position index, binary MessagePack + zstd)
- **VariantObservation** fields: position, ref_base, all_alleles (with counts), confidence, CIGAR, mapq, edit_distance, local_coverage (±50bp window), avg_base_quality, provenance
- **Coverage track**: quantized to nearest 5, RLE-compressed, full reference length
- **Mergeable**: combine multiple samples via position-centric merge (order-independent)

### `.phraya.queries` (query index, sidecar)
- For each query: list of alignment positions + scores above threshold (score_ratio ≥ 0.95)
- A query key appears **iff** it has ≥1 placement at score ≥ 0.95; `write_queries` drops reads whose filtered list is empty (so key count == placed-read count — do not count keys of reads with only sub-threshold hits as "aligned")
- Enables multi-mapping analysis: "exclude variants where >50% supporting reads multi-map"
- Separate file to keep merge fast (not needed for typical variant calling)
- Stored score is **absolute normalized identity** `(1 - edit/len)` for every placement (primary and alternates alike; `executor.rs` computes both this way) — not a primary-relative ratio. Absolute identity is invocation-independent, so it is directly comparable across reads and across reference spaces.

### `cross_space.phraya.queries` (cross-space sidecar, ADR-0011 / reference-palette mode)
- Written only by `phraya align --reference` mode. Type `CrossSpaceQueryIndex = query -> [(space, pos, identity)]` (`phraya-io/src/queries.rs`), spanning the whole presented reference palette rather than one coordinate space.
- Each placement carries absolute normalized identity `(1 - edit/len)`, so a cross-space margin (e.g. host 0.96 vs target 0.98) is recoverable straight off the sidecar. A read placed in N spaces appears once per space with directly comparable identities.
- Enacts **no** decision: no `--best`, no margin cutoff, no per-read classification — that is deferred to filter. Same 0.95 identity threshold and deterministic (BTreeMap + sorted) serialization as `.phraya.queries`.

### Reference-palette alignment mode (`phraya align --reference`, ADR-0011)
- New mode, mutually exclusive with traditional (`QUERY_ID`/`TARGET_ID`), `--worker`, and `--ensure`. `--output` is a **directory**.
- Aligns the plan's read pool against each repeatable `--reference` space independently → one `<output>/<label>.phraya` per space (label = palette name or content-hash prefix) plus one `cross_space.phraya.queries`.
- Each reference resolves by content hash against the plan's palette: **hit** reuses the planned sketch (`TargetContext::build_with_sketch`, no recompute); **miss** warns + sketches on the fly (tolerant default) or hard-errors under `--sealed` ("sealed" = fail-fast on unplanned reference; distinct from filter "strict").
- Read pool excludes every sequence whose hash is in the plan's palette (not just the presented subset), which makes the pool invocation-independent → composable: `align({A,B}) = align({A}) ∪ align({B})` (per-space `.phraya` byte-identical alone vs combined).

## K-mer Sketching

Sketching is implemented in `phraya-core/src/types.rs` using the `simd-minimizers` crate (Groot Koerkamp & Martayan, SEA 2025). Key types and functions:

```rust
// Type
pub struct MinimizerSketch { pub minimizers: Vec<(u64, u32)>, pub k: usize, pub w: usize }

// Functions (all in phraya_core::types)
pub fn sketch(sequence: &[u8], k: usize, w: usize) -> MinimizerSketch
pub fn sketch_sequence(seq: &Sequence, k: usize, w: usize) -> MinimizerSketch
pub fn sketch_sequence_default(seq: &Sequence) -> MinimizerSketch  // k=21, w=11
pub fn select_centroid(sketches: &[MinimizerSketch]) -> Option<usize>
pub fn compute_kmer_uniqueness(sketches: &[MinimizerSketch]) -> HashMap<u32, f64>

// Seeding (phraya_align::seeding)
pub struct Seed { pub query_pos: u32, pub target_pos: u32, pub minimizer: u64 }
pub fn find_seeds(q: &MinimizerSketch, t: &MinimizerSketch) -> Vec<Seed>
```

Parameters k=21, w=11 satisfy the simd-minimizers canonicality requirement (l = w+k-1 = 31, which is odd). These are reasonable defaults for bacterial genomics.

## Evidence Layer

### Plan-time (k-mer estimated)
- K-mer uniqueness: which k-mers appear in 1 vs many sequences
- Jaccard similarity: pairwise contig similarity from shared minimizers (used for centroid selection)
- Coverage estimate: k-mer depth across reference positions

### Alignment-time (ground truth)
- Per-position: multi_map_fraction, avg_score_ratio_gap to next-best alignment
- Polymorphic sites: all alleles with counts, reference base
- Invariant positions: where all samples match reference

## Filter Operations

### Supported Styles (all chainable)
1. **Threshold-based**: `--min-coverage 10 --min-mapq 30 --max-multi-map-fraction 0.3`
2. **Expression-based**: `--expr "coverage >= 10 && mapq > 30"` [Phase 2+]
3. **Named presets**: `--preset strict|tolerant`

### Output Formats
- Filtered `.phraya` (subset of records, same format for chaining)
- VCF (standard variant calling output)
- TSV/CSV (arbitrary column selection)

### Feature Space (available for filtering)
- Alignment: coverage, mapq, CIGAR complexity, edit_distance, multi_map_fraction, score_ratio_gap
- Context: edge_distance, local_gc, k-mer_uniqueness, in_homopolymer, in_tandem_repeat, snp_density (15bp/125bp/1000bp windows)
- Alleles: allele_frequency, ref_base, alt_bases
- Quality: avg_base_quality, confidence

## Phase 1 MVP (Shipped 2026-06-06)

**Complete:**
- Cases 2 (reads + ref), 3 (contigs + reads, auto-centroid), 4 (contigs only) working end-to-end
- Automatic centroid selection (Case 3): providing `--reference` overrides; omitting it triggers centroid selection. No separate flag needed.
- BAM/CRAM input via `noodles` (pure Rust, no htslib)
- `.phrayaplan` v2 with sketch reuse
- Real WFA O(s·n) alignment (wavefront-based, not diagonal DP)
- Real local coverage (±50bp window from alignment, not stubbed)
- mapq, avg_base_quality, confidence derived from input data (BAM records / alignment score)
- Tandem repeat detection wired end-to-end: annotation on variants, `exclude_tandem_repeats` filter option
- `phraya filter` threshold-based + named presets (strict / tolerant), outputs VCF/TSV/phraya
- `phraya-filter` crate library API with feature extractors (cigar_ops, allele_frequency, multi_map_fraction)

**Deferred to Phase 2+:**
- Expression-based filters (`--expr`)
- Variation hotspot estimation in plan
- Two-tier evidence (k-mer → alignment refinement)
- Python bindings / R integration
- Advanced coverage: per-base breakdown (not just ±50bp window)
- CRAM parsing with external reference (currently only reference-free/unmapped CRAM records)

## Dependencies

- **simd-minimizers** (`simd-minimizers` crate, Groot Koerkamp & Martayan, SEA 2025): canonical minimizer sketching via AVX2/NEON. Compile with `-C target-cpu=native` for full SIMD.
- **noodles** (`noodles-bam`, `noodles-cram`, `noodles-sam`): pure-Rust BAM/CRAM I/O
- **MessagePack** (`rmp-serde`): binary serialization for `.phraya`/`.phrayaplan` formats
- **Compression** (`zstd`): for plan/alignment file storage
- **Parallelism** (`rayon`): embarrassingly parallel tasks within single-node execution
  - `phraya align --ensure` processes missing chunks in parallel
  - Plan data shared via Arc (read-only after load)
  - `--threads N` overrides auto-detection
- Standard: `thiserror` (errors), `serde` (serialization), `clap` (CLI)

## Coding Conventions

- Stable Rust (not nightly for library crates)
- No unsafe unless necessary (none currently in WFA code)
- Tests alongside code (`#[cfg(test)] mod tests`)
- Benchmarks in `benches/` (criterion)
- Library-first: `phraya-filter` crate exposes public API, CLI is thin wrapper

## Implementation Notes

- **Score ratio hard-coded**: 0.95 threshold for multi-mapping is Phraya's opinion, not user-configurable
- **Centroid selection**: For case 3 (contigs + reads, no ref), select contig closest to k-mer space center (median Jaccard similarity) as reference coordinate space
- **MinimizerSketch positions**: `u32` (not `usize`) — matches simd-minimizers output directly
- **PHRAYAPLAN_VERSION = 2**: v1 files (Vec<MinimimizerSketch>) rejected; plan files are ephemeral (always regenerate with `phraya plan`)
- **`MinimimizerSketch` typo**: The old crate used `MinimimizerSketch` (extra 'i'). The current type is `MinimizerSketch`. Do not reintroduce the old name.
