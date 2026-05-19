# Phraya PRD
## Probabilistic Heuristic for Running Awesome Yield-Agnostic Alignments

**Version**: 0.3  
**Status**: Draft  
**Last Updated**: 2026-05-19

---

## Problem

SNP calling for bacterial genomics is fragmented by input type:

- **Read-based tools** (Snippy, CFSAN) use BWA + Freebayes. Rich evidence (depth, quality, strand bias) but slow and storage-heavy.
- **Assembly-based tools** (CSP2, Parsnp, MUMmer) align contigs. Fast and storage-efficient but all read-level uncertainty is lost.

This is a false dichotomy. Reads and assemblies represent the same biology at different certainty levels: reads are a probability distribution, assemblies are point estimates. The field treats them as categorically different when they're just different pre-processing stages of the same evidence.

**Consequences**:
- Mixed datasets (some reads, some assemblies, some both) have no single principled pipeline
- Tool-switching required for comparable analyses across labs
- Confidence reporting is inconsistent or absent
- **Deployment is a nightmare**: Installing BWA + MUMmer + samtools + their dependencies across platforms is fragile, especially in air-gapped/regulated/container environments
- **Performance is locked to decade-old tools**: BWA-MEM (2013), MUMmer 4 (2018) — neither optimized for 2026 hardware (AVX-512, ARM, GPU)

---

## Solution

**Phraya**: A unified SNP caller that accepts reads, assemblies, or any mixture, implements alignment algorithms natively in Rust, and produces variant calls with structured, source-aware confidence scores.

**Core principles**:

1. **Native implementations, not wrappers**: Every alignment algorithm is in Rust in the Phraya codebase. Zero runtime dependencies on BWA, MUMmer, Bowtie, samtools, or any external aligner. Single binary deployment.

2. **Speed-sensitivity continuum**: Multiple alignment strategies (`exact`, `balanced`, `fast`, `sketch`) in the same binary. Users choose accuracy/speed tradeoff with a flag, not by switching tools.

3. **Platform-native performance**: First-class optimization for 2026 hardware: high-end x86 (AVX2, AVX-512), Threadripper, ARM64/Apple Silicon, GPU acceleration. Same binary, runtime dispatch to optimal code path.

4. **Unified evidence model**: Structured `VariantEvidence` objects capture what kind of evidence supports each call (depth + quality, alignment identity, contig position, etc.). Full transparency for downstream filtering.

---

## Goals

### Must Have (v1.0)

1. Accept Illumina reads (PE/SE), assemblies (draft/complete), or both per sample
2. Zero binary dependencies — installed Phraya binary is sufficient for all alignment
3. Multiple alignment strategies on speed/sensitivity continuum
4. Structured evidence output (VCF + native format) with full transparency for filtering
5. Performance ≥ par with BWA-MEM2 (reads) and MUMmer 4 (assemblies)
6. Native ARM64 and Apple Silicon support with competitive performance
7. Correct for bacterial surveillance: haploid genomes, <1000 SNP divergence

### Should Have (v1.0 or v1.1)

8. Long read support (ONT/PacBio)
9. Hybrid mode: assembly scaffold + targeted read deployment at uncertain sites
10. GPU acceleration (CUDA) for large-scale pairwise workloads
11. Automatic reference selection via k-mer sketching
12. Caching of pairwise results for incremental cluster updates

### Won't Have (v1.0)

- Diploid/polyploid calling
- Structural variants beyond indels
- De novo assembly
- Functional annotation
- GUI

---

## Architecture Sketch

### Aligner Trait

```rust
pub trait Aligner: Send + Sync {
    fn align(
        &self,
        query: &SequenceSource,  // Reads | Assembly | Hybrid
        reference: &IndexedReference,
        params: &AlignParams,
    ) -> Result<AlignmentResult, AlignError>;
}
```

### Concrete Aligners

| Aligner | Input | Algorithm Family | Key Techniques |
|---|---|---|---|
| `ContigAligner` | Assembly | Seed-chain-extend | FM-index OR minimizer seeding → WFA extension |
| `ShortReadAligner` | Illumina reads | Seed-extend-rescue | FMD-index seeding → banded Smith-Waterman → paired-end rescue |
| `LongReadAligner` | ONT/PacBio | Minimizer-chain-align | Minimizer seeding → WFA with adaptive banding |
| `HybridAligner` | Reads + assembly | Targeted deployment | Contig alignment first, then targeted read mapping at uncertain sites only |

### Strategy Mapping

| Strategy | Speed | Sensitivity | Seeding | Extension | Use Case |
|---|---|---|---|---|---|
| `exact` | Slowest | Highest | FM-index, no threshold | Full WFA, wide band | Outbreak analysis, regulatory |
| `balanced` | Moderate | High | FM-index, min seed length | WFA adaptive band | General purpose (default) |
| `fast` | Fast | Good | Minimizer, sparse chain | Narrow band SW | Cluster screening, QC |
| `sketch` | Fastest | Approx | K-mer sketch only | None | Distance-only, triage |

Users select `--strategy`; code routes to appropriate parameters/algorithms.

### Platform Optimization

**SIMD** (automatic, runtime dispatch on x86):
- x86-64: SSE4.2 baseline, AVX2, AVX-512 compiled into same binary → CPUID dispatch
- ARM64: NEON (mandatory on aarch64, no dispatch needed)
- Kernels: WFA diagonal extension, Smith-Waterman cell fill, minimizer hash comparison

**Multi-core** (automatic):
- Rayon work-stealing for pairwise batch dispatch
- Threadripper 64–96 cores utilized transparently

**GPU** (opt-in, `--gpu` flag, requires `--features cuda` build):
- Pairwise `ContigAligner` batch offload
- WFA inner loop as CUDA kernel
- Target: ≥10× vs 64-core Threadripper for N≥100 assemblies

### Evidence Object

```rust
pub struct VariantEvidence {
    pub reference_pos: u64,
    pub ref_base: Base,
    pub alt_base: Base,
    pub call: VariantCall,
    pub filter_reasons: Vec<FilterReason>,
    pub confidence: f64,  // 0.0–1.0, calibrated
    pub evidence_source: EvidenceSource,  // Asm | Reads | Hybrid
    pub asm_evidence: Option<AsmEvidence>,
    pub read_evidence: Option<ReadEvidence>,
}

pub struct AsmEvidence {
    pub alignment_identity: f64,
    pub distance_to_contig_edge: u64,
    pub snp_density_windows: [f64; 3],  // 15bp, 125bp, 1000bp
    pub in_repeat_region: bool,
    // ...
}

pub struct ReadEvidence {
    pub depth: u32,
    pub alt_fraction: f64,
    pub mean_base_quality: f64,
    pub strand_bias_pvalue: f64,
    // ...
}
```

Every variant has full structured evidence. Users filter via `phraya filter` with explicit thresholds or presets.

---

## Output Formats

1. **`.phraya` files** (native): MessagePack binary, one per pairwise comparison. Full `VariantEvidence` for all sites (called + filtered). Caching unit for incremental runs.
2. **VCF**: Standard VCF 4.2 with Phraya-specific INFO/FORMAT tags for evidence fields.
3. **Core SNP alignment**: FASTA, one sequence per sample. Snippy-compatible for IQ-TREE/Gubbins.
4. **Distance matrix**: TSV pairwise distances (raw + filtered).
5. **QC summary**: Per-comparison audit trail (coverage, strategy used, filter stages).

---

## CLI Sketch

```bash
phraya snp --manifest samples.tsv --reference ref.fasta --strategy balanced --outdir results

phraya filter results/*.phraya --min-confidence 0.95 --min-depth 10 --preset strict-surveillance

phraya screen queries/ --reference db.fasta --strategy fast --threads 64 --gpu
```

Key flags:
- `--strategy {exact|balanced|fast|sketch}`
- `--gpu` (enable CUDA, requires build with `--features cuda`)
- `--preset {strict-surveillance|research|assembly-only|hybrid}` (for filtering)
- `--auto-ref <N>` (auto-select reference by sketch centrality)

---

## Technology Stack

**Language**: Rust (stable, not nightly for library crates)

**Key crates**:
- `bio` (rust-bio): BWT, FM-index, FMD-index scaffolds
- `wfa`: Pure Rust WFA implementation
- `needletail`: FASTQ/FASTA streaming
- `noodles`: VCF/BAM/GenBank I/O
- `rayon`: Parallel dispatch
- `multiversion`: x86 SIMD runtime dispatch
- `safe_arch` or `core::arch`: SIMD intrinsics
- `cudarc`: CUDA FFI (optional)

**Build**:
- CPU-only: zero external dependencies, single static binary
- CUDA: requires CUDA toolkit at build time, runtime detection graceful degrades if no GPU

---

## Validation

**Accuracy**:
- Simulated data with known SNPs across divergence/depth ranges
- Concordance with Snippy (reads) and CSP2 (assemblies) on real datasets
- Consistency checks: reads-only vs assembly-only vs hybrid for same sample
- Real outbreak validation (CFSAN/CDC datasets with known epidemiological linkage)

**Performance**:
- Benchmark N=10/50/100/500 samples, reads/assemblies/mixed, all strategies
- Hardware: x86 baseline (SSE), x86 AVX2, x86 AVX-512, Threadripper, ARM64/Graviton, Apple Silicon M4, NVIDIA A100
- Baselines: BWA-MEM2, MUMmer 4, CSP2, Snippy
- Targets: reads ≤2× BWA-MEM2, assemblies ≤1.5× MUMmer, GPU ≥10× CPU for N≥100

---

## Roadmap

**Phase 1** (8 weeks): Core + assembly alignment (`ContigAligner`, SSE/NEON baseline, `balanced` strategy only, assembly-only manifest)

**Phase 2** (10 weeks): Reads + SIMD hardening (`ShortReadAligner`, AVX2/AVX-512 dispatch, `exact`/`fast` strategies, `phraya filter`)

**Phase 3** (12 weeks): Hybrid + GPU (`HybridAligner`, CUDA batch kernel, mixed manifest, `sketch` strategy, long read scaffold)

**Phase 4** (14 weeks): Long reads + polish (`LongReadAligner` complete, docs, benchmarks, publication, v1.0)

**Total to v1.0**: ~44 weeks

---

## Open Questions

1. **Confidence scoring**: How to merge assembly identity + read depth into calibrated probability? Proposed: logistic regression on simulated ground truth.
2. **Pure Rust WFA vs WFA2-lib FFI**: Trade build simplicity for maturity. Start pure Rust, add FFI opt-in if needed.
3. **Indels**: Call and report in VCF, exclude from distance by default? Or include in `--strategy exact`? Needs user feedback.
4. **Repeat detection**: Auto-flag via self-alignment, or user-provided BED only? Proposed: auto-flag, filter only if user requests.
5. **Index versioning**: How to detect stale indexes? Proposed: version + parameter hash in index header, auto-rebuild on mismatch.
6. **Metal (Apple GPU)**: Post-v1.0. NEON CPU path is fast enough; Metal adds complexity without clear demand.
7. **Windows native**: WSL2 only for v1.0, native post-v1.0 if demand exists.

---

## Why "Phraya"?

The Chao Phraya river in Thailand is formed by the confluence of two rivers into a single navigable channel — apt for a tool that merges read-based and assembly-based evidence streams.

"Yield-Agnostic": works with 100× read depth or zero reads, treating sequencing yield as evidence to weight, not a prerequisite.
