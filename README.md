# Phraya

**Probabilistic Heuristic for Running Awesome Yield-Agnostic Alignments**

General-purpose pairwise sequence aligner for bacterial genomics. Produces rich alignment superpositions with deferred filtering for SNP calling, in-silico typing, classification, and other downstream analyses. Zero binary dependencies. Native Rust implementation with SIMD optimization (AVX2/NEON).

## Status

**Phase 1 MVP in development.** Architecture revision completed 2026-05-27. See [issue #58](https://github.com/CFSAN-Biostatistics/phraya/issues/58) for PRD.

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

Workspace with 6 crates:

- **phraya-core**: Core types (Sequence, VariantObservation, EvidenceLayer, CoverageTrack), errors
- **phraya-index**: K-mer sketching (wraps simd-minimizers), uniqueness, centroid selection, seeding
- **phraya-io**: FASTA/FASTQ parsing, `.phrayaplan`/`.phraya`/`.phraya.queries` formats (MessagePack + zstd)
- **phraya-align**: WFA extension, SIMD diagonal fill (SSE4.2/NEON)
- **phraya-filter**: Filtering library (threshold/expression/preset), output formatters (VCF/TSV/phraya)
- **phraya-cli**: Binary CLI (plan/plan-tasks/align/filter subcommands)

## File Formats

- **`.phrayaplan`**: Plan file (k-mer evidence + task list). Read-only during alignment. Binary MessagePack + zstd.
- **`.phraya`**: Position index (variant observations + coverage track). Mergeable. Binary MessagePack + zstd.
- **`.phraya.queries`**: Query index (multi-mapping alternatives per read). Sidecar file. Binary MessagePack + zstd.

## Building

```bash
cargo build --release
```

Requires Rust 1.75+. No external dependencies (BWA, minimap2, samtools).

## Design Decisions

- **Score ratio threshold**: 0.95 (hard-coded). Stores alternatives within 95% of best identity. Opinionated choice for storage efficiency.
- **K-mer parameters**: k=21, w=11 (defaults from simd-minimizers). Reasonable for bacterial genomes.
- **Coverage quantization**: Nearest 5. Enables RLE compression, negligible precision loss.
- **SIMD library**: Uses `simd-minimizers` crate (AVX2/NEON dispatch). State-of-the-art k-mer sketching.

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

## License

Unlicense. As a work product of the US Government (17 USC 105), Phraya is in the public domain.

## Contributing

See [issue #58](https://github.com/CFSAN-Biostatistics/phraya/issues/58) for Phase 1 PRD. Implementation contributions welcome.
