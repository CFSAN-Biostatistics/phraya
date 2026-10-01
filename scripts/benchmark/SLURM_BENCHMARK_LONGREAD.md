# Long-Read Alignment Benchmark

This document describes the long-read benchmarking extension for Phraya.

## Overview

Long reads (ONT/PacBio CLR, 10-50 kb) have fundamentally different error profiles than short reads: indels dominate, error rates are 10-15%, and traditional short-read aligners cannot handle them directly.

This benchmark compares:
- **Phraya `--strategy long-read`** — WFA-based seeded extension with automatic chunking support
- **minimap2 `-ax map-ont`** — ONT R10.4 preset (the de facto standard for ONT data)
- **minimap2 `-ax map-pb`** — PacBio CLR preset (for future comparison)

## Data Generation

`gen_synthetic_long.py` generates synthetic long reads with realistic error characteristics:

```python
# ONT R10.4-style errors (~12% total, 70% indels, 30% subs)
--read-error 0.12 --indel-frac 0.70
# Read length 10-30 kb
--min-len 10000 --max-len 30000
```

The truth sidecar records the pre-mutation fragment coordinates:
```
read_id  chrom  true_start  strand  true_end  n_subs  n_ins  n_del  indel_events
```

## Alignment Results

### Current Behavior: `WFA Memory Cap Abandonment`

**Critical finding**: `phraya align --strategy long-read` currently fails on realistic long-read error rates. The WFA is bounded by `default_max_s_cap`, which limits wavefront memory to 512 MB.

For a 20 kb read with 12% error rate:
- Expected edit distance: ~2,400
- `default_max_s_cap(20_000, 40_000) ≈ 768`
- WFA abandons at s ≈ 768, producing `AlignmentFailed` error

This is because `Strategy::LongRead` currently dispatches to the same `wfa_extend` path as other strategies, hitting the memory cap.

**Recommended fix**: Wire chunked WFA into `extend_chains` for `Strategy::LongRead`. Split reads into 5 kb chunks, extend each chunk via WFA independently (each chunk stays under the memory cap), and concatenate CIGARs.

## Running Locally

```bash
# Generate long-read test data
python3 scripts/benchmark/local/gen_synthetic_long.py \
    --genome-size 5000000 --num-reads 50 --min-len 8000 --max-len 20000 \
    --ref-out /tmp/long_ref.fa --reads-out /tmp/long_reads.fq

# Build Phraya with native CPU features
RUSTFLAGS="-C target-cpu=native" cargo build --release

# Run local benchmark
bash scripts/benchmark/local/run_local_bench_long.sh longread_test
```

## SLURM Benchmark Extension

### New Aligners
- `minimap2-ont.sh` — ONT R10.4 preset
- `minimap2-pb.sh` — PacBio CLR preset (add when needed)
- `phraya-longread.sh` — Phraya with `--strategy long-read`

### New Targets (T9a, T9b, T9c)
```
T9a  ecoli K12  | 4.6 Mb  | bacterial single-chromosome
T9b  ecoli  | 4.6 Mb  | same, for comparison
T9c  human chr21 | 48 Mb  | mammalian single chromosome
```

### New Accuracy Metric
`long_read_accuracy.py` computes placement accuracy from:
- phraya's `cross_space.phraya.queries` (msgpack+zstd)
- SAM/BAM output from minimap2

Uses the truth sidecar from `gen_synthetic_long.py`.

## Expected Results (Placeholder)

| Aligner | Wall (s) | RSS (GB) | PA | Notes |
|---------|----------|----------|----|-------|
| phraya-longread | TBD | TBD | TBD | Currently broken |
| minimap2-ont | TBD | TBD | TBD | baseline |

## Key Metrics

- **Wall time**: seconds per 20k reads
- **Peak RSS**: memory during alignment
- **PA (Placement Accuracy)**: fraction of reads placed within ±100 bp of true position
- **IEC (Indel Event Concordance)**: fraction of reads with correct indel count (future)