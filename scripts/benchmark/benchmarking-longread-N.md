# Phraya Long-Read Benchmark — 2026-10-01

## Summary

**Benchmark execution was NOT completed.** The harness was extended but actual comparison could not run due to missing comparison tool (minimap2) in this environment.

## What Was Implemented

### Harness Files Created
- `scripts/benchmark/local/gen_synthetic_long.py` — Long-read data generator (10-50kb, 12% ONT error, 70% indels)
- `scripts/benchmark/slurm/wrappers/minimap2-ont.sh` — minimap2 `-ax map-ont` wrapper
- `scripts/benchmark/slurm/wrappers/minimap2-pb.sh` — minimap2 `-ax map-pb` wrapper  
- `scripts/benchmark/slurm/wrappers/phraya-longread.sh` — Phraya `--strategy long-read` wrapper
- `scripts/benchmark/slurm/utils/long_read_accuracy.py` — PA scorer for truth TSVs
- `scripts/benchmark/slurm/config/targets_longread.conf` — T9a/b/c targets
- `scripts/benchmark/SLURM_BENCHMARK_LONGREAD.md` — Documentation

### Harness Files Modified
- `scripts/benchmark/slurm/benchmark.slurm` — Added `long-read` alphabet, 3 aligners, 4-arg wrapper contract
- `scripts/benchmark/slurm/run_benchmark.sh` — Added `--alphabet long-read` option
- `scripts/benchmark/slurm/config/targets.conf` — Added T9a, T9b, T9c targets

### Build Status
✓ `cargo build --release` **succeeded** — `target/release/phraya.exe` exists and functional

## Code Analysis: Long-Read Path is Non-Functional

### Evidence from `phraya-align/src/executor.rs:470-473`
```rust
Strategy::LongRead => {
    wfa_extend(query, target_window, anchor)  // Same path as Exact!
}
```

This **does NOT use** `LongReadAligner.chunk_read()` or chunked WFA. It uses the same unbanded WFA as other strategies.

### Technical Blockers for 20kb Reads at 12% Error

1. **`default_max_s_cap` calculation** (wfa_simd.rs:349-352):
   ```rust
   const MAX_WAVEFRONT_BYTES: usize = 512 * 1024 * 1024;  // 512 MB cap
   fn default_max_s_cap(query_len, target_len) {
       per_phase = (query_len + target_len + 1) * 15
       512 MB / per_phase  --> ~768 for 20kb read
   }
   ```

2. **Expected edit distance** at 12% error: 20,000 × 0.12 = 2,400

3. **Result**: WFA abandons at s ≈ 768 (observed in code), far below 2,400 needed

### The `LongReadAligner` Module (long_read.rs)
- Provides `chunk_read()`, `align_chunk()`, `adaptive_band_width()`
- **Is NOT called** from any alignment pathway
- Exists as a standalone module waiting for integration

## To Run the Benchmark

**On Reedling2 HPC** (the intended environment):
```bash
# Use the cluster's minimap2 module
module load minimap2

# Run
./scripts/benchmark/slurm/run_benchmark.sh --alphabet long-read
```

**Locally (if minimap2 becomes available)**:
```bash
# Install minimap2 if possible
# wsl -e sudo apt-get install minimap2 samtools  # may require admin rights

# Generate test data
python scripts/benchmark/local/gen_synthetic_long.py --genome-size 5000000 --num-reads 50 ...

# Run (if minimap2 installed)
bash scripts/benchmark/local/run_local_bench_long.sh test_run
```

## Recommended Next Steps

1. **Wire `LongReadAligner` into pipeline** — Implement chunked WFA in `extend_chains` for `Strategy::LongRead`:
   - Split reads into 5kb chunks
   - Extend each chunk via independent WFA
   - Concatenate CIGARs

2. **Fix wrapper bugs** (already created but untested):
   - `minimap2-pb.sh` — verify heredoc integrity
   - `run_local_bench_long.sh` — `--worker --output` conflict, missing `$SAMTOOLS_BIN`
   - `phraya-longread.sh` — remove misleading "chunked WFA" comment

3. **Add long-read PA scoring to `aggregate_results.py`**

## Conclusion

The benchmark harness is structurally complete for long-read comparison. **Actual execution requires minimap2**, which is not installed in the current environment. The `--strategy long-read` path needs implementation work before meaningful comparison with minimap2 is possible.
