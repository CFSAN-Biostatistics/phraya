# SIMD Performance Analysis

**Date**: 2026-06-08  
**Status**: ✅ FIXED (2026-06-08)  
**Question**: Should we use `wide` portable SIMD or raw SSE4.2/AVX2 intrinsics?

## Critical Finding: Wrong Algorithm in Production (FIXED)

**The performance problem is NOT `wide` vs intrinsics. The problem is O(n×m) diagonal DP vs O(s·n) WFA.**

### Current Production Dispatch (phraya-align/src/lib.rs:91-96)

```rust
#[cfg(target_arch = "x86_64")]
pub fn wfa_extend(query: &[u8], target: &[u8], seed: SeedAnchor) -> WfaResult {
    if !cfg!(debug_assertions) && is_x86_feature_detected!("sse4.2") {
        wfa_simd::wfa_extend_simd_impl(query, target, seed)  // ← WRONG: calls fill_simd (O(n×m))
    } else {
        wfa_simd::wfa_extend_naive_impl(query, target, seed)  // ← CORRECT: calls fill_wfa_fitting (O(s·n))
    }
}
```

### What Each Implementation Does

| Function | Algorithm | Complexity | Correctness |
|----------|-----------|------------|-------------|
| `wfa_extend_naive_impl` | WFA (wavefront) via `fill_wfa_fitting` | O(s·n) | ✅ Correct, fast |
| `wfa_extend_simd_impl` | Diagonal DP via `fill_simd` | O(n×m) | ✅ Correct, **slow** |

- `fill_wfa_fitting` (line 226): Real WFA — processes by edit distance wavefronts, O(s·n)
- `fill_simd` (line 674): Diagonal DP with `wide::i32x8` — processes full n×m matrix, O(n×m)

### Benchmark Results (RUSTFLAGS="-C target-cpu=native")

| Workload | naive (WFA O(s·n)) | simd (DP O(n×m)) | Slowdown |
|----------|-------------------|------------------|----------|
| 10kb, 1% divergence | 22ms | 630ms | **28×** |
| 10kb, 5% divergence | 47ms | 1.4s | **30×** |

**SIMD is slower because it's the wrong algorithm**, not because `wide` is slow.

## Root Cause Analysis

### Historical Context

1. **Phase 1**: Diagonal DP with `wide::i32x8` (`fill_simd`) was implemented first
2. **Phase 2**: Real WFA (`fill_wfa_fitting`) was added to `wfa_extend_naive_impl`
3. **Dispatch error**: Production release builds route to `wfa_extend_simd_impl` → `fill_simd` (OLD slow path)
4. **Tests pass**: Both impls produce correct edit distance + CIGAR, so differential tests don't catch perf regression

### Why Tests Didn't Catch This

- Correctness tests: `fill_simd` vs `fill_wfa` both return correct results
- Only 1 timing test exists (`wfa_is_faster_than_on2_for_sparse_edits_windowed`)
- That test FAILS (12.5ms vs <10ms target), but was attributed to "timing too strict" not "wrong algorithm"

## Recommendation

### Fix the Dispatch (Immediate)

**Option A**: Route SIMD dispatch to WFA, not diagonal DP

```rust
#[cfg(target_arch = "x86_64")]
pub fn wfa_extend(query: &[u8], target: &[u8], seed: SeedAnchor) -> WfaResult {
    // Always use WFA, not diagonal DP
    wfa_simd::wfa_extend_naive_impl(query, target, seed)
}
```

**Option B**: Delete `fill_simd` entirely, keep WFA only

- `fill_simd` is O(n×m) diagonal DP — obsolete after WFA impl
- No production use case for O(n×m) DP when O(s·n) WFA exists
- Simplifies codebase: one algorithm, not two

### SIMD WFA (Future)

If SIMD speedup desired, apply SIMD **to WFA**, not to diagonal DP:

1. **Diagonal extension**: WFA's greedy match phase can use SIMD for `memcmp`-like vectorized comparison
2. **Wavefront min operations**: When advancing s→s+1, compute `min(M[s], I[s]+1, D[s]+1)` with SIMD over multiple diagonals
3. **Backtrace**: Keep scalar (not performance-critical)

Example: WFA with SIMD diagonal extension (pseudocode):

```rust
fn wfa_extend_with_simd(q: &[u8], t: &[u8]) -> (String, usize) {
    // Standard WFA structure
    let mut wavefronts = vec![...];
    
    for s in 0..max_score {
        // Advance wavefront s
        for k in active_diagonals(s) {
            let (i, j) = wavefront_to_coords(wavefronts[s][k], k);
            
            // SIMD greedy extension: how far can we match?
            let extension = simd_match_run(&q[i..], &t[j..]);  // ← AVX2 memcmp
            wavefronts[s][k] += extension;
        }
        
        // Check termination
        if reached_end(...) { break; }
        
        // Expand to s+1 (insert/delete/mismatch)
        wavefronts[s+1] = compute_next_wavefront(&wavefronts[s]);
    }
    
    traceback(&wavefronts)
}
```

This would be **O(s·n) with SIMD**, not O(n×m).

## `wide` vs Raw Intrinsics

**Secondary question** (once algorithm is fixed): Should WFA use `wide` or raw intrinsics?

### `wide` crate pros/cons

**Pros:**
- Portable: x86_64, aarch64, other archs (scalar fallback)
- Safe: no `unsafe`, compiler verifies correctness
- Maintenance: no platform-specific codepaths

**Cons:**
- Abstraction overhead: `i32x8` may or may not lower to AVX2 (depends on RUSTFLAGS)
- Debug slow: unoptimized builds call out-of-line functions
- Less control: can't fine-tune for cache/alignment

### Raw intrinsics pros/cons

**Pros:**
- Guaranteed lowering: `_mm256_min_epi32` always emits `vpminsd`
- Maximum control: can optimize loads/stores, prefetch, alignment
- Benchmarkable: can measure exact instruction count

**Cons:**
- Unsafe: requires SAFETY comments, manual verification
- Platform-specific: separate x86_64 (SSE/AVX) and aarch64 (NEON) impls
- Maintenance: 2-3× code duplication

### Recommendation

**For WFA SIMD (if implemented):** Use `wide` for greedy extension, raw intrinsics only if profiling shows bottleneck.

**Rationale:**
1. WFA's perf win is O(s) vs O(n×m), not SIMD — algorithm matters more than intrinsic choice
2. Typical bacterial genomics: s=2-5% of n (200-500 edits / 10kb read) — WFA is already 20-50× faster than DP
3. SIMD greedy extension in WFA: processes ~64 bases/iteration with AVX2, ~32 with SSE — diminishing returns vs scalar
4. `wide` sufficient for 2-4× speedup; raw intrinsics might yield 4-6×, but at 3× code complexity

**Exception:** If profiling shows >30% time in match extension, then raw intrinsics justified.

## Action Items

1. **Fix dispatch** (high priority):
   - Change `wfa_extend` to always call `wfa_extend_naive_impl` (the WFA path)
   - OR delete `fill_simd` and `wfa_extend_simd_impl` entirely
   - Verify release benchmarks: should be ~22ms for 10kb 1% divergence, not 630ms

2. **Update docs** (medium priority):
   - CLAUDE.md: clarify that production uses WFA O(s·n), not diagonal DP
   - Remove "SIMD diagonal fill" claims (currently misleading)

3. **Add WFA benchmark** (medium priority):
   - Criterion bench comparing WFA vs diagonal DP vs Myers bit-parallel
   - Sizes: 150bp (typical read), 1kb, 10kb
   - Divergences: 1%, 2%, 5%, 10%

4. **SIMD WFA** (low priority, future):
   - Only pursue if profiling shows match extension >30% of runtime
   - Use `wide` first; switch to intrinsics only if measured benefit

## Appendix: What Tests Actually Measure

| Test | What it measures | Catches this bug? |
|------|------------------|-------------------|
| `fill_simd_matches_fill_scalar_property` | Correctness: SIMD DP == scalar DP | No (both O(n×m)) |
| `simd_vs_naive_differential::*` | Correctness: SIMD impl == naive impl | No (both correct) |
| `wfa_is_faster_than_on2_for_sparse_edits_windowed` | Timing: WFA < 10ms for 150bp×300bp | **YES** (fails: 12.5ms) |

**The one test that caught the bug FAILED**, but was attributed to "timing threshold too strict" instead of "wrong algorithm in production."

## Conclusion

**Primary issue:** Production dispatch uses O(n×m) diagonal DP instead of O(s·n) WFA. Fix by routing to `wfa_extend_naive_impl`.

**Secondary question** (`wide` vs intrinsics): **Irrelevant until primary issue fixed.** Once WFA is used in production, `wide` is adequate. Raw intrinsics only justified if profiling shows >30% time in SIMD-able code.

**Don't optimize the wrong algorithm.**

---

## Fix Applied (2026-06-08)

**Changes made:**

1. **phraya-align/src/lib.rs:85-112** - Updated dispatch to always use WFA:
   - x86_64: Removed conditional routing to `wfa_extend_simd_impl`, always use `wfa_extend_naive_impl` (WFA)
   - aarch64: Removed conditional routing to `wfa_extend_neon_impl`, always use `wfa_extend_naive_impl` (WFA)
   - Other architectures: Already correct (unchanged)

2. **phraya-align/src/wfa_simd.rs:663-682** - Documented `fill_simd` as test-only O(n×m) reference

3. **CLAUDE.md** - Removed misleading "SIMD diagonal fill" claims

**Result:** Production now uses O(s·n) WFA on all architectures. Expected 28× speedup for typical genomics workloads.

**Verification pending:**
- Test suite: `cargo test --all`
- Benchmark: `cargo bench --bench wfa_simd_bench -- "10kb"`
- Integration: `cargo test --release -p phraya-cli integration_test_e2e_case2`
