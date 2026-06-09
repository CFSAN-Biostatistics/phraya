# WFA SIMD Opportunities - Should We Use Platform Intrinsics?

**Date**: 2026-06-08  
**Context**: After fixing dispatch to use O(s·n) WFA instead of O(n×m) DP, would adding SIMD to WFA provide further speedup?

## Answer: Maybe - But Profile First

The current WFA implementation (phraya-align/src/wfa_simd.rs:226-322) has two SIMD opportunities:

### 1. Greedy Match Extension (`extend` closure, line 247-256)

**Current implementation:**
```rust
let extend = |q_pos: i32, k: i32| -> i32 {
    let mut i = q_pos;
    loop {
        let j = i - k;
        if i >= qn || j >= tn || j < 0 { break; }
        if q[i as usize] != t[j as usize] { break; }  // ← byte-by-byte comparison
        i += 1;
    }
    i
};
```

**SIMD opportunity:** Vectorized memcmp
- **AVX2**: Compare 32 bytes per instruction (vs 1 byte scalar)
- **SSE4.2**: Compare 16 bytes per instruction
- **NEON**: Compare 16 bytes per instruction

**Expected speedup:**
- Best case (long matches): 8-16× faster extension
- Typical bacterial genomics (95-99% identity): 2-4× faster overall if extension is >50% of runtime
- Worst case (many mismatches): ~1× (SIMD overhead dominates)

**Implementation complexity:**
- **Moderate**: ~50-80 lines of `unsafe` code per platform (AVX2, SSE4.2, NEON)
- Runtime dispatch via `is_x86_feature_detected!("avx2")`
- Must handle unaligned loads, partial chunks, edge cases

**Example (AVX2):**
```rust
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn extend_avx2(q: &[u8], t: &[u8], q_pos: usize, k: i32, qn: usize, tn: usize) -> usize {
    use std::arch::x86_64::*;
    let mut i = q_pos;
    
    loop {
        let j = i as i32 - k;
        if i >= qn || j < 0 || j as usize >= tn { break; }
        
        let j_usize = j as usize;
        let remaining = (qn - i).min(tn - j_usize);
        
        if remaining >= 32 {
            // SAFETY: bounds checked above, unaligned load OK
            let q_chunk = _mm256_loadu_si256(q.as_ptr().add(i) as *const __m256i);
            let t_chunk = _mm256_loadu_si256(t.as_ptr().add(j_usize) as *const __m256i);
            let cmp = _mm256_cmpeq_epi8(q_chunk, t_chunk);
            let mask = _mm256_movemask_epi8(cmp);
            
            if mask == -1 {
                // All 32 bytes match
                i += 32;
            } else {
                // Find first mismatch (count trailing ones in mask)
                i += mask.trailing_ones() as usize;
                break;
            }
        } else {
            // Scalar tail
            if q[i] != t[j_usize] { break; }
            i += 1;
        }
    }
    i
}
```

---

### 2. Wavefront Computation (line 279-309)

**Current implementation:**
```rust
for k in lo..=hi {
    // Compute best predecessor for diagonal k (scalar)
    let from_mm = prev[ki] + 1;
    let from_ins = prev[ki-1] + 1;
    let from_del = prev[ki+1];
    let best = max(from_mm, from_ins, from_del);
    
    wf_next[ki] = extend(best, k);
    ops_next[ki] = best_op;
}
```

**SIMD opportunity:** Process 8 diagonals in parallel
- Pack 8 consecutive diagonals into AVX2 registers
- Compute 8× `max(prev[ki], prev[ki-1]+1, prev[ki+1])` in parallel
- Requires gather/scatter for non-contiguous `prev[ki±1]` access

**Expected speedup:**
- Best case (many active diagonals): 4-6× faster wavefront computation
- Typical bacterial genomics (s=200-500 → ~10-50 active diagonals): 2-3× if wavefront is >30% of runtime
- Diminishing returns: WFA already O(s·n), not O(n×m)

**Implementation complexity:**
- **High**: ~150-200 lines, significant refactoring
- Must restructure loop to process 8 diagonals at once
- Gather/scatter instructions or manual packing
- More complex edge case handling

---

## Profiling Required

**Before implementing SIMD, measure:**

1. **What % of WFA time is in `extend` vs wavefront?**
   - If `extend` >50%: SIMD extension high ROI
   - If wavefront >30%: SIMD wavefront worth considering
   
2. **Average match run length**
   - If avg run <16 bytes: SIMD overhead dominates, not worth it
   - If avg run >32 bytes: SIMD big win

3. **Number of active diagonals per wavefront**
   - If <8: SIMD wavefront not worth it (underutilizes vector registers)
   - If >16: SIMD wavefront good ROI

**Profiling methods:**

```bash
# Method 1: cargo flamegraph (best)
cargo flamegraph --bench wfa_simd_bench -- "10kb"
# Look for hotspots in flamegraph.svg

# Method 2: perf (Linux)
RUSTFLAGS="-C target-cpu=native" cargo bench --bench wfa_simd_bench --no-run
perf record -g target/release/deps/wfa_simd_bench-* "10kb"
perf report

# Method 3: Instrumented code
# Add counters to fill_wfa_fitting:
// At top of function:
use std::sync::atomic::{AtomicUsize, Ordering};
static EXTEND_CALLS: AtomicUsize = AtomicUsize::new(0);
static EXTEND_BASES: AtomicUsize = AtomicUsize::new(0);

// In extend closure:
let start = i;
// ... loop ...
EXTEND_CALLS.fetch_add(1, Ordering::Relaxed);
EXTEND_BASES.fetch_add((i - start) as usize, Ordering::Relaxed);

// After benchmark:
println!("Avg match run: {} bases", EXTEND_BASES.load(Ordering::Relaxed) / EXTEND_CALLS.load(Ordering::Relaxed));
```

---

## Recommendation

### Step 1: Profile (1-2 hours)
Run flamegraph on realistic workload (150bp reads vs 10kb windows, 95% identity):
```bash
cargo install flamegraph
RUSTFLAGS="-C target-cpu=native" cargo flamegraph --release -p phraya-align --bench wfa_simd_bench -- "10kb"
```

### Step 2: If `extend` >50% of runtime → Implement SIMD Extension (1-2 days)
- Start with AVX2 (x86_64)
- Add SSE4.2 fallback
- Add NEON (aarch64)
- Runtime dispatch via `is_x86_feature_detected!`
- Benchmark to confirm 2-4× speedup

### Step 3: If wavefront >30% → Consider SIMD Wavefront (3-5 days)
- Higher complexity, more refactoring
- Only worth it if profiling shows clear bottleneck
- Expected speedup: 2-3×

### Step 4: Don't Implement Both Speculatively
- Profile → implement highest ROI → re-profile → repeat
- Diminishing returns: WFA already subquadratic

---

## When NOT to Use SIMD Intrinsics

**Don't implement SIMD if:**

1. **Profiling shows <20% time in targetable operations**
   - SIMD overhead (bounds checks, alignment, tail handling) dominates
   - 10% faster overall = weeks of `unsafe` code → not worth it

2. **`wide` crate already provides comparable speedup**
   - If `wide::i32x8` gets you 80% of intrinsic speedup with safe code → use `wide`
   - Only use raw intrinsics if measured gap >2×

3. **Algorithm is already fast enough**
   - Current WFA: ~100-200ms for 10kb workload
   - If that's acceptable for your use case → don't optimize further
   - Bacterial genomics: alignment not the bottleneck (I/O, variant calling are)

4. **Maintenance cost too high**
   - 3 platform variants (AVX2, SSE4.2, NEON) × testing/debugging
   - `unsafe` code requires careful review
   - More code to maintain for 2-3× speedup → may not be worth it

---

## Current Status

**No SIMD in WFA yet.** Current implementation is scalar only.

**`wide` crate was used in O(n×m) diagonal DP** (now deprecated for production):
- `wide::i32x8` for portable SIMD
- Worked, but wrong algorithm (28× slower than WFA)
- Not evidence that `wide` is bad — just that O(n×m) is bad

**For WFA, raw intrinsics likely better than `wide`:**
- `extend` needs byte-level `memcmp` → `_mm256_cmpeq_epi8` + `_mm256_movemask_epi8`
- `wide` doesn't provide byte-compare primitives (focused on i32/f32 arithmetic)
- Would need to drop down to `safe_arch` or raw intrinsics anyway

---

## Conclusion

**Answer to original question:** Yes, platform intrinsics (AVX2/NEON) *could* provide 2-4× speedup in WFA, but:

1. **Profile first** - Don't implement blind
2. **Target `extend` first** - Likely highest ROI
3. **Measure before/after** - Confirm speedup is real
4. **Consider maintenance cost** - 2× faster for 3× more code?

**Current priority:** WFA dispatch fix (done) → 28× speedup was the big win. SIMD on top of that is incremental (2-4×), not transformational.

**When to prioritize SIMD WFA:**
- If profiling shows `extend` >50% of runtime
- If avg match run >32 bytes
- If alignment is a measured bottleneck in real workflows

**When to skip SIMD WFA:**
- If current performance acceptable
- If other parts of pipeline are slower (I/O, variant calling)
- If maintenance cost outweighs 2-3× speedup
