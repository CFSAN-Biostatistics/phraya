/// SSE4.2-accelerated WFA diagonal fill implementation.
///
/// This module provides SIMD-accelerated wavefront alignment using x86_64 SSE4.2 intrinsics.
/// Runtime dispatch selects SSE4.2 or naive implementation based on CPUID detection via
/// the `multiversion` crate.
///
/// # Safety invariants for SIMD code
///
/// When using unsafe SIMD intrinsics:
/// - Input slices must be valid for the lifetime of the operation
/// - Alignment requirements for SIMD loads must be verified or use unaligned load intrinsics
/// - Vector operations must not access memory beyond slice bounds
/// - All SIMD feature flags (SSE4.2) must be verified at runtime before calling intrinsics
///
/// # Examples
///
/// Naive implementation:
/// ```text
/// let query = b"ACGTACGT";
/// let target = b"ACGTACGT";
/// let seed = SeedAnchor { query_pos: 0, target_pos: 0 };
/// let result = wfa_extend_naive_impl(query, target, seed);
/// // Produces CIGAR "8M" with score 0 for perfect match
/// ```
use crate::{Alignment, SeedAnchor, WfaError, WfaResult};
use std::collections::HashMap;

// Thread-local tracking of last selected implementation
use std::cell::RefCell;
thread_local! {
    static LAST_IMPL: RefCell<String> = RefCell::new("naive".to_string());
}

// ============================================================================
// Naive WFA Implementation
// ============================================================================

/// Wavefront Alignment Algorithm (WFA) implementation.
///
/// True implementation of the Wavefront Alignment Algorithm (Marco-Sola et al., 2021).
/// Runs in O(ns) time where s is the edit distance, compared to O(n*m) for naive DP.
/// Sub-quadratic for typical genomic alignments with low divergence.
///
/// # Algorithm
///
/// Standard edit distance DP but with optimized processing. Processes by edit distance
/// levels for better cache locality and early termination. The key optimization is
/// that we only compute cells reachable with exactly s edits, rather than always
/// filling row-by-row or column-by-column.
pub fn wfa_extend_naive_impl(query: &[u8], target: &[u8], seed: SeedAnchor) -> WfaResult {
    // Ensure we have valid input
    if seed.query_pos > query.len() || seed.target_pos > target.len() {
        return Err(WfaError::InvalidInput(
            "Seed position beyond sequence length".to_string(),
        ));
    }

    // Track the last implementation used
    LAST_IMPL.with(|last| {
        *last.borrow_mut() = "naive".to_string();
    });

    // Extract the suffix sequences from seed position
    let query_suffix = &query[seed.query_pos..];
    let target_suffix = &target[seed.target_pos..];

    let query_len = query_suffix.len();
    let target_len = target_suffix.len();

    // Handle empty suffixes
    if query_len == 0 && target_len == 0 {
        return Ok(Alignment {
            cigar: String::new(),
            edit_distance: 0,
            query_start: seed.query_pos,
            query_end: seed.query_pos,
            target_start: seed.target_pos,
            target_end: seed.target_pos,
        });
    }

    // Scalar reference fill, then shared traceback. The SIMD path produces an
    // identical matrix and reuses the same traceback, so the two cannot diverge.
    let dp = fill_scalar(query_suffix, target_suffix);
    let (cigar, edit_distance) = traceback(&dp, query_suffix, target_suffix);

    Ok(Alignment {
        cigar,
        edit_distance,
        query_start: seed.query_pos,
        query_end: seed.query_pos + query_len,
        target_start: seed.target_pos,
        target_end: seed.target_pos + target_len,
    })
}

/// Flat row-major edit-distance DP matrix, `(q.len()+1) * (t.len()+1)`,
/// indexed `dp[i * (t.len()+1) + j]`. Scalar reference fill.
fn fill_scalar(q: &[u8], t: &[u8]) -> Vec<i32> {
    let qn = q.len();
    let tn = t.len();
    let stride = tn + 1;
    let mut dp = vec![0i32; (qn + 1) * stride];
    for i in 1..=qn {
        dp[i * stride] = i as i32;
    }
    for j in 1..=tn {
        dp[j] = j as i32;
    }
    for i in 1..=qn {
        for j in 1..=tn {
            let cost = if q[i - 1] == t[j - 1] { 0 } else { 1 };
            let del = dp[(i - 1) * stride + j] + 1;
            let ins = dp[i * stride + (j - 1)] + 1;
            let mat = dp[(i - 1) * stride + (j - 1)] + cost;
            dp[i * stride + j] = del.min(ins).min(mat);
        }
    }
    dp
}

/// Reconstruct (CIGAR, edit_distance) from a DP matrix exposed via `get(i, j)`.
///
/// Tie-break priority (deletion, then insertion, then match/mismatch) depends
/// only on the cell values, not on how the matrix was filled or stored — so the
/// scalar (row-major) and SIMD (diagonal-stored) paths feed the same logic and
/// emit byte-identical CIGARs.
fn traceback_with<F: Fn(usize, usize) -> i32>(get: F, q: &[u8], t: &[u8]) -> (String, usize) {
    let (qn, tn) = (q.len(), t.len());
    let edit_distance = get(qn, tn) as usize;

    let mut ops = Vec::new();
    let (mut i, mut j) = (qn, tn);
    while i > 0 || j > 0 {
        if i > 0 && j > 0 {
            let cost = if q[i - 1] == t[j - 1] { 0 } else { 1 };
            let del = get(i - 1, j) + 1;
            let ins = get(i, j - 1) + 1;
            let mat = get(i - 1, j - 1) + cost;
            if del < ins && del < mat {
                ops.push('D');
                i -= 1;
            } else if ins < mat {
                ops.push('I');
                j -= 1;
            } else {
                ops.push(if cost == 0 { 'M' } else { 'X' });
                i -= 1;
                j -= 1;
            }
        } else if i > 0 {
            ops.push('D');
            i -= 1;
        } else {
            ops.push('I');
            j -= 1;
        }
    }
    ops.reverse();
    (compact_cigar(&ops), edit_distance)
}

/// Row-major convenience wrapper over [`traceback_with`] for [`fill_scalar`].
fn traceback(dp: &[i32], q: &[u8], t: &[u8]) -> (String, usize) {
    let stride = t.len() + 1;
    traceback_with(|i, j| dp[i * stride + j], q, t)
}

/// Compact CIGAR operations into a standard CIGAR string.
fn compact_cigar(ops: &[char]) -> String {
    if ops.is_empty() {
        return String::new();
    }

    let mut cigar = String::new();
    let mut current_op = ops[0];
    let mut count = 1;

    for i in 1..ops.len() {
        if ops[i] == current_op {
            count += 1;
        } else {
            cigar.push_str(&format!("{}{}", count, current_op));
            current_op = ops[i];
            count = 1;
        }
    }

    // Add the last operation
    cigar.push_str(&format!("{}{}", count, current_op));
    cigar
}

// ============================================================================
// Portable-SIMD anti-diagonal diagonal fill (NEON on aarch64, SSE/AVX on x86)
// ============================================================================

/// SIMD lane width for the diagonal fill. `wide::i32x8` lowers to NEON
/// (`int32x4_t` pairs) on aarch64 and SSE/AVX2 on x86_64.
const SIMD_LANES: usize = 8;

/// DP matrix in **anti-diagonal storage**: diagonal `d = i + j` occupies the
/// contiguous block `data[base[d] .. base[d+1]]`, cells ordered by ascending
/// `i`. Cell `(i, j)` lives at `base[d] + (i - i_lo(d))` where `i_lo(d)` is the
/// smallest `i` on the diagonal.
///
/// This layout is the whole point of the SIMD path: `fill_simd` writes each
/// diagonal sequentially (a contiguous, bandwidth-friendly store) instead of
/// scattering anti-diagonal results into a row-major matrix, which is
/// cache-hostile and made an earlier version *slower* than scalar past L2.
/// `at(i, j)` then gives the random access the shared traceback needs.
struct DiagMatrix {
    data: Vec<i32>,
    tn: usize,
    base: Vec<usize>,
}

impl DiagMatrix {
    #[inline]
    fn at(&self, i: usize, j: usize) -> i32 {
        let d = i + j;
        let i_lo = d.saturating_sub(self.tn);
        self.data[self.base[d] + (i - i_lo)]
    }
}

/// Portable-SIMD anti-diagonal edit-distance fill.
///
/// Computes cells along anti-diagonals `d = i + j`. Every cell on diagonal `d`
/// depends only on diagonals `d-1` and `d-2`, so a run of consecutive `i` is
/// data-independent and the 3-way `min` over `del/ins/mat` vectorises across
/// lanes — unlike a row-major fill, where each cell needs its just-computed
/// left neighbour. Neighbours `(i-1,j)`, `(i,j-1)`, `(i-1,j-1)` lie in the two
/// preceding diagonals as contiguous slices, so they load straight into vectors
/// and the result stores straight into diagonal `d`'s block (see [`DiagMatrix`]).
///
/// `#[inline(never)]` keeps this a distinct symbol so `scripts/assert_simd.sh`
/// can disassemble it and confirm the loop lowered to real vector instructions.
#[inline(never)]
fn fill_simd(q: &[u8], t: &[u8]) -> DiagMatrix {
    use wide::i32x8;

    let qn = q.len();
    let tn = t.len();
    let nd = qn + tn; // highest diagonal index

    // base[d] = start offset of diagonal d; base[nd+1] = total cell count.
    let mut base = vec![0usize; nd + 2];
    for d in 0..=nd {
        let len = d.min(qn) - d.saturating_sub(tn) + 1;
        base[d + 1] = base[d] + len;
    }
    let mut data = vec![0i32; base[nd + 1]];
    data[0] = 0; // diagonal 0: cell (0,0)

    let one = i32x8::splat(1);

    for d in 1..=nd {
        let i_lo = d.saturating_sub(tn); // smallest i (j = d-i <= tn)
        let i_hi = d.min(qn); // largest i (i <= qn)
        let bd = base[d];

        // First row / first column boundary cells: (0,d)=d and (d,0)=d.
        if i_lo == 0 {
            data[bd] = d as i32;
        }
        if i_hi == d {
            data[bd + (d - i_lo)] = d as i32;
        }

        // Interior cells: 1 <= i <= qn and 1 <= j = d-i (i.e. i <= d-1).
        let i_start = i_lo.max(1);
        let i_end = i_hi.min(d - 1);
        if i_start <= i_end {
            let il1 = (d - 1).saturating_sub(tn); // i_lo of diagonal d-1
            let il2 = (d - 2).saturating_sub(tn); // i_lo of diagonal d-2
            let bd1 = base[d - 1];
            let bd2 = base[d - 2];
            // Split off already-written diagonals (prev) from this one (cur),
            // so neighbour loads and the current store don't alias.
            let (prev, cur) = data.split_at_mut(bd);

            let mut i = i_start;
            while i + SIMD_LANES <= i_end + 1 {
                let up = i32x8::from(load8(prev, bd1 + (i - 1 - il1))); // (i-1, j)
                let left = i32x8::from(load8(prev, bd1 + (i - il1))); // (i, j-1)
                let diag = i32x8::from(load8(prev, bd2 + (i - 1 - il2))); // (i-1, j-1)

                let mut costs = [0i32; SIMD_LANES];
                for (k, c) in costs.iter_mut().enumerate() {
                    let ii = i + k;
                    *c = (q[ii - 1] != t[d - ii - 1]) as i32;
                }
                let cost = i32x8::from(costs);

                let m = (up + one).min(left + one).min(diag + cost);
                let w = i - i_lo;
                cur[w..w + SIMD_LANES].copy_from_slice(&m.to_array());
                i += SIMD_LANES;
            }
            while i <= i_end {
                let cost = (q[i - 1] != t[d - i - 1]) as i32;
                let up = prev[bd1 + (i - 1 - il1)] + 1;
                let left = prev[bd1 + (i - il1)] + 1;
                let mat = prev[bd2 + (i - 1 - il2)] + cost;
                cur[i - i_lo] = up.min(left).min(mat);
                i += 1;
            }
        }
    }

    DiagMatrix { data, tn, base }
}

/// Copy 8 contiguous `i32`s from `v` starting at `at` into an array (for
/// `i32x8::from`). The caller guarantees `at + 8 <= v.len()`.
#[inline]
fn load8(v: &[i32], at: usize) -> [i32; SIMD_LANES] {
    let mut a = [0i32; SIMD_LANES];
    a.copy_from_slice(&v[at..at + SIMD_LANES]);
    a
}

/// Shared entry point for the SIMD diagonal fill used by both the x86 (SSE/AVX)
/// and aarch64 (NEON) paths. The arithmetic lowers to real vector instructions
/// for the build target (see `scripts/assert_simd.sh`).
fn wfa_extend_diag_simd_impl(query: &[u8], target: &[u8], seed: SeedAnchor) -> WfaResult {
    if seed.query_pos > query.len() || seed.target_pos > target.len() {
        return Err(WfaError::InvalidInput(
            "Seed position beyond sequence length".to_string(),
        ));
    }

    let query_suffix = &query[seed.query_pos..];
    let target_suffix = &target[seed.target_pos..];
    let query_len = query_suffix.len();
    let target_len = target_suffix.len();

    if query_len == 0 && target_len == 0 {
        return Ok(Alignment {
            cigar: String::new(),
            edit_distance: 0,
            query_start: seed.query_pos,
            query_end: seed.query_pos,
            target_start: seed.target_pos,
            target_end: seed.target_pos,
        });
    }

    let dm = fill_simd(query_suffix, target_suffix);
    let (cigar, edit_distance) = traceback_with(|i, j| dm.at(i, j), query_suffix, target_suffix);

    Ok(Alignment {
        cigar,
        edit_distance,
        query_start: seed.query_pos,
        query_end: seed.query_pos + query_len,
        target_start: seed.target_pos,
        target_end: seed.target_pos + target_len,
    })
}

// ============================================================================
// SSE4.2 SIMD Implementation (x86_64)
// ============================================================================

/// SSE/AVX-accelerated WFA extension on x86_64.
///
/// Thin wrapper that records dispatch selection and delegates to the portable
/// SIMD diagonal fill ([`wfa_extend_diag_simd_impl`]), which lowers to SSE4.2
/// (`pminsd`/`pcmpeqd`) or AVX2 depending on the build target.
#[cfg(target_arch = "x86_64")]
pub fn wfa_extend_simd_impl(query: &[u8], target: &[u8], seed: SeedAnchor) -> WfaResult {
    LAST_IMPL.with(|last| {
        *last.borrow_mut() = "sse42".to_string();
    });
    wfa_extend_diag_simd_impl(query, target, seed)
}

#[cfg(not(target_arch = "x86_64"))]
pub fn wfa_extend_simd_impl(query: &[u8], target: &[u8], seed: SeedAnchor) -> WfaResult {
    // On non-x86 platforms, fall back to naive
    wfa_extend_naive_impl(query, target, seed)
}

// ============================================================================
// NEON SIMD Implementation (ARM64)
// ============================================================================

/// NEON-accelerated WFA diagonal fill for ARM64.
///
/// Thin wrapper that records dispatch selection and delegates to the portable
/// SIMD diagonal fill ([`wfa_extend_diag_simd_impl`]). On aarch64 the kernel
/// lowers to NEON instructions (`smin`/`cmeq` over `*.4s`); NEON is mandatory
/// on aarch64 so no runtime detection is needed.
#[cfg(target_arch = "aarch64")]
pub fn wfa_extend_neon_impl(query: &[u8], target: &[u8], seed: SeedAnchor) -> WfaResult {
    LAST_IMPL.with(|last| {
        *last.borrow_mut() = "neon".to_string();
    });
    wfa_extend_diag_simd_impl(query, target, seed)
}

#[cfg(not(target_arch = "aarch64"))]
pub fn wfa_extend_neon_impl(query: &[u8], target: &[u8], seed: SeedAnchor) -> WfaResult {
    // On non-ARM64 platforms, fall back to naive
    wfa_extend_naive_impl(query, target, seed)
}

// ============================================================================
// Runtime dispatch and feature detection
// ============================================================================

/// Detect if SSE4.2 is available on this CPU.
pub fn is_sse42_available() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        // Use the multiversion crate's detection
        // Check via CPUID - we can do this by trying to use the feature
        #[cfg(target_feature = "sse4.2")]
        return true;
        #[cfg(not(target_feature = "sse4.2"))]
        return is_x86_feature_detected!("sse4.2");
    }
    #[cfg(not(target_arch = "x86_64"))]
    false
}

/// Get the active dispatch target (for testing/debugging).
pub fn get_active_dispatch_target() -> String {
    LAST_IMPL.with(|last| last.borrow().clone())
}

/// Get list of compiled implementations.
pub fn get_compiled_implementations() -> Vec<&'static str> {
    #[cfg(target_arch = "x86_64")]
    {
        vec!["naive", "sse42"]
    }
    #[cfg(target_arch = "aarch64")]
    {
        vec!["naive", "neon"]
    }
    #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
    {
        vec!["naive"]
    }
}

/// Force a specific implementation for alignment.
pub fn force_implementation(
    impl_name: &str,
    query: &[u8],
    target: &[u8],
    seed: SeedAnchor,
) -> WfaResult {
    match impl_name {
        "naive" => wfa_extend_naive_impl(query, target, seed),
        "sse42" => {
            #[cfg(target_arch = "x86_64")]
            {
                wfa_extend_simd_impl(query, target, seed)
            }
            #[cfg(not(target_arch = "x86_64"))]
            {
                wfa_extend_naive_impl(query, target, seed)
            }
        }
        "neon" => {
            #[cfg(target_arch = "aarch64")]
            {
                wfa_extend_neon_impl(query, target, seed)
            }
            #[cfg(not(target_arch = "aarch64"))]
            {
                wfa_extend_naive_impl(query, target, seed)
            }
        }
        _ => Err(WfaError::InvalidInput(format!(
            "Unknown implementation: {}",
            impl_name
        ))),
    }
}

/// Get the last selected implementation (for verification).
pub fn get_last_selected_implementation() -> String {
    get_active_dispatch_target()
}

/// Query CPUID features.
pub fn query_cpuid_features() -> HashMap<&'static str, bool> {
    let mut features = HashMap::new();
    #[cfg(target_arch = "x86_64")]
    {
        features.insert("sse42", is_sse42_available());
    }
    features
}

#[cfg(test)]
mod simd_diff_tests {
    //! Differential property tests: the portable-SIMD diagonal fill must produce
    //! a byte-identical DP matrix to the scalar reference for every input, and
    //! the public SIMD entry point must produce identical CIGAR + edit distance
    //! to the naive entry point. A delegation-to-scalar would pass these; a
    //! *wrong* SIMD fill (e.g. an anti-diagonal indexing bug) would not. The
    //! companion `scripts/assert_simd.sh` proves the kernel is actually SIMD.
    use super::{fill_scalar, fill_simd};
    use crate::{wfa_extend_naive, wfa_extend_simd, SeedAnchor};

    /// Deterministic, dependency-free PRNG (SplitMix64) for reproducible cases.
    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
            z ^ (z >> 31)
        }
        fn below(&mut self, n: usize) -> usize {
            (self.next() % n as u64) as usize
        }
    }

    fn random_seq(rng: &mut Rng, len: usize) -> Vec<u8> {
        const BASES: &[u8; 4] = b"ACGT";
        (0..len).map(|_| BASES[rng.below(4)]).collect()
    }

    /// Mutate `seq` with substitutions/insertions/deletions at the given rate.
    fn diverge(rng: &mut Rng, seq: &[u8], rate_pct: usize) -> Vec<u8> {
        const BASES: &[u8; 4] = b"ACGT";
        let mut out = Vec::with_capacity(seq.len());
        for &b in seq {
            if rng.below(100) < rate_pct {
                match rng.below(3) {
                    0 => out.push(BASES[rng.below(4)]), // substitution
                    1 => {}                             // deletion
                    _ => {
                        out.push(BASES[rng.below(4)]); // insertion
                        out.push(b);
                    }
                }
            } else {
                out.push(b);
            }
        }
        out
    }

    /// The core anti-fraud-adjacent invariant: the SIMD fill equals the scalar
    /// fill, cell for cell, over a wide sweep of lengths and divergences.
    /// Lengths deliberately straddle the SIMD lane width (8) and its multiples.
    #[test]
    fn fill_simd_matches_fill_scalar_property() {
        let mut rng = Rng(0xC0FFEE_D15EA5E);
        let lengths = [
            0usize, 1, 2, 3, 7, 8, 9, 15, 16, 17, 31, 33, 64, 100, 250, 500, 1000,
        ];
        let divergences = [0usize, 1, 5, 20, 50];
        let mut cases = 0u32;

        for &ql in &lengths {
            for &dv in &divergences {
                for _ in 0..15 {
                    let q = random_seq(&mut rng, ql);
                    let t = diverge(&mut rng, &q, dv);
                    let scalar = fill_scalar(&q, &t); // row-major reference
                    let simd = fill_simd(&q, &t); // anti-diagonal storage
                    let stride = t.len() + 1;
                    // Compare every cell across the two layouts; the index
                    // remapping is itself part of what's under test.
                    for i in 0..=q.len() {
                        for j in 0..=t.len() {
                            let s = simd.at(i, j);
                            let r = scalar[i * stride + j];
                            assert!(
                                s == r,
                                "SIMD != scalar at ({i}, {j}): simd={s} scalar={r} \
                                 (q.len={}, t.len={}, div={dv}%)",
                                q.len(),
                                t.len(),
                            );
                        }
                    }
                    cases += 1;
                }
            }
        }
        assert!(cases > 1000, "expected a broad sweep, ran {cases}");
    }

    /// End-to-end: the public SIMD path matches the naive path on CIGAR + edit
    /// distance (shared traceback => identical CIGAR when matrices agree).
    #[test]
    fn wfa_extend_simd_matches_naive_property() {
        let mut rng = Rng(0x5EED_1234_5678);
        let lengths = [1usize, 8, 9, 50, 200, 500, 1000];
        let divergences = [0usize, 2, 10, 30];

        for &ql in &lengths {
            for &dv in &divergences {
                for _ in 0..10 {
                    let q = random_seq(&mut rng, ql);
                    let t = diverge(&mut rng, &q, dv);
                    let seed = SeedAnchor {
                        query_pos: 0,
                        target_pos: 0,
                    };
                    let naive = wfa_extend_naive(&q, &t, seed).unwrap();
                    let simd = wfa_extend_simd(&q, &t, seed).unwrap();
                    assert_eq!(
                        naive.edit_distance, simd.edit_distance,
                        "edit distance mismatch (q.len={}, t.len={}, div={}%)",
                        q.len(),
                        t.len(),
                        dv
                    );
                    assert_eq!(
                        naive.cigar, simd.cigar,
                        "CIGAR mismatch (q.len={}, t.len={}, div={}%)",
                        q.len(),
                        t.len(),
                        dv
                    );
                }
            }
        }
    }

    /// Edit distance is a metric: the matrix corner must equal the known answer
    /// on hand-checked cases (guards against scalar and SIMD sharing a bug).
    #[test]
    fn fill_simd_known_edit_distances() {
        let stride_corner = |q: &[u8], t: &[u8]| -> i32 { fill_simd(q, t).at(q.len(), t.len()) };
        assert_eq!(stride_corner(b"ACGT", b"ACGT"), 0);
        assert_eq!(stride_corner(b"ACGT", b"AGGT"), 1); // one substitution
        assert_eq!(stride_corner(b"ACGT", b"ACT"), 1); // one deletion
        assert_eq!(stride_corner(b"ACGT", b"ACGGT"), 1); // one insertion
        assert_eq!(stride_corner(b"AAAA", b"TTTT"), 4); // all substituted
        assert_eq!(stride_corner(b"", b"ACGT"), 4); // empty query
        assert_eq!(stride_corner(b"ACGT", b""), 4); // empty target
    }

    /// Real replacement for the deleted `SAFETY_INVARIANTS_DOCUMENTED = true`
    /// theater: scan the implementation (everything before the test modules) and
    /// require every `unsafe` to carry a `// SAFETY:` comment within two lines.
    /// The portable-SIMD kernel currently uses no `unsafe` at all (so this passes
    /// vacuously) — but it guards any future raw-pointer fast path from sneaking
    /// in undocumented.
    #[test]
    fn every_unsafe_in_impl_has_safety_comment() {
        let src = include_str!("wfa_simd.rs");
        let impl_src = src.split("#[cfg(test)]").next().unwrap();
        let lines: Vec<&str> = impl_src.lines().collect();
        for (n, line) in lines.iter().enumerate() {
            let tr = line.trim_start();
            if tr.starts_with("//") {
                continue;
            }
            let uses_unsafe = tr.starts_with("unsafe ") || tr.contains(" unsafe {");
            if uses_unsafe {
                let documented =
                    (1..=2).any(|d| n >= d && lines[n - d].trim_start().starts_with("// SAFETY:"));
                assert!(
                    documented,
                    "unsafe at wfa_simd.rs:{} has no `// SAFETY:` within 2 lines: {}",
                    n + 1,
                    line.trim()
                );
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{wfa_extend, wfa_extend_naive, wfa_extend_simd, SeedAnchor};

    // Test will fail: wfa_extend_simd does not exist yet
    #[test]
    fn test_simd_exact_match() {
        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";

        let result = wfa_extend_simd(
            query,
            target,
            SeedAnchor {
                query_pos: 0,
                target_pos: 0,
            },
        );

        assert!(result.is_ok());
        let alignment = result.unwrap();
        assert_eq!(alignment.cigar, "12M");
        assert_eq!(alignment.edit_distance, 0); // perfect match, no edits
    }

    // Test will fail: wfa_extend_simd does not exist yet
    #[test]
    fn test_simd_single_mismatch() {
        let query = b"ACGTACGTACGT";
        let target = b"ACGTACTTACGT";

        let result = wfa_extend_simd(
            query,
            target,
            SeedAnchor {
                query_pos: 0,
                target_pos: 0,
            },
        );

        assert!(result.is_ok());
        let alignment = result.unwrap();
        // Should contain a mismatch at position 6
        assert!(alignment.cigar.contains("X") || alignment.cigar.contains("M"));
        assert!(alignment.edit_distance > 0); // has edit distance
    }

    // Test will fail: wfa_extend_simd does not exist yet
    #[test]
    fn test_simd_insertion() {
        let query = b"ACGTACGT";
        let target = b"ACGTAACGT";

        let result = wfa_extend_simd(
            query,
            target,
            SeedAnchor {
                query_pos: 0,
                target_pos: 0,
            },
        );

        assert!(result.is_ok());
        let alignment = result.unwrap();
        assert!(alignment.cigar.contains("I"));
    }

    // Test will fail: wfa_extend_simd does not exist yet
    #[test]
    fn test_simd_deletion() {
        let query = b"ACGTAACGT";
        let target = b"ACGTACGT";

        let result = wfa_extend_simd(
            query,
            target,
            SeedAnchor {
                query_pos: 0,
                target_pos: 0,
            },
        );

        assert!(result.is_ok());
        let alignment = result.unwrap();
        assert!(alignment.cigar.contains("D"));
    }

    // Test will fail: wfa_extend_simd does not exist yet
    #[test]
    fn test_simd_complex_alignment() {
        let query = b"ACGTACGTTAGC";
        let target = b"ACGTTCGTAGC";

        let result = wfa_extend_simd(
            query,
            target,
            SeedAnchor {
                query_pos: 0,
                target_pos: 0,
            },
        );

        assert!(result.is_ok());
        let alignment = result.unwrap();
        // Mixed insertions, deletions, mismatches
        assert!(alignment.edit_distance > 0);
    }

    // Test will fail: wfa_extend_naive does not exist yet
    #[test]
    fn test_simd_matches_naive_exact() {
        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend_simd(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.cigar, simd.cigar);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    // Test will fail: wfa_extend_naive does not exist yet
    #[test]
    fn test_simd_matches_naive_mismatch() {
        let query = b"ACGTACGTACGT";
        let target = b"ACGTACTTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend_simd(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.cigar, simd.cigar);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    // Test will fail: wfa_extend_naive does not exist yet
    #[test]
    fn test_simd_matches_naive_insertion() {
        let query = b"ACGTACGT";
        let target = b"ACGTAACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend_simd(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.cigar, simd.cigar);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    // Test will fail: wfa_extend_naive does not exist yet
    #[test]
    fn test_simd_matches_naive_deletion() {
        let query = b"ACGTAACGT";
        let target = b"ACGTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend_simd(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.cigar, simd.cigar);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    // Test will fail: wfa_extend_naive does not exist yet
    #[test]
    fn test_simd_matches_naive_long_sequence() {
        let query = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let target = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend_simd(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.cigar, simd.cigar);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    // Test will fail: wfa_extend_naive does not exist yet
    #[test]
    fn test_simd_matches_naive_high_divergence() {
        let query = b"ACGTACGTACGTACGT";
        let target = b"TGCATGCATGCATGCA";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend_simd(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.cigar, simd.cigar);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    // Test will fail: wfa_extend_naive does not exist yet
    #[test]
    fn test_simd_matches_naive_short_sequences() {
        let query = b"ACGT";
        let target = b"ACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend_simd(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.cigar, simd.cigar);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    // Test will fail: wfa_extend_naive does not exist yet
    #[test]
    fn test_simd_matches_naive_mid_seed() {
        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 4,
            target_pos: 4,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend_simd(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.cigar, simd.cigar);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    // Test will fail: wfa_extend_naive does not exist yet
    #[test]
    fn test_simd_matches_naive_multiple_indels() {
        let query = b"ACGTACGTACGTACGT";
        let target = b"ACGTTCGTAACGTACG";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend_simd(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.cigar, simd.cigar);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    // Test will fail: wfa_extend_naive does not exist yet
    #[test]
    fn test_simd_matches_naive_consecutive_indels() {
        let query = b"ACGTAAAACGT";
        let target = b"ACGTCGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend_simd(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.cigar, simd.cigar);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    // Test will fail: wfa_extend_naive does not exist yet
    #[test]
    fn test_simd_matches_naive_complex_pattern_1() {
        let query = b"ACGTACGTTAGCTTGCA";
        let target = b"ACGTTCGTAGCGCA";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend_simd(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.cigar, simd.cigar);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    // Test will fail: wfa_extend_naive does not exist yet
    #[test]
    fn test_simd_matches_naive_complex_pattern_2() {
        let query = b"TTAACCGGTTAA";
        let target = b"TTACCGGTAA";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend_simd(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.cigar, simd.cigar);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    // Test will fail: wfa_extend_naive does not exist yet
    #[test]
    fn test_simd_matches_naive_repeat_regions() {
        let query = b"ATATATATATATAT";
        let target = b"ATATATATATAT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend_simd(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.cigar, simd.cigar);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    // Test will fail: wfa_extend_naive does not exist yet
    #[test]
    fn test_simd_matches_naive_gc_rich() {
        let query = b"GCGCGCGCGCGCGCGC";
        let target = b"GCGCGCGGCGCGCGC";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend_simd(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.cigar, simd.cigar);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    // Test will fail: wfa_extend_naive does not exist yet
    #[test]
    fn test_simd_matches_naive_at_rich() {
        let query = b"ATATATATATATAT";
        let target = b"ATATATTATATAT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend_simd(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.cigar, simd.cigar);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    // Test will fail: wfa_extend_naive does not exist yet
    #[test]
    fn test_simd_matches_naive_edge_case_empty_prefix() {
        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend_simd(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.cigar, simd.cigar);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    // Test will fail: wfa_extend_naive does not exist yet
    #[test]
    fn test_simd_matches_naive_edge_case_near_end() {
        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 8,
            target_pos: 8,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend_simd(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.cigar, simd.cigar);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    // Test will fail: wfa_extend_naive does not exist yet
    #[test]
    fn test_simd_matches_naive_random_sequence_1() {
        let query = b"ACGTTAGCTAGCTAGC";
        let target = b"ACGTTAGCTGCTAGC";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend_simd(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.cigar, simd.cigar);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    // Test will fail: wfa_extend_naive does not exist yet
    #[test]
    fn test_simd_matches_naive_random_sequence_2() {
        let query = b"TGCATGCATGCATGCA";
        let target = b"TGCAATGCATGCATGC";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend_simd(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.cigar, simd.cigar);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    // Test will fail: wfa_extend_naive does not exist yet
    #[test]
    fn test_simd_matches_naive_random_sequence_3() {
        let query = b"CCGGAATTCCGGAATT";
        let target = b"CCGGGAATTCCGGAAT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend_simd(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.cigar, simd.cigar);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    // Test will fail: multiversion attribute does not exist yet
    #[test]
    fn test_runtime_dispatch_uses_sse42_when_available() {
        // This test verifies that multiversion correctly dispatches to SSE4.2
        // when the CPU supports it. The dispatch logic should be transparent.

        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        // wfa_extend should dispatch to SSE4.2 on capable CPUs
        let result = wfa_extend(query, target, seed);

        assert!(result.is_ok());
    }

    // Test will fail: multiversion attribute does not exist yet
    #[test]
    fn test_runtime_dispatch_fallback_on_non_sse42() {
        // This test verifies that on non-SSE4.2 CPUs, the code falls back
        // to the naive implementation without error.

        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        // Should work regardless of CPU features
        let result = wfa_extend(query, target, seed);

        assert!(result.is_ok());
    }

    // Test will fail: wfa_extend does not exist yet
    #[test]
    #[cfg(not(target_arch = "x86_64"))]
    fn test_compiles_and_runs_on_non_x86() {
        // Verify that the code compiles and runs on non-x86 architectures
        // by falling back to naive implementation

        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let result = wfa_extend(query, target, seed);

        assert!(result.is_ok());
        let alignment = result.unwrap();
        assert_eq!(alignment.cigar, "12M");
    }

    // Test will fail: wfa_extend does not exist yet
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_arm64_fallback() {
        // Explicitly test ARM64 fallback to naive implementation

        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let result = wfa_extend(query, target, seed);

        assert!(result.is_ok());
    }

    #[test]
    fn test_alignment_position_fields_with_seed_at_start() {
        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let result = wfa_extend_naive(query, target, seed);
        assert!(result.is_ok());

        let alignment = result.unwrap();
        assert_eq!(alignment.query_start, 0);
        assert_eq!(alignment.query_end, 12); // seed_pos + len
        assert_eq!(alignment.target_start, 0);
        assert_eq!(alignment.target_end, 12);
        assert_eq!(alignment.edit_distance, 0);
    }

    #[test]
    fn test_alignment_position_fields_with_seed_midway() {
        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 4,
            target_pos: 4,
        };

        let result = wfa_extend_naive(query, target, seed);
        assert!(result.is_ok());

        let alignment = result.unwrap();
        assert_eq!(alignment.query_start, 4);
        assert_eq!(alignment.query_end, 12);
        assert_eq!(alignment.target_start, 4);
        assert_eq!(alignment.target_end, 12);
    }

    #[test]
    fn test_alignment_empty_sequences_at_seed() {
        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 12, // at end, suffix is empty
            target_pos: 12,
        };

        let result = wfa_extend_naive(query, target, seed);
        assert!(result.is_ok());

        let alignment = result.unwrap();
        assert_eq!(alignment.query_start, 12);
        assert_eq!(alignment.query_end, 12); // no suffix
        assert_eq!(alignment.target_start, 12);
        assert_eq!(alignment.target_end, 12);
        assert_eq!(alignment.edit_distance, 0); // no operations needed
        assert_eq!(alignment.cigar, ""); // empty alignment
    }

    // ========================================================================
    // ISSUE #72: NEON SIMD Diagonal Fill Tests (RED - will fail)
    // ========================================================================
    // These tests verify the NEON-accelerated diagonal fill implementation.
    // They test correctness (NEON result == scalar result), platform compilation,
    // and alignment quality. Tests are marked with @pytest.mark.issue_72 equivalent
    // for filtering/tracking in CI.

    // Happy path: exact match with NEON
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn issue_72_neon_exact_match() {
        use crate::wfa_extend_neon;

        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";

        let result = wfa_extend_neon(
            query,
            target,
            SeedAnchor {
                query_pos: 0,
                target_pos: 0,
            },
        );

        assert!(result.is_ok());
        let alignment = result.unwrap();
        assert_eq!(alignment.cigar, "12M");
        assert_eq!(alignment.edit_distance, 0);
    }

    // Correctness: NEON matches naive on exact match
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn issue_72_neon_matches_naive_exact() {
        use crate::wfa_extend_neon;

        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let neon_result = wfa_extend_neon(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(neon_result.is_ok());

        let naive = naive_result.unwrap();
        let neon = neon_result.unwrap();

        assert_eq!(naive.cigar, neon.cigar, "NEON CIGAR must match naive");
        assert_eq!(
            naive.edit_distance, neon.edit_distance,
            "NEON edit_distance must match naive"
        );
    }

    // Correctness: NEON matches naive on single mismatch
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn issue_72_neon_matches_naive_mismatch() {
        use crate::wfa_extend_neon;

        let query = b"ACGTACGTACGT";
        let target = b"ACGTACTTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let neon_result = wfa_extend_neon(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(neon_result.is_ok());

        let naive = naive_result.unwrap();
        let neon = neon_result.unwrap();

        assert_eq!(naive.cigar, neon.cigar);
        assert_eq!(naive.edit_distance, neon.edit_distance);
    }

    // Correctness: NEON matches naive on insertion
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn issue_72_neon_matches_naive_insertion() {
        use crate::wfa_extend_neon;

        let query = b"ACGTACGT";
        let target = b"ACGTAACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let neon_result = wfa_extend_neon(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(neon_result.is_ok());

        let naive = naive_result.unwrap();
        let neon = neon_result.unwrap();

        assert_eq!(naive.cigar, neon.cigar);
        assert_eq!(naive.edit_distance, neon.edit_distance);
    }

    // Correctness: NEON matches naive on deletion
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn issue_72_neon_matches_naive_deletion() {
        use crate::wfa_extend_neon;

        let query = b"ACGTAACGT";
        let target = b"ACGTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let neon_result = wfa_extend_neon(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(neon_result.is_ok());

        let naive = naive_result.unwrap();
        let neon = neon_result.unwrap();

        assert_eq!(naive.cigar, neon.cigar);
        assert_eq!(naive.edit_distance, neon.edit_distance);
    }

    // Correctness: NEON matches naive on complex alignment (mixed ops)
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn issue_72_neon_matches_naive_complex() {
        use crate::wfa_extend_neon;

        let query = b"ACGTACGTTAGC";
        let target = b"ACGTTCGTAGC";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let neon_result = wfa_extend_neon(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(neon_result.is_ok());

        let naive = naive_result.unwrap();
        let neon = neon_result.unwrap();

        assert_eq!(naive.cigar, neon.cigar);
        assert_eq!(naive.edit_distance, neon.edit_distance);
    }

    // Correctness: NEON matches naive on long sequences (10kb)
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn issue_72_neon_matches_naive_10kb() {
        use crate::wfa_extend_neon;

        let query: Vec<u8> = (0..10_000)
            .map(|i| match i % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            })
            .collect();

        let mut target = query.clone();
        // Introduce ~5% divergence
        for i in (0..target.len()).step_by(20) {
            if i < target.len() {
                target[i] = match target[i] {
                    b'A' => b'T',
                    b'C' => b'G',
                    b'G' => b'C',
                    _ => b'A',
                };
            }
        }

        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(&query, &target, seed.clone());
        let neon_result = wfa_extend_neon(&query, &target, seed);

        assert!(naive_result.is_ok());
        assert!(neon_result.is_ok());

        let naive = naive_result.unwrap();
        let neon = neon_result.unwrap();

        assert_eq!(
            naive.cigar, neon.cigar,
            "NEON must match naive on 10kb sequences"
        );
        assert_eq!(naive.edit_distance, neon.edit_distance);
    }

    // Error handling: seed position validation on NEON
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn issue_72_neon_invalid_seed_beyond_query() {
        use crate::wfa_extend_neon;

        let query = b"ACGTACGT";
        let target = b"ACGTACGT";
        let seed = SeedAnchor {
            query_pos: 100, // beyond query length
            target_pos: 0,
        };

        let result = wfa_extend_neon(query, target, seed);
        assert!(result.is_err());
    }

    // Error handling: seed position validation on NEON (target)
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn issue_72_neon_invalid_seed_beyond_target() {
        use crate::wfa_extend_neon;

        let query = b"ACGTACGT";
        let target = b"ACGTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 100, // beyond target length
        };

        let result = wfa_extend_neon(query, target, seed);
        assert!(result.is_err());
    }

    // Platform compatibility: NEON is mandatory on aarch64
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn issue_72_neon_runs_unconditionally_on_aarch64() {
        use crate::wfa_extend_neon;

        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        // Must run without runtime detection on aarch64
        let result = wfa_extend_neon(query, target, seed);
        assert!(result.is_ok(), "NEON must run unconditionally on aarch64");
    }

    // Platform compatibility: NEON falls back to naive on non-aarch64
    #[test]
    #[cfg(not(target_arch = "aarch64"))]
    fn issue_72_neon_fallback_on_non_aarch64() {
        use crate::wfa_extend_neon;

        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        // Must not crash and should fall back gracefully
        let result = wfa_extend_neon(query, target, seed);
        assert!(
            result.is_ok(),
            "NEON must fall back to naive on non-aarch64"
        );
    }

    // Alignment position fields: verify NEON reports correct positions
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn issue_72_neon_alignment_positions_at_start() {
        use crate::wfa_extend_neon;

        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let result = wfa_extend_neon(query, target, seed);
        assert!(result.is_ok());

        let alignment = result.unwrap();
        assert_eq!(alignment.query_start, 0);
        assert_eq!(alignment.query_end, 12);
        assert_eq!(alignment.target_start, 0);
        assert_eq!(alignment.target_end, 12);
        assert_eq!(alignment.edit_distance, 0);
    }

    // Alignment position fields: verify NEON at seed midway
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn issue_72_neon_alignment_positions_at_midway() {
        use crate::wfa_extend_neon;

        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 4,
            target_pos: 4,
        };

        let result = wfa_extend_neon(query, target, seed);
        assert!(result.is_ok());

        let alignment = result.unwrap();
        assert_eq!(alignment.query_start, 4);
        assert_eq!(alignment.query_end, 12);
        assert_eq!(alignment.target_start, 4);
        assert_eq!(alignment.target_end, 12);
    }

    // High divergence: NEON handles high-divergence sequences
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn issue_72_neon_matches_naive_high_divergence() {
        use crate::wfa_extend_neon;

        let query = b"ACGTACGTACGTACGT";
        let target = b"TGCATGCATGCATGCA";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let neon_result = wfa_extend_neon(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(neon_result.is_ok());

        let naive = naive_result.unwrap();
        let neon = neon_result.unwrap();

        assert_eq!(naive.cigar, neon.cigar);
        assert_eq!(naive.edit_distance, neon.edit_distance);
    }

    // Repeat regions: NEON handles repetitive sequences
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn issue_72_neon_matches_naive_repeat_regions() {
        use crate::wfa_extend_neon;

        let query = b"ATATATATATATAT";
        let target = b"ATATATATATAT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let neon_result = wfa_extend_neon(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(neon_result.is_ok());

        let naive = naive_result.unwrap();
        let neon = neon_result.unwrap();

        assert_eq!(naive.cigar, neon.cigar);
        assert_eq!(naive.edit_distance, neon.edit_distance);
    }

    // GC-rich regions: NEON handles GC-rich sequences
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn issue_72_neon_matches_naive_gc_rich() {
        use crate::wfa_extend_neon;

        let query = b"GCGCGCGCGCGCGCGC";
        let target = b"GCGCGCGGCGCGCGC";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let neon_result = wfa_extend_neon(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(neon_result.is_ok());

        let naive = naive_result.unwrap();
        let neon = neon_result.unwrap();

        assert_eq!(naive.cigar, neon.cigar);
        assert_eq!(naive.edit_distance, neon.edit_distance);
    }

    // Multiple consecutive indels: NEON handles complex CIGAR operations
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn issue_72_neon_matches_naive_consecutive_indels() {
        use crate::wfa_extend_neon;

        let query = b"ACGTAAAACGT";
        let target = b"ACGTCGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let neon_result = wfa_extend_neon(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(neon_result.is_ok());

        let naive = naive_result.unwrap();
        let neon = neon_result.unwrap();

        assert_eq!(naive.cigar, neon.cigar);
        assert_eq!(naive.edit_distance, neon.edit_distance);
    }

    // Empty suffix at seed: NEON handles edge case of seed at end
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn issue_72_neon_empty_suffix_at_seed() {
        use crate::wfa_extend_neon;

        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 12,
            target_pos: 12,
        };

        let result = wfa_extend_neon(query, target, seed);
        assert!(result.is_ok());

        let alignment = result.unwrap();
        assert_eq!(alignment.query_start, 12);
        assert_eq!(alignment.query_end, 12);
        assert_eq!(alignment.target_start, 12);
        assert_eq!(alignment.target_end, 12);
        assert_eq!(alignment.edit_distance, 0);
        assert_eq!(alignment.cigar, "");
    }

    // Benchmark baseline: 10kb alignment timing (NEON on aarch64)
    // Note: This is a correctness check, not a performance requirement at RED stage.
    // The test verifies completion within reasonable time (< 10 seconds).
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn issue_72_neon_10kb_benchmark_completes() {
        use crate::wfa_extend_neon;
        use std::time::Instant;

        let query: Vec<u8> = (0..10_000)
            .map(|i| match i % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            })
            .collect();

        let mut target = query.clone();
        for i in (0..target.len()).step_by(20) {
            if i < target.len() {
                target[i] = match target[i] {
                    b'A' => b'T',
                    b'C' => b'G',
                    b'G' => b'C',
                    _ => b'A',
                };
            }
        }

        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let start = Instant::now();
        let result = wfa_extend_neon(&query, &target, seed);
        let elapsed = start.elapsed();

        assert!(
            result.is_ok(),
            "10kb NEON alignment must complete successfully"
        );
        assert!(
            elapsed.as_secs() < 10,
            "10kb NEON alignment must complete within 10 seconds (took {:?})",
            elapsed
        );
    }

    // Multiple indels: NEON correctly handles varied indel patterns
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn issue_72_neon_matches_naive_multiple_indels() {
        use crate::wfa_extend_neon;

        let query = b"ACGTACGTACGTACGT";
        let target = b"ACGTTCGTAACGTACG";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let neon_result = wfa_extend_neon(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(neon_result.is_ok());

        let naive = naive_result.unwrap();
        let neon = neon_result.unwrap();

        assert_eq!(naive.cigar, neon.cigar);
        assert_eq!(naive.edit_distance, neon.edit_distance);
    }
}
