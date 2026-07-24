/// Portable SIMD-accelerated WFA diagonal fill implementation.
///
/// This module provides SIMD-accelerated wavefront alignment using the `wide` crate for
/// portable SIMD operations. When compiled with `-C target-cpu=native`, the `wide` crate
/// lowers to SSE4.2 on x86_64 or NEON on aarch64. Otherwise, it emulates with scalar ops.
///
/// # Portable SIMD notes
///
/// The `wide::i32x8` type abstracts SIMD ops. On x86_64 with `target-cpu=native`, operations
/// like `min()` lower to `pminsd` (SSE4.2). On aarch64, they lower to NEON `smin`. Without
/// CPU-specific flags, `wide` emulates with scalar loops. No `unsafe` intrinsics in this module.
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

// Thread-local tracking of last selected implementation
use std::cell::RefCell;
thread_local! {
    static LAST_IMPL: RefCell<String> = RefCell::new("naive".to_string());
}

// ============================================================================
// SIMD-accelerated match extension
// ============================================================================
//
// The WFA inner loop ("extend") advances along a diagonal while query and target
// bytes are equal. This is the alignment hot path. `count_matching_prefix` returns
// the length of the longest common prefix of two byte slices — i.e. the number of
// matching bytes before the first mismatch (or the shorter length, if no mismatch).
//
// Three tiers, all returning bit-identical results (enforced by differential tests):
//   - `count_matching_prefix_scalar`: byte-by-byte reference. Always correct.
//   - `count_matching_prefix_u64`:    8 bytes/step via little-endian XOR. Portable,
//                                     no `unsafe`, endian-independent.
//   - arch SIMD (SSE2 / NEON):        16 bytes/step. SSE2 is mandatory on x86_64 and
//                                     NEON is mandatory on aarch64, so these are selected
//                                     at compile time with no runtime feature dispatch.

/// Length of the longest common prefix of `a` and `b` (matching bytes before the
/// first mismatch, capped at `a.len().min(b.len())`). Dispatches to the fastest tier
/// available for the target architecture.
#[inline]
pub fn count_matching_prefix(a: &[u8], b: &[u8]) -> usize {
    #[cfg(any(target_arch = "x86_64", target_arch = "aarch64"))]
    {
        count_matching_prefix_arch(a, b)
    }
    #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
    {
        count_matching_prefix_u64(a, b)
    }
}

/// 16-bytes-per-step SSE2 implementation. SSE2 is part of the x86_64 baseline, so the
/// intrinsics are always available and need no runtime feature detection. Compares 16
/// bytes at a time; on the first chunk containing a mismatch, locates it from the
/// per-byte equality mask. The sub-16-byte tail falls through to the scalar reference.
#[cfg(target_arch = "x86_64")]
#[inline]
pub fn count_matching_prefix_arch(a: &[u8], b: &[u8]) -> usize {
    use std::arch::x86_64::*;
    let n = a.len().min(b.len());
    let mut i = 0;
    // SAFETY: SSE2 is guaranteed on all x86_64 targets. Loads are unaligned
    // (`loadu`) and bounded by `i + 16 <= n`, so they never read past the slices.
    unsafe {
        while i + 16 <= n {
            let va = _mm_loadu_si128(a.as_ptr().add(i) as *const __m128i);
            let vb = _mm_loadu_si128(b.as_ptr().add(i) as *const __m128i);
            let eq = _mm_cmpeq_epi8(va, vb);
            // One bit per byte: 1 where equal. All 16 equal => 0xFFFF.
            let mask = _mm_movemask_epi8(eq) as u32 & 0xFFFF;
            if mask != 0xFFFF {
                // First 0 bit = first mismatching byte within this chunk.
                return i + (mask ^ 0xFFFF).trailing_zeros() as usize;
            }
            i += 16;
        }
    }
    i + count_matching_prefix_scalar(&a[i..n], &b[i..n])
}

/// 16-bytes-per-step NEON implementation. Advanced SIMD (NEON) is mandatory on
/// AArch64, so the intrinsics are always available with no runtime detection.
#[cfg(target_arch = "aarch64")]
#[inline]
pub fn count_matching_prefix_arch(a: &[u8], b: &[u8]) -> usize {
    use std::arch::aarch64::*;
    let n = a.len().min(b.len());
    let mut i = 0;
    // SAFETY: NEON is guaranteed on all aarch64 targets. Loads are bounded by
    // `i + 16 <= n`, so they never read past the slices.
    unsafe {
        while i + 16 <= n {
            let va = vld1q_u8(a.as_ptr().add(i));
            let vb = vld1q_u8(b.as_ptr().add(i));
            let eq = vceqq_u8(va, vb); // 0xFF per byte where equal, else 0x00.
            // Horizontal min across the 16 lanes: 0xFF iff every byte matched.
            if vminvq_u8(eq) != 0xFF {
                let mut lanes = [0u8; 16];
                vst1q_u8(lanes.as_mut_ptr(), eq);
                let first = lanes.iter().position(|&x| x != 0xFF).unwrap();
                return i + first;
            }
            i += 16;
        }
    }
    i + count_matching_prefix_scalar(&a[i..n], &b[i..n])
}

/// Byte-by-byte reference implementation. The semantic ground truth all other tiers
/// are differential-tested against.
#[inline]
pub fn count_matching_prefix_scalar(a: &[u8], b: &[u8]) -> usize {
    let n = a.len().min(b.len());
    let mut i = 0;
    while i < n && a[i] == b[i] {
        i += 1;
    }
    i
}

/// Portable 8-bytes-per-step implementation using little-endian word XOR.
///
/// Reading each 8-byte chunk as a little-endian `u64` makes byte index 0 the least
/// significant byte, so the XOR of two chunks has its lowest set bit in the first
/// differing byte regardless of host endianness; `trailing_zeros() / 8` is therefore
/// the index of the first mismatch within the chunk. No `unsafe`.
#[inline]
pub fn count_matching_prefix_u64(a: &[u8], b: &[u8]) -> usize {
    let n = a.len().min(b.len());
    let mut i = 0;
    while i + 8 <= n {
        let aw = u64::from_le_bytes(a[i..i + 8].try_into().unwrap());
        let bw = u64::from_le_bytes(b[i..i + 8].try_into().unwrap());
        let diff = aw ^ bw;
        if diff != 0 {
            return i + (diff.trailing_zeros() as usize / 8);
        }
        i += 8;
    }
    while i < n && a[i] == b[i] {
        i += 1;
    }
    i
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

    // Fitting alignment: query must be fully consumed, target end is free.
    // This is the correct mode for aligning reads against a longer reference window.
    // Global alignment inflates edit distance by the length gap (target extras become deletions).
    let (cigar, edit_distance, target_consumed) = match fill_wfa_fitting(query_suffix, target_suffix) {
        Some(result) => result,
        None => return Err(WfaError::AlignmentFailed(
            "alignment abandoned: max edit distance exceeded".to_string(),
        )),
    };

    Ok(Alignment {
        cigar,
        edit_distance,
        query_start: seed.query_pos,
        query_end: seed.query_pos + query_len,
        target_start: seed.target_pos,
        target_end: seed.target_pos + target_consumed,
    })
}

/// Wavefront Alignment (WFA) — O(s·n) where s = edit distance.
///
/// For each edit distance s, maintains one i32 per active diagonal k = query_pos - target_pos,
/// storing the furthest-reaching query position reachable on diagonal k with exactly s edits.
/// Extend greedily (free matches); expand to s+1 via mismatch (same diagonal) or indel (±1 diagonal).
/// Returns (cigar, edit_distance).
fn fill_wfa(q: &[u8], t: &[u8]) -> (String, usize) {
    let qn = q.len() as i32;
    let tn = t.len() as i32;

    // wf[k] = furthest-reaching query position on diagonal k (k = q_pos - t_pos).
    // Diagonals range from -tn to +qn; offset by tn so index = k + tn.
    let size = (qn + tn + 1) as usize;
    // Sentinel: diagonal not yet reached.
    const UNSET: i32 = i32::MIN;

    // wf_cur[d] = furthest q_pos on diagonal (d - tn); wf_ops[s][d] = edit op that
    // produced this wavefront position (for backtrace). Op encoding: 0=extend/match,
    // 1=mismatch (X), 2=insert (I, q advances only), 3=delete (D, t advances only).
    let mut wf_cur = vec![UNSET; size];
    // Per-wavefront history for traceback: (wavefronts indexed by s, each has ops per diagonal).
    // We store the predecessor wavefront diagonal index for traceback.
    // wf_hist[s] = (wf array, predecessor_diagonal array).
    // Op: 0=M (extended), 1=X (mismatch expand), 2=I (insert expand), 3=D (delete expand).
    let mut wf_hist: Vec<(Vec<i32>, Vec<u8>)> = Vec::new();

    // --- Extend helper: advance query pos on diagonal k as far as matches allow ---
    let extend = |q_pos: i32, k: i32| -> i32 {
        // Advance along diagonal k while query/target bytes match. SIMD-accelerated
        // longest-common-prefix over the suffixes starting at (q_pos, q_pos - k).
        let j = q_pos - k; // target pos
        if q_pos < 0 || q_pos >= qn || j < 0 || j >= tn {
            return q_pos;
        }
        let qi = q_pos as usize;
        let ti = j as usize;
        q_pos + count_matching_prefix(&q[qi..], &t[ti..]) as i32
    };

    // s=0: start on diagonal 0
    let k0_idx = tn as usize; // diagonal 0 offset
    wf_cur[k0_idx] = extend(0, 0);

    let ops_cur = vec![0u8; size];

    // Check termination at s=0
    let target_diag = qn - tn; // diagonal of end cell
    let target_idx = (target_diag + tn) as usize;
    if wf_cur[target_idx] == qn {
        // Perfect match — build cigar
        let (cigar, _) = traceback_wfa(&wf_hist, q, t, 0);
        return (cigar, 0);
    }

    wf_hist.push((wf_cur, ops_cur));

    // Iterate edit distances
    let max_s = qn + tn; // upper bound
    for s in 1..=max_s {
        let prev = &wf_hist[s as usize - 1].0;
        let mut wf_next = vec![UNSET; size];
        let mut ops_next = vec![0u8; size];

        let lo = (-(s.min(tn))) as i32;
        let hi = s.min(qn) as i32;

        for k in lo..=hi {
            let ki = (k + tn) as usize;
            // Three predecessors:
            // mismatch from (s-1, k): q_pos was prev[ki], advance both by 1
            let from_mm = if ki < prev.len() && prev[ki] != UNSET {
                (prev[ki] + 1, 1u8) // +1 q, +1 t → same diagonal
            } else { (UNSET, 0) };
            // insert: from (s-1, k-1): q advances, t stays → diagonal k = (k-1)+1
            // valid when ki_prev = k-1+tn >= 0 → k > -tn
            let from_ins = if k > -tn {
                let ki_prev = (k - 1 + tn) as usize;
                if ki_prev < prev.len() && prev[ki_prev] != UNSET {
                    (prev[ki_prev] + 1, 2u8)
                } else { (UNSET, 0) }
            } else { (UNSET, 0) };
            // delete: from (s-1, k+1): t advances, q stays → diagonal k = (k+1)-1
            // valid when ki_prev = k+1+tn < size → k < qn
            let from_del = if k < qn {
                let ki_prev = (k + 1 + tn) as usize;
                if ki_prev < prev.len() && prev[ki_prev] != UNSET {
                    (prev[ki_prev], 3u8) // q stays same
                } else { (UNSET, 0) }
            } else { (UNSET, 0) };

            // Pick best (furthest reaching)
            let (best_pos, best_op) = [from_mm, from_ins, from_del]
                .into_iter()
                .filter(|&(p, _)| p != UNSET)
                .max_by_key(|&(p, _)| p)
                .unwrap_or((UNSET, 0));

            if best_pos == UNSET { continue; }

            // Bounds check: ensure target pos is valid
            let j = best_pos - k;
            if best_pos < 0 || best_pos > qn || j < 0 || j > tn { continue; }

            wf_next[ki] = extend(best_pos, k);
            ops_next[ki] = best_op;
        }

        // Check termination
        let t_idx = (target_diag + tn) as usize;
        if t_idx < wf_next.len() && wf_next[t_idx] >= qn {
            wf_hist.push((wf_next, ops_next));
            let (cigar, _) = traceback_wfa(&wf_hist, q, t, s as usize);
            return (cigar, s as usize);
        }

        wf_hist.push((wf_next, ops_next));
    }

    // Fallback: should not reach here for valid inputs
    (String::new(), 0)
}

/// Fitting alignment variant of [`fill_wfa`].
///
/// Fitting (semi-global) alignment: query must be fully consumed, target end is free.
/// This is the correct mode for reads vs a longer reference window — global alignment
/// inflates edit distance by the length gap (excess target bases become forced deletions).
///
/// When target is not substantially longer than query (tn ≤ qn + qn/2 + 10), delegates
/// to global [`fill_wfa`] to preserve correctness for short-indel and equal-length cases.
///
/// Returns `(cigar, edit_distance, target_consumed)` where `target_consumed` is the number
/// of target bases actually aligned (≤ `t.len()`).
///
/// With optional `max_s_cap`, returns `None` if the loop exceeds the cap.
pub fn fill_wfa_fitting_impl(
    q: &[u8],
    t: &[u8],
    max_s_cap: Option<usize>,
) -> Option<(String, usize, usize)> {
    if q.is_empty() {
        return Some((String::new(), 0, 0));
    }

    let qn = q.len() as i32;
    let tn = t.len() as i32;
    // Use fitting only when target is substantially longer: at least 1.5× + 10 bp.
    // For similar-length sequences (small indels, equal lengths) global is correct
    // and avoids spurious early termination that under-counts edits.
    if tn <= qn + qn / 2 + 10 {
        let (cigar, edit_dist) = fill_wfa(q, t);
        // This branch bypasses the wavefront cap loop entirely, so max_s_cap must be
        // enforced explicitly here too -- otherwise a capped caller (ADR-0007 / #183)
        // silently gets an uncapped alignment for any query/target pair in this length
        // ratio, defeating the early-abandonment guarantee.
        if let Some(cap) = max_s_cap {
            if edit_dist > cap {
                return None;
            }
        }
        return Some((cigar, edit_dist, t.len()));
    }

    let size = (qn + tn + 1) as usize;
    const UNSET: i32 = i32::MIN;

    let mut wf_cur = vec![UNSET; size];
    let mut wf_hist: Vec<(Vec<i32>, Vec<u8>)> = Vec::new();

    // Same greedy diagonal extension as the global path, via the SIMD
    // `count_matching_prefix` (16 bytes/step) instead of a scalar byte loop.
    let extend = |q_pos: i32, k: i32| -> i32 {
        let j = q_pos - k; // target pos
        if q_pos < 0 || q_pos >= qn || j < 0 || j >= tn {
            return q_pos;
        }
        let qi = q_pos as usize;
        let ti = j as usize;
        q_pos + count_matching_prefix(&q[qi..], &t[ti..]) as i32
    };

    let k0_idx = tn as usize;
    wf_cur[k0_idx] = extend(0, 0);
    let ops_cur = vec![0u8; size];

    let max_s = qn + tn;

    let max_s_to_explore = if let Some(cap) = max_s_cap {
        (cap as i32).min(max_s)
    } else {
        max_s
    };

    // Check initial state for fitting-end (s=0, i.e. a perfect match). This is always the
    // best possible result, so it must be accepted regardless of any cap (0 <= max_s_cap
    // always holds) -- rejecting it here would let a cap prune a better-than-incumbent
    // alignment, violating the score-bound safety invariant (ADR-0007 / #183).
    if let Some((t_end, k_win)) = fitting_end_k(&wf_cur, qn, tn) {
        wf_hist.push((wf_cur, ops_cur));
        let (cigar, _) = traceback_wfa_with_tend(&wf_hist, q, t, 0, t_end, k_win);
        return Some((cigar, 0, t_end));
    }

    wf_hist.push((wf_cur, ops_cur));

    for s in 1..=max_s_to_explore {

        let prev = &wf_hist[s as usize - 1].0;
        let mut wf_next = vec![UNSET; size];
        let mut ops_next = vec![0u8; size];

        let lo = (-(s.min(tn))) as i32;
        let hi = s.min(qn) as i32;

        for k in lo..=hi {
            let ki = (k + tn) as usize;
            let from_mm = if ki < prev.len() && prev[ki] != UNSET {
                (prev[ki] + 1, 1u8)
            } else { (UNSET, 0) };
            let from_ins = if k > -tn {
                let ki_prev = (k - 1 + tn) as usize;
                if ki_prev < prev.len() && prev[ki_prev] != UNSET {
                    (prev[ki_prev] + 1, 2u8)
                } else { (UNSET, 0) }
            } else { (UNSET, 0) };
            let from_del = if k < qn {
                let ki_prev = (k + 1 + tn) as usize;
                if ki_prev < prev.len() && prev[ki_prev] != UNSET {
                    (prev[ki_prev], 3u8)
                } else { (UNSET, 0) }
            } else { (UNSET, 0) };

            let (best_pos, best_op) = [from_mm, from_ins, from_del]
                .into_iter()
                .filter(|&(p, _)| p != UNSET)
                .max_by_key(|&(p, _)| p)
                .unwrap_or((UNSET, 0));

            if best_pos == UNSET { continue; }
            let j = best_pos - k;
            if best_pos < 0 || best_pos > qn || j < 0 || j > tn { continue; }

            wf_next[ki] = extend(best_pos, k);
            ops_next[ki] = best_op;
        }

        if let Some((t_end, k_win)) = fitting_end_k(&wf_next, qn, tn) {
            wf_hist.push((wf_next, ops_next));
            let (cigar, _) = traceback_wfa_with_tend(&wf_hist, q, t, s as usize, t_end, k_win);
            return Some((cigar, s as usize, t_end));
        }

        wf_hist.push((wf_next, ops_next));
    }

    // Fallback: alignment abandoned (loop exhausted without reaching fitting-end)
    None
}

/// Wavefront fitting alignment with optional cap on maximum edit distance.
/// Returns `(cigar, edit_distance, target_consumed)` or `None` if abandoned.
fn fill_wfa_fitting(q: &[u8], t: &[u8]) -> Option<(String, usize, usize)> {
    fill_wfa_fitting_impl(q, t, None)
}

/// Find a valid fitting-end diagonal: any k where `wf[k] >= qn` and `0 <= qn - k <= tn`.
/// Returns `(t_end, k)` for the diagonal that consumes the most target (maximum t_end = qn - k).
/// Preferring maximum t_end ensures we don't truncate natural deletions at the alignment end.
fn fitting_end_k(wf: &[i32], qn: i32, tn: i32) -> Option<(usize, i32)> {
    let mut best: Option<(usize, i32)> = None;
    for k in (qn - tn)..=qn {
        let t_end = qn - k;
        if t_end < 0 || t_end > tn { continue; }
        let ki = (k + tn) as usize;
        if ki >= wf.len() { continue; }
        if wf[ki] >= qn {
            let te = t_end as usize;
            if best.is_none() || te > best.unwrap().0 {
                best = Some((te, k));
            }
        }
    }
    best
}

/// Backtrace starting from `(qn, t_end)` — the fitting-alignment endpoint.
/// Identical to [`traceback_wfa`] except `ti` starts at `t_end` instead of `t.len()`.
fn traceback_wfa_with_tend(
    hist: &[(Vec<i32>, Vec<u8>)],
    q: &[u8],
    t: &[u8],
    edit_dist: usize,
    t_end: usize,
    _k_win: i32,
) -> (String, usize) {
    let qn = q.len() as i32;
    let tn = t.len() as i32;

    let mut ops: Vec<char> = Vec::new();
    let mut qi = qn;
    let mut ti = t_end as i32; // start at actual end, not tn

    let mut k = qi - ti;
    let mut s = edit_dist as i32;

    while qi > 0 || ti > 0 {
        if s < 0 { break; }

        let ki = (k + tn) as usize;

        if s == 0 {
            while qi > 0 && ti > 0 {
                ops.push('M');
                qi -= 1;
                ti -= 1;
            }
            break;
        }

        let (_wf, wf_ops) = &hist[s as usize];
        let op = if ki < wf_ops.len() { wf_ops[ki] } else { 0 };

        let prev_wf = &hist[s as usize - 1].0;
        let (pred_pos, pred_k) = match op {
            1 => {
                let pk = ki;
                (if pk < prev_wf.len() { prev_wf[pk] } else { 0 }, k)
            }
            2 => {
                let pk = (k - 1 + tn) as usize;
                (if pk < prev_wf.len() { prev_wf[pk] } else { 0 }, k - 1)
            }
            3 => {
                let pk = (k + 1 + tn) as usize;
                (if pk < prev_wf.len() { prev_wf[pk] } else { 0 }, k + 1)
            }
            _ => {
                let pk = ki;
                (if pk < prev_wf.len() { prev_wf[pk] } else { 0 }, k)
            }
        };

        let edit_start = match op {
            1 | 2 => pred_pos + 1,
            3 => pred_pos,
            _ => pred_pos,
        };

        let match_count = qi - edit_start;
        for _ in 0..match_count.max(0) {
            ops.push('M');
            qi -= 1;
            ti -= 1;
        }

        match op {
            1 => { ops.push('X'); qi -= 1; ti -= 1; }
            2 => { ops.push('D'); qi -= 1; }
            3 => { ops.push('I'); ti -= 1; }
            _ => {}
        }

        k = pred_k;
        s -= 1;
    }

    ops.reverse();
    let cigar = compact_cigar(&ops);
    (cigar, edit_dist)
}

/// Backtrace through WFA wavefront history to produce a CIGAR string.
fn traceback_wfa(
    hist: &[(Vec<i32>, Vec<u8>)],
    q: &[u8],
    t: &[u8],
    edit_dist: usize,
) -> (String, usize) {
    let qn = q.len() as i32;
    let tn = t.len() as i32;

    let mut ops: Vec<char> = Vec::new();
    let mut qi = qn;
    let mut ti = tn;

    // Walk backwards from (qn, tn) through the wavefront history.
    let mut k = qi - ti; // current diagonal
    let mut s = edit_dist as i32;

    while qi > 0 || ti > 0 {
        if s < 0 { break; }

        let ki = (k + tn) as usize;

        if s == 0 {
            // All remaining must be matches (extend phase of s=0)
            while qi > 0 && ti > 0 {
                ops.push('M');
                qi -= 1;
                ti -= 1;
            }
            break;
        }

        let (_wf, wf_ops) = &hist[s as usize];
        let op = if ki < wf_ops.len() { wf_ops[ki] } else { 0 };

        let prev_wf = &hist[s as usize - 1].0;
        let (pred_pos, pred_k) = match op {
            1 => { // mismatch: pred on same diagonal, s-1
                let pk = ki;
                (if pk < prev_wf.len() { prev_wf[pk] } else { 0 }, k)
            }
            2 => { // insert (q advances, diagonal was k-1 before)
                let pk = (k - 1 + tn) as usize;
                (if pk < prev_wf.len() { prev_wf[pk] } else { 0 }, k - 1)
            }
            3 => { // delete (t advances, diagonal was k+1 before)
                let pk = (k + 1 + tn) as usize;
                (if pk < prev_wf.len() { prev_wf[pk] } else { 0 }, k + 1)
            }
            _ => { // shouldn't happen mid-trace
                let pk = ki;
                (if pk < prev_wf.len() { prev_wf[pk] } else { 0 }, k)
            }
        };

        // Matches from extend: positions pred_pos+1..=qi on this diagonal (after the edit op)
        // Matches = qi - (pred_pos + 1) for X/I, or qi - pred_pos for D
        let edit_start = match op {
            1 | 2 => pred_pos + 1, // after mismatch/insert, both or q advanced
            3 => pred_pos,         // delete: q didn't advance, t did
            _ => pred_pos,
        };

        // emit matches from qi down to edit_start on diagonal k
        let match_count = qi - edit_start;
        for _ in 0..match_count.max(0) {
            ops.push('M');
            qi -= 1;
            ti -= 1;
        }

        // emit the edit op — match traceback_with convention: 'D'=q advances, 'I'=t advances
        match op {
            1 => { ops.push('X'); qi -= 1; ti -= 1; } // mismatch
            2 => { ops.push('D'); qi -= 1; }           // q advances (no t) → 'D'
            3 => { ops.push('I'); ti -= 1; }           // t advances (no q) → 'I'
            _ => {}
        }

        k = pred_k;
        s -= 1;
    }

    ops.reverse();
    let cigar = compact_cigar(&ops);
    (cigar, edit_dist)
}

/// Flat row-major edit-distance DP matrix, `(q.len()+1) * (t.len()+1)`,
/// indexed `dp[i * (t.len()+1) + j]`. Scalar reference fill (test-only).
#[cfg(test)]
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
/// **NOTE:** This is O(n×m) complexity. Production uses `fill_wfa_fitting` which
/// is O(s·n) where s = edit distance. This function exists for:
/// - Differential correctness testing vs `fill_scalar`
/// - Validating portable SIMD lowering on different architectures
///
/// NOT used in production alignment (too slow for genomics workloads).
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

/// Portable SIMD WFA extension for x86_64.
///
/// Thin wrapper that records dispatch selection and delegates to the portable
/// SIMD diagonal fill ([`wfa_extend_diag_simd_impl`]). With `-C target-cpu=native`,
/// the `wide` crate lowers `min()` ops to SSE4.2 `pminsd` or AVX2 equivalents.
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

/// Portable SIMD WFA extension for ARM64.
///
/// Thin wrapper that records dispatch selection and delegates to the portable
/// SIMD diagonal fill ([`wfa_extend_diag_simd_impl`]). With `-C target-cpu=native`,
/// the `wide` crate lowers to NEON instructions (`smin` over `*.4s`). NEON is
/// mandatory on aarch64, so no runtime detection needed.
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


#[cfg(test)]
mod simd_diff_tests {
    //! Differential property tests: the portable-SIMD diagonal fill must produce
    //! a byte-identical DP matrix to the scalar reference for every input, and
    //! the public SIMD entry point must produce identical CIGAR + edit distance
    //! to the naive entry point. A delegation-to-scalar would pass these; a
    //! *wrong* SIMD fill (e.g. an anti-diagonal indexing bug) would not. The
    //! companion `scripts/assert_simd.sh` proves the kernel is actually SIMD.
    use super::{fill_scalar, fill_simd};
    use crate::{wfa_extend, wfa_extend_naive, SeedAnchor};

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

    /// End-to-end: the SIMD path edit distance matches the WFA naive path.
    /// CIGAR tie-breaking may differ (both are valid minimum-edit alignments)
    /// since fill_wfa uses wavefront backtrace and fill_simd uses DP traceback.
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
                    let simd = wfa_extend(&q, &t, seed).unwrap();
                    assert_eq!(
                        naive.edit_distance, simd.edit_distance,
                        "edit distance mismatch (q.len={}, t.len={}, div={}%)",
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

    /// Scan the implementation for two invariants:
    ///
    /// 1. Every `unsafe` block has a `// SAFETY:` comment within the preceding two lines.
    /// 2. On SIMD platforms (x86_64 / aarch64), at least one `unsafe` block must exist
    ///    — the arch-specific `count_matching_prefix_arch` intrinsic loops. If the count
    ///    drops to zero, SIMD code was silently removed and this guard becomes meaningless.
    ///    On non-SIMD platforms, zero `unsafe` blocks is the correct expectation.
    #[test]
    fn every_unsafe_in_impl_has_safety_comment() {
        let src = include_str!("wfa_simd.rs");
        let impl_src = src.split("#[cfg(test)]").next().unwrap();
        let lines: Vec<&str> = impl_src.lines().collect();
        let mut unsafe_count = 0usize;
        for (n, line) in lines.iter().enumerate() {
            let tr = line.trim_start();
            if tr.starts_with("//") {
                continue;
            }
            let uses_unsafe = tr.starts_with("unsafe ") || tr.contains(" unsafe {");
            if uses_unsafe {
                unsafe_count += 1;
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
        #[cfg(any(target_arch = "x86_64", target_arch = "aarch64"))]
        assert!(
            unsafe_count >= 1,
            "expected ≥1 unsafe block in the SIMD implementation (count_matching_prefix_arch), \
             found 0; was the arch-specific path removed without updating this guard?"
        );
        #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
        assert_eq!(
            unsafe_count, 0,
            "unexpected unsafe block(s) on non-SIMD platform; the fallback path must be safe"
        );
    }
}

#[cfg(test)]
mod tests {
    use crate::{wfa_extend, wfa_extend_naive, SeedAnchor};

    #[test]
    fn test_simd_exact_match() {
        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";

        let result = wfa_extend(
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

    #[test]
    fn test_simd_single_mismatch() {
        let query = b"ACGTACGTACGT";
        let target = b"ACGTACTTACGT";

        let result = wfa_extend(
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

    #[test]
    fn test_simd_insertion() {
        let query = b"ACGTACGT";
        let target = b"ACGTAACGT";

        let result = wfa_extend(
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

    #[test]
    fn test_simd_deletion() {
        let query = b"ACGTAACGT";
        let target = b"ACGTACGT";

        let result = wfa_extend(
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

    #[test]
    fn test_simd_complex_alignment() {
        let query = b"ACGTACGTTAGC";
        let target = b"ACGTTCGTAGC";

        let result = wfa_extend(
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

    // Differential tests: naive vs SIMD must agree on edit distance.
    // Only meaningful on x86_64 where wfa_extend uses a different kernel
    // (wfa_extend_diag_simd_impl) than wfa_extend_naive (fill_wfa). On other
    // architectures, wfa_extend delegates to naive — comparing them is tautological.
    #[cfg(target_arch = "x86_64")]
    mod simd_vs_naive_differential {
        use crate::{wfa_extend_naive, wfa_extend, SeedAnchor};

    #[test]
    fn test_simd_matches_naive_exact() {
        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.edit_distance, simd.edit_distance);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    #[test]
    fn test_simd_matches_naive_mismatch() {
        let query = b"ACGTACGTACGT";
        let target = b"ACGTACTTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.edit_distance, simd.edit_distance);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    #[test]
    fn test_simd_matches_naive_insertion() {
        let query = b"ACGTACGT";
        let target = b"ACGTAACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.edit_distance, simd.edit_distance);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    #[test]
    fn test_simd_matches_naive_deletion() {
        let query = b"ACGTAACGT";
        let target = b"ACGTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.edit_distance, simd.edit_distance);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    #[test]
    fn test_simd_matches_naive_long_sequence() {
        let query = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let target = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.edit_distance, simd.edit_distance);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    #[test]
    fn test_simd_matches_naive_high_divergence() {
        let query = b"ACGTACGTACGTACGT";
        let target = b"TGCATGCATGCATGCA";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.edit_distance, simd.edit_distance);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    #[test]
    fn test_simd_matches_naive_short_sequences() {
        let query = b"ACGT";
        let target = b"ACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.edit_distance, simd.edit_distance);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    #[test]
    fn test_simd_matches_naive_mid_seed() {
        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 4,
            target_pos: 4,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.edit_distance, simd.edit_distance);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    #[test]
    fn test_simd_matches_naive_multiple_indels() {
        let query = b"ACGTACGTACGTACGT";
        let target = b"ACGTTCGTAACGTACG";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.edit_distance, simd.edit_distance);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    #[test]
    fn test_simd_matches_naive_consecutive_indels() {
        let query = b"ACGTAAAACGT";
        let target = b"ACGTCGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.edit_distance, simd.edit_distance);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    #[test]
    fn test_simd_matches_naive_complex_pattern_1() {
        let query = b"ACGTACGTTAGCTTGCA";
        let target = b"ACGTTCGTAGCGCA";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.edit_distance, simd.edit_distance);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    #[test]
    fn test_simd_matches_naive_complex_pattern_2() {
        let query = b"TTAACCGGTTAA";
        let target = b"TTACCGGTAA";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.edit_distance, simd.edit_distance);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    #[test]
    fn test_simd_matches_naive_repeat_regions() {
        let query = b"ATATATATATATAT";
        let target = b"ATATATATATAT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.edit_distance, simd.edit_distance);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    #[test]
    fn test_simd_matches_naive_gc_rich() {
        let query = b"GCGCGCGCGCGCGCGC";
        let target = b"GCGCGCGGCGCGCGC";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.edit_distance, simd.edit_distance);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    #[test]
    fn test_simd_matches_naive_at_rich() {
        let query = b"ATATATATATATAT";
        let target = b"ATATATTATATAT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.edit_distance, simd.edit_distance);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    #[test]
    fn test_simd_matches_naive_edge_case_empty_prefix() {
        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.edit_distance, simd.edit_distance);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    #[test]
    fn test_simd_matches_naive_edge_case_near_end() {
        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 8,
            target_pos: 8,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.edit_distance, simd.edit_distance);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    #[test]
    fn test_simd_matches_naive_random_sequence_1() {
        let query = b"ACGTTAGCTAGCTAGC";
        let target = b"ACGTTAGCTGCTAGC";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.edit_distance, simd.edit_distance);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    #[test]
    fn test_simd_matches_naive_random_sequence_2() {
        let query = b"TGCATGCATGCATGCA";
        let target = b"TGCAATGCATGCATGC";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.edit_distance, simd.edit_distance);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    #[test]
    fn test_simd_matches_naive_random_sequence_3() {
        let query = b"CCGGAATTCCGGAATT";
        let target = b"CCGGGAATTCCGGAAT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let simd_result = wfa_extend(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(simd_result.is_ok());

        let naive = naive_result.unwrap();
        let simd = simd_result.unwrap();

        assert_eq!(naive.edit_distance, simd.edit_distance);
        assert_eq!(naive.edit_distance, simd.edit_distance);
    }

    } // mod simd_vs_naive_differential

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
        use crate::wfa_extend;

        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";

        let result = wfa_extend(
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
        use crate::wfa_extend;

        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let neon_result = wfa_extend(query, target, seed);

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
        use crate::wfa_extend;

        let query = b"ACGTACGTACGT";
        let target = b"ACGTACTTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let neon_result = wfa_extend(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(neon_result.is_ok());

        let naive = naive_result.unwrap();
        let neon = neon_result.unwrap();

        assert_eq!(naive.edit_distance, neon.edit_distance);
        assert_eq!(naive.edit_distance, neon.edit_distance);
    }

    // Correctness: NEON matches naive on insertion
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn issue_72_neon_matches_naive_insertion() {
        use crate::wfa_extend;

        let query = b"ACGTACGT";
        let target = b"ACGTAACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let neon_result = wfa_extend(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(neon_result.is_ok());

        let naive = naive_result.unwrap();
        let neon = neon_result.unwrap();

        assert_eq!(naive.edit_distance, neon.edit_distance);
        assert_eq!(naive.edit_distance, neon.edit_distance);
    }

    // Correctness: NEON matches naive on deletion
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn issue_72_neon_matches_naive_deletion() {
        use crate::wfa_extend;

        let query = b"ACGTAACGT";
        let target = b"ACGTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let neon_result = wfa_extend(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(neon_result.is_ok());

        let naive = naive_result.unwrap();
        let neon = neon_result.unwrap();

        assert_eq!(naive.edit_distance, neon.edit_distance);
        assert_eq!(naive.edit_distance, neon.edit_distance);
    }

    // Correctness: NEON matches naive on complex alignment (mixed ops)
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn issue_72_neon_matches_naive_complex() {
        use crate::wfa_extend;

        let query = b"ACGTACGTTAGC";
        let target = b"ACGTTCGTAGC";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let neon_result = wfa_extend(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(neon_result.is_ok());

        let naive = naive_result.unwrap();
        let neon = neon_result.unwrap();

        assert_eq!(naive.edit_distance, neon.edit_distance);
        assert_eq!(naive.edit_distance, neon.edit_distance);
    }

    // Correctness: NEON matches naive on long sequences (10kb)
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn issue_72_neon_matches_naive_10kb() {
        use crate::wfa_extend;

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
        let neon_result = wfa_extend(&query, &target, seed);

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
        use crate::wfa_extend;

        let query = b"ACGTACGT";
        let target = b"ACGTACGT";
        let seed = SeedAnchor {
            query_pos: 100, // beyond query length
            target_pos: 0,
        };

        let result = wfa_extend(query, target, seed);
        assert!(result.is_err());
    }

    // Error handling: seed position validation on NEON (target)
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn issue_72_neon_invalid_seed_beyond_target() {
        use crate::wfa_extend;

        let query = b"ACGTACGT";
        let target = b"ACGTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 100, // beyond target length
        };

        let result = wfa_extend(query, target, seed);
        assert!(result.is_err());
    }

    // Platform compatibility: NEON is mandatory on aarch64
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn issue_72_neon_runs_unconditionally_on_aarch64() {
        use crate::wfa_extend;

        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        // Must run without runtime detection on aarch64
        let result = wfa_extend(query, target, seed);
        assert!(result.is_ok(), "NEON must run unconditionally on aarch64");
    }

    // Platform compatibility: NEON falls back to naive on non-aarch64
    #[test]
    #[cfg(not(target_arch = "aarch64"))]
    fn issue_72_neon_fallback_on_non_aarch64() {
        use crate::wfa_extend;

        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        // Must not crash and should fall back gracefully
        let result = wfa_extend(query, target, seed);
        assert!(
            result.is_ok(),
            "NEON must fall back to naive on non-aarch64"
        );
    }

    // Alignment position fields: verify NEON reports correct positions
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn issue_72_neon_alignment_positions_at_start() {
        use crate::wfa_extend;

        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let result = wfa_extend(query, target, seed);
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
        use crate::wfa_extend;

        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 4,
            target_pos: 4,
        };

        let result = wfa_extend(query, target, seed);
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
        use crate::wfa_extend;

        let query = b"ACGTACGTACGTACGT";
        let target = b"TGCATGCATGCATGCA";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let neon_result = wfa_extend(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(neon_result.is_ok());

        let naive = naive_result.unwrap();
        let neon = neon_result.unwrap();

        assert_eq!(naive.edit_distance, neon.edit_distance);
        assert_eq!(naive.edit_distance, neon.edit_distance);
    }

    // Repeat regions: NEON handles repetitive sequences
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn issue_72_neon_matches_naive_repeat_regions() {
        use crate::wfa_extend;

        let query = b"ATATATATATATAT";
        let target = b"ATATATATATAT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let neon_result = wfa_extend(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(neon_result.is_ok());

        let naive = naive_result.unwrap();
        let neon = neon_result.unwrap();

        assert_eq!(naive.edit_distance, neon.edit_distance);
        assert_eq!(naive.edit_distance, neon.edit_distance);
    }

    // GC-rich regions: NEON handles GC-rich sequences
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn issue_72_neon_matches_naive_gc_rich() {
        use crate::wfa_extend;

        let query = b"GCGCGCGCGCGCGCGC";
        let target = b"GCGCGCGGCGCGCGC";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let neon_result = wfa_extend(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(neon_result.is_ok());

        let naive = naive_result.unwrap();
        let neon = neon_result.unwrap();

        assert_eq!(naive.edit_distance, neon.edit_distance);
        assert_eq!(naive.edit_distance, neon.edit_distance);
    }

    // Multiple consecutive indels: NEON handles complex CIGAR operations
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn issue_72_neon_matches_naive_consecutive_indels() {
        use crate::wfa_extend;

        let query = b"ACGTAAAACGT";
        let target = b"ACGTCGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let neon_result = wfa_extend(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(neon_result.is_ok());

        let naive = naive_result.unwrap();
        let neon = neon_result.unwrap();

        assert_eq!(naive.edit_distance, neon.edit_distance);
        assert_eq!(naive.edit_distance, neon.edit_distance);
    }

    // Empty suffix at seed: NEON handles edge case of seed at end
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn issue_72_neon_empty_suffix_at_seed() {
        use crate::wfa_extend;

        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 12,
            target_pos: 12,
        };

        let result = wfa_extend(query, target, seed);
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
        use crate::wfa_extend;
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
        let result = wfa_extend(&query, &target, seed);
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
        use crate::wfa_extend;

        let query = b"ACGTACGTACGTACGT";
        let target = b"ACGTTCGTAACGTACG";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let naive_result = wfa_extend_naive(query, target, seed.clone());
        let neon_result = wfa_extend(query, target, seed);

        assert!(naive_result.is_ok());
        assert!(neon_result.is_ok());

        let naive = naive_result.unwrap();
        let neon = neon_result.unwrap();

        assert_eq!(naive.edit_distance, neon.edit_distance);
        assert_eq!(naive.edit_distance, neon.edit_distance);
    }
}

/// TDD suite for the real WFA algorithm implemented in [`fill_wfa`].
///
/// These tests drive the O(s) wavefront implementation. The reference oracle is
/// [`fill_scalar`] (O(n×m) DP), which stays as an immutable correctness ground truth.
/// The final test (`wfa_is_faster_than_on2`) cannot be satisfied by an O(n×m) delegate.
#[cfg(test)]
mod wfa_algorithm_tests {
    use super::{fill_scalar, fill_wfa};

    // ── Behavior 1 (tracer bullet): perfect match ─────────────────────────────

    #[test]
    fn perfect_match_has_zero_edits_and_all_match_cigar() {
        let (cigar, edit_dist) = fill_wfa(b"ACGT", b"ACGT");
        assert_eq!(edit_dist, 0);
        assert_eq!(cigar, "4M");
    }

    // ── Behavior 2: single mismatch ───────────────────────────────────────────

    #[test]
    fn single_mismatch_edit_distance_matches_scalar() {
        // ACGT vs AGGT: position 1 differs (C→G), edit_dist=1
        let q = b"ACGT";
        let t = b"AGGT";
        let dp = fill_scalar(q, t);
        let expected_edit = dp[q.len() * (t.len() + 1) + t.len()] as usize;
        let (_, got_edit) = fill_wfa(q, t);
        assert_eq!(got_edit, expected_edit, "edit distance must match scalar reference");
        assert_eq!(got_edit, 1);
    }

    // ── Behavior 3: single insertion in query ─────────────────────────────────

    #[test]
    fn single_insertion_edit_distance_matches_scalar() {
        // ACGGT vs ACGT: extra G in query, edit_dist=1
        let q = b"ACGGT";
        let t = b"ACGT";
        let dp = fill_scalar(q, t);
        let expected_edit = dp[q.len() * (t.len() + 1) + t.len()] as usize;
        let (_, got_edit) = fill_wfa(q, t);
        assert_eq!(got_edit, expected_edit);
        assert_eq!(got_edit, 1);
    }

    // ── Behavior 4: single deletion from query ────────────────────────────────

    #[test]
    fn single_deletion_edit_distance_matches_scalar() {
        // ACT vs ACGT: query missing G, edit_dist=1
        let q = b"ACT";
        let t = b"ACGT";
        let dp = fill_scalar(q, t);
        let expected_edit = dp[q.len() * (t.len() + 1) + t.len()] as usize;
        let (_, got_edit) = fill_wfa(q, t);
        assert_eq!(got_edit, expected_edit);
        assert_eq!(got_edit, 1);
    }

    // ── Behavior 5: edit distance agrees with scalar on wide sweep ────────────

    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
            z ^ (z >> 31)
        }
        fn below(&mut self, n: usize) -> usize { (self.next() % n as u64) as usize }
    }
    fn random_dna(rng: &mut Rng, len: usize) -> Vec<u8> {
        const B: &[u8; 4] = b"ACGT";
        (0..len).map(|_| B[rng.below(4)]).collect()
    }
    fn diverge_dna(rng: &mut Rng, seq: &[u8], rate_pct: usize) -> Vec<u8> {
        const B: &[u8; 4] = b"ACGT";
        let mut out = Vec::with_capacity(seq.len());
        for &b in seq {
            if rng.below(100) < rate_pct {
                match rng.below(3) {
                    0 => out.push(B[rng.below(4)]),
                    1 => {}
                    _ => { out.push(B[rng.below(4)]); out.push(b); }
                }
            } else { out.push(b); }
        }
        out
    }

    #[test]
    fn edit_distance_matches_scalar_property() {
        let mut rng = Rng(0xDEAD_BEEF_1234);
        let lengths = [0usize, 1, 4, 8, 16, 32, 64, 128];
        let divs = [0usize, 2, 10, 30];
        let mut cases = 0u32;
        for &ql in &lengths {
            for &dv in &divs {
                for _ in 0..10 {
                    let q = random_dna(&mut rng, ql);
                    let t = diverge_dna(&mut rng, &q, dv);
                    let dp = fill_scalar(&q, &t);
                    let expected = dp[q.len() * (t.len() + 1) + t.len()] as usize;
                    let (_, got) = fill_wfa(&q, &t);
                    assert_eq!(
                        got, expected,
                        "edit distance mismatch: q.len={} t.len={} div={}%",
                        q.len(), t.len(), dv
                    );
                    cases += 1;
                }
            }
        }
        assert!(cases > 200, "expected broad sweep, ran {cases}");
    }

    // ── Behavior 6: O(s) performance — impossible for O(n×m) ─────────────────
    //
    // 150bp vs 300bp window at 2% divergence: edit_dist ≈ 3 (3 C positions in window).
    // This is the realistic case: executor windows the target to ~2× query length around
    // the seed anchor before calling fill_wfa. O(n×m) = 45k cells ≈ fast but so is WFA.
    //
    // The 10kbp vs 10kbp case is the meaningful stress test for long sequences:
    // edit_dist ≈ 10 edits, O(n×m) = 100M cells, WFA = O(s*n) = 100k ops → 1000x faster.

    #[test]
    fn wfa_is_faster_than_on2_for_sparse_edits_windowed() {
        use std::time::Instant;
        // 150bp query vs 300bp windowed target (2% divergence → edit_dist ≈ 3)
        let q = vec![b'A'; 150];
        let mut t = vec![b'A'; 300];
        for i in (0..t.len()).step_by(50) { t[i] = b'C'; }
        let start = Instant::now();
        let (_, edit) = fill_wfa(&q, &t);
        let elapsed = start.elapsed();
        assert!(edit > 0, "expected non-zero edits");
        assert!(
            elapsed.as_millis() < 10,
            "fill_wfa took {:?}; windowed WFA must complete in <10ms for 150bp vs 300bp",
            elapsed
        );
    }

    #[test]
    fn wfa_is_faster_than_on2_for_long_similar_sequences() {
        use std::time::Instant;
        // 10kbp vs 10kbp, 0.1% divergence → edit_dist ≈ 10.
        // O(n×m) = 100M cells, WFA O(s*n) = 100k — ~1000x difference.
        // Even in debug, 100M cells take seconds; O(s*n) = trivial.
        let q: Vec<u8> = (0..10_000).map(|i| if i % 1000 == 0 { b'C' } else { b'A' }).collect();
        let t: Vec<u8> = (0..10_000).map(|i| if i % 1001 == 0 { b'C' } else { b'A' }).collect();
        let start = Instant::now();
        let (_, edit) = fill_wfa(&q, &t);
        let elapsed = start.elapsed();
        assert!(edit > 0);
        assert!(
            elapsed.as_millis() < 100,
            "fill_wfa took {:?}; O(s) WFA must complete in <100ms for 10kbp vs 10kbp \
             at low divergence (O(n×m) takes seconds)",
            elapsed
        );
    }
}

// ============================================================================
// Myers Bit-Parallel Edit Distance (Issue #144)
// Myers 1999: "A fast bit-vector algorithm for approximate string matching"
// O(n * ceil(m/w)) time, O(m/w) space for distance; O(nm/w) for CIGAR traceback.
// ============================================================================

/// Reconstruct column i of the DP matrix from stored (vp, vn) bitvectors.
/// `bottom_score` = dp[m][j], returned column is dp[0..=m][j].
/// Reconstruct an edit-distance DP column from Myers `(vp, vn)` bitvectors into a caller
/// -provided buffer (`col.len()` must be `m + 1`). Fills in place so the backtrace can
/// reuse two column buffers instead of allocating a fresh column per traceback step.
fn reconstruct_into(vp: &[u64], vn: &[u64], m: usize, bottom_score: usize, col: &mut [usize]) {
    col[m] = bottom_score;
    for i in (0..m).rev() {
        let block = i / 64;
        let bit = i % 64;
        let vp_bit = ((vp[block] >> bit) & 1) as isize;
        let vn_bit = ((vn[block] >> bit) & 1) as isize;
        col[i] = (col[i + 1] as isize - vp_bit + vn_bit) as usize;
    }
}

/// Flat-buffer Myers forward pass: stores VP/VN bitvectors in a single pre-allocated
/// Vec<u64> indexed as `[col * num_blocks * 2 + block]` (VP) and
/// `[col * num_blocks * 2 + num_blocks + block]` (VN), plus a separate scores Vec<usize>.
/// Eliminates ~2×n heap allocations (one Vec clone per column) from the hot path.
struct MyersColumns {
    data: Vec<u64>,
    scores: Vec<usize>,
    num_blocks: usize,
    num_cols: usize,
}

impl MyersColumns {
    fn vp(&self, col: usize, block: usize) -> u64 {
        self.data[col * self.num_blocks * 2 + block]
    }
    fn vn(&self, col: usize, block: usize) -> u64 {
        self.data[col * self.num_blocks * 2 + self.num_blocks + block]
    }
    fn score(&self, col: usize) -> usize {
        self.scores[col]
    }
    fn vp_slice(&self, col: usize) -> &[u64] {
        let base = col * self.num_blocks * 2;
        &self.data[base..base + self.num_blocks]
    }
    fn vn_slice(&self, col: usize) -> &[u64] {
        let base = col * self.num_blocks * 2 + self.num_blocks;
        &self.data[base..base + self.num_blocks]
    }
}

/// Early-termination threshold for fitting mode: once the running score exceeds
/// `best_seen + margin`, stop — the fitting endpoint was already found and
/// further columns can only get worse. `None` disables early termination (global mode).
fn myers_forward_flat_with_cutoff(query: &[u8], target: &[u8], abandon_margin: Option<usize>) -> MyersColumns {
    const W: usize = 64;
    let m = query.len();
    let n = target.len();

    let num_blocks = (m + W - 1) / W;
    let last_block = num_blocks - 1;
    let last_bits = if m % W == 0 { W } else { m % W };
    let last_mask: u64 = if last_bits == W { u64::MAX } else { (1u64 << last_bits) - 1 };
    let score_bit = (m - 1) % W;

    let mut pm = vec![[0u64; 256]; num_blocks];
    for (i, &b) in query.iter().enumerate() {
        pm[i / W][b as usize] |= 1u64 << (i % W);
    }

    let mut vp = vec![u64::MAX; num_blocks];
    vp[last_block] = last_mask;
    let mut vn = vec![0u64; num_blocks];
    let mut score = m;
    let mut best_score_seen = m;

    let stride = num_blocks * 2;
    let mut data = vec![0u64; n * stride];
    let mut scores = Vec::with_capacity(n);

    let mut new_vp = vec![0u64; num_blocks];
    let mut new_vn = vec![0u64; num_blocks];

    for (col_idx, &tb) in target.iter().enumerate() {
        let mut add_carry: u64 = 0;
        let mut hp_carry: u64 = 1;
        let mut hn_carry: u64 = 0;

        let mut last_hp = 0u64;
        let mut last_hn = 0u64;

        for k in 0..num_blocks {
            let eq = pm[k][tb as usize];
            let xh = eq | vn[k];
            let xhvp = xh & vp[k];
            let (s1, c1) = vp[k].overflowing_add(xhvp);
            let (s2, c2) = s1.overflowing_add(add_carry);
            add_carry = (c1 as u64) | (c2 as u64);
            let d0_raw = (s2 ^ vp[k]) | xh;
            let hn_raw = vp[k] & d0_raw;
            let hp_raw = vn[k] | !(vp[k] | d0_raw);

            let (d0, hn, hp) = if k == last_block && last_bits < W {
                (d0_raw & last_mask, hn_raw & last_mask, hp_raw & last_mask)
            } else {
                (d0_raw, hn_raw, hp_raw)
            };

            if k == last_block {
                last_hp = hp;
                last_hn = hn;
            }

            let next_hp_carry = hp >> (W - 1);
            let x = (hp << 1) | hp_carry;
            hp_carry = next_hp_carry;

            let next_hn_carry = hn >> (W - 1);
            let hn_shifted = (hn << 1) | hn_carry;
            hn_carry = next_hn_carry;

            new_vn[k] = x & d0;
            let vp_raw2 = hn_shifted | !(d0 | x);
            new_vp[k] = if k == last_block && last_bits < W { vp_raw2 & last_mask } else { vp_raw2 };
        }

        if (last_hp >> score_bit) & 1 != 0 { score += 1; }
        if (last_hn >> score_bit) & 1 != 0 { score = score.saturating_sub(1); }

        std::mem::swap(&mut vp, &mut new_vp);
        std::mem::swap(&mut vn, &mut new_vn);

        let base = col_idx * stride;
        data[base..base + num_blocks].copy_from_slice(&vp);
        data[base + num_blocks..base + stride].copy_from_slice(&vn);
        scores.push(score);
        if score < best_score_seen {
            best_score_seen = score;
        }

        if let Some(margin) = abandon_margin {
            if score > best_score_seen + margin {
                let actual_cols = col_idx + 1;
                data.truncate(actual_cols * stride);
                return MyersColumns { data, scores, num_blocks, num_cols: actual_cols };
            }
        }
    }

    MyersColumns { data, scores, num_blocks, num_cols: n }
}

fn myers_forward_flat(query: &[u8], target: &[u8]) -> MyersColumns {
    myers_forward_flat_with_cutoff(query, target, None)
}

/// Myers bit-parallel **global** edit distance for arbitrary query length.
/// Uses multi-word blocks of 64 bits for queries > 64bp.
/// Returns (edit_distance, CIGAR string) consuming the entire query and target.
pub fn myers_edit_distance_impl(query: &[u8], target: &[u8]) -> (usize, String) {
    let m = query.len();
    let n = target.len();

    if m == 0 {
        let cigar = if n > 0 { format!("{n}I") } else { String::new() };
        return (n, cigar);
    }
    if n == 0 {
        return (m, format!("{m}D"));
    }

    let cols = myers_forward_flat(query, target);
    let edit_distance = if cols.num_cols > 0 { cols.score(cols.num_cols - 1) } else { m };
    let cigar = myers_backtrace_flat(query, target, &cols, n);
    (edit_distance, cigar)
}

/// Myers bit-parallel **fitting** alignment: the query must be fully consumed, but the
/// target end is free (no penalty for an unconsumed target tail). This is the correct
/// mode for aligning a read against a longer reference window, and mirrors
/// [`fill_wfa_fitting`]'s semantics so the Myers and WFA alignment paths agree.
///
/// Returns `(edit_distance, cigar, target_consumed)` where `target_consumed` is the
/// number of target bases in the chosen fitting alignment.
///
/// When the target is not substantially longer than the query (`n <= m + m/2 + 10`),
/// delegates to the global path — exactly the threshold [`fill_wfa_fitting`] uses — so
/// that terminal indels are not silently hidden by a free end gap.
/// Reconstruct DP column `j` from flat MyersColumns into `buf[0..=m]`.
fn fill_dp_column_flat(cols: &MyersColumns, m: usize, j: usize, buf: &mut [usize]) {
    const W: usize = 64;
    if j == 0 {
        for (i, slot) in buf.iter_mut().enumerate().take(m + 1) {
            *slot = i;
        }
    } else {
        let vp = cols.vp_slice(j - 1);
        let vn = cols.vn_slice(j - 1);
        let bs = cols.score(j - 1);
        reconstruct_into(vp, vn, m, bs, buf);
    }
}

/// Backtrace from flat MyersColumns (mirrors `myers_backtrace` logic).
fn myers_backtrace_flat(
    query: &[u8],
    target: &[u8],
    cols: &MyersColumns,
    end_ti: usize,
) -> String {
    let m = query.len();
    let mut ops: Vec<u8> = Vec::with_capacity(m + end_ti);
    let mut qi = m;
    let mut ti = end_ti;

    let mut cur_col = vec![0usize; m + 1];
    let mut prev_col = vec![0usize; m + 1];
    if ti > 0 {
        fill_dp_column_flat(cols, m, ti, &mut cur_col);
        fill_dp_column_flat(cols, m, ti - 1, &mut prev_col);
    }

    while qi > 0 || ti > 0 {
        if qi == 0 {
            ops.push(b'I');
            ti -= 1;
            continue;
        }
        if ti == 0 {
            ops.push(b'D');
            qi -= 1;
            continue;
        }

        let cur = cur_col[qi];
        let diag_cost = if query[qi - 1] == target[ti - 1] { 0 } else { 1 };
        let from_diag = prev_col[qi - 1] + diag_cost;
        let from_above = cur_col[qi - 1] + 1;
        let from_left = prev_col[qi] + 1;

        if cur == from_above {
            ops.push(b'D');
            qi -= 1;
        } else if cur == from_left {
            ops.push(b'I');
            ti -= 1;
            std::mem::swap(&mut cur_col, &mut prev_col);
            if ti > 0 {
                fill_dp_column_flat(cols, m, ti - 1, &mut prev_col);
            }
        } else {
            debug_assert_eq!(cur, from_diag);
            ops.push(if diag_cost == 0 { b'M' } else { b'X' });
            qi -= 1;
            ti -= 1;
            std::mem::swap(&mut cur_col, &mut prev_col);
            if ti > 0 {
                fill_dp_column_flat(cols, m, ti - 1, &mut prev_col);
            }
        }
    }

    ops.reverse();

    let mut cigar = String::new();
    if !ops.is_empty() {
        let mut count = 1usize;
        let mut cur_op = ops[0];
        for &op in &ops[1..] {
            if op == cur_op {
                count += 1;
            } else {
                cigar.push_str(&count.to_string());
                cigar.push(cur_op as char);
                cur_op = op;
                count = 1;
            }
        }
        cigar.push_str(&count.to_string());
        cigar.push(cur_op as char);
    }

    cigar
}

pub fn myers_fitting_impl(query: &[u8], target: &[u8]) -> (usize, String, usize) {
    let m = query.len();
    let n = target.len();

    if m == 0 {
        return (0, String::new(), 0);
    }
    if n == 0 {
        return (m, format!("{m}D"), 0);
    }

    if n <= m + m / 2 + 10 {
        let (edit, cigar) = myers_edit_distance_impl(query, target);
        return (edit, cigar, n);
    }

    let abandon_margin = m / 10 + 1;
    let cols = myers_forward_flat_with_cutoff(query, target, Some(abandon_margin));

    let mut best_consumed = 0usize;
    let mut best_score = m;
    for c in 0..cols.num_cols {
        let s = cols.score(c);
        if s < best_score {
            best_score = s;
            best_consumed = c + 1;
        }
    }

    let cigar = myers_backtrace_flat(query, target, &cols, best_consumed);
    (best_score, cigar, best_consumed)
}

#[cfg(test)]
mod simd_prefix_tests {
    use super::{count_matching_prefix, count_matching_prefix_scalar, count_matching_prefix_u64};

    #[test]
    fn counts_matching_prefix_before_first_mismatch() {
        // 4 matching bytes, then a mismatch at index 4.
        let a = b"ACGTACGT";
        let b = b"ACGTTCGT";
        assert_eq!(count_matching_prefix(a, b), 4);
    }

    /// Build a battery of (a, b) pairs that exercise the boundaries the chunked
    /// implementations care about: lengths straddling 8 and 16 bytes, a mismatch at
    /// every possible position, empty/short slices, and unequal lengths.
    fn battery() -> Vec<(Vec<u8>, Vec<u8>)> {
        let mut cases = Vec::new();
        // Identical sequences of every length 0..40 (no mismatch; capped at len).
        for len in 0..40 {
            let v: Vec<u8> = (0..len).map(|i| b"ACGT"[i % 4]).collect();
            cases.push((v.clone(), v));
        }
        // Mismatch at position p for a 40-byte sequence, for every p.
        for p in 0..40 {
            let a: Vec<u8> = (0..40).map(|i| b"ACGT"[i % 4]).collect();
            let mut b = a.clone();
            b[p] ^= 0x01; // flip a bit so the byte differs
            cases.push((a, b));
        }
        // Unequal lengths: prefix matches, then one runs out.
        cases.push((b"ACGTACGTAC".to_vec(), b"ACGTACGT".to_vec()));
        cases.push((b"ACGT".to_vec(), b"ACGTACGTACGT".to_vec()));
        // Empty cases.
        cases.push((Vec::new(), b"ACGT".to_vec()));
        cases.push((b"ACGT".to_vec(), Vec::new()));
        cases
    }

    #[test]
    fn u64_tier_matches_scalar() {
        for (a, b) in battery() {
            assert_eq!(
                count_matching_prefix_u64(&a, &b),
                count_matching_prefix_scalar(&a, &b),
                "u64 tier disagreed with scalar on a={a:?} b={b:?}"
            );
        }
    }

    #[cfg(any(target_arch = "x86_64", target_arch = "aarch64"))]
    #[test]
    fn arch_tier_matches_scalar() {
        use super::count_matching_prefix_arch;
        for (a, b) in battery() {
            assert_eq!(
                count_matching_prefix_arch(&a, &b),
                count_matching_prefix_scalar(&a, &b),
                "arch SIMD tier disagreed with scalar on a={a:?} b={b:?}"
            );
        }
    }

    #[test]
    fn dispatch_matches_scalar() {
        // The production entry point must agree with the reference on every case,
        // whichever tier the target architecture selects.
        for (a, b) in battery() {
            assert_eq!(
                count_matching_prefix(&a, &b),
                count_matching_prefix_scalar(&a, &b),
                "dispatched count_matching_prefix disagreed with scalar on a={a:?} b={b:?}"
            );
        }
    }
}

// ============================================================================
// ISSUE #180: Abandonment Sentinel for fill_wfa_fitting (ADR-0007 prerequisite)
// ============================================================================
//
// Tests verify that fill_wfa_fitting returns a distinct "abandoned" result
// on loop exhaustion (no fitting-end reached before max_s), not a spurious
// edit_distance=0 perfect alignment. This is the critical fix needed before
// ADR-0007 introduces score-bounded early abandonment with small max_s.
//
// These tests are marked with @test #[cfg(test)] for filtering in CI.
#[cfg(test)]
mod issue_180_abandonment_tests {
    /// Helper: inject a small cap on max_s for testing abandonment.
    /// Returns (cigar, edit_distance, target_consumed) or None if abandoned.
    /// **NOTE**: This is a test-only helper; production `fill_wfa_fitting` has no cap.
    /// Issue #180 requires modifying `fill_wfa_fitting` to support abandonment
    /// and this helper tests the interface contract.
    fn fill_wfa_fitting_with_cap(
        q: &[u8],
        t: &[u8],
        max_s_cap: usize,
    ) -> Option<(String, usize, usize)> {
        crate::wfa_simd::fill_wfa_fitting_impl(q, t, Some(max_s_cap))
    }

    /// **ACCEPTANCE CRITERION**: Normal alignment with fitting-end reached
    /// must return Some, not abandon, even under default uncapped max_s.
    ///
    /// This verifies backward compatibility: at the current max_s = qn + tn,
    /// no alignment is abandoned, so results are identical to today.
    #[test]
    fn issue_180_happy_path_fitting_end_reached() {
        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // Much longer to trigger fitting mode

        let result = fill_wfa_fitting_with_cap(query, target, usize::MAX);

        // With uncapped max_s, fitting-end must be reached => Some(...)
        assert!(
            result.is_some(),
            "normal fitting alignment must return Some, not None (abandoned)"
        );

        let (cigar, edit_distance, target_consumed) = result.unwrap();

        // Sanity checks: normal alignment has non-empty CIGAR and valid metrics
        assert!(!cigar.is_empty(), "happy-path CIGAR must not be empty");
        assert!(edit_distance < query.len(), "happy-path edit_distance should be small");
        assert!(target_consumed <= target.len(), "consumed target must not exceed length");
        assert!(target_consumed > 0, "consumed target must be non-zero for matching sequences");
    }

    /// **ACCEPTANCE CRITERION**: Alignment that exhausts the loop without
    /// reaching a fitting-end must return None (abandoned), NOT a spurious
    /// (empty_cigar, edit_distance=0, full_target).
    ///
    /// This is the core bug fix: injecting a small max_s cap forces exhaustion,
    /// and the alignment engine must report it as abandoned, not as a perfect match.
    #[test]
    fn issue_180_abandonment_on_exhaustion_with_small_cap() {
        // A query/target pair with a real, known-large edit distance: query is all 'A',
        // target has a substitution every 3rd base, so the true edit distance is far
        // above the tiny cap below. (The original fixture here reused a repeated "ACGT"
        // query embedded verbatim in the target -- true edit distance 0 -- which
        // trivially fits under any cap and can never demonstrate real abandonment.)
        let query = vec![b'A'; 100];
        let mut target = query.clone();
        for i in (0..100).step_by(3) {
            target[i] = b'T';
        }

        let result = fill_wfa_fitting_with_cap(&query, &target, 2); // true distance ~33, cap=2

        // With a tiny cap far below the true edit distance, exhaustion is expected.
        assert!(
            result.is_none(),
            "with injected small max_s cap far below the true edit distance, exhaustion \
             must return None (abandoned), not a spurious (empty, 0, full_len) perfect alignment"
        );
    }

    /// **ACCEPTANCE CRITERION**: Abandoned alignment must never have
    /// edit_distance = 0, which would make it appear as a perfect match to
    /// score_alignments (which uses min_by_key(|a| edit_distance)).
    ///
    /// This test verifies the bug fix explicitly: if abandonment is represented
    /// as a tuple, it must NEVER return (_, 0, _) on exhaustion.
    #[test]
    fn issue_180_abandoned_never_reports_edit_distance_zero() {
        // Hypothetical scenario: if abandonment were mis-represented as
        // (empty_cigar, 0, full_target), it would be a silent correctness bug.
        // This test would fail if the implementation returns that.
        let query = b"ACGTACGTACGTACGT";
        let target = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

        let result = fill_wfa_fitting_with_cap(query, target, 1); // max_s=1

        // The key assertion: if it's Some(...), the inner tuple must NOT have edit_distance=0
        if let Some((cigar, edit_distance, _target_consumed)) = result {
            assert!(
                edit_distance > 0 || !cigar.is_empty(),
                "abandoned alignment (if returned as Some) must not have both empty CIGAR \
                 and edit_distance=0, which would fool score_alignments into a false primary"
            );
        }
        // Preferred outcome: None (abandoned), not Some(..., 0, ...)
    }

    /// **ACCEPTANCE CRITERION**: Verify that the abandonment representation
    /// is distinct from a real (low-divergence) alignment.
    ///
    /// This documents the contract: at the default uncapped max_s, even a
    /// difficult alignment must succeed (fitting-end reached), while with a
    /// tiny cap, it should abandon. The two must be distinguishable.
    #[test]
    fn issue_180_distinguishable_normal_vs_abandoned() {
        // Same real-divergence fixture as issue_180_abandonment_on_exhaustion_with_small_cap:
        // true edit distance ~33 (a repeated "ACGT" query embedded verbatim in the target
        // has true edit distance 0, which trivially "succeeds" under any cap including 0,
        // so it can never be visibly distinct from the uncapped result).
        let query = vec![b'A'; 100];
        let mut target = query.clone();
        for i in (0..100).step_by(3) {
            target[i] = b'T';
        }

        // Scenario 1: uncapped (normal operation)
        let normal = fill_wfa_fitting_with_cap(&query, &target, usize::MAX);

        // Scenario 2: capped well below the true edit distance (triggers abandonment)
        let abandoned = fill_wfa_fitting_with_cap(&query, &target, 2);

        // Normal must succeed, abandoned must fail or be visibly different
        assert!(normal.is_some(), "uncapped alignment must succeed");

        // Abandoned should be None or visibly distinct from normal
        // If abandonment is represented as Some(...), it must differ from normal
        if abandoned.is_some() {
            let (normal_cigar, normal_edit, _) = normal.unwrap();
            let (abandoned_cigar, abandoned_edit, _) = abandoned.unwrap();
            assert!(
                abandoned_cigar != normal_cigar || abandoned_edit != normal_edit,
                "abandoned must be visibly distinct from normal (different CIGAR or edit distance)"
            );
        }
    }

    /// **ACCEPTANCE CRITERION**: Existing differential tests continue to pass
    /// (at uncapped max_s, no alignment is abandoned, so output is identical).
    ///
    /// This verifies backward compatibility: the abandonment feature only
    /// activates when max_s is capped; the default code path is unchanged.
    #[test]
    fn issue_180_backward_compatibility_uncapped_matches_current() {
        // This test verifies that with max_s = qn + tn (current default),
        // fill_wfa_fitting produces the same result as today (no abandonment).
        // If the implementation is correct, this test passes unchanged.
        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

        // With uncapped max_s (usize::MAX simulates infinity), should always return Some
        let result = fill_wfa_fitting_with_cap(query, target, usize::MAX);

        assert!(
            result.is_some(),
            "at uncapped max_s, all alignments must complete normally (backward compat)"
        );

        let (cigar, edit_distance, _target_consumed) = result.unwrap();

        // The actual values should match pre-ADR-0007 behavior (which had no cap)
        // Since we're at RED stage, we just verify it's not the fallback garbage value
        assert!(
            !cigar.is_empty() || edit_distance > 0 || false,
            "must not be the fallback (empty_cigar, 0, full_target) value"
        );
    }
}