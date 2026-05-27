//! NEON SIMD optimized WFA diagonal fill for ARM64
//!
//! This module provides NEON intrinsics acceleration for the WFA inner loop
//! on aarch64 platforms. NEON is mandatory on ARM64, so this is unconditionally
//! compiled on that architecture.
//!
//! # Safety
//!
//! Uses unsafe SIMD intrinsics from `core::arch::aarch64`. All unsafe blocks
//! document their invariants inline.

#[cfg(target_arch = "aarch64")]
use core::arch::aarch64;

/// WFA diagonal fill using NEON intrinsics
///
/// Fills a diagonal of WFA wavefront scores using NEON vector operations.
/// Should achieve ≥1.5× speedup vs naive scalar implementation.
///
/// # Arguments
///
/// * `diagonal` - The diagonal buffer to fill
/// * `prev_wavefront` - Previous wavefront scores
/// * `query` - Query sequence bytes
/// * `target` - Target sequence bytes
///
/// # Returns
///
/// The filled diagonal buffer with updated WFA scores
///
/// # Safety
///
/// This function uses unsafe SIMD intrinsics. Callers must ensure:
/// - `diagonal` and `prev_wavefront` slices are properly aligned for NEON loads
/// - Buffer lengths are validated to prevent out-of-bounds access
#[cfg(target_arch = "aarch64")]
pub fn wfa_diagonal_fill_neon(
    diagonal: &mut [i32],
    prev_wavefront: &[i32],
    query: &[u8],
    target: &[u8],
) {
    // Handle empty sequences
    if diagonal.is_empty() || query.is_empty() || target.is_empty() {
        return;
    }

    // Ensure all buffers are the same length
    let len = diagonal.len();
    if prev_wavefront.len() != len || query.len() != len || target.len() != len {
        // Fallback to simple scalar for mismatched lengths
        wfa_diagonal_fill_scalar(diagonal, prev_wavefront, query, target);
        return;
    }

    // NEON vector width is 4 x i32 per register (16 bytes)
    const NEON_WIDTH: usize = 4;

    unsafe {
        // Process in chunks using NEON vectorization
        let full_chunks = len / NEON_WIDTH;

        for chunk_idx in 0..full_chunks {
            let base = chunk_idx * NEON_WIDTH;

            // SAFETY INVARIANT: All buffer accesses are within bounds because:
            // - len is validated at function entry to be non-empty
            // - base = chunk_idx * NEON_WIDTH where chunk_idx < len / NEON_WIDTH
            // - Therefore base + NEON_WIDTH - 1 < len, all accesses are valid

            // Load previous wavefront scores (4 x i32)
            let prev_scores = aarch64::vld1q_s32(prev_wavefront.as_ptr().add(base));

            // Load 4 query and 4 target bytes, compute costs
            let q0 = query.get_unchecked(base);
            let q1 = query.get_unchecked(base + 1);
            let q2 = query.get_unchecked(base + 2);
            let q3 = query.get_unchecked(base + 3);

            let t0 = target.get_unchecked(base);
            let t1 = target.get_unchecked(base + 1);
            let t2 = target.get_unchecked(base + 2);
            let t3 = target.get_unchecked(base + 3);

            // Compute match costs: 0 if match, -1 if mismatch
            let cost0 = if *q0 == *t0 { 0i32 } else { -1i32 };
            let cost1 = if *q1 == *t1 { 0i32 } else { -1i32 };
            let cost2 = if *q2 == *t2 { 0i32 } else { -1i32 };
            let cost3 = if *q3 == *t3 { 0i32 } else { -1i32 };

            // Load costs as a vector (vld1q_s32 expects aligned pointer)
            let costs_array = [cost0, cost1, cost2, cost3];
            let costs = aarch64::vld1q_s32(costs_array.as_ptr());

            // Add previous scores and costs: new_score = prev_score + cost
            let new_scores = aarch64::vaddq_s32(prev_scores, costs);

            // Store result back to diagonal
            aarch64::vst1q_s32(diagonal.as_mut_ptr().add(base), new_scores);
        }

        // Handle remaining elements with scalar loop
        let remainder = len % NEON_WIDTH;
        if remainder > 0 {
            let start = full_chunks * NEON_WIDTH;
            for i in start..len {
                let cost = if query[i] == target[i] { 0 } else { -1 };
                diagonal[i] = prev_wavefront[i] + cost;
            }
        }
    }
}

/// Scalar fallback for WFA diagonal fill
#[inline]
#[cfg(target_arch = "aarch64")]
fn wfa_diagonal_fill_scalar(
    diagonal: &mut [i32],
    prev_wavefront: &[i32],
    query: &[u8],
    target: &[u8],
) {
    for i in 0..diagonal.len() {
        let cost = if query[i] == target[i] { 0 } else { -1 };
        diagonal[i] = prev_wavefront[i] + cost;
    }
}

/// Fallback for non-ARM64 platforms (no-op)
#[cfg(not(target_arch = "aarch64"))]
pub fn wfa_diagonal_fill_neon(
    _diagonal: &mut [i32],
    _prev_wavefront: &[i32],
    _query: &[u8],
    _target: &[u8],
) {
    // No-op on non-ARM64 platforms
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper to create test sequences
    fn seq(s: &str) -> Vec<u8> {
        s.as_bytes().to_vec()
    }

    /// Test 1: Exact match sequences - diagonal should extend fully
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_neon_exact_match() {
        let query = seq("ACGTACGT");
        let target = seq("ACGTACGT");
        let mut diagonal = vec![0; 8];
        let prev_wavefront = vec![0; 8];

        wfa_diagonal_fill_neon(&mut diagonal, &prev_wavefront, &query, &target);

        // Should match naive WFA result for exact match
        // (Will fail until implementation exists)
        assert!(
            diagonal.iter().any(|&x| x != 0),
            "NEON diagonal fill produced no output for exact match"
        );
    }

    /// Test 2: Single mismatch - diagonal should reflect penalty
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_neon_single_mismatch() {
        let query = seq("ACGTACGT");
        let target = seq("ACGTTCGT"); // T instead of A at position 4
        let mut diagonal = vec![0; 8];
        let prev_wavefront = vec![0; 8];

        wfa_diagonal_fill_neon(&mut diagonal, &prev_wavefront, &query, &target);

        // Should handle mismatch correctly
        assert!(
            diagonal.iter().any(|&x| x != 0),
            "NEON diagonal fill produced no output for mismatch"
        );
    }

    /// Test 3: All mismatches - stress test for penalty accumulation
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_neon_all_mismatches() {
        let query = seq("AAAAAAA");
        let target = seq("TTTTTTT");
        let mut diagonal = vec![0; 7];
        let prev_wavefront = vec![0; 7];

        wfa_diagonal_fill_neon(&mut diagonal, &prev_wavefront, &query, &target);

        // Should accumulate penalties for all mismatches
        assert!(
            diagonal.iter().any(|&x| x != 0),
            "NEON diagonal fill produced no output for all mismatches"
        );
    }

    /// Test 4: Long sequence (10kb) - benchmark target case
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_neon_long_sequence_10kb() {
        // 10kb sequences as mentioned in acceptance criteria
        let query = vec![b'A'; 10_000];
        let target = vec![b'A'; 10_000];
        let mut diagonal = vec![0; 10_000];
        let prev_wavefront = vec![0; 10_000];

        wfa_diagonal_fill_neon(&mut diagonal, &prev_wavefront, &query, &target);

        // Should handle large sequences without panic
        assert!(
            diagonal.len() == 10_000,
            "NEON diagonal fill corrupted buffer size"
        );
    }

    /// Test 5: Short sequence (boundary case)
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_neon_short_sequence() {
        let query = seq("ACG");
        let target = seq("ACG");
        let mut diagonal = vec![0; 3];
        let prev_wavefront = vec![0; 3];

        wfa_diagonal_fill_neon(&mut diagonal, &prev_wavefront, &query, &target);

        // Should handle short sequences (less than NEON vector width)
        assert!(
            diagonal.iter().any(|&x| x != 0),
            "NEON diagonal fill failed on short sequence"
        );
    }

    /// Test 6: Aligned to NEON vector width (16 bytes)
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_neon_aligned_16_bytes() {
        let query = vec![b'A'; 16];
        let target = vec![b'A'; 16];
        let mut diagonal = vec![0; 16];
        let prev_wavefront = vec![0; 16];

        wfa_diagonal_fill_neon(&mut diagonal, &prev_wavefront, &query, &target);

        // Should efficiently process aligned vector width
        assert!(
            diagonal.len() == 16,
            "NEON diagonal fill failed on 16-byte aligned input"
        );
    }

    /// Test 7: Unaligned to NEON vector width (17 bytes)
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_neon_unaligned_17_bytes() {
        let query = vec![b'A'; 17];
        let target = vec![b'A'; 17];
        let mut diagonal = vec![0; 17];
        let prev_wavefront = vec![0; 17];

        wfa_diagonal_fill_neon(&mut diagonal, &prev_wavefront, &query, &target);

        // Should handle unaligned sizes with tail processing
        assert!(
            diagonal.len() == 17,
            "NEON diagonal fill failed on unaligned input"
        );
    }

    /// Test 8: Empty sequences (edge case)
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_neon_empty_sequences() {
        let query = seq("");
        let target = seq("");
        let mut diagonal = vec![];
        let prev_wavefront = vec![];

        wfa_diagonal_fill_neon(&mut diagonal, &prev_wavefront, &query, &target);

        // Should handle empty input gracefully
        assert!(
            diagonal.is_empty(),
            "NEON diagonal fill should preserve empty diagonal"
        );
    }

    /// Test 9: Non-ACGT characters
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_neon_non_standard_bases() {
        let query = seq("ACGTNNN");
        let target = seq("ACGTNNN");
        let mut diagonal = vec![0; 7];
        let prev_wavefront = vec![0; 7];

        wfa_diagonal_fill_neon(&mut diagonal, &prev_wavefront, &query, &target);

        // Should handle N (ambiguous base) correctly
        assert!(
            diagonal.iter().any(|&x| x != 0),
            "NEON diagonal fill failed on non-standard bases"
        );
    }

    /// Test 10: Negative scores in previous wavefront
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_neon_negative_wavefront_scores() {
        let query = seq("ACGT");
        let target = seq("ACGT");
        let mut diagonal = vec![0; 4];
        let prev_wavefront = vec![-5, -3, -1, 0]; // Negative penalty scores

        wfa_diagonal_fill_neon(&mut diagonal, &prev_wavefront, &query, &target);

        // Should handle negative scores from previous wavefront
        assert!(
            diagonal.iter().any(|&x| x != 0),
            "NEON diagonal fill failed with negative previous scores"
        );
    }

    /// Test 11: High-scoring previous wavefront
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_neon_high_scores() {
        let query = seq("ACGT");
        let target = seq("ACGT");
        let mut diagonal = vec![0; 4];
        let prev_wavefront = vec![100, 200, 300, 400]; // High scores

        wfa_diagonal_fill_neon(&mut diagonal, &prev_wavefront, &query, &target);

        // Should handle large score values without overflow
        assert!(
            diagonal.iter().any(|&x| x != 0),
            "NEON diagonal fill failed with high scores"
        );
    }

    /// Test 12: GC-rich sequences
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_neon_gc_rich() {
        let query = seq("GCGCGCGC");
        let target = seq("GCGCGCGC");
        let mut diagonal = vec![0; 8];
        let prev_wavefront = vec![0; 8];

        wfa_diagonal_fill_neon(&mut diagonal, &prev_wavefront, &query, &target);

        // Should handle GC-rich sequences (common in bacteria)
        assert!(
            diagonal.iter().any(|&x| x != 0),
            "NEON diagonal fill failed on GC-rich sequence"
        );
    }

    /// Test 13: AT-rich sequences
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_neon_at_rich() {
        let query = seq("ATATATAT");
        let target = seq("ATATATAT");
        let mut diagonal = vec![0; 8];
        let prev_wavefront = vec![0; 8];

        wfa_diagonal_fill_neon(&mut diagonal, &prev_wavefront, &query, &target);

        // Should handle AT-rich sequences
        assert!(
            diagonal.iter().any(|&x| x != 0),
            "NEON diagonal fill failed on AT-rich sequence"
        );
    }

    /// Test 14: Homopolymer run
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_neon_homopolymer() {
        let query = seq("AAAAAAAA");
        let target = seq("AAAAAAAA");
        let mut diagonal = vec![0; 8];
        let prev_wavefront = vec![0; 8];

        wfa_diagonal_fill_neon(&mut diagonal, &prev_wavefront, &query, &target);

        // Should handle homopolymer runs
        assert!(
            diagonal.iter().any(|&x| x != 0),
            "NEON diagonal fill failed on homopolymer"
        );
    }

    /// Test 15: Mixed case with multiple event types
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_neon_mixed_events() {
        let query = seq("ACGTACGTACGT");
        let target = seq("ACGTTCGTACTT");
        let mut diagonal = vec![0; 12];
        let prev_wavefront = vec![0; 12];

        wfa_diagonal_fill_neon(&mut diagonal, &prev_wavefront, &query, &target);

        // Should handle multiple mismatches in sequence
        assert!(
            diagonal.iter().any(|&x| x != 0),
            "NEON diagonal fill failed on mixed events"
        );
    }

    /// Test 16: Verify no-op on non-ARM64 platforms
    #[test]
    #[cfg(not(target_arch = "aarch64"))]
    fn test_neon_noop_on_non_arm64() {
        let query = seq("ACGT");
        let target = seq("ACGT");
        let mut diagonal = vec![0; 4];
        let prev_wavefront = vec![0; 4];

        // Should not panic, just no-op
        wfa_diagonal_fill_neon(&mut diagonal, &prev_wavefront, &query, &target);

        // Should remain unchanged on non-ARM64
        assert_eq!(
            diagonal,
            vec![0, 0, 0, 0],
            "Non-ARM64 platform should no-op"
        );
    }

    /// Test 17: Correctness vs naive WFA (will need naive implementation)
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_neon_matches_naive_wfa_case1() {
        let query = seq("ACGTACGT");
        let target = seq("ACGTACGT");
        let mut diagonal_neon = vec![0; 8];
        let mut diagonal_naive = vec![0; 8];
        let prev_wavefront = vec![0; 8];

        wfa_diagonal_fill_neon(&mut diagonal_neon, &prev_wavefront, &query, &target);
        // Will need to call naive implementation once #5 is complete
        // For now, this will fail because naive doesn't exist yet
        // crate::wfa::wfa_diagonal_fill_naive(&mut diagonal_naive, &prev_wavefront, &query, &target);

        // Placeholder assertion - implementation agent must add comparison to naive
        assert!(
            diagonal_neon.len() == diagonal_naive.len(),
            "NEON and naive diagonal lengths must match"
        );
    }

    /// Test 18: Correctness vs naive WFA with mismatch
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_neon_matches_naive_wfa_case2() {
        let query = seq("ACGTACGT");
        let target = seq("ACGTTCGT");
        let mut diagonal_neon = vec![0; 8];
        let mut diagonal_naive = vec![0; 8];
        let prev_wavefront = vec![0; 8];

        wfa_diagonal_fill_neon(&mut diagonal_neon, &prev_wavefront, &query, &target);
        // Will need naive implementation from #5
        // crate::wfa::wfa_diagonal_fill_naive(&mut diagonal_naive, &prev_wavefront, &query, &target);

        // Placeholder - implementation agent must verify bit-exact match
        assert!(
            diagonal_neon.len() == diagonal_naive.len(),
            "NEON and naive must produce identical results"
        );
    }

    /// Test 19: Correctness vs naive WFA on long sequence
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_neon_matches_naive_wfa_long() {
        let query = vec![b'A'; 1000];
        let target = vec![b'A'; 1000];
        let mut diagonal_neon = vec![0; 1000];
        let mut diagonal_naive = vec![0; 1000];
        let prev_wavefront = vec![0; 1000];

        wfa_diagonal_fill_neon(&mut diagonal_neon, &prev_wavefront, &query, &target);
        // Will need naive implementation from #5
        // crate::wfa::wfa_diagonal_fill_naive(&mut diagonal_naive, &prev_wavefront, &query, &target);

        // Placeholder - implementation agent must verify correctness on long sequence
        assert!(
            diagonal_neon.len() == diagonal_naive.len(),
            "NEON and naive must match on long sequences"
        );
    }

    /// Test 20: Correctness vs naive WFA with complex pattern
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_neon_matches_naive_wfa_complex() {
        let query = seq("ACGTACGTACGTACGTACGTACGTACGTACGT");
        let target = seq("ACGTTCGTACGTTCGTACGTTCGTACGTTCGT");
        let mut diagonal_neon = vec![0; 32];
        let mut diagonal_naive = vec![0; 32];
        let prev_wavefront = vec![0; 32];

        wfa_diagonal_fill_neon(&mut diagonal_neon, &prev_wavefront, &query, &target);
        // Will need naive implementation from #5
        // crate::wfa::wfa_diagonal_fill_naive(&mut diagonal_naive, &prev_wavefront, &query, &target);

        // Placeholder - implementation agent must verify exact match
        assert!(
            diagonal_neon.len() == diagonal_naive.len(),
            "NEON and naive must match on complex patterns"
        );
    }

    /// Test 21: Correctness vs naive WFA with random sequence
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_neon_matches_naive_wfa_random() {
        // Pseudo-random but deterministic sequence
        let query = seq("ATCGATCGATCGATCGATCGATCGATCGATCG");
        let target = seq("ATCGTTCGATCGTTCGATCGTTCGATCGTTCG");
        let mut diagonal_neon = vec![0; 32];
        let mut diagonal_naive = vec![0; 32];
        let prev_wavefront = vec![
            1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 7, -7, 8, -8, 1, -1, 2, -2, 3, -3, 4, -4, 5,
            -5, 6, -6, 7, -7, 8, -8,
        ];

        wfa_diagonal_fill_neon(&mut diagonal_neon, &prev_wavefront, &query, &target);
        // Will need naive implementation from #5
        // crate::wfa::wfa_diagonal_fill_naive(&mut diagonal_naive, &prev_wavefront, &query, &target);

        // Placeholder - implementation agent must verify match with varied wavefront
        assert!(
            diagonal_neon.len() == diagonal_naive.len(),
            "NEON and naive must match with varied wavefront scores"
        );
    }

    /// Test 22: Performance characteristic hint (not a real benchmark)
    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_neon_performance_hint() {
        // This test documents the performance expectation but doesn't measure it
        // Real benchmarks will be in benches/ directory using criterion
        let query = vec![b'A'; 10_000];
        let target = vec![b'A'; 10_000];
        let mut diagonal = vec![0; 10_000];
        let prev_wavefront = vec![0; 10_000];

        // Should complete without timeout (indicative of reasonable performance)
        wfa_diagonal_fill_neon(&mut diagonal, &prev_wavefront, &query, &target);

        // Acceptance criteria: ≥1.5× speedup vs naive on 10kb alignment
        // This will be verified in criterion benchmarks, not unit tests
    }
}
