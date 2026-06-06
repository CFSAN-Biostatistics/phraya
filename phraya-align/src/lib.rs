pub mod executor;
pub mod seeding;
pub mod wfa_simd;
pub mod wfa_simd_dispatch;

pub use executor::{align_task_with_config, AlignConfig, Strategy};
pub use seeding::{find_seeds, Seed};

#[cfg(test)]
mod local_coverage_tests;

/// Seed anchor position for WFA extension.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SeedAnchor {
    pub query_pos: usize,
    pub target_pos: usize,
}

/// Alignment result from WFA.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alignment {
    pub cigar: String,
    pub edit_distance: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub target_start: usize,
    pub target_end: usize,
}

/// Result type for WFA operations.
pub type WfaResult = Result<Alignment, WfaError>;

/// Errors that can occur during WFA alignment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WfaError {
    InvalidInput(String),
    AlignmentFailed(String),
}

/// Scored alignments with primary and alternatives.
#[derive(Debug, Clone, PartialEq)]
pub struct ScoredAlignments {
    pub primary: Alignment,
    pub alternatives: Vec<Alignment>,
}

// Core WFA implementation - naive baseline
/// Naive scalar WFA extension implementation.
///
/// This is the baseline implementation that all optimized versions match against.
/// Uses simple dynamic programming without SIMD optimization.
pub fn wfa_extend_naive(query: &[u8], target: &[u8], seed: SeedAnchor) -> WfaResult {
    wfa_simd::wfa_extend_naive_impl(query, target, seed)
}

/// SSE4.2-accelerated WFA extension.
///
/// Uses SSE4.2 intrinsics for diagonal fill operations.
/// Falls back to naive implementation if SSE4.2 is not available.
pub fn wfa_extend_simd(query: &[u8], target: &[u8], seed: SeedAnchor) -> WfaResult {
    #[cfg(target_arch = "x86_64")]
    {
        wfa_simd::wfa_extend_simd_impl(query, target, seed)
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        // On non-x86 platforms, fall back to naive
        wfa_simd::wfa_extend_naive_impl(query, target, seed)
    }
}

/// NEON-accelerated WFA extension for ARM64.
///
/// Uses NEON intrinsics for diagonal fill operations on aarch64.
/// NEON is mandatory on ARM64, so no runtime detection is needed.
/// Falls back to naive implementation on non-aarch64 platforms.
pub fn wfa_extend_neon(query: &[u8], target: &[u8], seed: SeedAnchor) -> WfaResult {
    wfa_simd::wfa_extend_neon_impl(query, target, seed)
}

/// Runtime-dispatched WFA extension using multiversion.
///
/// Automatically selects SSE4.2 (if available) or naive implementation
/// based on runtime CPU feature detection via the multiversion crate.
// The SIMD diagonal fill is a *release-time* win: unoptimised `wide` lowers to
// out-of-line calls and is slower than the scalar fill, so debug builds (dev /
// `cargo test`) use scalar and only optimised builds — the only ones used for
// real workloads — pay for SIMD. Correctness of the SIMD kernel is covered in
// debug by the direct `fill_simd` vs `fill_scalar` differential test.
#[cfg(target_arch = "x86_64")]
pub fn wfa_extend(query: &[u8], target: &[u8], seed: SeedAnchor) -> WfaResult {
    if !cfg!(debug_assertions) && is_x86_feature_detected!("sse4.2") {
        wfa_simd::wfa_extend_simd_impl(query, target, seed)
    } else {
        wfa_simd::wfa_extend_naive_impl(query, target, seed)
    }
}

#[cfg(target_arch = "aarch64")]
pub fn wfa_extend(query: &[u8], target: &[u8], seed: SeedAnchor) -> WfaResult {
    // NEON is mandatory on aarch64; use it in release, scalar in debug.
    if cfg!(debug_assertions) {
        wfa_simd::wfa_extend_naive_impl(query, target, seed)
    } else {
        wfa_simd::wfa_extend_neon_impl(query, target, seed)
    }
}

#[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
pub fn wfa_extend(query: &[u8], target: &[u8], seed: SeedAnchor) -> WfaResult {
    // No SIMD path for other architectures.
    wfa_simd::wfa_extend_naive_impl(query, target, seed)
}

/// Score alignments by normalized edit distance and filter alternatives.
pub fn score_alignments(alignments: &[Alignment], query_len: usize) -> ScoredAlignments {
    if alignments.is_empty() {
        panic!("At least one alignment required");
    }

    // Find primary (best) alignment
    let primary = alignments
        .iter()
        .min_by_key(|a| a.edit_distance)
        .unwrap()
        .clone();

    // Filter alternatives with score_ratio >= 0.95
    let alternatives = alignments
        .iter()
        .filter(|a| *a != &primary)
        .filter(|alt| {
            let primary_norm = 1.0 - (primary.edit_distance as f64 / query_len as f64);
            let alt_norm = 1.0 - (alt.edit_distance as f64 / query_len as f64);
            let score_ratio = alt_norm / primary_norm;
            score_ratio >= 0.95
        })
        .cloned()
        .collect();

    ScoredAlignments {
        primary,
        alternatives,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_single_alignment() {
        let alignment = Alignment {
            cigar: "100M".to_string(),
            edit_distance: 0,
            query_start: 0,
            query_end: 100,
            target_start: 0,
            target_end: 100,
        };

        let scored = score_alignments(&[alignment.clone()], 100);
        assert_eq!(scored.primary, alignment);
        assert!(scored.alternatives.is_empty());
    }

    #[test]
    fn test_similar_alignments_both_stored() {
        let primary = Alignment {
            cigar: "100M".to_string(),
            edit_distance: 2,
            query_start: 0,
            query_end: 100,
            target_start: 0,
            target_end: 100,
        };
        let alternative = Alignment {
            cigar: "98M1I1D".to_string(),
            edit_distance: 3, // normalized: (1 - 3/100) / (1 - 2/100) ≈ 0.9898 > 0.95
            query_start: 0,
            query_end: 100,
            target_start: 0,
            target_end: 100,
        };

        let scored = score_alignments(&[primary.clone(), alternative.clone()], 100);
        assert_eq!(scored.primary, primary);
        assert_eq!(scored.alternatives.len(), 1);
        assert_eq!(scored.alternatives[0], alternative);
    }

    #[test]
    fn test_distant_alignments_only_primary() {
        let primary = Alignment {
            cigar: "100M".to_string(),
            edit_distance: 2,
            query_start: 0,
            query_end: 100,
            target_start: 0,
            target_end: 100,
        };
        let alternative = Alignment {
            cigar: "80M10I10D".to_string(),
            edit_distance: 20, // normalized: (1 - 20/100) / (1 - 2/100) ≈ 0.816 < 0.95
            query_start: 0,
            query_end: 100,
            target_start: 0,
            target_end: 100,
        };

        let scored = score_alignments(&[primary.clone(), alternative], 100);
        assert_eq!(scored.primary, primary);
        assert!(scored.alternatives.is_empty());
    }

    #[test]
    fn test_query_len_normalization() {
        // With longer query, same edit distance has better (higher) normalized score
        let primary_short = Alignment {
            cigar: "10M".to_string(),
            edit_distance: 1,
            query_start: 0,
            query_end: 10,
            target_start: 0,
            target_end: 10,
        };
        let alt_short = Alignment {
            cigar: "9M1D".to_string(),
            edit_distance: 1,
            query_start: 0,
            query_end: 10,
            target_start: 0,
            target_end: 10,
        };

        // Same edit distance, same query length → score_ratio = 1.0
        let scored_short = score_alignments(&[primary_short, alt_short.clone()], 10);
        assert_eq!(scored_short.alternatives.len(), 1);

        // Now with longer query
        let primary_long = Alignment {
            cigar: "100M".to_string(),
            edit_distance: 1,
            query_start: 0,
            query_end: 100,
            target_start: 0,
            target_end: 100,
        };
        let alt_long = Alignment {
            cigar: "99M1D".to_string(),
            edit_distance: 1,
            query_start: 0,
            query_end: 100,
            target_start: 0,
            target_end: 100,
        };

        let scored_long = score_alignments(&[primary_long, alt_long], 100);
        assert_eq!(scored_long.alternatives.len(), 1);
    }
}
