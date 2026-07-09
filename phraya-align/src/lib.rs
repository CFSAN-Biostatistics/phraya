pub mod executor;
pub mod seeding;
pub mod wfa_simd;
pub mod wfa_simd_dispatch;

pub use executor::{
    align_read, align_task_with_config, extend_alternates_bounded, score_bound_max_s,
    wfa_extend_capped, AlignConfig, AlignStats, Strategy, TargetContext, WindowedCoverage,
};
pub use seeding::{
    build_minimizer_index, find_seeds, find_seeds_indexed, find_seeds_indexed_capped,
    seed_occurrence_cap, MinimizerIndex, Seed,
};

#[cfg(test)]
mod local_coverage_tests;

#[cfg(test)]
mod issue_185_tests;

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

/// WFA O(s·n) extension — the single production entry point.
///
/// Aligns `query[seed.query_pos..]` against `target[seed.target_pos..]`.
/// SIMD acceleration runs inside `count_matching_prefix` (SSE2 on x86_64,
/// NEON on aarch64, u64-XOR elsewhere); there is no separate dispatch at
/// this level. Use `wfa_extend_naive` in tests that want the deterministic
/// scalar implementation explicitly.
pub fn wfa_extend(query: &[u8], target: &[u8], seed: SeedAnchor) -> WfaResult {
    wfa_simd::wfa_extend_naive_impl(query, target, seed)
}

/// Myers' bit-parallel edit distance for short-read throughput optimization.
///
/// Issue #144: Implements Myers (1999) bit-parallel algorithm for edit distance computation.
/// Packs DP into bitvectors, advancing ~64 cells per machine word with bitwise operations.
/// Produces identical edit distance and CIGAR output to scalar WFA implementation.
///
/// # Arguments
/// - `query`: The query sequence
/// - `target`: The target sequence
///
/// # Returns
/// A tuple of (edit_distance, cigar_string) matching scalar WFA semantics.
///
/// # Performance
/// Target: ~10x faster than WFA on short reads (≤500bp), measurably faster than
/// portable-SIMD path in release builds.
pub fn myers_edit_distance(query: &[u8], target: &[u8]) -> (usize, String) {
    wfa_simd::myers_edit_distance_impl(query, target)
}

/// Myers fitting extension: the alignment-path counterpart to [`wfa_extend`].
///
/// Aligns `query[seed.query_pos..]` against `target[seed.target_pos..]` in *fitting*
/// mode — the query is fully consumed but the target end is free, so a read windowed
/// against a longer reference is not penalised for the unconsumed tail. This mirrors
/// [`wfa_extend`]'s semantics, producing the same edit distance, CIGAR (M/X/I/D), and
/// consumed target span, so the Myers and WFA strategies are interchangeable for
/// variant calling. Faster than WFA for short reads and higher-divergence alignments.
pub fn myers_extend(query: &[u8], target: &[u8], seed: SeedAnchor) -> WfaResult {
    if seed.query_pos > query.len() || seed.target_pos > target.len() {
        return Err(WfaError::InvalidInput(
            "Seed position beyond sequence length".to_string(),
        ));
    }

    let query_suffix = &query[seed.query_pos..];
    let target_suffix = &target[seed.target_pos..];

    let (edit_distance, cigar, target_consumed) =
        wfa_simd::myers_fitting_impl(query_suffix, target_suffix);

    Ok(Alignment {
        cigar,
        edit_distance,
        query_start: seed.query_pos,
        query_end: seed.query_pos + query_suffix.len(),
        target_start: seed.target_pos,
        target_end: seed.target_pos + target_consumed,
    })
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

    #[test]
    #[should_panic(expected = "At least one alignment required")]
    fn score_alignments_panics_on_empty_slice() {
        score_alignments(&[], 100);
    }

    #[test]
    fn myers_extend_rejects_seed_beyond_query_length() {
        let query = b"ACGT";
        let target = b"ACGTACGT";
        let seed = SeedAnchor {
            query_pos: 5, // beyond query.len() == 4
            target_pos: 0,
        };

        let result = myers_extend(query, target, seed);
        assert!(matches!(result, Err(WfaError::InvalidInput(_))));
    }

    #[test]
    fn myers_extend_rejects_seed_beyond_target_length() {
        let query = b"ACGT";
        let target = b"ACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 10, // beyond target.len() == 4
        };

        let result = myers_extend(query, target, seed);
        assert!(matches!(result, Err(WfaError::InvalidInput(_))));
    }
}
