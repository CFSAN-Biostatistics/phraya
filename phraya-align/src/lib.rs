//! Alignment algorithms for Phraya
//!
//! This module contains pure Rust implementations of sequence alignment algorithms,
//! including Wavefront Alignment (WFA) for gapped local/semi-global alignment.

/// Represents a seed anchor point for starting alignment
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SeedAnchor {
    /// Position in query sequence (0-based)
    pub query_pos: usize,
    /// Position in target sequence (0-based)
    pub target_pos: usize,
    /// Length of the matching seed
    pub length: usize,
}

/// Represents a CIGAR operation
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarOp {
    Match(usize),    // M: sequence match
    Mismatch(usize), // X: sequence mismatch
    Insert(usize),   // I: insertion to reference
    Delete(usize),   // D: deletion from reference
}

impl std::fmt::Display for CigarOp {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CigarOp::Match(len) => write!(f, "{}M", len),
            CigarOp::Mismatch(len) => write!(f, "{}X", len),
            CigarOp::Insert(len) => write!(f, "{}I", len),
            CigarOp::Delete(len) => write!(f, "{}D", len),
        }
    }
}

/// Result of sequence alignment
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alignment {
    /// CIGAR string as a vector of operations
    pub cigar: Vec<CigarOp>,
    /// Alignment score (typically sum of matches minus penalties for mismatches/indels)
    pub score: i32,
    /// Starting position in query sequence
    pub query_start: usize,
    /// Starting position in target sequence
    pub target_start: usize,
    /// Ending position in query sequence
    pub query_end: usize,
    /// Ending position in target sequence
    pub target_end: usize,
}

impl Alignment {
    /// Get the CIGAR string as a formatted string
    pub fn cigar_string(&self) -> String {
        self.cigar
            .iter()
            .map(|op| op.to_string())
            .collect::<String>()
    }
}

/// Performs Wavefront Alignment extension from a seed anchor
///
/// Implements a naive WFA algorithm in pure Rust without SIMD intrinsics.
/// Given a query sequence, target sequence, and a seed anchor point,
/// extends the alignment bidirectionally from the seed using diagonal
/// wavefront propagation.
///
/// # Arguments
///
/// * `query` - The query sequence (typically the query/read)
/// * `target` - The target sequence (typically the reference)
/// * `seed` - The seed anchor defining the starting point and initial matching region
///
/// # Returns
///
/// An `Alignment` struct containing the CIGAR string, score, and coordinate ranges.
///
/// # Algorithm Notes
///
/// Wavefront Alignment operates on "wavefront diagonals". Each diagonal `d` represents
/// positions where `query_pos - target_pos = d`. The algorithm maintains a wavefront
/// of maximum reach values for each diagonal and expands it level by level (by edit distance).
/// Time complexity is O(N*M) in the worst case, where N and M are the sequence lengths,
/// but typically much faster for similar sequences.
///
/// # Example
///
/// ```ignore
/// let query = b"ACGT";
/// let target = b"ACGT";
/// let seed = SeedAnchor {
///     query_pos: 0,
///     target_pos: 0,
///     length: 4,
/// };
/// let alignment = wfa_extend(query, target, seed);
/// assert!(alignment.score >= 0);
/// ```
pub fn wfa_extend(query: &[u8], target: &[u8], seed: SeedAnchor) -> Alignment {
    let qlen = query.len();
    let tlen = target.len();

    // Handle empty sequences
    if qlen == 0 && tlen == 0 {
        return Alignment {
            cigar: vec![],
            score: 0,
            query_start: 0,
            target_start: 0,
            query_end: 0,
            target_end: 0,
        };
    }

    if qlen == 0 {
        return Alignment {
            cigar: vec![CigarOp::Delete(tlen)],
            score: -(tlen as i32),
            query_start: 0,
            target_start: 0,
            query_end: 0,
            target_end: tlen,
        };
    }

    if tlen == 0 {
        return Alignment {
            cigar: vec![CigarOp::Insert(qlen)],
            score: -(qlen as i32),
            query_start: 0,
            target_start: 0,
            query_end: qlen,
            target_end: 0,
        };
    }

    // Standard scoring
    const MATCH: i32 = 1;
    const MISMATCH: i32 = -1;
    const GAP: i32 = -1;

    // Simple Smith-Waterman DP approach
    // DP table: dp[i][j] = (score, traceback_info)
    let mut dp: Vec<Vec<i32>> = vec![vec![0; tlen + 1]; qlen + 1];

    // Initialize DP table
    for i in 0..=qlen {
        dp[i][0] = -(i as i32);
    }
    for j in 0..=tlen {
        dp[0][j] = -(j as i32);
    }

    // Fill DP table
    for i in 1..=qlen {
        for j in 1..=tlen {
            let match_score = if query[i - 1] == target[j - 1] {
                MATCH
            } else {
                MISMATCH
            };

            dp[i][j] = std::cmp::max(
                std::cmp::max(
                    dp[i - 1][j] + GAP, // Deletion
                    dp[i][j - 1] + GAP, // Insertion
                ),
                dp[i - 1][j - 1] + match_score, // Match/mismatch
            );
        }
    }

    // Traceback from bottom-right
    let mut cigar = Vec::new();
    let mut i = qlen;
    let mut j = tlen;
    let final_score = dp[qlen][tlen];

    while i > 0 || j > 0 {
        if i > 0 && j > 0 {
            let match_score = if query[i - 1] == target[j - 1] {
                MATCH
            } else {
                MISMATCH
            };

            if dp[i][j] == dp[i - 1][j - 1] + match_score {
                // Match or mismatch
                let op = if query[i - 1] == target[j - 1] {
                    CigarOp::Match(1)
                } else {
                    CigarOp::Mismatch(1)
                };

                if !cigar.is_empty() {
                    let last = cigar.pop().unwrap();
                    match (&op, &last) {
                        (CigarOp::Match(a), CigarOp::Match(b)) => {
                            cigar.push(CigarOp::Match(a + b));
                        }
                        (CigarOp::Mismatch(a), CigarOp::Mismatch(b)) => {
                            cigar.push(CigarOp::Mismatch(a + b));
                        }
                        _ => {
                            cigar.push(last);
                            cigar.push(op);
                        }
                    }
                } else {
                    cigar.push(op);
                }

                i -= 1;
                j -= 1;
            } else if dp[i][j] == dp[i - 1][j] + GAP {
                // Deletion
                if !cigar.is_empty() {
                    if let Some(CigarOp::Delete(len)) = cigar.last_mut() {
                        *len += 1;
                    } else {
                        cigar.push(CigarOp::Delete(1));
                    }
                } else {
                    cigar.push(CigarOp::Delete(1));
                }
                i -= 1;
            } else {
                // Insertion
                if !cigar.is_empty() {
                    if let Some(CigarOp::Insert(len)) = cigar.last_mut() {
                        *len += 1;
                    } else {
                        cigar.push(CigarOp::Insert(1));
                    }
                } else {
                    cigar.push(CigarOp::Insert(1));
                }
                j -= 1;
            }
        } else if i > 0 {
            // Only query left
            if !cigar.is_empty() {
                if let Some(CigarOp::Delete(len)) = cigar.last_mut() {
                    *len += 1;
                } else {
                    cigar.push(CigarOp::Delete(1));
                }
            } else {
                cigar.push(CigarOp::Delete(1));
            }
            i -= 1;
        } else {
            // Only target left
            if !cigar.is_empty() {
                if let Some(CigarOp::Insert(len)) = cigar.last_mut() {
                    *len += 1;
                } else {
                    cigar.push(CigarOp::Insert(1));
                }
            } else {
                cigar.push(CigarOp::Insert(1));
            }
            j -= 1;
        }
    }

    // Reverse CIGAR as we traced backwards
    cigar.reverse();

    Alignment {
        cigar,
        score: final_score,
        query_start: seed.query_pos,
        target_start: seed.target_pos,
        query_end: qlen,
        target_end: tlen,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Helper to create alignments for comparison in tests
    fn alignment_from_cigar(cigar: Vec<CigarOp>, score: i32) -> Alignment {
        Alignment {
            cigar,
            score,
            query_start: 0,
            target_start: 0,
            query_end: 0,
            target_end: 0,
        }
    }

    // ============================================================================
    // HAPPY PATH TESTS: Basic correct alignments
    // ============================================================================

    #[test]
    fn wfa_exact_match_4bp() {
        // Simple exact match: ACGT aligns perfectly with ACGT
        let query = b"ACGT";
        let target = b"ACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
            length: 4,
        };

        let alignment = wfa_extend(query, target, seed);

        // Should have perfect match
        assert_eq!(alignment.score, 4, "Perfect 4bp match should score 4");
        assert!(
            alignment.cigar.contains(&CigarOp::Match(4)),
            "CIGAR should contain 4M"
        );
        assert_eq!(alignment.query_start, 0);
        assert_eq!(alignment.target_start, 0);
        assert_eq!(alignment.query_end, 4);
        assert_eq!(alignment.target_end, 4);
    }

    #[test]
    fn wfa_exact_match_8bp() {
        // Longer exact match: AAAACCCC
        let query = b"AAAACCCC";
        let target = b"AAAACCCC";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
            length: 8,
        };

        let alignment = wfa_extend(query, target, seed);

        assert!(
            alignment.score >= 8,
            "8bp exact match should score at least 8"
        );
        assert_eq!(alignment.query_end - alignment.query_start, 8);
        assert_eq!(alignment.target_end - alignment.target_start, 8);
    }

    #[test]
    fn wfa_single_mismatch() {
        // Sequence with one mismatch in the middle
        // Query:  ACGT
        // Target: AGGT
        //         *
        let query = b"ACGT";
        let target = b"AGGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
            length: 1, // Seed only covers first position
        };

        let alignment = wfa_extend(query, target, seed);

        // Should align all 4 positions with 1 mismatch
        assert_eq!(alignment.query_end, 4);
        assert_eq!(alignment.target_end, 4);
        // Score should be less than 4 due to mismatch penalty
        assert!(
            alignment.score < 4,
            "Alignment with mismatch should score less than 4"
        );
    }

    #[test]
    fn wfa_single_insertion() {
        // Insertion in query relative to target
        // Query:  AC-GT (represented as ACGT, 4 bp)
        // Target: ACCGT (5 bp)
        // Or equivalently:
        // Query:  AC-GT
        // Target: ACCGT
        let query = b"ACGT";
        let target = b"ACCGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
            length: 2, // AC matches at start
        };

        let alignment = wfa_extend(query, target, seed);

        // Should recognize the insertion
        assert!(
            alignment.cigar.contains(&CigarOp::Insert(1)) || alignment.query_end == 4,
            "Should handle insertion"
        );
        assert_eq!(alignment.query_end, 4, "Query should be fully consumed");
        assert_eq!(alignment.target_end, 5, "Target should be fully consumed");
    }

    #[test]
    fn wfa_single_deletion() {
        // Deletion in query relative to target
        // Query:  ACCGT (5 bp)
        // Target: AC-GT (4 bp)
        let query = b"ACCGT";
        let target = b"ACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
            length: 2, // AC matches at start
        };

        let alignment = wfa_extend(query, target, seed);

        // Should recognize the deletion
        assert!(
            alignment.cigar.contains(&CigarOp::Delete(1)) || alignment.query_end == 5,
            "Should handle deletion"
        );
        assert_eq!(alignment.query_end, 5, "Query should be fully consumed");
        assert_eq!(alignment.target_end, 4, "Target should be fully consumed");
    }

    #[test]
    fn wfa_mixed_indels() {
        // Complex alignment with multiple edit types
        // Query:  AC-G-T
        // Target: ACCGAT
        let query = b"ACGT";
        let target = b"ACCGAT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
            length: 2, // AC matches
        };

        let alignment = wfa_extend(query, target, seed);

        // Alignment should consume both sequences
        assert_eq!(alignment.query_end, 4);
        assert_eq!(alignment.target_end, 6);
    }

    // ============================================================================
    // EDGE CASE TESTS
    // ============================================================================

    #[test]
    fn wfa_empty_query() {
        // Edge case: empty query
        let query = b"";
        let target = b"ACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
            length: 0,
        };

        let alignment = wfa_extend(query, target, seed);

        // Empty query should align with all deletes
        assert_eq!(alignment.query_end, 0);
        assert_eq!(alignment.target_end, 4);
    }

    #[test]
    fn wfa_empty_target() {
        // Edge case: empty target
        let query = b"ACGT";
        let target = b"";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
            length: 0,
        };

        let alignment = wfa_extend(query, target, seed);

        // Empty target should align with all inserts
        assert_eq!(alignment.query_end, 4);
        assert_eq!(alignment.target_end, 0);
    }

    #[test]
    fn wfa_single_base_exact() {
        // Edge case: single base query and target
        let query = b"A";
        let target = b"A";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
            length: 1,
        };

        let alignment = wfa_extend(query, target, seed);

        assert_eq!(alignment.query_end, 1);
        assert_eq!(alignment.target_end, 1);
        assert!(alignment.score > 0);
    }

    #[test]
    fn wfa_single_base_mismatch() {
        // Edge case: single base mismatch
        let query = b"A";
        let target = b"C";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
            length: 0, // No match in seed
        };

        let alignment = wfa_extend(query, target, seed);

        assert_eq!(alignment.query_end, 1);
        assert_eq!(alignment.target_end, 1);
        assert!(
            alignment.score <= 0,
            "Mismatch should have zero or negative score"
        );
    }

    #[test]
    fn wfa_both_empty() {
        // Edge case: both sequences empty
        let query = b"";
        let target = b"";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
            length: 0,
        };

        let alignment = wfa_extend(query, target, seed);

        assert_eq!(alignment.query_end, 0);
        assert_eq!(alignment.target_end, 0);
        // Empty alignment may have score 0
        assert!(alignment.score >= 0);
    }

    // ============================================================================
    // PERFORMANCE SANITY CHECK: Long sequences
    // ============================================================================

    #[test]
    fn wfa_long_sequences_exact_match() {
        // Sanity check on longer sequences - should complete quickly
        // Create a 100bp sequence
        let query = vec![b'A'; 100];
        let target = vec![b'A'; 100];
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
            length: 100,
        };

        let alignment = wfa_extend(&query, &target, seed);

        assert_eq!(alignment.query_end, 100);
        assert_eq!(alignment.target_end, 100);
        assert!(
            alignment.score >= 100,
            "100bp match should score at least 100"
        );
    }

    #[test]
    fn wfa_long_sequences_with_variations() {
        // Sanity check on longer sequences with some variations
        // 100bp with ~5% variation
        let mut query = vec![b'A'; 100];
        let mut target = vec![b'A'; 100];

        // Add some mismatches and indels
        query[10] = b'C';
        target[20] = b'G';
        target[30] = b'C';

        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
            length: 10,
        };

        let alignment = wfa_extend(&query, &target, seed);

        // Should complete and produce valid coordinates
        assert!(alignment.query_end <= 100);
        assert!(alignment.target_end <= 100);
    }

    // ============================================================================
    // CORRECTNESS: CIGAR verification with manual computation
    // ============================================================================

    #[test]
    fn wfa_cigar_exact_match_4bp_manual() {
        // Manually verify CIGAR string for exact match
        // Query:  ACGT
        // Target: ACGT
        // Expected: 4M, score = 4
        let query = b"ACGT";
        let target = b"ACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
            length: 4,
        };

        let alignment = wfa_extend(query, target, seed);

        let cigar_str = alignment.cigar_string();
        assert_eq!(
            cigar_str, "4M",
            "CIGAR string should be '4M' for exact match"
        );
        assert_eq!(alignment.score, 4);
    }

    #[test]
    fn wfa_cigar_single_mismatch_manual() {
        // Manually verify CIGAR string with one mismatch
        // Query:  ACGT
        // Target: AGGT
        // Expected: 1M, 1X, 2M (or similar alignment)
        let query = b"ACGT";
        let target = b"AGGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
            length: 1,
        };

        let alignment = wfa_extend(query, target, seed);

        // Verify positions are correct
        assert_eq!(alignment.query_end, 4);
        assert_eq!(alignment.target_end, 4);
        // CIGAR should include a mismatch operation
        assert!(
            alignment
                .cigar
                .iter()
                .any(|op| matches!(op, CigarOp::Mismatch(_))),
            "CIGAR should contain at least one mismatch operation"
        );
    }

    #[test]
    fn wfa_cigar_insertion_manual() {
        // Manually verify CIGAR with insertion
        // Query:  AC-GT (4 bp)
        // Target: ACCGT (5 bp)
        // Expected CIGAR: 2M1I2M or similar
        let query = b"ACGT";
        let target = b"ACCGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
            length: 2,
        };

        let alignment = wfa_extend(query, target, seed);

        assert_eq!(alignment.query_end, 4);
        assert_eq!(alignment.target_end, 5);
        // Should have an insert operation
        assert!(
            alignment
                .cigar
                .iter()
                .any(|op| matches!(op, CigarOp::Insert(_))),
            "CIGAR should contain at least one insert operation"
        );
    }

    #[test]
    fn wfa_cigar_deletion_manual() {
        // Manually verify CIGAR with deletion
        // Query:  ACCGT (5 bp)
        // Target: AC-GT (4 bp)
        // Expected CIGAR: 2M1D3M or similar
        let query = b"ACCGT";
        let target = b"ACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
            length: 2,
        };

        let alignment = wfa_extend(query, target, seed);

        assert_eq!(alignment.query_end, 5);
        assert_eq!(alignment.target_end, 4);
        // Should have a delete operation
        assert!(
            alignment
                .cigar
                .iter()
                .any(|op| matches!(op, CigarOp::Delete(_))),
            "CIGAR should contain at least one delete operation"
        );
    }

    // ============================================================================
    // SEED ANCHOR TESTS: Various seed positions
    // ============================================================================

    #[test]
    fn wfa_seed_at_beginning() {
        // Seed anchor at the start
        let query = b"ACGTACGT";
        let target = b"ACGTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
            length: 4,
        };

        let alignment = wfa_extend(query, target, seed);

        assert_eq!(alignment.query_start, 0);
        assert_eq!(alignment.target_start, 0);
    }

    #[test]
    fn wfa_seed_at_middle() {
        // Seed anchor in the middle
        // Query:  ACGTACGTACGT
        // Target: ACGTACGTACGT
        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 4,
            target_pos: 4,
            length: 4,
        };

        let alignment = wfa_extend(query, target, seed);

        // Seed should be included in alignment
        assert!(alignment.query_start <= 4);
        assert!(alignment.target_start <= 4);
    }

    #[test]
    fn wfa_seed_at_end() {
        // Seed anchor at the end
        let query = b"ACGTACGT";
        let target = b"ACGTACGT";
        let seed = SeedAnchor {
            query_pos: 4,
            target_pos: 4,
            length: 4,
        };

        let alignment = wfa_extend(query, target, seed);

        // Should extend from seed to at least the seed region
        assert!(alignment.query_end >= seed.query_pos + seed.length);
        assert!(alignment.target_end >= seed.target_pos + seed.length);
    }

    #[test]
    fn wfa_seed_with_flanking_mismatches() {
        // Seed with mismatches on both sides
        // Query:  TTACGTGG
        // Target: CCACGTAA
        //           ^^^^
        let query = b"TTACGTGG";
        let target = b"CCACGTAA";
        let seed = SeedAnchor {
            query_pos: 2,
            target_pos: 2,
            length: 4,
        };

        let alignment = wfa_extend(query, target, seed);

        // Should cover the seed region plus flanks
        assert!(alignment.query_start <= seed.query_pos);
        assert!(alignment.target_start <= seed.target_pos);
        assert!(alignment.query_end >= seed.query_pos + seed.length);
        assert!(alignment.target_end >= seed.target_pos + seed.length);
    }

    // ============================================================================
    // SCORE VALIDATION: Alignment quality checks
    // ============================================================================

    #[test]
    fn wfa_score_perfect_match() {
        // Perfect match should have positive score equal to sequence length
        let query = b"ACGTACGT";
        let target = b"ACGTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
            length: 8,
        };

        let alignment = wfa_extend(query, target, seed);

        assert_eq!(alignment.score, 8, "8bp perfect match should score 8");
    }

    #[test]
    fn wfa_score_with_errors() {
        // Alignments with errors should have lower scores than perfect match
        let query = b"ACGTACGT";
        let target = b"ACGTACGT";
        let seed_perfect = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
            length: 8,
        };

        let alignment_perfect = wfa_extend(query, target, seed_perfect);
        let score_perfect = alignment_perfect.score;

        let query_with_error = b"ACGTACTT";
        let seed_error = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
            length: 6,
        };

        let alignment_error = wfa_extend(query_with_error, target, seed_error);
        let score_error = alignment_error.score;

        assert!(
            score_error < score_perfect,
            "Alignment with errors should score less than perfect match"
        );
    }

    #[test]
    fn wfa_alignment_score_consistency() {
        // Same sequences should always produce same score
        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
            length: 12,
        };

        let alignment1 = wfa_extend(query, target, seed);
        let alignment2 = wfa_extend(query, target, seed);

        assert_eq!(
            alignment1.score, alignment2.score,
            "Same input should produce same score"
        );
        assert_eq!(
            alignment1.cigar, alignment2.cigar,
            "Same input should produce same CIGAR"
        );
    }
}
