//! Phraya Index: FM-index and k-mer search structures
//!
//! This module provides efficient indexing structures for reference sequences,
//! including FM-index for seed-based alignment and k-mer queries for confidence scoring.

/// FM-index data structure for rapid sequence search
///
/// Stores Burrows-Wheeler Transform (BWT), sampled suffix array, and occurrence table
/// to support O(m) k-mer uniqueness queries without rebuilding the index.
#[derive(Debug)]
pub struct FmIndex {
    /// Reference sequence (stored for validation and diagnostics)
    pub reference: Vec<u8>,
    /// Burrows-Wheeler Transform of the reference
    pub bwt: Vec<u8>,
    /// Sampled suffix array (stores every `sample_rate`-th suffix array entry)
    pub sa_samples: Vec<usize>,
    /// Sample rate for suffix array (e.g., every 128 bases)
    pub sample_rate: usize,
    /// Occurrence table: cumulative count of each base up to each BWT position
    /// Structure: For each base in [A,C,G,T], stores cumulative counts
    pub occ_table: Vec<[usize; 4]>,
}

impl FmIndex {
    /// Construct an FM-index from a reference sequence
    ///
    /// # Arguments
    /// * `reference` - DNA sequence (ACGT only, case-insensitive)
    /// * `sample_rate` - Suffix array sample rate (typically 128)
    ///
    /// # Returns
    /// * `Result<FmIndex>` - Constructed index
    ///
    /// # Panics
    /// Panics if reference contains non-DNA characters
    pub fn build(reference: &[u8], sample_rate: usize) -> Result<Self, FmIndexError> {
        if reference.is_empty() {
            return Err(FmIndexError::EmptyReference);
        }
        if sample_rate == 0 {
            return Err(FmIndexError::InvalidSampleRate);
        }

        todo!("FM-index construction not implemented")
    }

    /// Query whether a k-mer appears exactly once in the reference
    ///
    /// # Arguments
    /// * `kmer` - The k-mer sequence (ACGT only)
    /// * `k` - Length of the k-mer (default 31)
    ///
    /// # Returns
    /// * `bool` - True if k-mer appears exactly once, false otherwise
    ///
    /// # Panics
    /// Panics if k-mer contains non-DNA characters or if k > reference length
    pub fn is_kmer_unique(&self, kmer: &[u8], k: usize) -> Result<bool, FmIndexError> {
        if kmer.len() != k {
            return Err(FmIndexError::KmerLengthMismatch { expected: k, got: kmer.len() });
        }
        if k == 0 {
            return Err(FmIndexError::InvalidKmerLength(0));
        }
        if k > self.reference.len() {
            return Err(FmIndexError::KmerLongerThanReference { k, ref_len: self.reference.len() });
        }

        todo!("K-mer uniqueness query not implemented")
    }

    /// Count k-mer occurrences in the reference
    ///
    /// # Arguments
    /// * `kmer` - The k-mer sequence (ACGT only)
    /// * `k` - Length of the k-mer
    ///
    /// # Returns
    /// * Count of occurrences (0, 1, or >1)
    pub fn count_kmer_occurrences(&self, kmer: &[u8], k: usize) -> Result<usize, FmIndexError> {
        if kmer.len() != k {
            return Err(FmIndexError::KmerLengthMismatch { expected: k, got: kmer.len() });
        }

        todo!("K-mer counting not implemented")
    }
}

/// Errors in FM-index construction and querying
#[derive(Debug, thiserror::Error, PartialEq, Eq, Clone)]
pub enum FmIndexError {
    #[error("Empty reference sequence")]
    EmptyReference,

    #[error("Invalid sample rate (must be > 0)")]
    InvalidSampleRate,

    #[error("Invalid k-mer length")]
    InvalidKmerLength(usize),

    #[error("K-mer length mismatch: expected {expected}, got {got}")]
    KmerLengthMismatch { expected: usize, got: usize },

    #[error("K-mer longer than reference: k={k}, ref_len={ref_len}")]
    KmerLongerThanReference { k: usize, ref_len: usize },

    #[error("Invalid character in sequence")]
    InvalidCharacter(u8),

    #[error("Invalid BWT")]
    InvalidBwt,

    #[error("Suffix array construction failed")]
    SuffixArrayError,
}

#[cfg(test)]
mod tests {
    use super::*;

    // ========== HAPPY PATH TESTS ==========

    /// Test: Construct index from small valid sequence
    #[test]
    fn test_build_small_sequence() {
        let reference = b"ACGTACGTACGT";
        let index = FmIndex::build(reference, 4);
        assert!(index.is_ok(), "Should successfully build index from valid sequence");
    }

    /// Test: Query unique k-mer in small sequence
    #[test]
    fn test_query_unique_kmer() {
        let reference = b"ACGTACGTAAAA";
        let index = FmIndex::build(reference, 4).expect("Failed to build index");

        // "AAAA" appears once at the end
        let result = index.is_kmer_unique(b"AAAA", 4);
        assert!(result.is_ok());
        assert!(result.unwrap(), "AAAA should be unique in ACGTACGTAAAA");
    }

    /// Test: Query non-unique k-mer
    #[test]
    fn test_query_non_unique_kmer() {
        let reference = b"ACGTACGTACGT";
        let index = FmIndex::build(reference, 4).expect("Failed to build index");

        // "ACGT" appears twice
        let result = index.is_kmer_unique(b"ACGT", 4);
        assert!(result.is_ok());
        assert!(!result.unwrap(), "ACGT should not be unique (appears twice)");
    }

    /// Test: Query k-mer that appears 3+ times
    #[test]
    fn test_query_highly_repetitive_kmer() {
        let reference = b"ATATATATATAT";
        let index = FmIndex::build(reference, 4).expect("Failed to build index");

        // "AT" appears 6 times
        let result = index.is_kmer_unique(b"AT", 2);
        assert!(result.is_ok());
        assert!(!result.unwrap(), "AT should not be unique in ATATATATATAT");
    }

    /// Test: Default k-mer length (k=31)
    #[test]
    fn test_default_kmer_length() {
        // Create a 50bp reference with one 31-mer that's unique
        let reference = b"ACGTACGTACGTACGTACGTACGTACGTACGTTTTTTTTTTTTTTTTTTTT";
        assert_eq!(reference.len(), 50);

        let index = FmIndex::build(reference, 8).expect("Failed to build index");

        // First 31 bases should be unique (they don't repeat)
        let kmer = &reference[0..31];
        let result = index.is_kmer_unique(kmer, 31);
        assert!(result.is_ok());
        assert!(result.unwrap(), "First 31 bases should be unique");
    }

    // ========== EDGE CASE TESTS ==========

    /// Test: Single base sequence
    #[test]
    fn test_single_base_sequence() {
        let reference = b"A";
        let index = FmIndex::build(reference, 1);
        assert!(index.is_ok(), "Should handle single base");

        let index = index.unwrap();
        let result = index.is_kmer_unique(b"A", 1);
        assert!(result.is_ok());
        assert!(result.unwrap(), "Single A should be unique");
    }

    /// Test: All same base
    #[test]
    fn test_all_same_base() {
        let reference = b"AAAAAAAAAA";
        let index = FmIndex::build(reference, 4).expect("Failed to build index");

        // Any k-mer of As is not unique
        let result = index.is_kmer_unique(b"AAA", 3);
        assert!(result.is_ok());
        assert!(!result.unwrap(), "AAA should not be unique in all-A sequence");

        // Single A is not unique (appears 10 times)
        let result = index.is_kmer_unique(b"A", 1);
        assert!(result.is_ok());
        assert!(!result.unwrap(), "Single A should not be unique");
    }

    /// Test: K-mer length equals reference length
    #[test]
    fn test_kmer_equals_reference_length() {
        let reference = b"ACGTACGT";
        let index = FmIndex::build(reference, 4).expect("Failed to build index");

        // Query entire sequence - must be unique
        let result = index.is_kmer_unique(b"ACGTACGT", 8);
        assert!(result.is_ok());
        assert!(result.unwrap(), "Entire reference should be unique");
    }

    /// Test: Very short reference with long k-mer request
    #[test]
    fn test_short_reference_long_kmer_request() {
        let reference = b"ACGT";
        let index = FmIndex::build(reference, 2).expect("Failed to build index");

        // Requesting a 5-mer from a 4bp reference should error
        let result = index.is_kmer_unique(b"ACGTA", 5);
        assert!(result.is_err());
        assert!(matches!(
            result.unwrap_err(),
            FmIndexError::KmerLongerThanReference { .. }
        ));
    }

    // ========== ERROR HANDLING TESTS ==========

    /// Test: Empty reference
    #[test]
    fn test_empty_reference() {
        let reference = b"";
        let result = FmIndex::build(reference, 4);
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), FmIndexError::EmptyReference));
    }

    /// Test: Invalid sample rate (zero)
    #[test]
    fn test_invalid_sample_rate_zero() {
        let reference = b"ACGTACGT";
        let result = FmIndex::build(reference, 0);
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), FmIndexError::InvalidSampleRate));
    }

    /// Test: K-mer length mismatch
    #[test]
    fn test_kmer_length_mismatch() {
        let reference = b"ACGTACGTACGT";
        let index = FmIndex::build(reference, 4).expect("Failed to build index");

        // Requesting 4-mer with 5-byte buffer
        let result = index.is_kmer_unique(b"ACGTA", 4);
        assert!(result.is_err());
        assert!(matches!(
            result.unwrap_err(),
            FmIndexError::KmerLengthMismatch { expected: 4, got: 5 }
        ));
    }

    /// Test: Zero k-mer length
    #[test]
    fn test_zero_kmer_length() {
        let reference = b"ACGTACGTACGT";
        let index = FmIndex::build(reference, 4).expect("Failed to build index");

        let result = index.is_kmer_unique(b"", 0);
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), FmIndexError::InvalidKmerLength(0)));
    }

    // ========== CORRECTNESS TESTS ==========

    /// Test: Known unique k-mer identified correctly
    #[test]
    fn test_known_unique_kmer_correctness() {
        // Sequence with a known unique motif at the end
        let reference = b"ACGTACGTACGTAAAAAAAAAAA";
        let index = FmIndex::build(reference, 4).expect("Failed to build index");

        // AAAAAA (6 As) is unique
        let result = index.is_kmer_unique(b"AAAAAA", 6);
        assert!(result.is_ok());
        assert!(result.unwrap(), "AAAAAA should be unique at end");

        // But AAAAA (5 As) might not be (depends on boundary)
        // Just verify it computes without panic
        let result2 = index.is_kmer_unique(b"AAAAA", 5);
        assert!(result2.is_ok());
    }

    /// Test: K-mer at sequence boundary
    #[test]
    fn test_kmer_at_start_boundary() {
        let reference = b"CCCCGGGGAAAA";
        let index = FmIndex::build(reference, 4).expect("Failed to build index");

        // First 4 bases CCCC
        let result = index.is_kmer_unique(b"CCCC", 4);
        assert!(result.is_ok());
        assert!(result.unwrap(), "CCCC at start should be unique");
    }

    /// Test: K-mer at sequence end boundary
    #[test]
    fn test_kmer_at_end_boundary() {
        let reference = b"CCCCGGGGAAAA";
        let index = FmIndex::build(reference, 4).expect("Failed to build index");

        // Last 4 bases AAAA
        let result = index.is_kmer_unique(b"AAAA", 4);
        assert!(result.is_ok());
        assert!(result.unwrap(), "AAAA at end should be unique");
    }

    /// Test: Overlapping k-mers
    #[test]
    fn test_overlapping_kmers() {
        // Reference: ACGTAGCGTAGCG
        // 4-mers: ACGT, CGTA, GTAG, TAGC, AGCG, GCGT, CGTA (duplicate)
        let reference = b"ACGTAGCGTAGCG";
        let index = FmIndex::build(reference, 4).expect("Failed to build index");

        // CGTA appears twice (positions 1 and 8)
        let result = index.is_kmer_unique(b"CGTA", 4);
        assert!(result.is_ok());
        assert!(!result.unwrap(), "CGTA should appear twice");

        // AGCG appears once
        let result = index.is_kmer_unique(b"AGCG", 4);
        assert!(result.is_ok());
        assert!(result.unwrap(), "AGCG should be unique");
    }

    /// Test: Case insensitivity
    #[test]
    fn test_case_insensitivity() {
        let reference = b"AcGtAcGt";
        let index = FmIndex::build(reference, 4).expect("Should handle mixed case");

        // Query should work with uppercase
        let result = index.is_kmer_unique(b"ACGT", 4);
        assert!(result.is_ok());
    }

    /// Test: Complex repetitive pattern
    #[test]
    fn test_complex_repetitive_pattern() {
        // Reference with tandem repeats: ACGTACGTACGTACGT
        let reference = b"ACGTACGTACGTACGT";
        let index = FmIndex::build(reference, 4).expect("Failed to build index");

        // ACGT repeats 4 times
        let result = index.is_kmer_unique(b"ACGT", 4);
        assert!(result.is_ok());
        assert!(!result.unwrap(), "ACGT should repeat 4 times");

        // ACGTACGT repeats 2 times
        let result = index.is_kmer_unique(b"ACGTACGT", 8);
        assert!(result.is_ok());
        assert!(!result.unwrap(), "ACGTACGT should repeat 2 times");
    }

    // ========== PERFORMANCE TESTS ==========

    /// Test: Index construction for 5Mbp reference completes in <10 seconds
    #[test]
    #[ignore] // Only run manually or in benchmark suite
    fn perf_test_5mbp_construction() {
        use std::time::Instant;

        // Generate a 5Mbp reference with some structure
        let mut reference = Vec::with_capacity(5_000_000);
        let bases = [b'A', b'C', b'G', b'T'];
        for i in 0..5_000_000 {
            reference.push(bases[i % 4]);
        }

        let start = Instant::now();
        let result = FmIndex::build(&reference, 128);
        let elapsed = start.elapsed();

        assert!(result.is_ok(), "Should successfully build 5Mbp index");
        assert!(
            elapsed.as_secs() < 10,
            "Index construction should complete in < 10 seconds, took: {:?}",
            elapsed
        );
    }

    /// Test: K-mer queries scale to large references
    #[test]
    #[ignore] // Only run manually or in benchmark suite
    fn perf_test_kmer_query_large_reference() {
        use std::time::Instant;

        // 1Mbp reference
        let mut reference = Vec::with_capacity(1_000_000);
        let bases = [b'A', b'C', b'G', b'T'];
        for i in 0..1_000_000 {
            reference.push(bases[i % 4]);
        }

        let index = FmIndex::build(&reference, 128).expect("Failed to build index");

        // Query 1000 different k-mers
        let kmer = b"ACGTACGTACGTACGTACGTACGTACGTACGT"; // 31-mer
        let start = Instant::now();
        for _ in 0..1000 {
            let _ = index.is_kmer_unique(kmer, 31);
        }
        let elapsed = start.elapsed();

        // Should complete reasonably fast (O(31) per query)
        println!("1000 queries on 1Mbp index: {:?}", elapsed);
    }

    /// Test: Index building with various sample rates
    #[test]
    fn test_index_construction_various_sample_rates() {
        let reference = b"ACGTACGTACGTACGTACGTACGTACGTACGT";

        // Small sample rate
        let result1 = FmIndex::build(reference, 1);
        assert!(result1.is_ok(), "Should accept sample_rate=1");

        // Medium sample rate
        let result2 = FmIndex::build(reference, 8);
        assert!(result2.is_ok(), "Should accept sample_rate=8");

        // Large sample rate (close to reference length)
        let result3 = FmIndex::build(reference, 64);
        assert!(result3.is_ok(), "Should accept large sample_rate");
    }

    /// Test: Stress test with many k-mers of varying length
    #[test]
    fn test_stress_test_many_kmers() {
        let reference = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let index = FmIndex::build(reference, 8).expect("Failed to build index");

        // Query k-mers of different lengths
        for k in 1..=31 {
            if k <= reference.len() {
                let kmer = &reference[0..k];
                let result = index.is_kmer_unique(kmer, k);
                assert!(result.is_ok(), "Should handle k-mer of length {}", k);
            }
        }
    }

    /// Test: K-mer count method
    #[test]
    fn test_kmer_count_occurrences() {
        let reference = b"ACGTACGTACGTACGT";
        let index = FmIndex::build(reference, 4).expect("Failed to build index");

        // ACGT appears 4 times
        let count = index.count_kmer_occurrences(b"ACGT", 4);
        assert!(count.is_ok());
        assert_eq!(count.unwrap(), 4, "ACGT should appear 4 times");

        // Entire sequence appears once
        let count = index.count_kmer_occurrences(b"ACGTACGTACGTACGT", 16);
        assert!(count.is_ok());
        assert_eq!(count.unwrap(), 1, "Full sequence should appear once");
    }
}
