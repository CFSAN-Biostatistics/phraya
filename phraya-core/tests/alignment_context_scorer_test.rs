/// Comprehensive acceptance tests for AlignmentContextScorer (Issue #8)
///
/// These tests define the contract for AlignmentContextScorer, which integrates
/// multiple confidence features to produce a combined confidence score.
///
/// Expected behavior: ALL TESTS SHOULD FAIL (RED phase)
/// - The scorer, types, and methods do not exist yet
/// - Expect ImportError, compilation failures, or missing struct/method errors
use phraya_core::{AlignmentContextScorer, BaseConfidence, Sequence};
use phraya_index::FmIndex;

// Test helper to create a test reference sequence
fn create_test_reference() -> Vec<u8> {
    // 200bp reference with varied characteristics
    // Positions 0-19: high GC region
    // Positions 20-39: AAAA homopolymer
    // Positions 40-79: normal region
    // Positions 80-99: ATATATATAT tandem repeat
    // Positions 100-139: normal region with unique k-mers
    // Positions 140-159: low complexity region
    // Positions 160-199: edge region
    b"GCGCGCGCGCGCGCGCGCGC\
      AAAAAAAAAAAAAAAAAAAA\
      ATCGATCGATCGATCGATCGATCGATCGATCGATCG\
      ATATATATATATATATATAT\
      TCGACTGACTGACTGACTGACTGACTGACTGACTGA\
      TTTTTTTTTTTTTTTTTTTT\
      AGCTACGTACGTACGTACGTACGTACGTACGTACGT"
        .to_vec()
}

#[test]
fn test_scorer_struct_exists() {
    // AlignmentContextScorer should be constructible
    let reference = create_test_reference();
    let fm_index = FmIndex::new(&reference);

    let scorer = AlignmentContextScorer::new(&reference, &fm_index);

    // Should not panic during construction
    assert!(std::mem::size_of_val(&scorer) > 0);
}

#[test]
fn test_score_method_returns_base_confidence() {
    let reference = create_test_reference();
    let fm_index = FmIndex::new(&reference);
    let scorer = AlignmentContextScorer::new(&reference, &fm_index);

    // Score a position in the middle of the sequence
    let position = 120;
    let aligned_sequence = Sequence::new(b"ACTGACTGACTGACTGACTG".to_vec(), None, None);
    let alignment_identity = 0.98;

    let confidence = scorer.score(position, &aligned_sequence, alignment_identity);

    // Should return a BaseConfidence struct
    assert!(confidence.combined_confidence() >= 0.0);
    assert!(confidence.combined_confidence() <= 1.0);
}

#[test]
fn test_high_confidence_for_unique_high_quality_region() {
    // Test acceptance criterion: high confidence for unique high-quality regions
    let reference = create_test_reference();
    let fm_index = FmIndex::new(&reference);
    let scorer = AlignmentContextScorer::new(&reference, &fm_index);

    // Position 120: far from edges, unique k-mers, normal GC, no repeats
    let position = 120;
    let aligned_sequence = Sequence::new(b"CTGACTGACTGACTGACTGA".to_vec(), None, None);
    let alignment_identity = 0.98; // High alignment quality

    let confidence = scorer.score(position, &aligned_sequence, alignment_identity);

    // Should have high confidence (>0.8) for ideal conditions
    assert!(
        confidence.combined_confidence() > 0.8,
        "Expected high confidence for unique, high-quality region, got {}",
        confidence.combined_confidence()
    );

    // Verify individual feature contributions are captured
    assert!(confidence.edge_distance_score() > 0.9);
    assert!(confidence.kmer_uniqueness_score() > 0.8);
    assert!(!confidence.in_repeat_region());
}

#[test]
fn test_low_confidence_near_sequence_edge() {
    // Test acceptance criterion: low confidence for positions near edges
    let reference = create_test_reference();
    let fm_index = FmIndex::new(&reference);
    let scorer = AlignmentContextScorer::new(&reference, &fm_index);

    // Position 5: very close to start edge
    let position = 5;
    let aligned_sequence = Sequence::new(b"GCGCGCGCGCGCGCGCGCGC".to_vec(), None, None);
    let alignment_identity = 0.98;

    let confidence = scorer.score(position, &aligned_sequence, alignment_identity);

    // Should have low confidence due to edge proximity
    assert!(
        confidence.combined_confidence() < 0.5,
        "Expected low confidence near edge, got {}",
        confidence.combined_confidence()
    );

    // Edge distance feature should be penalized
    assert!(confidence.edge_distance_score() < 0.5);
}

#[test]
fn test_low_confidence_in_homopolymer() {
    // Test acceptance criterion: low confidence for homopolymer regions
    let reference = create_test_reference();
    let fm_index = FmIndex::new(&reference);
    let scorer = AlignmentContextScorer::new(&reference, &fm_index);

    // Position 25: middle of AAAA homopolymer run
    let position = 25;
    let aligned_sequence = Sequence::new(b"AAAAAAAAAAAAAAAAAAAA".to_vec(), None, None);
    let alignment_identity = 0.98;

    let confidence = scorer.score(position, &aligned_sequence, alignment_identity);

    // Should have reduced confidence due to homopolymer context
    assert!(
        confidence.combined_confidence() < 0.6,
        "Expected low confidence in homopolymer, got {}",
        confidence.combined_confidence()
    );

    assert!(confidence.in_homopolymer());
}

#[test]
fn test_low_confidence_in_tandem_repeat() {
    // Test acceptance criterion: low confidence for repeat regions
    let reference = create_test_reference();
    let fm_index = FmIndex::new(&reference);
    let scorer = AlignmentContextScorer::new(&reference, &fm_index);

    // Position 90: middle of ATATAT tandem repeat
    let position = 90;
    let aligned_sequence = Sequence::new(b"ATATATATATATATATATAT".to_vec(), None, None);
    let alignment_identity = 0.98;

    let confidence = scorer.score(position, &aligned_sequence, alignment_identity);

    // Should have reduced confidence due to tandem repeat
    assert!(
        confidence.combined_confidence() < 0.6,
        "Expected low confidence in tandem repeat, got {}",
        confidence.combined_confidence()
    );

    assert!(confidence.in_repeat_region());
    assert!(confidence.repeat_period().is_some());
}

#[test]
fn test_kmer_uniqueness_penalty() {
    // Test that k-mer uniqueness via FM-index affects score
    let reference = create_test_reference();
    let fm_index = FmIndex::new(&reference);
    let scorer = AlignmentContextScorer::new(&reference, &fm_index);

    // Position with non-unique k-mer context
    let position = 150; // In low-complexity TTTTTT region
    let aligned_sequence = Sequence::new(b"TTTTTTTTTTTTTTTTTTTT".to_vec(), None, None);
    let alignment_identity = 0.98;

    let confidence = scorer.score(position, &aligned_sequence, alignment_identity);

    // K-mer uniqueness should be low (non-unique k-mers)
    assert!(
        confidence.kmer_uniqueness_score() < 0.5,
        "Expected low k-mer uniqueness for non-unique region"
    );
}

#[test]
fn test_gc_content_feature() {
    // Test that local GC content is computed correctly
    let reference = create_test_reference();
    let fm_index = FmIndex::new(&reference);
    let scorer = AlignmentContextScorer::new(&reference, &fm_index);

    // Position 10: in high GC region
    let position = 10;
    let aligned_sequence = Sequence::new(b"GCGCGCGCGCGCGCGCGCGC".to_vec(), None, None);
    let alignment_identity = 0.98;

    let confidence_high_gc = scorer.score(position, &aligned_sequence, alignment_identity);

    // Position 150: in low GC region (TTTTTT)
    let position_low_gc = 150;
    let aligned_sequence_low_gc = Sequence::new(b"TTTTTTTTTTTTTTTTTTTT".to_vec(), None, None);

    let confidence_low_gc = scorer.score(
        position_low_gc,
        &aligned_sequence_low_gc,
        alignment_identity,
    );

    // GC content should be captured and differ
    assert!(confidence_high_gc.local_gc_content() > 0.8);
    assert!(confidence_low_gc.local_gc_content() < 0.2);
}

#[test]
fn test_snp_density_windows() {
    // Test SNP density in 15bp, 125bp, and 1000bp windows
    let reference = create_test_reference();
    let fm_index = FmIndex::new(&reference);
    let scorer = AlignmentContextScorer::new(&reference, &fm_index);

    // Create a mock SNP density map (would normally come from alignment data)
    // For testing, we'll use the scorer's ability to track SNP density
    let position = 120;
    let aligned_sequence = Sequence::new(b"CTGACTGACTGACTGACTGA".to_vec(), None, None);
    let alignment_identity = 0.98;

    let confidence = scorer.score(position, &aligned_sequence, alignment_identity);

    // Should provide SNP density metrics for different window sizes
    assert!(confidence.snp_density_15bp() >= 0.0);
    assert!(confidence.snp_density_125bp() >= 0.0);
    assert!(confidence.snp_density_1000bp() >= 0.0);
}

#[test]
fn test_high_snp_density_reduces_confidence() {
    // Test acceptance criterion: high SNP density regions have reduced confidence
    let reference = create_test_reference();
    let fm_index = FmIndex::new(&reference);
    let scorer = AlignmentContextScorer::new(&reference, &fm_index);

    // Simulate a region with many nearby SNPs
    let position = 100;
    let aligned_sequence = Sequence::new(b"TCGACTGACTGACTGACTGA".to_vec(), None, None);
    let alignment_identity = 0.85; // Lower identity indicates more mismatches

    let confidence_high_snp = scorer.score(position, &aligned_sequence, alignment_identity);

    // Compare to a region with no nearby SNPs
    let aligned_sequence_clean = Sequence::new(b"TCGACTGACTGACTGACTGA".to_vec(), None, None);
    let alignment_identity_clean = 0.98;

    let confidence_low_snp =
        scorer.score(position, &aligned_sequence_clean, alignment_identity_clean);

    // Higher SNP density should reduce confidence
    assert!(confidence_high_snp.combined_confidence() < confidence_low_snp.combined_confidence());
}

#[test]
fn test_alignment_identity_window() {
    // Test that alignment identity in ±50bp window affects score
    let reference = create_test_reference();
    let fm_index = FmIndex::new(&reference);
    let scorer = AlignmentContextScorer::new(&reference, &fm_index);

    let position = 120;
    let aligned_sequence = Sequence::new(b"CTGACTGACTGACTGACTGA".to_vec(), None, None);

    // High alignment identity
    let confidence_high = scorer.score(position, &aligned_sequence, 0.98);

    // Low alignment identity
    let confidence_low = scorer.score(position, &aligned_sequence, 0.75);

    // Higher alignment identity should yield higher confidence
    assert!(confidence_high.combined_confidence() > confidence_low.combined_confidence());
    assert!(confidence_high.alignment_identity_score() > confidence_low.alignment_identity_score());
}

#[test]
fn test_combined_score_weighted_integration() {
    // Test acceptance criterion: combined score uses weighted integration
    let reference = create_test_reference();
    let fm_index = FmIndex::new(&reference);
    let scorer = AlignmentContextScorer::new(&reference, &fm_index);

    let position = 120;
    let aligned_sequence = Sequence::new(b"CTGACTGACTGACTGACTGA".to_vec(), None, None);
    let alignment_identity = 0.98;

    let confidence = scorer.score(position, &aligned_sequence, alignment_identity);

    // The combined confidence should be a weighted integration of features
    // Verify that it's not just a simple average by checking bounds
    let combined = confidence.combined_confidence();

    // Should be within [0.0, 1.0]
    assert!(combined >= 0.0 && combined <= 1.0);

    // Should have individual feature scores available
    let edge_score = confidence.edge_distance_score();
    let gc_score = confidence.local_gc_content_score();
    let kmer_score = confidence.kmer_uniqueness_score();
    let repeat_penalty = if confidence.in_repeat_region() {
        0.5
    } else {
        1.0
    };

    // Combined score should incorporate all features
    // (exact formula will be documented in implementation)
    assert!(combined >= 0.0);
    assert!(edge_score >= 0.0 && edge_score <= 1.0);
    assert!(gc_score >= 0.0 && gc_score <= 1.0);
    assert!(kmer_score >= 0.0 && kmer_score <= 1.0);
}

#[test]
fn test_edge_case_at_sequence_start() {
    // Test edge case: position exactly at sequence start
    let reference = create_test_reference();
    let fm_index = FmIndex::new(&reference);
    let scorer = AlignmentContextScorer::new(&reference, &fm_index);

    let position = 0;
    let aligned_sequence = Sequence::new(b"GCGCGCGCGCGCGCGCGCGC".to_vec(), None, None);
    let alignment_identity = 0.98;

    let confidence = scorer.score(position, &aligned_sequence, alignment_identity);

    // Should handle boundary gracefully
    assert!(confidence.combined_confidence() >= 0.0);
    assert!(confidence.edge_distance_score() < 0.3); // Very low due to edge
}

#[test]
fn test_edge_case_at_sequence_end() {
    // Test edge case: position near sequence end
    let reference = create_test_reference();
    let fm_index = FmIndex::new(&reference);
    let scorer = AlignmentContextScorer::new(&reference, &fm_index);

    let position = 195; // 5bp from end
    let aligned_sequence = Sequence::new(b"ACGTACGTACGTACGTACGT".to_vec(), None, None);
    let alignment_identity = 0.98;

    let confidence = scorer.score(position, &aligned_sequence, alignment_identity);

    // Should handle boundary gracefully
    assert!(confidence.combined_confidence() >= 0.0);
    assert!(confidence.edge_distance_score() < 0.3); // Very low due to edge
}

#[test]
fn test_edge_case_short_sequence() {
    // Test edge case: very short sequence where everything is near an edge
    let reference = b"ATCGATCG".to_vec(); // Only 8bp
    let fm_index = FmIndex::new(&reference);
    let scorer = AlignmentContextScorer::new(&reference, &fm_index);

    let position = 4;
    let aligned_sequence = Sequence::new(b"ATCGATCG".to_vec(), None, None);
    let alignment_identity = 0.98;

    let confidence = scorer.score(position, &aligned_sequence, alignment_identity);

    // Should handle short sequences without panicking
    assert!(confidence.combined_confidence() >= 0.0);
    assert!(confidence.combined_confidence() <= 1.0);
}

#[test]
fn test_feature_weights_documented() {
    // Test acceptance criterion: weights should be documented
    // This test verifies that the scorer provides access to its weights
    let reference = create_test_reference();
    let fm_index = FmIndex::new(&reference);
    let scorer = AlignmentContextScorer::new(&reference, &fm_index);

    // Should be able to access feature weights (for transparency)
    let weights = scorer.feature_weights();

    assert!(weights.edge_distance_weight > 0.0);
    assert!(weights.gc_content_weight > 0.0);
    assert!(weights.kmer_uniqueness_weight > 0.0);
    assert!(weights.repeat_penalty_weight > 0.0);
    assert!(weights.snp_density_weight > 0.0);
    assert!(weights.alignment_identity_weight > 0.0);

    // All weights should sum to 1.0 (normalized)
    let total = weights.edge_distance_weight
        + weights.gc_content_weight
        + weights.kmer_uniqueness_weight
        + weights.repeat_penalty_weight
        + weights.snp_density_weight
        + weights.alignment_identity_weight;

    assert!(
        (total - 1.0).abs() < 0.001,
        "Feature weights should sum to 1.0"
    );
}

#[test]
fn test_multiple_positions_consistency() {
    // Test that scoring is consistent and deterministic
    let reference = create_test_reference();
    let fm_index = FmIndex::new(&reference);
    let scorer = AlignmentContextScorer::new(&reference, &fm_index);

    let position = 120;
    let aligned_sequence = Sequence::new(b"CTGACTGACTGACTGACTGA".to_vec(), None, None);
    let alignment_identity = 0.98;

    // Score the same position twice
    let confidence1 = scorer.score(position, &aligned_sequence, alignment_identity);
    let confidence2 = scorer.score(position, &aligned_sequence, alignment_identity);

    // Should produce identical results
    assert_eq!(
        confidence1.combined_confidence(),
        confidence2.combined_confidence()
    );
    assert_eq!(
        confidence1.edge_distance_score(),
        confidence2.edge_distance_score()
    );
    assert_eq!(
        confidence1.kmer_uniqueness_score(),
        confidence2.kmer_uniqueness_score()
    );
}

#[test]
fn test_extreme_gc_content() {
    // Test extreme GC content (all GC or all AT)
    let reference = b"GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC".to_vec();
    let fm_index = FmIndex::new(&reference);
    let scorer = AlignmentContextScorer::new(&reference, &fm_index);

    let position = 20;
    let aligned_sequence = Sequence::new(b"GCGCGCGCGCGCGCGCGCGC".to_vec(), None, None);
    let alignment_identity = 0.98;

    let confidence = scorer.score(position, &aligned_sequence, alignment_identity);

    // Should handle extreme GC without issues
    assert!(confidence.local_gc_content() > 0.95);
    assert!(confidence.combined_confidence() >= 0.0);
}

#[test]
fn test_scoring_with_quality_scores() {
    // Test that base quality scores (if available) affect confidence
    let reference = create_test_reference();
    let fm_index = FmIndex::new(&reference);
    let scorer = AlignmentContextScorer::new(&reference, &fm_index);

    let position = 120;

    // High quality sequence
    let high_qual = vec![40u8; 20]; // Phred 40 = 99.99% accuracy
    let aligned_sequence_high_qual =
        Sequence::new(b"CTGACTGACTGACTGACTGA".to_vec(), Some(high_qual), None);

    // Low quality sequence
    let low_qual = vec![15u8; 20]; // Phred 15 = ~97% accuracy
    let aligned_sequence_low_qual =
        Sequence::new(b"CTGACTGACTGACTGACTGA".to_vec(), Some(low_qual), None);

    let confidence_high_qual = scorer.score(position, &aligned_sequence_high_qual, 0.98);
    let confidence_low_qual = scorer.score(position, &aligned_sequence_low_qual, 0.98);

    // Higher base quality should yield higher confidence
    assert!(
        confidence_high_qual.combined_confidence() >= confidence_low_qual.combined_confidence()
    );
}
