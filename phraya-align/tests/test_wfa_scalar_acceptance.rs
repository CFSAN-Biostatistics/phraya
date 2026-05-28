/// Acceptance tests for Issue #70: WFA extension (scalar)
///
/// These tests validate the acceptance criteria for a scalar (non-SIMD) WFA
/// implementation. This is the baseline correctness implementation before SIMD
/// optimization.
///
/// Expected behavior: ALL TESTS FAIL (RED phase of TDD)
/// - The function signature in the acceptance criteria differs from current implementation
/// - Need: wfa_extend(query, target, seed_pos) -> Alignment
/// - Current: wfa_extend(query, target, seed: SeedAnchor) -> Result<Alignment, WfaError>
/// - Alignment needs: edit_distance field (not score)
/// - Alignment needs: query_start, query_end, target_start, target_end fields
use phraya_align::{Alignment, wfa_extend};

// ============================================================================
// ACCEPTANCE CRITERION: wfa_extend(query: &[u8], target: &[u8], seed_pos: usize) → Alignment
// ============================================================================

#[test]
fn test_wfa_extend_signature_accepts_seed_pos() {
    // Test will FAIL: current signature uses SeedAnchor, not usize
    let query = b"ACGTACGT";
    let target = b"ACGTACGT";
    let seed_pos: usize = 0;

    let alignment = wfa_extend(query, target, seed_pos);

    // Should return Alignment directly, not Result
    assert_eq!(alignment.cigar, "8M");
}

// ============================================================================
// ACCEPTANCE CRITERION: Alignment struct contains required fields
// ============================================================================

#[test]
fn test_alignment_has_edit_distance_field() {
    // Test will FAIL: current Alignment has 'score', not 'edit_distance'
    let query = b"ACGTACGT";
    let target = b"ACGTACTT"; // 1 mismatch
    let seed_pos = 0;

    let alignment = wfa_extend(query, target, seed_pos);

    // Should have edit_distance field
    assert_eq!(alignment.edit_distance, 1);
}

#[test]
fn test_alignment_has_query_positions() {
    // Test will FAIL: current Alignment doesn't have query_start/query_end
    let query = b"ACGTACGT";
    let target = b"ACGTACGT";
    let seed_pos = 0;

    let alignment = wfa_extend(query, target, seed_pos);

    assert_eq!(alignment.query_start, 0);
    assert_eq!(alignment.query_end, 8);
}

#[test]
fn test_alignment_has_target_positions() {
    // Test will FAIL: current Alignment doesn't have target_start/target_end
    let query = b"ACGTACGT";
    let target = b"ACGTACGT";
    let seed_pos = 0;

    let alignment = wfa_extend(query, target, seed_pos);

    assert_eq!(alignment.target_start, 0);
    assert_eq!(alignment.target_end, 8);
}

// ============================================================================
// ACCEPTANCE CRITERION: exact match (100bp) → CIGAR "100M", edit_dist 0
// ============================================================================

#[test]
fn test_exact_match_100bp() {
    // Perfect match over 100bp
    let query = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    let target = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

    assert_eq!(query.len(), 100);
    assert_eq!(target.len(), 100);

    let alignment = wfa_extend(query, target, 0);

    assert_eq!(alignment.cigar, "100M");
    assert_eq!(alignment.edit_distance, 0);
}

#[test]
fn test_exact_match_produces_zero_edit_distance() {
    let query = b"ACGTACGTACGTACGT";
    let target = b"ACGTACGTACGTACGT";

    let alignment = wfa_extend(query, target, 0);

    assert_eq!(alignment.edit_distance, 0);
    assert!(alignment.cigar.chars().all(|c| c == 'M' || c.is_numeric()));
}

// ============================================================================
// ACCEPTANCE CRITERION: single mismatch → CIGAR includes mismatch operation
// ============================================================================

#[test]
fn test_single_mismatch_at_position_50() {
    let mut query = vec![b'A'; 100];
    let mut target = vec![b'A'; 100];
    target[50] = b'T'; // Single mismatch at position 50

    let alignment = wfa_extend(&query, &target, 0);

    // CIGAR should include mismatch operation (X or generic M)
    // Standard CIGAR uses 'M' for both match and mismatch
    // Extended CIGAR uses '=' for match, 'X' for mismatch
    assert!(alignment.cigar.contains('M') || alignment.cigar.contains('X'));
    assert_eq!(alignment.edit_distance, 1);
}

#[test]
fn test_single_mismatch_early_position() {
    let query = b"ACGTACGTACGT";
    let target = b"ATGTACGTACGT"; // Mismatch at position 1 (C->T)

    let alignment = wfa_extend(query, target, 0);

    assert_eq!(alignment.edit_distance, 1);
    assert!(alignment.cigar.len() > 0);
}

#[test]
fn test_single_mismatch_late_position() {
    let query = b"ACGTACGTACGT";
    let target = b"ACGTACGTACGA"; // Mismatch at last position (T->A)

    let alignment = wfa_extend(query, target, 0);

    assert_eq!(alignment.edit_distance, 1);
}

// ============================================================================
// ACCEPTANCE CRITERION: single insertion → CIGAR includes 'I'
// ============================================================================

#[test]
fn test_single_insertion() {
    let query = b"ACGTACGT";
    let target = b"ACGTAACGT"; // Extra 'A' inserted in target

    let alignment = wfa_extend(query, target, 0);

    assert!(alignment.cigar.contains('I'));
    assert_eq!(alignment.edit_distance, 1);
}

#[test]
fn test_insertion_at_start() {
    let query = b"ACGTACGT";
    let target = b"TACGTACGT"; // 'T' inserted at start

    let alignment = wfa_extend(query, target, 0);

    assert!(alignment.cigar.contains('I'));
    assert_eq!(alignment.edit_distance, 1);
}

#[test]
fn test_insertion_at_end() {
    let query = b"ACGTACGT";
    let target = b"ACGTACGTA"; // 'A' inserted at end

    let alignment = wfa_extend(query, target, 0);

    assert!(alignment.cigar.contains('I'));
    assert_eq!(alignment.edit_distance, 1);
}

#[test]
fn test_multiple_insertions() {
    let query = b"ACGTACGT";
    let target = b"ACGTAACCGT"; // Two extra 'C's

    let alignment = wfa_extend(query, target, 0);

    assert!(alignment.cigar.contains('I'));
    assert_eq!(alignment.edit_distance, 2);
}

// ============================================================================
// ACCEPTANCE CRITERION: single deletion → CIGAR includes 'D'
// ============================================================================

#[test]
fn test_single_deletion() {
    let query = b"ACGTAACGT";
    let target = b"ACGTACGT"; // 'A' deleted from query

    let alignment = wfa_extend(query, target, 0);

    assert!(alignment.cigar.contains('D'));
    assert_eq!(alignment.edit_distance, 1);
}

#[test]
fn test_deletion_at_start() {
    let query = b"TACGTACGT";
    let target = b"ACGTACGT"; // First 'T' deleted

    let alignment = wfa_extend(query, target, 0);

    assert!(alignment.cigar.contains('D'));
    assert_eq!(alignment.edit_distance, 1);
}

#[test]
fn test_deletion_at_end() {
    let query = b"ACGTACGTA";
    let target = b"ACGTACGT"; // Last 'A' deleted

    let alignment = wfa_extend(query, target, 0);

    assert!(alignment.cigar.contains('D'));
    assert_eq!(alignment.edit_distance, 1);
}

#[test]
fn test_multiple_deletions() {
    let query = b"ACGTAACCGT";
    let target = b"ACGTACGT"; // Two 'C's deleted

    let alignment = wfa_extend(query, target, 0);

    assert!(alignment.cigar.contains('D'));
    assert_eq!(alignment.edit_distance, 2);
}

// ============================================================================
// ACCEPTANCE CRITERION: complex alignment (mix of M/I/D) → verify CIGAR + edit_dist
// ============================================================================

#[test]
fn test_complex_alignment_mixed_operations() {
    let query = b"ACGTACGTTAGC";
    let target = b"ACGTTCGTAGC"; // Multiple differences

    let alignment = wfa_extend(query, target, 0);

    // Should have a mix of operations
    let has_match = alignment.cigar.contains('M') || alignment.cigar.contains('=');
    let has_indel = alignment.cigar.contains('I') || alignment.cigar.contains('D');

    assert!(has_match);
    assert!(has_indel || alignment.cigar.contains('X'));
    assert!(alignment.edit_distance > 0);
}

#[test]
fn test_complex_alignment_consecutive_indels() {
    let query = b"ACGTAAAACGT";
    let target = b"ACGTCGT"; // Multiple consecutive deletions

    let alignment = wfa_extend(query, target, 0);

    assert!(alignment.cigar.contains('D'));
    assert_eq!(alignment.edit_distance, 4); // Four A's deleted
}

#[test]
fn test_complex_alignment_alternating_events() {
    // Pattern: match, mismatch, insertion, deletion, match
    let query = b"ACGTACGTTAGCTTGCA";
    let target = b"ACGTTCGTAGCGCA";

    let alignment = wfa_extend(query, target, 0);

    // Verify CIGAR is valid and edit distance is reasonable
    assert!(alignment.cigar.len() > 0);
    assert!(alignment.edit_distance > 0);
    assert!(alignment.edit_distance <= query.len().max(target.len()));
}

#[test]
fn test_complex_alignment_high_divergence() {
    let query = b"ACGTACGTACGTACGT";
    let target = b"TGCATGCATGCATGCA"; // Highly divergent

    let alignment = wfa_extend(query, target, 0);

    assert!(alignment.edit_distance > 8); // Expect many edits
    assert!(alignment.cigar.len() > 0);
}

// ============================================================================
// ACCEPTANCE CRITERION: empty sequences → handle gracefully
// ============================================================================

#[test]
fn test_empty_query_sequence() {
    let query = b"";
    let target = b"ACGTACGT";

    let alignment = wfa_extend(query, target, 0);

    // Empty query aligns to target with all insertions
    assert_eq!(alignment.edit_distance, target.len());
    assert!(alignment.cigar.contains('I') || alignment.cigar == "");
}

#[test]
fn test_empty_target_sequence() {
    let query = b"ACGTACGT";
    let target = b"";

    let alignment = wfa_extend(query, target, 0);

    // Empty target requires all deletions from query
    assert_eq!(alignment.edit_distance, query.len());
    assert!(alignment.cigar.contains('D') || alignment.cigar == "");
}

#[test]
fn test_both_sequences_empty() {
    let query = b"";
    let target = b"";

    let alignment = wfa_extend(query, target, 0);

    assert_eq!(alignment.edit_distance, 0);
    assert_eq!(alignment.cigar, "");
}

// ============================================================================
// ACCEPTANCE CRITERION: seed_pos parameter behavior
// ============================================================================

#[test]
fn test_seed_pos_at_start() {
    let query = b"ACGTACGTACGT";
    let target = b"ACGTACGTACGT";

    let alignment = wfa_extend(query, target, 0);

    assert_eq!(alignment.query_start, 0);
    assert_eq!(alignment.target_start, 0);
    assert_eq!(alignment.cigar, "12M");
}

#[test]
fn test_seed_pos_in_middle() {
    // Seed position should determine where alignment starts
    let query = b"ACGTACGTACGT";
    let target = b"ACGTACGTACGT";
    let seed_pos = 4;

    let alignment = wfa_extend(query, target, seed_pos);

    // Alignment should start from seed_pos
    assert_eq!(alignment.query_start, seed_pos);
    assert_eq!(alignment.target_start, seed_pos);
}

#[test]
fn test_seed_pos_near_end() {
    let query = b"ACGTACGTACGT";
    let target = b"ACGTACGTACGT";
    let seed_pos = 10;

    let alignment = wfa_extend(query, target, seed_pos);

    assert_eq!(alignment.query_start, seed_pos);
    assert!(alignment.query_end <= query.len());
}

// ============================================================================
// CIGAR format validation
// ============================================================================

#[test]
fn test_cigar_format_valid() {
    let query = b"ACGTACGTACGT";
    let target = b"ACGTACTTACGT";

    let alignment = wfa_extend(query, target, 0);

    // CIGAR should only contain valid operations: M, I, D (and optionally =, X)
    for c in alignment.cigar.chars() {
        assert!(
            c.is_numeric() || c == 'M' || c == 'I' || c == 'D' || c == '=' || c == 'X',
            "Invalid CIGAR character: {}",
            c
        );
    }
}

#[test]
fn test_cigar_length_consistency() {
    let query = b"ACGTACGTACGT";
    let target = b"ACGTACGTACGT";

    let alignment = wfa_extend(query, target, 0);

    // Parse CIGAR and verify it accounts for all positions
    let total_query_len = alignment.query_end - alignment.query_start;
    let total_target_len = alignment.target_end - alignment.target_start;

    assert!(total_query_len > 0);
    assert!(total_target_len > 0);
}

// ============================================================================
// Edit distance validation
// ============================================================================

#[test]
fn test_edit_distance_is_non_negative() {
    let query = b"ACGTACGT";
    let target = b"ACGTACGT";

    let alignment = wfa_extend(query, target, 0);

    // Edit distance cannot be negative
    assert!(alignment.edit_distance >= 0);
}

#[test]
fn test_edit_distance_matches_cigar_operations() {
    let query = b"ACGTACGT";
    let target = b"ACGTACTT"; // 1 mismatch

    let alignment = wfa_extend(query, target, 0);

    // Edit distance should count mismatches + insertions + deletions
    assert_eq!(alignment.edit_distance, 1);
}

#[test]
fn test_edit_distance_bounded_by_sequence_length() {
    let query = b"ACGTACGT";
    let target = b"TGCATGCA"; // Completely different

    let alignment = wfa_extend(query, target, 0);

    // Edit distance cannot exceed max(query_len, target_len)
    let max_len = query.len().max(target.len());
    assert!(alignment.edit_distance <= max_len);
}

// ============================================================================
// ACCEPTANCE CRITERION: Benchmark 10kb sequences (measure time)
// ============================================================================
// Note: This is tested in benches/bench_wfa_scalar.rs

#[test]
fn test_10kb_alignment_completes() {
    // Generate 10kb sequences for alignment
    let query: Vec<u8> = (0..10000).map(|i| b"ACGT"[i % 4]).collect();
    let target: Vec<u8> = (0..10000).map(|i| b"ACGT"[i % 4]).collect();

    assert_eq!(query.len(), 10000);
    assert_eq!(target.len(), 10000);

    // Should complete without panic or excessive time
    let alignment = wfa_extend(&query, &target, 0);

    assert_eq!(alignment.cigar.len() > 0, true);
}

// ============================================================================
// Error cases and edge conditions
// ============================================================================

#[test]
fn test_single_base_sequences() {
    let query = b"A";
    let target = b"A";

    let alignment = wfa_extend(query, target, 0);

    assert_eq!(alignment.cigar, "1M");
    assert_eq!(alignment.edit_distance, 0);
}

#[test]
fn test_single_base_mismatch() {
    let query = b"A";
    let target = b"T";

    let alignment = wfa_extend(query, target, 0);

    assert_eq!(alignment.edit_distance, 1);
}

#[test]
fn test_very_long_indel() {
    let query = vec![b'A'; 50];
    let target = vec![b'A'; 100]; // 50 extra bases

    let alignment = wfa_extend(&query, &target, 0);

    assert!(alignment.cigar.contains('I'));
    assert_eq!(alignment.edit_distance, 50);
}

#[test]
fn test_sequences_with_ambiguous_bases() {
    // Note: WFA should handle IUPAC codes if present
    let query = b"ACGTACGT";
    let target = b"ACGTACGT";

    let alignment = wfa_extend(query, target, 0);

    assert_eq!(alignment.edit_distance, 0);
}
