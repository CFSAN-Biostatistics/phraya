//! Acceptance tests for k-mer sketching integration (issue #63)
//!
//! This test suite covers the integration of simd-minimizers crate with
//! phraya-index for k-mer sketching of Sequence objects.
//!
//! RED PHASE: These tests are expected to FAIL until implementation is complete.

use phraya_core::types::Sequence;
use phraya_index::{sketch, sketch_default, Sketch};

// ============================================================================
// HAPPY PATH: Basic sketching functionality
// ============================================================================

#[test]
fn test_sketch_basic_sequence() {
    // Should be able to sketch a simple bacterial sequence with custom parameters
    let bases = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let seq = Sequence::new(bases, None, "test_seq".to_string(), None);

    let sketch = sketch(&seq, 21, 11);

    // Verify sketch has expected properties
    assert_eq!(sketch.k(), 21);
    assert_eq!(sketch.w(), 11);
    assert!(!sketch.is_empty());
}

#[test]
fn test_sketch_default_parameters() {
    // Should sketch with default k=21, w=11
    let bases = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let seq = Sequence::new(bases, None, "default_test".to_string(), None);

    let sketch = sketch_default(&seq);

    assert_eq!(sketch.k(), 21, "Default k should be 21");
    assert_eq!(sketch.w(), 11, "Default w should be 11");
}

#[test]
fn test_sketch_deterministic() {
    // Same sequence should produce identical sketches every time
    let bases = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let seq = Sequence::new(bases, None, "deterministic_test".to_string(), None);

    let sketches: Vec<Sketch> = (0..10)
        .map(|_| sketch_default(&seq))
        .collect();

    // All 10 sketches should be identical
    for i in 1..sketches.len() {
        assert_eq!(
            sketches[0], sketches[i],
            "Sketch {} differs from sketch 0",
            i
        );
    }
}

#[test]
fn test_sketch_different_sequences_differ() {
    // Different sequences should produce different sketches
    let bases1 = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let bases2 = b"TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA".to_vec();

    let seq1 = Sequence::new(bases1, None, "seq1".to_string(), None);
    let seq2 = Sequence::new(bases2, None, "seq2".to_string(), None);

    let sketch1 = sketch_default(&seq1);
    let sketch2 = sketch_default(&seq2);

    assert_ne!(sketch1, sketch2, "Different sequences should produce different sketches");
}

#[test]
fn test_sketch_with_quality_scores() {
    // Sketching should work with sequences that have quality scores
    let bases = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let quality = vec![30u8; bases.len()];
    let seq = Sequence::new(bases, Some(quality), "qual_test".to_string(), None);

    let sketch = sketch_default(&seq);

    // Quality scores should not affect the sketch (only bases matter)
    assert!(!sketch.is_empty());
}

#[test]
fn test_sketch_bacterial_genome_size() {
    // Test with a realistic bacterial genome fragment (100kb)
    let size = 100_000;
    let bases: Vec<u8> = (0..size)
        .map(|i| match i % 4 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        })
        .collect();

    let seq = Sequence::new(bases, None, "bacterial_fragment".to_string(), None);
    let sketch = sketch_default(&seq);

    // Should successfully sketch a large sequence
    assert!(!sketch.is_empty());
    // Sketch should be smaller than the original sequence (compression)
    assert!(sketch.len() < size);
}

// ============================================================================
// EDGE CASES: Empty and short sequences
// ============================================================================

#[test]
fn test_sketch_empty_sequence() {
    // Empty sequence should return empty sketch
    let bases = Vec::new();
    let seq = Sequence::new(bases, None, "empty".to_string(), None);

    let sketch = sketch_default(&seq);

    assert!(sketch.is_empty(), "Empty sequence should produce empty sketch");
    assert_eq!(sketch.len(), 0);
}

#[test]
fn test_sketch_sequence_shorter_than_k() {
    // Sequence shorter than k-mer size should handle gracefully
    let bases = b"ACGT".to_vec(); // Only 4 bases, k=21
    let seq = Sequence::new(bases, None, "short".to_string(), None);

    let sketch = sketch_default(&seq);

    // Should return empty or handle gracefully (not panic)
    assert!(sketch.len() == 0 || sketch.len() <= 1);
}

#[test]
fn test_sketch_sequence_exactly_k_bases() {
    // Sequence exactly k bases long
    let bases = b"ACGTACGTACGTACGTACGTA".to_vec(); // Exactly 21 bases
    let seq = Sequence::new(bases, None, "exactly_k".to_string(), None);

    let sketch = sketch_default(&seq);

    // Should successfully sketch (single k-mer)
    assert!(sketch.len() <= 1);
}

#[test]
fn test_sketch_sequence_just_over_k() {
    // Sequence just slightly longer than k
    let bases = b"ACGTACGTACGTACGTACGTAAA".to_vec(); // 24 bases (k+3)
    let seq = Sequence::new(bases, None, "k_plus_3".to_string(), None);

    let sketch = sketch_default(&seq);

    assert!(!sketch.is_empty());
}

#[test]
fn test_sketch_single_base_sequence() {
    // Single base sequence
    let bases = b"A".to_vec();
    let seq = Sequence::new(bases, None, "single".to_string(), None);

    let sketch = sketch_default(&seq);

    assert!(sketch.is_empty());
}

// ============================================================================
// PARAMETER VARIATIONS: Different k and w values
// ============================================================================

#[test]
fn test_sketch_small_k_and_w() {
    // Test with smaller k and w values
    let bases = b"ACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let seq = Sequence::new(bases, None, "small_params".to_string(), None);

    let sketch = sketch(&seq, 7, 3);

    assert_eq!(sketch.k(), 7);
    assert_eq!(sketch.w(), 3);
    assert!(!sketch.is_empty());
}

#[test]
fn test_sketch_large_k_and_w() {
    // Test with larger k and w values
    let bases = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let seq = Sequence::new(bases, None, "large_params".to_string(), None);

    let sketch = sketch(&seq, 31, 15);

    assert_eq!(sketch.k(), 31);
    assert_eq!(sketch.w(), 15);
}

#[test]
fn test_sketch_k_equals_w() {
    // Window size equal to k-mer size (minimum valid case)
    let bases = b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let seq = Sequence::new(bases, None, "k_eq_w".to_string(), None);

    let sketch = sketch(&seq, 11, 11);

    assert_eq!(sketch.k(), 11);
    assert_eq!(sketch.w(), 11);
}

#[test]
fn test_sketch_standard_bacterial_params() {
    // Standard parameters for bacterial genomics (k=21, w=11)
    let bases = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let seq = Sequence::new(bases, None, "bacterial_standard".to_string(), None);

    let sketch1 = sketch(&seq, 21, 11);
    let sketch2 = sketch_default(&seq);

    assert_eq!(sketch1, sketch2, "sketch(seq, 21, 11) should equal sketch_default(seq)");
}

// ============================================================================
// SEQUENCE PATTERNS: Homopolymers, repeats, complex
// ============================================================================

#[test]
fn test_sketch_homopolymer_sequence() {
    // All identical bases
    let bases = vec![b'A'; 1000];
    let seq = Sequence::new(bases, None, "homopolymer".to_string(), None);

    let sketch = sketch_default(&seq);

    // Homopolymer should produce a very small sketch (all k-mers identical)
    assert!(!sketch.is_empty());
    // Should be heavily compressed
    assert!(sketch.len() < 100);
}

#[test]
fn test_sketch_alternating_pattern() {
    // Alternating bases pattern
    let bases: Vec<u8> = (0..1000)
        .map(|i| if i % 2 == 0 { b'A' } else { b'T' })
        .collect();
    let seq = Sequence::new(bases, None, "alternating".to_string(), None);

    let sketch = sketch_default(&seq);

    assert!(!sketch.is_empty());
}

#[test]
fn test_sketch_tandem_repeat() {
    // Tandem repeat sequence (same pattern repeated)
    let pattern = b"ACGT";
    let bases: Vec<u8> = pattern.iter().cycle().take(1000).copied().collect();
    let seq = Sequence::new(bases, None, "tandem_repeat".to_string(), None);

    let sketch = sketch_default(&seq);

    // Tandem repeats should produce moderately sized sketch
    assert!(!sketch.is_empty());
}

#[test]
fn test_sketch_random_like_sequence() {
    // Pseudo-random bacterial-like sequence
    let bases: Vec<u8> = (0..10000)
        .map(|i: u32| {
            let hash = i.wrapping_mul(2654435761);
            match hash % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            }
        })
        .collect();
    let seq = Sequence::new(bases, None, "random_like".to_string(), None);

    let sketch = sketch_default(&seq);

    // Random sequence should produce sketch with good coverage
    assert!(!sketch.is_empty());
    assert!(sketch.len() > 100);
}

#[test]
fn test_sketch_gc_rich_sequence() {
    // GC-rich sequence (typical for some bacterial genomes)
    let bases: Vec<u8> = (0..5000)
        .map(|i| if i % 2 == 0 { b'G' } else { b'C' })
        .collect();
    let seq = Sequence::new(bases, None, "gc_rich".to_string(), None);

    let sketch = sketch_default(&seq);

    assert!(!sketch.is_empty());
}

#[test]
fn test_sketch_at_rich_sequence() {
    // AT-rich sequence
    let bases: Vec<u8> = (0..5000)
        .map(|i| if i % 2 == 0 { b'A' } else { b'T' })
        .collect();
    let seq = Sequence::new(bases, None, "at_rich".to_string(), None);

    let sketch = sketch_default(&seq);

    assert!(!sketch.is_empty());
}

// ============================================================================
// SKETCH PROPERTIES: API and structure
// ============================================================================

#[test]
fn test_sketch_has_k_method() {
    // Sketch should expose k parameter
    let bases = b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let seq = Sequence::new(bases, None, "k_method_test".to_string(), None);

    let sketch = sketch(&seq, 15, 7);

    assert_eq!(sketch.k(), 15);
}

#[test]
fn test_sketch_has_w_method() {
    // Sketch should expose w parameter
    let bases = b"ACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let seq = Sequence::new(bases, None, "w_method_test".to_string(), None);

    let sketch = sketch(&seq, 13, 5);

    assert_eq!(sketch.w(), 5);
}

#[test]
fn test_sketch_has_len_method() {
    // Sketch should have len() method
    let bases = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let seq = Sequence::new(bases, None, "len_method_test".to_string(), None);

    let sketch = sketch_default(&seq);

    let length = sketch.len();
    assert!(length > 0);
}

#[test]
fn test_sketch_has_is_empty_method() {
    // Sketch should have is_empty() method
    let bases = Vec::new();
    let seq = Sequence::new(bases, None, "empty_method_test".to_string(), None);

    let sketch = sketch_default(&seq);

    assert!(sketch.is_empty());
}

#[test]
fn test_sketch_implements_clone() {
    // Sketch should be cloneable
    let bases = b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let seq = Sequence::new(bases, None, "clone_test".to_string(), None);

    let sketch1 = sketch_default(&seq);
    let sketch2 = sketch1.clone();

    assert_eq!(sketch1, sketch2);
}

#[test]
fn test_sketch_implements_debug() {
    // Sketch should implement Debug for logging
    let bases = b"ACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let seq = Sequence::new(bases, None, "debug_test".to_string(), None);

    let sketch = sketch_default(&seq);

    let debug_str = format!("{:?}", sketch);
    assert!(!debug_str.is_empty());
}

#[test]
fn test_sketch_implements_eq() {
    // Sketch should implement equality comparison
    let bases = b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let seq = Sequence::new(bases, None, "eq_test".to_string(), None);

    let sketch1 = sketch_default(&seq);
    let sketch2 = sketch_default(&seq);

    assert_eq!(sketch1, sketch2);
}

// ============================================================================
// IDENTICAL SEQUENCES: Various representations
// ============================================================================

#[test]
fn test_identical_sequences_identical_sketches() {
    // Two sequences with identical bases should produce identical sketches
    let bases1 = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let bases2 = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();

    let seq1 = Sequence::new(bases1, None, "seq1".to_string(), None);
    let seq2 = Sequence::new(bases2, None, "seq2".to_string(), None);

    let sketch1 = sketch_default(&seq1);
    let sketch2 = sketch_default(&seq2);

    assert_eq!(sketch1, sketch2);
}

#[test]
fn test_identical_with_different_metadata() {
    // Same bases, different IDs/descriptions should produce same sketch
    let bases = b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();

    let seq1 = Sequence::new(
        bases.clone(),
        None,
        "seq_id_1".to_string(),
        Some("Description 1".to_string()),
    );
    let seq2 = Sequence::new(
        bases.clone(),
        None,
        "seq_id_2".to_string(),
        Some("Description 2".to_string()),
    );

    let sketch1 = sketch_default(&seq1);
    let sketch2 = sketch_default(&seq2);

    assert_eq!(sketch1, sketch2, "Sketches should be identical regardless of metadata");
}

#[test]
fn test_identical_with_different_quality() {
    // Same bases, different quality scores should produce same sketch
    let bases = b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let quality1 = vec![20u8; bases.len()];
    let quality2 = vec![40u8; bases.len()];

    let seq1 = Sequence::new(bases.clone(), Some(quality1), "seq1".to_string(), None);
    let seq2 = Sequence::new(bases.clone(), Some(quality2), "seq2".to_string(), None);

    let sketch1 = sketch_default(&seq1);
    let sketch2 = sketch_default(&seq2);

    assert_eq!(sketch1, sketch2, "Quality scores should not affect sketch");
}

// ============================================================================
// PERFORMANCE: Benchmark tests
// ============================================================================

#[test]
fn test_performance_5mbp_ecoli_genome() {
    use std::time::Instant;

    // Sketch a 5Mbp E. coli-like genome
    let size = 5_000_000;
    let bases: Vec<u8> = (0..size)
        .map(|i: u32| {
            // Use a pseudo-random pattern to simulate bacterial genome
            let hash = i.wrapping_mul(2654435761);
            match hash % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            }
        })
        .collect();

    let seq = Sequence::new(bases, None, "ecoli_5mbp".to_string(), None);

    let start = Instant::now();
    let sketch = sketch_default(&seq);
    let elapsed = start.elapsed();

    // Should complete in reasonable time (< 5 seconds for 5Mbp)
    assert!(
        elapsed.as_secs() < 5,
        "Sketching 5Mbp took {:?}, expected < 5s",
        elapsed
    );

    // Sanity check: sketch should be non-empty and compressed
    assert!(!sketch.is_empty());
    assert!(sketch.len() < size / 10, "Sketch should be compressed");
}

#[test]
fn test_performance_1mbp_sequence() {
    use std::time::Instant;

    // Benchmark 1Mbp sequence
    let size = 1_000_000;
    let bases: Vec<u8> = (0..size)
        .map(|i| match i % 4 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        })
        .collect();

    let seq = Sequence::new(bases, None, "1mbp".to_string(), None);

    let start = Instant::now();
    let sketch = sketch_default(&seq);
    let elapsed = start.elapsed();

    // Should be very fast (< 1 second)
    assert!(
        elapsed.as_secs() < 1,
        "Sketching 1Mbp took {:?}, expected < 1s",
        elapsed
    );

    assert!(!sketch.is_empty());
}

#[test]
fn test_performance_multiple_sketches() {
    use std::time::Instant;

    // Test performance of sketching multiple sequences
    let size = 100_000;
    let bases: Vec<u8> = (0..size).map(|i| match i % 4 {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        _ => b'T',
    }).collect();

    let sequences: Vec<Sequence> = (0..10)
        .map(|i| Sequence::new(bases.clone(), None, format!("seq_{}", i), None))
        .collect();

    let start = Instant::now();
    for seq in &sequences {
        let _sketch = sketch_default(seq);
    }
    let elapsed = start.elapsed();

    // Should sketch 10x100kb sequences in < 1 second total
    assert!(
        elapsed.as_secs() < 1,
        "Sketching 10x100kb took {:?}, expected < 1s",
        elapsed
    );
}

// ============================================================================
// SIMD PLATFORM DETECTION: Verify automatic dispatch
// ============================================================================

#[test]
fn test_simd_dispatch_works() {
    // This test verifies that simd-minimizers automatic dispatch works
    // by ensuring sketching succeeds on any platform (AVX2/SSE4.2/NEON/scalar)
    let bases = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let seq = Sequence::new(bases, None, "simd_test".to_string(), None);

    let sketch = sketch_default(&seq);

    // Should work regardless of platform
    assert!(!sketch.is_empty());
}

#[test]
fn test_simd_determinism_across_calls() {
    // SIMD implementations should give deterministic results
    let bases = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let seq = Sequence::new(bases, None, "simd_determinism".to_string(), None);

    // Sketch multiple times
    let sketches: Vec<Sketch> = (0..5)
        .map(|_| sketch_default(&seq))
        .collect();

    // All should be identical
    for i in 1..sketches.len() {
        assert_eq!(sketches[0], sketches[i]);
    }
}

// ============================================================================
// ERROR HANDLING: Graceful degradation
// ============================================================================

#[test]
fn test_sketch_handles_ambiguous_bases() {
    // Sequence with ambiguous bases (N, R, Y, etc.)
    // Should handle gracefully or skip ambiguous bases
    let bases = b"ACGTNNNACGTRYACGT".to_vec();
    let seq = Sequence::new(bases, None, "ambiguous".to_string(), None);

    // Should not panic, might produce smaller sketch
    let sketch = sketch_default(&seq);

    // May be empty or have fewer minimizers, but should not crash
    assert!(sketch.len() >= 0);
}

#[test]
fn test_sketch_handles_lowercase_bases() {
    // Lowercase bases should be handled (DNA is case-insensitive)
    let bases = b"acgtacgtacgtacgtacgtacgtacgtacgtacgt".to_vec();
    let seq = Sequence::new(bases, None, "lowercase".to_string(), None);

    let sketch = sketch_default(&seq);

    assert!(!sketch.is_empty());
}

#[test]
fn test_sketch_mixed_case_equals_uppercase() {
    // Mixed case should produce same sketch as uppercase
    let bases_upper = b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let bases_mixed = b"AcGtAcGtAcGtAcGtAcGtAcGtAcGtAcGt".to_vec();

    let seq_upper = Sequence::new(bases_upper, None, "upper".to_string(), None);
    let seq_mixed = Sequence::new(bases_mixed, None, "mixed".to_string(), None);

    let sketch_upper = sketch_default(&seq_upper);
    let sketch_mixed = sketch_default(&seq_mixed);

    assert_eq!(sketch_upper, sketch_mixed, "Case should not affect sketch");
}

// ============================================================================
// SKETCH COMPARISON: Finding shared minimizers
// ============================================================================

#[test]
fn test_shared_minimizers_identical_sequences() {
    // Two identical sequences should share all minimizers
    let bases = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let seq1 = Sequence::new(bases.clone(), None, "seq1".to_string(), None);
    let seq2 = Sequence::new(bases.clone(), None, "seq2".to_string(), None);

    let sketch1 = sketch_default(&seq1);
    let sketch2 = sketch_default(&seq2);

    // If Sketch has a method to find shared minimizers, test it
    // For now, just verify sketches are equal
    assert_eq!(sketch1, sketch2);
}

#[test]
fn test_shared_minimizers_partial_overlap() {
    // Two sequences with partial overlap
    let bases1 = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let bases2 = b"ACGTACGTACGTACGTACGTTTTTTTTTTTTTTTTTTTTTTTTTT".to_vec();

    let seq1 = Sequence::new(bases1, None, "overlap1".to_string(), None);
    let seq2 = Sequence::new(bases2, None, "overlap2".to_string(), None);

    let sketch1 = sketch_default(&seq1);
    let sketch2 = sketch_default(&seq2);

    // Sketches should differ (partial overlap)
    assert_ne!(sketch1, sketch2);
}

#[test]
fn test_shared_minimizers_no_overlap() {
    // Completely different sequences
    let bases1 = vec![b'A'; 1000];
    let bases2 = vec![b'C'; 1000];

    let seq1 = Sequence::new(bases1, None, "all_a".to_string(), None);
    let seq2 = Sequence::new(bases2, None, "all_c".to_string(), None);

    let sketch1 = sketch_default(&seq1);
    let sketch2 = sketch_default(&seq2);

    assert_ne!(sketch1, sketch2);
}
