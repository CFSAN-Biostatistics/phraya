// Confidence feature computation for variant calling
// Functions to compute per-position features used for confidence scoring

#[cfg(test)]
mod tests {
    use super::*;

    // Tests for distance_to_edge
    #[test]
    fn test_distance_to_edge_at_start() {
        assert_eq!(distance_to_edge(0, 100), 0);
    }

    #[test]
    fn test_distance_to_edge_at_end() {
        assert_eq!(distance_to_edge(99, 100), 0);
    }

    #[test]
    fn test_distance_to_edge_at_middle() {
        assert_eq!(distance_to_edge(50, 100), 50);
    }

    #[test]
    fn test_distance_to_edge_closer_to_start() {
        assert_eq!(distance_to_edge(10, 100), 10);
    }

    #[test]
    fn test_distance_to_edge_closer_to_end() {
        assert_eq!(distance_to_edge(90, 100), 9);
    }

    #[test]
    fn test_distance_to_edge_single_base_sequence() {
        assert_eq!(distance_to_edge(0, 1), 0);
    }

    #[test]
    fn test_distance_to_edge_two_base_sequence() {
        assert_eq!(distance_to_edge(0, 2), 0);
        assert_eq!(distance_to_edge(1, 2), 0);
    }

    #[test]
    fn test_distance_to_edge_short_sequence_middle() {
        // In a 5bp sequence, position 2 (0-indexed) is equidistant: 2 from start, 2 from end
        assert_eq!(distance_to_edge(2, 5), 2);
    }

    // Tests for local_gc_content
    #[test]
    fn test_local_gc_content_all_gc() {
        let seq = b"GCGCGCGCGCGCGCGCGCGC";
        assert_eq!(local_gc_content(seq, 10, 10), 1.0);
    }

    #[test]
    fn test_local_gc_content_all_at() {
        let seq = b"ATATATATATATATATATATAT";
        assert_eq!(local_gc_content(seq, 10, 10), 0.0);
    }

    #[test]
    fn test_local_gc_content_half_gc() {
        let seq = b"ATATATGCGCGCATATGCGCGC";
        // Window of 20bp centered at position 10
        let gc = local_gc_content(seq, 10, 10);
        assert!((gc - 0.5).abs() < 0.01);
    }

    #[test]
    fn test_local_gc_content_at_start_boundary() {
        let seq = b"GCGCATATATAT";
        // Position 0, window_size 10 - should only count available bases
        let gc = local_gc_content(seq, 0, 10);
        // First 11 bases (pos 0 +/- 10, but only goes to position 10): GCGCATATATA
        // 4 GC out of 11 = ~0.36
        assert!((gc - 4.0 / 11.0).abs() < 0.01);
    }

    #[test]
    fn test_local_gc_content_at_end_boundary() {
        let seq = b"ATATATATATGCGC";
        let pos = seq.len() - 1;
        // Position 13 (last), window_size 10 - should only count available bases
        let gc = local_gc_content(seq, pos, 10);
        // Positions 3-13 (11 bases): TATATATGCGC = 4 GC out of 11
        assert!((gc - 4.0 / 11.0).abs() < 0.01);
    }

    #[test]
    fn test_local_gc_content_short_sequence() {
        let seq = b"GCGC";
        assert_eq!(local_gc_content(seq, 1, 20), 1.0);
    }

    #[test]
    fn test_local_gc_content_single_base() {
        let seq = b"G";
        assert_eq!(local_gc_content(seq, 0, 5), 1.0);
    }

    #[test]
    fn test_local_gc_content_window_size_zero() {
        let seq = b"ATGC";
        // Window size 0 means only the base at position
        assert_eq!(local_gc_content(seq, 2, 0), 1.0); // G
        assert_eq!(local_gc_content(seq, 0, 0), 0.0); // A
    }

    #[test]
    fn test_local_gc_content_default_20bp_window() {
        // 10 AT + 10 GC + 10 GC + 10 AT = 40 bases
        let seq = b"ATATATATATGCGCGCGCGCGCGCGCGCGCATATATATAT";
        // Position 20 (middle), window ±20 means positions 0-40 (41 bases total, but seq is 40)
        // 20 GC bases in the middle section
        let gc = local_gc_content(seq, 20, 20);
        assert!((gc - 0.5).abs() < 0.01);
    }

    #[test]
    fn test_local_gc_content_lowercase() {
        let seq = b"atgcgcgcgcgc";
        // Should handle lowercase letters
        let gc = local_gc_content(seq, 6, 5);
        assert!(gc > 0.5);
    }

    // Tests for in_homopolymer
    #[test]
    fn test_in_homopolymer_a_run_length_4() {
        let seq = b"ATCGAAAA";
        assert!(in_homopolymer(seq, 4, 4)); // First A of run
        assert!(in_homopolymer(seq, 5, 4)); // Second A
        assert!(in_homopolymer(seq, 6, 4)); // Third A
        assert!(in_homopolymer(seq, 7, 4)); // Fourth A
    }

    #[test]
    fn test_in_homopolymer_t_run_length_5() {
        let seq = b"GCATTTTTGC";
        assert!(in_homopolymer(seq, 3, 4)); // First T of 5-T run
        assert!(in_homopolymer(seq, 4, 4));
        assert!(in_homopolymer(seq, 5, 4));
        assert!(in_homopolymer(seq, 6, 4));
        assert!(in_homopolymer(seq, 7, 4)); // Fifth T
    }

    #[test]
    fn test_in_homopolymer_g_run_length_6() {
        let seq = b"ATGGGGGGCA";
        assert!(in_homopolymer(seq, 2, 4)); // First G of 6-G run
        assert!(in_homopolymer(seq, 7, 4)); // Last G of run
    }

    #[test]
    fn test_in_homopolymer_c_run_length_4() {
        let seq = b"ATCCCCGAT";
        assert!(in_homopolymer(seq, 2, 4)); // First C
        assert!(in_homopolymer(seq, 5, 4)); // Last C
    }

    #[test]
    fn test_not_in_homopolymer_length_3() {
        let seq = b"ATCGAAACG";
        assert!(!in_homopolymer(seq, 4, 4)); // Only 3 As
        assert!(!in_homopolymer(seq, 5, 4));
        assert!(!in_homopolymer(seq, 6, 4));
    }

    #[test]
    fn test_not_in_homopolymer_single_base() {
        let seq = b"ATCGATCG";
        assert!(!in_homopolymer(seq, 0, 4));
        assert!(!in_homopolymer(seq, 4, 4));
    }

    #[test]
    fn test_not_in_homopolymer_alternating() {
        let seq = b"ATATATAT";
        for pos in 0..seq.len() {
            assert!(!in_homopolymer(seq, pos, 4));
        }
    }

    #[test]
    fn test_in_homopolymer_at_sequence_start() {
        let seq = b"AAAAATGC";
        assert!(in_homopolymer(seq, 0, 4));
        assert!(in_homopolymer(seq, 4, 4)); // Last A of run
    }

    #[test]
    fn test_in_homopolymer_at_sequence_end() {
        let seq = b"ATGCGGGG";
        assert!(in_homopolymer(seq, 4, 4));
        assert!(in_homopolymer(seq, 7, 4)); // Last G
    }

    #[test]
    fn test_in_homopolymer_entire_sequence() {
        let seq = b"AAAAAAA";
        for pos in 0..seq.len() {
            assert!(in_homopolymer(seq, pos, 4));
        }
    }

    #[test]
    fn test_in_homopolymer_min_length_5() {
        let seq = b"ATCGAAAAAGC";
        assert!(!in_homopolymer(seq, 4, 5)); // Run of 5, not 6
        assert!(in_homopolymer(seq, 4, 5)); // Should be true for exactly 5
    }

    #[test]
    fn test_in_homopolymer_min_length_6() {
        let seq = b"ATCGAAAAAAGC";
        // 6 As
        assert!(in_homopolymer(seq, 4, 6));
        assert!(in_homopolymer(seq, 9, 6));
    }

    #[test]
    fn test_in_homopolymer_min_length_1() {
        let seq = b"ATCG";
        // Every position is in a "homopolymer" of length >= 1
        assert!(in_homopolymer(seq, 0, 1));
        assert!(in_homopolymer(seq, 1, 1));
        assert!(in_homopolymer(seq, 2, 1));
        assert!(in_homopolymer(seq, 3, 1));
    }

    #[test]
    fn test_in_homopolymer_lowercase() {
        let seq = b"atcgaaaaa";
        // Should handle lowercase
        assert!(in_homopolymer(seq, 4, 4));
    }

    #[test]
    fn test_in_homopolymer_mixed_case() {
        let seq = b"ATCGaAaAa";
        // Mixed case should still count as homopolymer if same letter
        assert!(in_homopolymer(seq, 4, 4));
    }

    #[test]
    fn test_in_homopolymer_adjacent_runs() {
        let seq = b"AAAAGGGG";
        // Two adjacent runs of different bases
        assert!(in_homopolymer(seq, 0, 4)); // In A run
        assert!(in_homopolymer(seq, 3, 4)); // Last A
        assert!(!in_homopolymer(seq, 3, 5)); // A run is only 4
        assert!(in_homopolymer(seq, 4, 4)); // First G
        assert!(in_homopolymer(seq, 7, 4)); // Last G
    }

    #[test]
    fn test_in_homopolymer_short_sequence() {
        let seq = b"AAA";
        assert!(!in_homopolymer(seq, 0, 4)); // Sequence too short
        assert!(!in_homopolymer(seq, 2, 4));
    }

    // Integration/manual verification tests
    #[test]
    fn test_manual_gc_calculation_verification() {
        // Manual verification: ATGCGCATATGC
        // G=3, C=3, total=12, GC=6/12=0.5
        let seq = b"ATGCGCATATGC";
        let gc = local_gc_content(seq, 6, 6); // Should cover whole sequence
        assert!((gc - 0.5).abs() < 0.001);
    }

    #[test]
    fn test_manual_gc_with_exact_window() {
        // AAAAGCGCAAAA (12 bases)
        // Window at position 6 (G), window_size=2 means positions 4-8: AGCGC
        // GC: 4 out of 5 = 0.8
        let seq = b"AAAAGCGCAAAA";
        let gc = local_gc_content(seq, 6, 2);
        assert!((gc - 0.8).abs() < 0.001);
    }

    #[test]
    fn test_edge_distance_and_homopolymer_correlation() {
        // Real-world scenario: homopolymer at sequence edge
        let seq = b"AAAAAATGCATGCGGGGGG";
        // Position 0-5: in AAAAAA homopolymer, distance_to_edge = 0
        assert!(in_homopolymer(seq, 0, 4));
        assert_eq!(distance_to_edge(0, seq.len()), 0);

        // Position 13-18: in GGGGGG homopolymer, distance_to_edge close to 0
        assert!(in_homopolymer(seq, seq.len() - 1, 4));
        assert_eq!(distance_to_edge(seq.len() - 1, seq.len()), 0);
    }
}
