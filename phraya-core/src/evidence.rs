/// Evidence layer extraction and serialization tests
///
/// These tests define the API contract for population evidence extraction.
/// All tests should FAIL until the implementation is complete.

#[cfg(test)]
mod tests {
    use crate::{EvidenceLayer, Sequence, extract_evidence};
    use std::collections::HashMap;

    /// Helper to create a minimal test sequence
    fn make_sequence(id: &str, bases: &str) -> Sequence {
        Sequence {
            id: id.to_string(),
            bases: bases.as_bytes().to_vec(),
            quality: None,
        }
    }

    #[test]
    fn test_evidence_layer_struct_exists() {
        // This test verifies that EvidenceLayer struct exists with expected fields
        let evidence = EvidenceLayer {
            kmer_uniqueness: HashMap::new(),
            invariant_positions: vec![],
            polymorphic_sites: HashMap::new(),
            sample_count: 0,
            reference_length: 0,
        };

        assert_eq!(evidence.sample_count, 0);
        assert_eq!(evidence.reference_length, 0);
    }

    #[test]
    fn test_extract_evidence_with_reference() {
        // Test extraction with a reference sequence provided
        let reference = make_sequence("ref", "ACGTACGT");
        let seq1 = make_sequence("sample1", "ACGTACGT");
        let seq2 = make_sequence("sample2", "ACGTTCGT");
        let seq3 = make_sequence("sample3", "ACGTTCGT");

        let sequences = vec![seq1, seq2, seq3];
        let evidence = extract_evidence(&sequences, Some(&reference));

        // Should detect position 4 as polymorphic (A in reference, T in 2/3 samples)
        assert_eq!(evidence.sample_count, 3);
        assert_eq!(evidence.reference_length, 8);
        assert!(evidence.polymorphic_sites.contains_key(&4));

        // Position 4 should show A->T variant with frequency 2/3
        let site = &evidence.polymorphic_sites[&4];
        assert_eq!(site.reference_base, b'A');
        assert_eq!(site.allele_counts.get(&b'T'), Some(&2));
        assert_eq!(site.allele_counts.get(&b'A'), Some(&1));
    }

    #[test]
    fn test_extract_evidence_without_reference_msa_mode() {
        // Test extraction in MSA mode (no reference)
        let seq1 = make_sequence("sample1", "ACGTACGT");
        let seq2 = make_sequence("sample2", "ACGTTCGT");
        let seq3 = make_sequence("sample3", "ACGTTCGT");

        let sequences = vec![seq1, seq2, seq3];
        let evidence = extract_evidence(&sequences, None);

        // In MSA mode, should still detect position 4 as polymorphic
        assert_eq!(evidence.sample_count, 3);
        assert!(evidence.polymorphic_sites.contains_key(&4));

        // Should detect both A and T alleles at position 4
        let site = &evidence.polymorphic_sites[&4];
        assert_eq!(site.allele_counts.len(), 2);
        assert!(site.allele_counts.contains_key(&b'A'));
        assert!(site.allele_counts.contains_key(&b'T'));
    }

    #[test]
    fn test_invariant_positions_detection() {
        // All samples match reference at these positions
        let reference = make_sequence("ref", "ACGTACGT");
        let seq1 = make_sequence("sample1", "ACGTACGT");
        let seq2 = make_sequence("sample2", "ACGTACGT");
        let seq3 = make_sequence("sample3", "ACGTACGT");

        let sequences = vec![seq1, seq2, seq3];
        let evidence = extract_evidence(&sequences, Some(&reference));

        // All 8 positions should be invariant
        assert_eq!(evidence.invariant_positions.len(), 8);
        assert_eq!(evidence.polymorphic_sites.len(), 0);
    }

    #[test]
    fn test_kmer_uniqueness_map() {
        // Test that k-mer uniqueness is computed from reference
        let reference = make_sequence("ref", "ACGTACGTACGTACGT");
        let sequences = vec![make_sequence("sample1", "ACGTACGTACGTACGT")];

        let evidence = extract_evidence(&sequences, Some(&reference));

        // Should have uniqueness scores for positions
        assert!(!evidence.kmer_uniqueness.is_empty());

        // Positions with repeated k-mers should have lower uniqueness scores
        // (actual values depend on k-mer size, but there should be variation)
        let uniqueness_values: Vec<f64> = evidence.kmer_uniqueness.values().copied().collect();
        let min_score = uniqueness_values
            .iter()
            .fold(f64::INFINITY, |a, &b| a.min(b));
        let max_score = uniqueness_values
            .iter()
            .fold(f64::NEG_INFINITY, |a, &b| a.max(b));

        // Repeated sequence should show variation in uniqueness
        assert!(max_score > min_score);
    }

    #[test]
    fn test_serialization_to_json() {
        // Test that evidence layer can be serialized to JSON
        let reference = make_sequence("ref", "ACGT");
        let seq1 = make_sequence("sample1", "ACGT");
        let seq2 = make_sequence("sample2", "ACTT");

        let sequences = vec![seq1, seq2];
        let evidence = extract_evidence(&sequences, Some(&reference));

        let json = serde_json::to_string(&evidence).expect("serialization should succeed");

        // Should contain key fields
        assert!(json.contains("sample_count"));
        assert!(json.contains("reference_length"));
        assert!(json.contains("polymorphic_sites"));
        assert!(json.contains("invariant_positions"));
        assert!(json.contains("kmer_uniqueness"));
    }

    #[test]
    fn test_deserialization_from_json() {
        // Test round-trip serialization
        let reference = make_sequence("ref", "ACGT");
        let seq1 = make_sequence("sample1", "ACGT");
        let seq2 = make_sequence("sample2", "ACTT");

        let sequences = vec![seq1, seq2];
        let evidence = extract_evidence(&sequences, Some(&reference));

        let json = serde_json::to_string(&evidence).expect("serialization should succeed");
        let deserialized: EvidenceLayer =
            serde_json::from_str(&json).expect("deserialization should succeed");

        // Verify fields match
        assert_eq!(deserialized.sample_count, evidence.sample_count);
        assert_eq!(deserialized.reference_length, evidence.reference_length);
        assert_eq!(
            deserialized.polymorphic_sites.len(),
            evidence.polymorphic_sites.len()
        );
        assert_eq!(
            deserialized.invariant_positions.len(),
            evidence.invariant_positions.len()
        );
    }

    #[test]
    fn test_deterministic_extraction() {
        // Same inputs should produce identical evidence layers
        let reference = make_sequence("ref", "ACGTACGT");
        let seq1 = make_sequence("sample1", "ACGTTCGT");
        let seq2 = make_sequence("sample2", "ACGTACGT");

        let sequences1 = vec![seq1.clone(), seq2.clone()];
        let sequences2 = vec![seq1, seq2];

        let evidence1 = extract_evidence(&sequences1, Some(&reference));
        let evidence2 = extract_evidence(&sequences2, Some(&reference));

        // Serialize to compare
        let json1 = serde_json::to_string(&evidence1).unwrap();
        let json2 = serde_json::to_string(&evidence2).unwrap();

        assert_eq!(json1, json2);
    }

    #[test]
    fn test_evidence_layer_size_is_linear() {
        // Evidence layer size should be O(positions), not O(sequences * positions)
        let reference = make_sequence("ref", "ACGTACGTACGTACGT");

        // Create many samples with variations
        let mut sequences = vec![];
        for i in 0..100 {
            let mut bases = "ACGTACGTACGTACGT".to_string();
            if i % 2 == 0 {
                bases.replace_range(4..5, "T");
            }
            sequences.push(make_sequence(&format!("sample{}", i), &bases));
        }

        let evidence = extract_evidence(&sequences, Some(&reference));

        // Verify that evidence layer is compact
        assert_eq!(evidence.sample_count, 100);

        // Polymorphic sites should only have one entry per position
        assert!(evidence.polymorphic_sites.len() <= reference.bases.len());

        // Each polymorphic site should contain counts, not individual observations
        for (_pos, site) in &evidence.polymorphic_sites {
            let total_count: usize = site.allele_counts.values().sum();
            assert_eq!(
                total_count, 100,
                "should aggregate counts, not store individual obs"
            );
        }
    }

    #[test]
    fn test_mixed_variation_patterns() {
        // Test with mix of invariant, biallelic, and triallelic sites
        let reference = make_sequence("ref", "ACGTACGT");
        let sequences = vec![
            make_sequence("s1", "ACGTACGT"), // matches reference
            make_sequence("s2", "TCGTACGT"), // pos 0: A->T
            make_sequence("s3", "GCGTACGT"), // pos 0: A->G (triallelic at pos 0)
            make_sequence("s4", "ACGTACGT"), // matches reference
            make_sequence("s5", "ACTTACGT"), // pos 3: T->T (already T in ref, so no change)
            make_sequence("s6", "ACGTCCGT"), // pos 4: A->C, pos 5: C->C (no change)
        ];

        let evidence = extract_evidence(&sequences, Some(&reference));

        // Position 0 should be triallelic (A, T, G)
        assert!(evidence.polymorphic_sites.contains_key(&0));
        let site0 = &evidence.polymorphic_sites[&0];
        assert_eq!(site0.allele_counts.len(), 3);

        // Position 4 should have variation
        assert!(evidence.polymorphic_sites.contains_key(&4));
    }

    #[test]
    fn test_empty_sequences() {
        // Edge case: empty sequence list
        let sequences: Vec<Sequence> = vec![];
        let evidence = extract_evidence(&sequences, None);

        assert_eq!(evidence.sample_count, 0);
        assert_eq!(evidence.reference_length, 0);
        assert!(evidence.polymorphic_sites.is_empty());
        assert!(evidence.invariant_positions.is_empty());
    }

    #[test]
    fn test_single_sequence() {
        // Edge case: single sequence (no variation possible)
        let reference = make_sequence("ref", "ACGT");
        let seq1 = make_sequence("sample1", "ACGT");

        let sequences = vec![seq1];
        let evidence = extract_evidence(&sequences, Some(&reference));

        assert_eq!(evidence.sample_count, 1);
        assert_eq!(evidence.polymorphic_sites.len(), 0);
        assert_eq!(evidence.invariant_positions.len(), 4);
    }

    #[test]
    fn test_unequal_length_sequences() {
        // Edge case: sequences of different lengths
        let reference = make_sequence("ref", "ACGTACGT");
        let seq1 = make_sequence("sample1", "ACGT"); // shorter
        let seq2 = make_sequence("sample2", "ACGTACGTAC"); // longer
        let seq3 = make_sequence("sample3", "ACGTACGT"); // matches ref

        let sequences = vec![seq1, seq2, seq3];
        let evidence = extract_evidence(&sequences, Some(&reference));

        // Should handle length differences gracefully
        assert_eq!(evidence.sample_count, 3);
        assert_eq!(evidence.reference_length, 8);
    }

    #[test]
    fn test_allele_frequency_calculation() {
        // Verify allele frequencies are accurately computed
        let reference = make_sequence("ref", "AAAA");
        let sequences = vec![
            make_sequence("s1", "TAAA"), // T at position 0
            make_sequence("s2", "TAAA"), // T at position 0
            make_sequence("s3", "TAAA"), // T at position 0
            make_sequence("s4", "AAAA"), // A at position 0
        ];

        let evidence = extract_evidence(&sequences, Some(&reference));

        let site0 = &evidence.polymorphic_sites[&0];
        assert_eq!(site0.allele_counts.get(&b'T'), Some(&3));
        assert_eq!(site0.allele_counts.get(&b'A'), Some(&1));

        // Could also compute frequency as 3/4 = 0.75 for T
        let total: usize = site0.allele_counts.values().sum();
        let t_freq = *site0.allele_counts.get(&b'T').unwrap() as f64 / total as f64;
        assert!((t_freq - 0.75).abs() < 0.01);
    }

    #[test]
    fn test_json_human_readable() {
        // JSON should be human-readable (pretty-printed or at least not binary)
        let reference = make_sequence("ref", "ACGT");
        let seq1 = make_sequence("sample1", "TCGT");

        let sequences = vec![seq1];
        let evidence = extract_evidence(&sequences, Some(&reference));

        let json =
            serde_json::to_string_pretty(&evidence).expect("pretty serialization should succeed");

        // Should be multi-line and readable
        assert!(json.contains('\n'));
        assert!(json.len() > 50); // Not just "{}"
    }

    #[test]
    fn test_performance_target_50_samples() {
        // Performance test: 50 samples @ 5Mbp each in <30 seconds
        // This test defines the performance contract but will fail until optimized
        use std::time::Instant;

        let size = 5_000_000;
        let sample_count = 50;

        // Generate reference
        let reference_bases: String = (0..size)
            .map(|i| match i % 4 {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                _ => 'T',
            })
            .collect();
        let reference = make_sequence("ref", &reference_bases);

        // Generate samples with occasional variations
        let mut sequences = vec![];
        for sample_idx in 0..sample_count {
            let mut bases = reference_bases.clone();
            // Add ~100 SNPs per sample at deterministic positions
            for snp_idx in 0..100 {
                let pos = (sample_idx * 1000 + snp_idx * 5000) % size;
                bases.replace_range(pos..pos + 1, "T");
            }
            sequences.push(make_sequence(&format!("sample{}", sample_idx), &bases));
        }

        let start = Instant::now();
        let evidence = extract_evidence(&sequences, Some(&reference));
        let elapsed = start.elapsed();

        // Verify correctness
        assert_eq!(evidence.sample_count, sample_count);
        assert_eq!(evidence.reference_length, size);

        // Performance target: <30 seconds
        assert!(
            elapsed.as_secs() < 30,
            "extraction took {:?}, expected <30s",
            elapsed
        );
    }
}
