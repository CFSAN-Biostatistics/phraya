// This module will contain core types for Phraya.
// Tests are written first (TDD RED phase).

#[cfg(test)]
mod tests {
    use super::*;

    // ===== Sequence type tests =====

    #[test]
    fn sequence_creation_without_quality() {
        let seq = Sequence::new(
            b"ACGT".to_vec(),
            None,
            "seq1".to_string(),
            Some("test sequence".to_string()),
        );
        assert_eq!(seq.len(), 4);
        assert_eq!(seq.id(), "seq1");
        assert_eq!(seq.description(), Some("test sequence"));
        assert!(seq.quality_scores().is_none());
    }

    #[test]
    fn sequence_creation_with_quality() {
        let seq = Sequence::new(
            b"ACGT".to_vec(),
            Some(vec![30, 35, 40, 38]),
            "seq2".to_string(),
            None,
        );
        assert_eq!(seq.len(), 4);
        assert_eq!(seq.quality_at(0), Some(30));
        assert_eq!(seq.quality_at(1), Some(35));
        assert_eq!(seq.quality_at(2), Some(40));
        assert_eq!(seq.quality_at(3), Some(38));
        assert_eq!(seq.quality_at(4), None);
    }

    #[test]
    fn sequence_avg_quality_calculation() {
        let seq = Sequence::new(
            b"ACGT".to_vec(),
            Some(vec![20, 30, 40, 30]),
            "seq3".to_string(),
            None,
        );
        assert_eq!(seq.avg_quality(), Some(30.0));
    }

    #[test]
    fn sequence_avg_quality_without_scores() {
        let seq = Sequence::new(b"ACGT".to_vec(), None, "seq4".to_string(), None);
        assert_eq!(seq.avg_quality(), None);
    }

    #[test]
    fn sequence_empty() {
        let seq = Sequence::new(b"".to_vec(), None, "empty".to_string(), None);
        assert_eq!(seq.len(), 0);
        assert_eq!(seq.avg_quality(), None);
    }

    #[test]
    #[should_panic(expected = "quality scores length must match sequence length")]
    fn sequence_quality_length_mismatch_panics() {
        let _seq = Sequence::new(
            b"ACGT".to_vec(),
            Some(vec![30, 35]), // Only 2 quality scores for 4 bases
            "bad".to_string(),
            None,
        );
    }

    #[test]
    fn sequence_serialization() {
        let seq = Sequence::new(
            b"ACGT".to_vec(),
            Some(vec![30, 35, 40, 38]),
            "seq5".to_string(),
            Some("description".to_string()),
        );

        let json = serde_json::to_string(&seq).expect("serialization failed");
        let deserialized: Sequence = serde_json::from_str(&json).expect("deserialization failed");

        assert_eq!(deserialized.len(), 4);
        assert_eq!(deserialized.id(), "seq5");
        assert_eq!(deserialized.quality_at(0), Some(30));
    }

    // ===== VariantObservation type tests =====

    #[test]
    fn variant_observation_creation() {
        let mut alleles = std::collections::HashMap::new();
        alleles.insert(b'A', 10);
        alleles.insert(b'T', 5);

        let obs = VariantObservation::new(
            100,
            b'A',
            alleles,
            0.95,
            "10M".to_string(),
            60,
            2,
            vec![10, 12, 15, 18, 20],
            35.5,
            "sample1:read42".to_string(),
        );

        assert_eq!(obs.position(), 100);
        assert_eq!(obs.ref_base(), b'A');
        assert_eq!(obs.confidence(), 0.95);
        assert_eq!(obs.mapq(), 60);
        assert_eq!(obs.edit_distance(), 2);
        assert_eq!(obs.avg_base_quality(), 35.5);
        assert_eq!(obs.provenance(), "sample1:read42");
    }

    #[test]
    fn variant_observation_allele_counts() {
        let mut alleles = std::collections::HashMap::new();
        alleles.insert(b'A', 10);
        alleles.insert(b'T', 5);
        alleles.insert(b'G', 2);

        let obs = VariantObservation::new(
            200,
            b'A',
            alleles.clone(),
            0.98,
            "5M1I4M".to_string(),
            50,
            1,
            vec![20, 22, 25],
            40.0,
            "sample2:read99".to_string(),
        );

        let all_alleles = obs.all_alleles();
        assert_eq!(all_alleles.get(&b'A'), Some(&10));
        assert_eq!(all_alleles.get(&b'T'), Some(&5));
        assert_eq!(all_alleles.get(&b'G'), Some(&2));

        let total: u32 = all_alleles.values().sum();
        assert_eq!(total, 17);
    }

    #[test]
    fn variant_observation_local_coverage() {
        let obs = VariantObservation::new(
            150,
            b'C',
            [(b'C', 8), (b'T', 2)].into_iter().collect(),
            0.90,
            "20M".to_string(),
            55,
            0,
            vec![8, 9, 10, 10, 10, 12, 15, 18, 20, 22],
            38.0,
            "sample3:read1".to_string(),
        );

        let coverage = obs.local_coverage();
        assert_eq!(coverage.len(), 10);
        assert_eq!(coverage[0], 8);
        assert_eq!(coverage[9], 22);
    }

    #[test]
    fn variant_observation_serialization() {
        let obs = VariantObservation::new(
            300,
            b'G',
            [(b'G', 15), (b'A', 3)].into_iter().collect(),
            0.99,
            "25M".to_string(),
            60,
            0,
            vec![15, 16, 18, 20],
            42.0,
            "sample4:read5".to_string(),
        );

        let json = serde_json::to_string(&obs).expect("serialization failed");
        let deserialized: VariantObservation =
            serde_json::from_str(&json).expect("deserialization failed");

        assert_eq!(deserialized.position(), 300);
        assert_eq!(deserialized.ref_base(), b'G');
        assert_eq!(deserialized.mapq(), 60);
    }

    // ===== EvidenceLayer type tests =====

    #[test]
    fn evidence_layer_creation() {
        let mut kmer_uniqueness = std::collections::HashMap::new();
        kmer_uniqueness.insert(100, 1.0);
        kmer_uniqueness.insert(200, 0.5);

        let mut polymorphic_sites = std::collections::HashMap::new();
        polymorphic_sites.insert(150, vec![b'A', b'T']);

        let mut invariant_positions = std::collections::HashSet::new();
        invariant_positions.insert(50);
        invariant_positions.insert(51);

        let mut multi_map_fraction = std::collections::HashMap::new();
        multi_map_fraction.insert(100, 0.2);

        let mut avg_score_ratio_gap = std::collections::HashMap::new();
        avg_score_ratio_gap.insert(100, 0.15);

        let evidence = EvidenceLayer::new(
            kmer_uniqueness,
            polymorphic_sites,
            invariant_positions,
            multi_map_fraction,
            avg_score_ratio_gap,
        );

        assert_eq!(evidence.kmer_uniqueness().get(&100), Some(&1.0));
        assert_eq!(evidence.kmer_uniqueness().get(&200), Some(&0.5));
        assert_eq!(
            evidence.polymorphic_sites().get(&150),
            Some(&vec![b'A', b'T'])
        );
        assert!(evidence.invariant_positions().contains(&50));
        assert!(evidence.invariant_positions().contains(&51));
        assert_eq!(evidence.multi_map_fraction().get(&100), Some(&0.2));
        assert_eq!(evidence.avg_score_ratio_gap().get(&100), Some(&0.15));
    }

    #[test]
    fn evidence_layer_empty() {
        let evidence = EvidenceLayer::new(
            std::collections::HashMap::new(),
            std::collections::HashMap::new(),
            std::collections::HashSet::new(),
            std::collections::HashMap::new(),
            std::collections::HashMap::new(),
        );

        assert!(evidence.kmer_uniqueness().is_empty());
        assert!(evidence.polymorphic_sites().is_empty());
        assert!(evidence.invariant_positions().is_empty());
        assert!(evidence.multi_map_fraction().is_empty());
        assert!(evidence.avg_score_ratio_gap().is_empty());
    }

    #[test]
    fn evidence_layer_serialization() {
        let mut kmer_uniqueness = std::collections::HashMap::new();
        kmer_uniqueness.insert(100, 1.0);

        let evidence = EvidenceLayer::new(
            kmer_uniqueness,
            std::collections::HashMap::new(),
            std::collections::HashSet::new(),
            std::collections::HashMap::new(),
            std::collections::HashMap::new(),
        );

        let json = serde_json::to_string(&evidence).expect("serialization failed");
        let deserialized: EvidenceLayer =
            serde_json::from_str(&json).expect("deserialization failed");

        assert_eq!(deserialized.kmer_uniqueness().get(&100), Some(&1.0));
    }

    // ===== CoverageTrack type stub tests =====

    #[test]
    fn coverage_track_stub_exists() {
        // Just verify the type exists - implementation is for a separate slice
        let _track: CoverageTrack;
    }

    // ===== Error type tests =====

    #[test]
    fn parse_error_creation() {
        let err = ParseError::InvalidFormat("bad FASTA".to_string());
        assert_eq!(format!("{}", err), "invalid format: bad FASTA");
    }

    #[test]
    fn io_error_creation() {
        let err = IoError::FileNotFound("/path/to/file".to_string());
        assert_eq!(format!("{}", err), "file not found: /path/to/file");
    }

    #[test]
    fn alignment_error_creation() {
        let err = AlignmentError::NoSeeds("seq1".to_string());
        assert_eq!(format!("{}", err), "no seeds found for sequence: seq1");
    }

    #[test]
    fn filter_error_creation() {
        let err = FilterError::InvalidThreshold("coverage".to_string(), -5.0);
        assert_eq!(format!("{}", err), "invalid threshold for coverage: -5");
    }

    #[test]
    fn error_serialization() {
        // Errors should be serializable for structured logging
        let err = ParseError::InvalidFormat("test".to_string());
        let json = serde_json::to_string(&err).expect("serialization failed");
        let deserialized: ParseError = serde_json::from_str(&json).expect("deserialization failed");
        assert_eq!(format!("{}", deserialized), format!("{}", err));
    }

    // ===== Property tests =====

    #[test]
    fn property_quality_scores_match_sequence_length() {
        // Property: if quality scores are present, they must match sequence length
        let sequences = vec![
            (b"ACGT".to_vec(), vec![30, 35, 40, 38]),
            (b"A".to_vec(), vec![40]),
            (b"ACGTACGT".to_vec(), vec![30, 30, 30, 30, 40, 40, 40, 40]),
        ];

        for (bases, quality) in sequences {
            let seq = Sequence::new(
                bases.clone(),
                Some(quality.clone()),
                "test".to_string(),
                None,
            );
            assert_eq!(bases.len(), quality.len());
            assert_eq!(seq.len(), quality.len());
        }
    }

    #[test]
    fn property_allele_counts_are_positive() {
        // Property: all allele counts must be positive
        let mut alleles = std::collections::HashMap::new();
        alleles.insert(b'A', 10);
        alleles.insert(b'T', 5);
        alleles.insert(b'G', 1);

        let obs = VariantObservation::new(
            100,
            b'A',
            alleles.clone(),
            0.95,
            "10M".to_string(),
            60,
            0,
            vec![10],
            35.0,
            "test".to_string(),
        );

        for count in obs.all_alleles().values() {
            assert!(*count > 0, "all allele counts must be positive");
        }
    }

    #[test]
    fn property_allele_counts_sum_correctly() {
        // Property: sum of allele counts should equal total coverage at position
        let mut alleles = std::collections::HashMap::new();
        alleles.insert(b'A', 10);
        alleles.insert(b'T', 5);
        alleles.insert(b'G', 2);
        alleles.insert(b'C', 3);

        let obs = VariantObservation::new(
            100,
            b'A',
            alleles.clone(),
            0.95,
            "10M".to_string(),
            60,
            0,
            vec![20],
            35.0,
            "test".to_string(),
        );

        let sum: u32 = obs.all_alleles().values().sum();
        assert_eq!(sum, 20, "sum of allele counts should match total coverage");
    }

    #[test]
    fn property_ref_base_should_be_in_alleles() {
        // Property: reference base should typically be present in alleles map
        let mut alleles = std::collections::HashMap::new();
        alleles.insert(b'A', 10);
        alleles.insert(b'T', 5);

        let obs = VariantObservation::new(
            100,
            b'A',
            alleles.clone(),
            0.95,
            "10M".to_string(),
            60,
            0,
            vec![15],
            35.0,
            "test".to_string(),
        );

        assert!(
            obs.all_alleles().contains_key(&obs.ref_base()),
            "reference base should be in alleles map"
        );
    }

    #[test]
    fn property_confidence_in_valid_range() {
        // Property: confidence should be between 0.0 and 1.0
        let obs = VariantObservation::new(
            100,
            b'A',
            [(b'A', 10)].into_iter().collect(),
            0.95,
            "10M".to_string(),
            60,
            0,
            vec![10],
            35.0,
            "test".to_string(),
        );

        let conf = obs.confidence();
        assert!(
            conf >= 0.0 && conf <= 1.0,
            "confidence must be in range [0.0, 1.0]"
        );
    }

    #[test]
    fn property_mapq_in_valid_range() {
        // Property: MAPQ should be in valid range [0, 60]
        let obs = VariantObservation::new(
            100,
            b'A',
            [(b'A', 10)].into_iter().collect(),
            0.95,
            "10M".to_string(),
            60,
            0,
            vec![10],
            35.0,
            "test".to_string(),
        );

        let mapq = obs.mapq();
        assert!(mapq <= 60, "MAPQ should be in valid range [0, 60]");
    }

    #[test]
    fn property_kmer_uniqueness_in_valid_range() {
        // Property: k-mer uniqueness values should be between 0.0 and 1.0
        let mut kmer_uniqueness = std::collections::HashMap::new();
        kmer_uniqueness.insert(100, 1.0);
        kmer_uniqueness.insert(200, 0.5);
        kmer_uniqueness.insert(300, 0.0);

        let evidence = EvidenceLayer::new(
            kmer_uniqueness.clone(),
            std::collections::HashMap::new(),
            std::collections::HashSet::new(),
            std::collections::HashMap::new(),
            std::collections::HashMap::new(),
        );

        for value in evidence.kmer_uniqueness().values() {
            assert!(
                *value >= 0.0 && *value <= 1.0,
                "k-mer uniqueness must be in range [0.0, 1.0]"
            );
        }
    }

    #[test]
    fn property_multi_map_fraction_in_valid_range() {
        // Property: multi-mapping fraction should be between 0.0 and 1.0
        let mut multi_map_fraction = std::collections::HashMap::new();
        multi_map_fraction.insert(100, 0.0);
        multi_map_fraction.insert(200, 0.5);
        multi_map_fraction.insert(300, 1.0);

        let evidence = EvidenceLayer::new(
            std::collections::HashMap::new(),
            std::collections::HashMap::new(),
            std::collections::HashSet::new(),
            multi_map_fraction.clone(),
            std::collections::HashMap::new(),
        );

        for value in evidence.multi_map_fraction().values() {
            assert!(
                *value >= 0.0 && *value <= 1.0,
                "multi-map fraction must be in range [0.0, 1.0]"
            );
        }
    }
}
