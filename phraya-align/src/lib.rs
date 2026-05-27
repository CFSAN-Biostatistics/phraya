// Placeholder for the unified aligner implementation
// This will be implemented when #9 is tackled

#[cfg(test)]
mod tests {
    use phraya_core::{BaseConfidence, EvidenceLayer, Sequence, VariantObservation};

    // Test helper to create a simple DNA sequence
    fn make_sequence(data: &[u8]) -> Sequence {
        Sequence::new(data.to_vec(), None)
    }

    // Test helper to create a sequence with quality scores
    fn make_sequence_with_quality(data: &[u8], quality: &[u8]) -> Sequence {
        Sequence::new(data.to_vec(), Some(quality.to_vec()))
    }

    #[test]
    fn test_aligner_exists() {
        // RED: Aligner struct should exist
        let _aligner = super::Aligner::new();
    }

    #[test]
    fn test_aligner_basic_alignment() {
        // RED: Test basic alignment of identical sequences
        let aligner = super::Aligner::new();
        let query = make_sequence(b"ACGTACGTACGT");
        let target = make_sequence(b"ACGTACGTACGT");
        let evidence = EvidenceLayer::empty();

        let observations = aligner.align(&query, &target, &evidence);

        // Should return empty variant observations for identical sequences
        assert_eq!(observations.len(), 0);
    }

    #[test]
    fn test_aligner_snp_detection() {
        // RED: Detect single SNP
        let aligner = super::Aligner::new();
        let query = make_sequence(b"ACGTACGTACGT");
        let target = make_sequence(b"ACGTCCGTACGT");
        //                                     ^
        let evidence = EvidenceLayer::empty();

        let observations = aligner.align(&query, &target, &evidence);

        // Should detect one variant at position 5
        assert_eq!(observations.len(), 1);
        assert_eq!(observations[0].position, 5);
        assert_eq!(observations[0].ref_base, b'A');
        assert_eq!(observations[0].alt_base, b'C');
    }

    #[test]
    fn test_aligner_multiple_snps() {
        // RED: Detect multiple SNPs
        let aligner = super::Aligner::new();
        let query = make_sequence(b"ACGTACGTACGT");
        let target = make_sequence(b"TCGTACGTACCG");
        //                               ^       ^
        let evidence = EvidenceLayer::empty();

        let observations = aligner.align(&query, &target, &evidence);

        // Should detect two variants
        assert_eq!(observations.len(), 2);
        assert_eq!(observations[0].position, 0);
        assert_eq!(observations[1].position, 10);
    }

    #[test]
    fn test_aligner_insertion() {
        // RED: Detect insertion
        let aligner = super::Aligner::new();
        let query = make_sequence(b"ACGTACGT");
        let target = make_sequence(b"ACGTTACGT");
        //                                  ^
        let evidence = EvidenceLayer::empty();

        let observations = aligner.align(&query, &target, &evidence);

        // Should detect insertion around position 4
        assert!(observations.len() > 0);
        assert!(observations.iter().any(|obs| obs.is_insertion()));
    }

    #[test]
    fn test_aligner_deletion() {
        // RED: Detect deletion
        let aligner = super::Aligner::new();
        let query = make_sequence(b"ACGTTACGT");
        let target = make_sequence(b"ACGTACGT");
        let evidence = EvidenceLayer::empty();

        let observations = aligner.align(&query, &target, &evidence);

        // Should detect deletion around position 4
        assert!(observations.len() > 0);
        assert!(observations.iter().any(|obs| obs.is_deletion()));
    }

    #[test]
    fn test_aligner_assembly_length_sequences() {
        // RED: Test with 10kb sequences (lower end of assembly range)
        let aligner = super::Aligner::new();
        let query_data = vec![b'A'; 10_000];
        let target_data = vec![b'A'; 10_000];
        let query = make_sequence(&query_data);
        let target = make_sequence(&target_data);
        let evidence = EvidenceLayer::empty();

        let observations = aligner.align(&query, &target, &evidence);

        // Should successfully align without errors
        assert_eq!(observations.len(), 0);
    }

    #[test]
    fn test_aligner_large_assembly_sequences() {
        // RED: Test with 500kb sequences (upper end of assembly range)
        let aligner = super::Aligner::new();
        let mut query_data = vec![b'A'; 500_000];
        let mut target_data = vec![b'A'; 500_000];
        // Add a few SNPs
        query_data[1000] = b'T';
        target_data[1000] = b'C';
        query_data[250_000] = b'G';
        target_data[250_000] = b'T';

        let query = make_sequence(&query_data);
        let target = make_sequence(&target_data);
        let evidence = EvidenceLayer::empty();

        let observations = aligner.align(&query, &target, &evidence);

        // Should detect the two SNPs
        assert_eq!(observations.len(), 2);
    }

    #[test]
    fn test_aligner_read_length_sequences() {
        // RED: Test with 100bp sequences (short read length)
        let aligner = super::Aligner::new();
        let query = make_sequence(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        let target = make_sequence(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        let evidence = EvidenceLayer::empty();

        let observations = aligner.align(&query, &target, &evidence);

        assert_eq!(observations.len(), 0);
    }

    #[test]
    fn test_aligner_1kb_sequences() {
        // RED: Test with 1kb sequences (long read length)
        let aligner = super::Aligner::new();
        let query_data = vec![b'C'; 1000];
        let target_data = vec![b'C'; 1000];
        let query = make_sequence(&query_data);
        let target = make_sequence(&target_data);
        let evidence = EvidenceLayer::empty();

        let observations = aligner.align(&query, &target, &evidence);

        assert_eq!(observations.len(), 0);
    }

    #[test]
    fn test_aligner_adapts_seed_density_for_long_sequences() {
        // RED: Verify that seed density is adapted for long sequences
        // This is an indirect test - we verify behavior, not implementation
        let aligner = super::Aligner::new();
        let short_query = make_sequence(&vec![b'G'; 1000]);
        let long_query = make_sequence(&vec![b'G'; 100_000]);
        let target = make_sequence(&vec![b'G'; 100_000]);
        let evidence = EvidenceLayer::empty();

        // Both should successfully align
        let short_obs = aligner.align(&short_query, &target, &evidence);
        let long_obs = aligner.align(&long_query, &target, &evidence);

        assert_eq!(short_obs.len(), 0);
        assert_eq!(long_obs.len(), 0);
    }

    #[test]
    fn test_aligner_confidence_metadata() {
        // RED: Verify that variant observations include confidence metadata
        let aligner = super::Aligner::new();
        let query = make_sequence(b"ACGTACGTACGT");
        let target = make_sequence(b"TCGTACGTACGT");
        let evidence = EvidenceLayer::empty();

        let observations = aligner.align(&query, &target, &evidence);

        assert_eq!(observations.len(), 1);
        let obs = &observations[0];

        // Should have confidence metadata
        assert!(obs.confidence.combined_score >= 0.0 && obs.confidence.combined_score <= 1.0);
        assert!(obs.confidence.alignment_quality.is_some());
    }

    #[test]
    fn test_aligner_evidence_layer_improves_confidence() {
        // RED: Verify that evidence layer affects confidence scoring
        let aligner = super::Aligner::new();
        let query = make_sequence(b"ACGTACGTACGT");
        let target = make_sequence(b"TCGTACGTACGT");

        // First alignment with empty evidence
        let empty_evidence = EvidenceLayer::empty();
        let obs_no_evidence = aligner.align(&query, &target, &empty_evidence);

        // Create mock evidence that says position 0 is polymorphic with high frequency
        let mut mock_evidence = EvidenceLayer::empty();
        mock_evidence.mark_polymorphic_site(0, b'A', b'T', 0.45);

        let obs_with_evidence = aligner.align(&query, &target, &mock_evidence);

        // Confidence should be higher when evidence supports the variant
        assert!(
            obs_with_evidence[0].confidence.combined_score
                > obs_no_evidence[0].confidence.combined_score
        );
    }

    #[test]
    fn test_aligner_evidence_layer_penalizes_novel_variants() {
        // RED: Verify that evidence layer penalizes variants not seen in population
        let aligner = super::Aligner::new();
        let query = make_sequence(b"ACGTACGTACGT");
        let target = make_sequence(b"TCGTACGTACGT");

        // Empty evidence means this variant is novel
        let empty_evidence = EvidenceLayer::empty();
        let obs_novel = aligner.align(&query, &target, &empty_evidence);

        // Evidence marking position 0 as invariant (everyone has 'A')
        let mut evidence_invariant = EvidenceLayer::empty();
        evidence_invariant.mark_invariant_position(0, b'A');

        let obs_invariant = aligner.align(&query, &target, &evidence_invariant);

        // Novel variant against invariant position should have lower confidence
        assert!(
            obs_invariant[0].confidence.combined_score < obs_novel[0].confidence.combined_score
        );
    }

    #[test]
    fn test_aligner_integration_ecoli_contig() {
        // RED: Integration test with realistic E. coli sequence
        // This uses a small synthetic fragment that mimics E. coli characteristics
        let aligner = super::Aligner::new();

        // 1kb synthetic sequence with ~50% GC (typical E. coli)
        let reference_seq = b"ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\
                              ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\
                              GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG\
                              GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG\
                              ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\
                              ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\
                              GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG\
                              GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG\
                              ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\
                              ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\
                              GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG\
                              GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG\
                              ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\
                              ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\
                              GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG\
                              GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG";

        // Query with a few SNPs
        let mut query_seq = reference_seq.to_vec();
        query_seq[100] = b'T'; // A->T
        query_seq[500] = b'G'; // C->G

        let reference = make_sequence(reference_seq);
        let query = make_sequence(&query_seq);

        // Minimal evidence layer
        let evidence = EvidenceLayer::empty();

        let observations = aligner.align(&query, &reference, &evidence);

        // Should detect the two SNPs
        assert_eq!(observations.len(), 2);
        assert_eq!(observations[0].position, 100);
        assert_eq!(observations[1].position, 500);

        // Both should have reasonable confidence
        for obs in &observations {
            assert!(obs.confidence.combined_score > 0.0);
            assert!(obs.confidence.combined_score <= 1.0);
        }
    }

    #[test]
    fn test_aligner_quality_scores_affect_confidence() {
        // RED: Verify quality scores influence confidence
        let aligner = super::Aligner::new();

        // Same variant in two sequences with different qualities
        let query_high_qual = make_sequence_with_quality(
            b"TCGTACGTACGT",
            &vec![40u8; 12], // High quality (Q40)
        );
        let query_low_qual = make_sequence_with_quality(
            b"TCGTACGTACGT",
            &vec![10u8; 12], // Low quality (Q10)
        );
        let target = make_sequence(b"ACGTACGTACGT");
        let evidence = EvidenceLayer::empty();

        let obs_high_qual = aligner.align(&query_high_qual, &target, &evidence);
        let obs_low_qual = aligner.align(&query_low_qual, &target, &evidence);

        // High quality should result in higher confidence
        assert!(
            obs_high_qual[0].confidence.combined_score > obs_low_qual[0].confidence.combined_score
        );
    }

    #[test]
    fn test_aligner_pipeline_stages() {
        // RED: Verify the full pipeline executes all stages
        // This is an indirect test that the pipeline includes:
        // 1. Sketch seeding
        // 2. WFA extension
        // 3. Confidence scoring with evidence
        // 4. Variant observation emission

        let aligner = super::Aligner::new();
        let query = make_sequence(b"ACGTACGTACGTACGTACGT");
        let target = make_sequence(b"ACGTACGTCCGTACGTACGT");
        //                                     ^^
        let evidence = EvidenceLayer::empty();

        let observations = aligner.align(&query, &target, &evidence);

        // Should successfully complete all pipeline stages
        assert_eq!(observations.len(), 2);

        // Each observation should have:
        // - Position from CIGAR parsing
        assert!(observations[0].position < 20);
        // - Confidence from scorer
        assert!(observations[0].confidence.combined_score >= 0.0);
        // - Ref/alt bases from alignment
        assert!(observations[0].ref_base != observations[0].alt_base);
    }

    #[test]
    fn test_aligner_handles_ambiguous_bases() {
        // RED: Test handling of N (ambiguous) bases
        let aligner = super::Aligner::new();
        let query = make_sequence(b"ACGTNACGTACGT");
        let target = make_sequence(b"ACGTAACGTACGT");
        let evidence = EvidenceLayer::empty();

        let observations = aligner.align(&query, &target, &evidence);

        // Should handle N appropriately (may or may not report as variant)
        // Key requirement: should not panic
        assert!(observations.len() <= 1);
    }

    #[test]
    fn test_aligner_empty_sequences() {
        // RED: Edge case - empty sequences
        let aligner = super::Aligner::new();
        let query = make_sequence(b"");
        let target = make_sequence(b"");
        let evidence = EvidenceLayer::empty();

        let observations = aligner.align(&query, &target, &evidence);

        // Should handle gracefully
        assert_eq!(observations.len(), 0);
    }

    #[test]
    fn test_aligner_very_divergent_sequences() {
        // RED: Test with highly divergent sequences (50% difference)
        let aligner = super::Aligner::new();
        let query = make_sequence(b"AAAAAAAAAAAAAAAAAAAA");
        let target = make_sequence(b"TTTTTTTTTTTTTTTTTTTT");
        let evidence = EvidenceLayer::empty();

        let observations = aligner.align(&query, &target, &evidence);

        // Should detect all 20 mismatches
        assert_eq!(observations.len(), 20);

        // All should be SNPs
        for obs in &observations {
            assert_eq!(obs.ref_base, b'A');
            assert_eq!(obs.alt_base, b'T');
        }
    }

    #[test]
    fn test_aligner_complex_indel_region() {
        // RED: Test complex region with nearby indels
        let aligner = super::Aligner::new();
        let query = make_sequence(b"ACGTACGTAAACGTACGT");
        let target = make_sequence(b"ACGTACGTACGT");
        //                               Remove AAA
        let evidence = EvidenceLayer::empty();

        let observations = aligner.align(&query, &target, &evidence);

        // Should detect the deletion
        assert!(observations.len() > 0);
        assert!(observations.iter().any(|obs| obs.is_deletion()));
    }
}
