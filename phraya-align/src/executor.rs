use phraya_core::types::{Sequence, VariantObservation};
use phraya_io::plan::PhrayaPlan;
use std::collections::HashMap;

/// Result of a single alignment task.
#[derive(Debug, Clone)]
pub struct AlignmentResult {
    /// Variant observations at polymorphic sites
    pub variants: Vec<VariantObservation>,
    /// Coverage track (position → count), quantized to nearest 5
    pub coverage_track: Vec<u32>,
    /// Query index: list of all alignment positions and scores
    pub query_positions: Vec<u32>,
}

/// Execute a single alignment task: query vs target.
pub fn align_task(
    query: &Sequence,
    target: &Sequence,
    _plan: &PhrayaPlan,
) -> Option<AlignmentResult> {
    // For now, only handle sequences of equal length (simple alignment)
    if query.len() != target.len() {
        return None;
    }

    let mut variants = Vec::new();
    let coverage_track = vec![1; query.len()];
    let query_positions = (0..query.len() as u32).collect();

    // Extract variants where query != target
    for pos in 0..query.len() {
        let query_base = query.bases()[pos];
        let target_base = target.bases()[pos];

        if query_base != target_base {
            // Create a variant observation
            let mut alleles = HashMap::new();
            alleles.insert(query_base, 1);

            let variant = VariantObservation::new(
                pos as u32,
                target_base,
                alleles,
                1.0,                 // confidence
                "1M".to_string(),    // simple CIGAR
                60,                  // perfect MAPQ
                0,                   // edit_distance (will be updated with real alignment)
                vec![1],             // local coverage
                60.0,                // avg quality
                "query".to_string(), // provenance
            );
            variants.push(variant);
        }
    }

    Some(AlignmentResult {
        variants,
        coverage_track,
        query_positions,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_perfect_match_no_variants() {
        let query = Sequence::new(b"ACGTACGT".to_vec(), None, "query1".to_string(), None);
        let target = Sequence::new(b"ACGTACGT".to_vec(), None, "ref".to_string(), None);
        let plan = PhrayaPlan::new(
            phraya_io::plan::UseCase::ReadsWithRef,
            vec!["test".to_string()],
            "2026-06-01T00:00:00Z".to_string(),
            vec![],
            HashMap::new(),
            vec![],
        );

        let result = align_task(&query, &target, &plan);
        assert!(result.is_some());

        let result = result.unwrap();
        assert!(
            result.variants.is_empty(),
            "Perfect match should have no variants"
        );
        assert_eq!(
            result.coverage_track.len(),
            query.len(),
            "Coverage track should match query length"
        );
        assert_eq!(
            result.coverage_track[0], 1,
            "All positions should have coverage 1"
        );
    }

    #[test]
    fn test_single_snp_creates_variant() {
        // Query has T at position 2, target has C (SNP)
        let query = Sequence::new(b"ACTACGT".to_vec(), None, "query1".to_string(), None);
        let target = Sequence::new(b"ACCACGT".to_vec(), None, "ref".to_string(), None);
        let plan = PhrayaPlan::new(
            phraya_io::plan::UseCase::ReadsWithRef,
            vec!["test".to_string()],
            "2026-06-01T00:00:00Z".to_string(),
            vec![],
            HashMap::new(),
            vec![],
        );

        let result = align_task(&query, &target, &plan);
        assert!(result.is_some());

        let result = result.unwrap();
        assert_eq!(
            result.variants.len(),
            1,
            "One SNP should produce one variant"
        );

        let var = &result.variants[0];
        assert_eq!(
            var.position(),
            2,
            "Variant should be at position 2 (0-indexed)"
        );
        assert_eq!(var.ref_base(), b'C', "Reference base should be C");
        assert!(
            var.all_alleles().contains_key(&b'T'),
            "Allele T should be present"
        );
    }
}
