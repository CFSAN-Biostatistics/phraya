use crate::{score_alignments, wfa_extend, SeedAnchor};
use phraya_core::types::{Sequence, VariantObservation};
use phraya_index::{find_seeds, sketch_sequence_default};
use phraya_io::plan::PhrayaPlan;
use std::collections::HashMap;

/// Result of a single alignment task.
#[derive(Debug, Clone)]
pub struct AlignmentResult {
    /// Variant observations at polymorphic sites
    pub variants: Vec<VariantObservation>,
    /// Coverage track (position → count), quantized to nearest 5
    pub coverage_track: Vec<u32>,
    /// Query index: (target_position, normalized_score) for primary + alternatives
    pub query_positions: Vec<(u32, f64)>,
}

/// Execute a single alignment task: query vs target.
pub fn align_task(
    query: &Sequence,
    target: &Sequence,
    _plan: &PhrayaPlan,
) -> Option<AlignmentResult> {
    let query_sketch = sketch_sequence_default(query);
    let target_sketch = sketch_sequence_default(target);
    let seeds = find_seeds(&query_sketch, &target_sketch);

    // Collect WFA alignments from each seed; fall back to anchor (0,0) for short sequences.
    let mut alignments = Vec::new();

    let anchors: Vec<SeedAnchor> = if seeds.is_empty() {
        vec![SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        }]
    } else {
        seeds
            .iter()
            .map(|s| SeedAnchor {
                query_pos: s.query_pos as usize,
                target_pos: s.target_pos as usize,
            })
            .collect()
    };

    for anchor in anchors {
        match wfa_extend(query.bases(), target.bases(), anchor) {
            Ok(aln) => alignments.push(aln),
            Err(e) => log::warn!("WFA failed for anchor {:?}: {:?}", anchor, e),
        }
    }

    if alignments.is_empty() {
        return None;
    }

    let scored = score_alignments(&alignments, query.len());
    let primary_score =
        1.0 - (scored.primary.edit_distance as f64 / query.len().max(1) as f64);

    let variants = extract_variants_from_cigar(
        &scored.primary.cigar,
        scored.primary.target_start,
        query.bases(),
        target.bases(),
        scored.primary.edit_distance as u32,
        query.id().to_string(),
    );

    let coverage_track = compute_coverage_track(&scored, target.len());

    let mut query_positions = vec![(scored.primary.target_start as u32, primary_score)];
    for alt in &scored.alternatives {
        let alt_score = 1.0 - (alt.edit_distance as f64 / query.len().max(1) as f64);
        query_positions.push((alt.target_start as u32, alt_score));
    }

    Some(AlignmentResult {
        variants,
        coverage_track,
        query_positions,
    })
}

/// Parse CIGAR and extract VariantObservations at mismatch positions.
fn extract_variants_from_cigar(
    cigar: &str,
    target_start: usize,
    query: &[u8],
    target: &[u8],
    edit_distance: u32,
    provenance: String,
) -> Vec<VariantObservation> {
    let mut variants = Vec::new();
    let mut q_pos = 0usize;
    let mut t_pos = target_start;

    let ops = parse_cigar(cigar);
    for (count, op) in ops {
        match op {
            'M' => {
                q_pos += count;
                t_pos += count;
            }
            'X' => {
                // Mismatch: one VariantObservation per position
                for i in 0..count {
                    let qp = q_pos + i;
                    let tp = t_pos + i;
                    if qp < query.len() && tp < target.len() {
                        let alt_base = query[qp];
                        let ref_base = target[tp];
                        let mut alleles = HashMap::new();
                        alleles.insert(alt_base, 1u32);
                        variants.push(VariantObservation::new(
                            tp as u32,
                            ref_base,
                            alleles,
                            1.0,
                            cigar.to_string(),
                            60,
                            edit_distance,
                            vec![1],
                            60.0,
                            provenance.clone(),
                        ));
                    }
                }
                q_pos += count;
                t_pos += count;
            }
            'I' => {
                q_pos += count;
            }
            'D' => {
                t_pos += count;
            }
            _ => {}
        }
    }

    variants
}

fn parse_cigar(cigar: &str) -> Vec<(usize, char)> {
    let mut ops = Vec::new();
    let mut count_str = String::new();
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            count_str.push(ch);
        } else {
            let count: usize = count_str.parse().unwrap_or(1);
            count_str.clear();
            ops.push((count, ch));
        }
    }
    ops
}

fn compute_coverage_track(
    scored: &crate::ScoredAlignments,
    target_len: usize,
) -> Vec<u32> {
    let mut track = vec![0u32; target_len];

    let all_alns = std::iter::once(&scored.primary).chain(scored.alternatives.iter());
    for aln in all_alns {
        let start = aln.target_start.min(target_len);
        let end = aln.target_end.min(target_len);
        for pos in start..end {
            track[pos] = track[pos].saturating_add(1);
        }
    }

    // Quantize to nearest 5
    for v in track.iter_mut() {
        *v = (((*v as usize + 2) / 5) * 5) as u32;
    }

    track
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_plan() -> PhrayaPlan {
        PhrayaPlan::new(
            phraya_io::plan::UseCase::ReadsWithRef,
            vec!["test".to_string()],
            "2026-06-01T00:00:00Z".to_string(),
            vec![],
            HashMap::new(),
            vec![],
        )
    }

    #[test]
    fn test_align_task_handles_indel() {
        // Query has a deletion relative to target: target has 'T' at position 4 that query lacks.
        // Currently returns None due to equal-length guard — must use WFA instead.
        let query = Sequence::new(b"ACGACGT".to_vec(), None, "query_del".to_string(), None);
        let target = Sequence::new(b"ACGTACGT".to_vec(), None, "ref".to_string(), None);
        let plan = make_plan();

        let result = align_task(&query, &target, &plan);
        assert!(
            result.is_some(),
            "align_task must handle different-length sequences via WFA"
        );
    }

    #[test]
    fn test_perfect_match_no_variants() {
        let query = Sequence::new(b"ACGTACGT".to_vec(), None, "query1".to_string(), None);
        let target = Sequence::new(b"ACGTACGT".to_vec(), None, "ref".to_string(), None);
        let plan = make_plan();

        let result = align_task(&query, &target, &plan);
        assert!(result.is_some());

        let result = result.unwrap();
        assert!(
            result.variants.is_empty(),
            "Perfect match should have no variants"
        );
        assert_eq!(
            result.coverage_track.len(),
            target.len(),
            "Coverage track should match target length"
        );
    }

    #[test]
    fn test_query_positions_carry_scores() {
        let query = Sequence::new(b"ACGTACGT".to_vec(), None, "q".to_string(), None);
        let target = Sequence::new(b"ACGTACGT".to_vec(), None, "ref".to_string(), None);
        let plan = make_plan();

        let result = align_task(&query, &target, &plan).expect("alignment should succeed");

        assert!(!result.query_positions.is_empty(), "should have at least one position");
        let (_pos, score) = result.query_positions[0];
        assert!(
            score > 0.0 && score <= 1.0,
            "score must be in (0.0, 1.0], got {score}"
        );
    }

    #[test]
    fn test_variant_cigar_reflects_wfa_not_stub() {
        // SNP: T at position 2, otherwise identical 7-base sequences
        let query = Sequence::new(b"ACTACGT".to_vec(), None, "q".to_string(), None);
        let target = Sequence::new(b"ACCACGT".to_vec(), None, "ref".to_string(), None);
        let plan = make_plan();

        let result = align_task(&query, &target, &plan).expect("alignment should succeed");
        assert_eq!(result.variants.len(), 1);

        let cigar = result.variants[0].cigar();
        assert_ne!(cigar, "1M", "CIGAR must come from WFA, not stub");
        // WFA over 7 equal-length bases with 1 mismatch produces something like "2M1X4M"
        assert!(
            cigar.contains('X') || cigar.contains('M'),
            "CIGAR should contain M or X ops, got: {cigar}"
        );
        assert!(cigar.len() > 2, "CIGAR should represent the full alignment, got: {cigar}");
    }

    #[test]
    fn test_single_snp_creates_variant() {
        // Query has T at position 2, target has C (SNP)
        let query = Sequence::new(b"ACTACGT".to_vec(), None, "query1".to_string(), None);
        let target = Sequence::new(b"ACCACGT".to_vec(), None, "ref".to_string(), None);
        let plan = make_plan();

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
