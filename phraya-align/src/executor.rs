use crate::seeding::find_seeds;
use crate::{score_alignments, wfa_extend, SeedAnchor};
use phraya_core::types::{sketch_sequence_default, Sequence, VariantObservation};
use phraya_io::plan::PhrayaPlan;
use std::collections::{HashMap, HashSet};

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
    plan: &PhrayaPlan,
) -> Option<AlignmentResult> {
    // Reuse pre-computed sketches from plan if available; fall back to recomputing
    let query_sketch = plan
        .get_sketch(query.id())
        .cloned()
        .unwrap_or_else(|| sketch_sequence_default(query));
    let target_sketch = plan
        .get_sketch(target.id())
        .cloned()
        .unwrap_or_else(|| sketch_sequence_default(target));
    let seeds = find_seeds(&query_sketch, &target_sketch);

    // Convert seeds to full-query anchors (query_pos=0, target_pos=target-query offset).
    // Seeds mid-query would miss variants before the seed; aligning from query position 0
    // ensures the full query is aligned. Deduplicate by target_start to avoid redundant calls.
    let mut alignments = Vec::new();

    let anchors: Vec<SeedAnchor> = if seeds.is_empty() {
        vec![SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        }]
    } else {
        let mut seen = HashSet::new();
        let mut result = Vec::new();
        for s in &seeds {
            let target_start = (s.target_pos as i64 - s.query_pos as i64).max(0) as usize;
            if seen.insert(target_start) {
                result.push(SeedAnchor {
                    query_pos: 0,
                    target_pos: target_start,
                });
            }
        }
        result
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
    let primary_score = 1.0 - (scored.primary.edit_distance as f64 / query.len().max(1) as f64);

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
            // WFA convention: 'I' = target has extra bases (standard 'D'); 'D' = query has extra.
            'I' => {
                t_pos += count;
            }
            'D' => {
                q_pos += count;
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

fn compute_coverage_track(scored: &crate::ScoredAlignments, target_len: usize) -> Vec<u32> {
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
            HashMap::new(),
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

        assert!(
            !result.query_positions.is_empty(),
            "should have at least one position"
        );
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
        assert!(
            cigar.len() > 2,
            "CIGAR should represent the full alignment, got: {cigar}"
        );
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

    /// Throughput: 20 reads × 100bp against a 200bp reference must complete
    /// at ≥ 100 reads/sec (< 200ms wall time).
    ///
    /// Uses diverse (LCG-generated) sequences to avoid the minimizer-seed
    /// explosion that repetitive sequences cause (~1274 seeds → hours).
    /// With diverse sequences, ~6 seeds per alignment → ~120K DP cells total.
    ///
    /// CURRENTLY FAILS in debug builds: the naive O(n×m) DP is too slow for
    /// even short sequences. Production target requires the true WFA wavefront
    /// algorithm (sub-quadratic time/space) — the current DP is a placeholder.
    #[test]
    fn issue_88_throughput_100_reads_per_sec() {
        fn diverse_dna(len: usize, seed: u64) -> Vec<u8> {
            let mut x = seed;
            (0..len)
                .map(|_| {
                    x = x
                        .wrapping_mul(6364136223846793005)
                        .wrapping_add(1442695040888963407);
                    b"ACGT"[((x >> 33) & 3) as usize]
                })
                .collect()
        }

        let ref_seq = diverse_dna(200, 42);
        let read_seq: Vec<u8> = ref_seq[..100].to_vec();

        let target = Sequence::new(ref_seq, None, "ref".to_string(), None);
        let plan = make_plan();

        let start = std::time::Instant::now();
        for i in 0..20 {
            let query = Sequence::new(read_seq.clone(), None, format!("read{i}"), None);
            let _ = align_task(&query, &target, &plan);
        }
        let elapsed = start.elapsed();

        assert!(
            elapsed.as_millis() < 200,
            "20 alignments (150bp vs 1000bp) took {}ms — below 100 reads/sec target.\n\
             The naive O(n×m) DP must be replaced with true WFA wavefront algorithm.",
            elapsed.as_millis()
        );
    }
}
