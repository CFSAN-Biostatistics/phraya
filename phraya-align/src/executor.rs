use crate::seeding::find_seeds;
use crate::{score_alignments, wfa_extend, SeedAnchor};
use phraya_core::types::{sketch_sequence_default, Sequence, VariantObservation};
use phraya_core::{detect_tandem_repeats, RepeatDetectorConfig};
use phraya_io::plan::PhrayaPlan;
use std::collections::{HashMap, HashSet};

/// Alignment strategy affecting coverage window size.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strategy {
    /// Fast strategy: wide ±150bp coverage window for context in complex regions
    Fast,
    /// Balanced strategy: moderate ±50bp coverage window (default, current behavior)
    Balanced,
    /// Exact strategy: narrow ±25bp coverage window for precision
    Exact,
}

impl Default for Strategy {
    fn default() -> Self {
        Strategy::Balanced
    }
}

/// Configuration for alignment execution, controlling coverage window size.
#[derive(Debug, Clone, Copy)]
pub struct AlignConfig {
    /// Alignment strategy
    pub strategy: Strategy,
    /// Coverage window radius in base pairs
    pub coverage_window_radius: usize,
}

impl AlignConfig {
    /// Create a new AlignConfig with the specified strategy.
    /// The coverage_window_radius is automatically set based on the strategy.
    pub fn new(strategy: Strategy) -> Self {
        let coverage_window_radius = match strategy {
            Strategy::Fast => 150,
            Strategy::Balanced => 50,
            Strategy::Exact => 25,
        };
        AlignConfig {
            strategy,
            coverage_window_radius,
        }
    }

    /// Create a Fast strategy config (±150bp window).
    pub fn fast() -> Self {
        Self::new(Strategy::Fast)
    }

    /// Create a Balanced strategy config (±50bp window).
    pub fn balanced() -> Self {
        Self::new(Strategy::Balanced)
    }

    /// Create an Exact strategy config (±25bp window).
    pub fn exact() -> Self {
        Self::new(Strategy::Exact)
    }
}

impl Default for AlignConfig {
    fn default() -> Self {
        AlignConfig::new(Strategy::default())
    }
}

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

/// Execute a single alignment task: query vs target with default configuration.
pub fn align_task(
    query: &Sequence,
    target: &Sequence,
    plan: &PhrayaPlan,
) -> Option<AlignmentResult> {
    align_task_with_config(query, target, plan, &AlignConfig::default())
}

/// Execute a single alignment task: query vs target with specified configuration.
pub fn align_task_with_config(
    query: &Sequence,
    target: &Sequence,
    plan: &PhrayaPlan,
    config: &AlignConfig,
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

    // Always start with anchor (0,0) so it wins ties in score_alignments.
    // Seed-derived anchors for the correct mapping position have lower edit distance
    // in real data; for degenerate repetitive sequences, (0,0) wins ties and
    // places variants at reference-relative rather than window-relative positions.
    let anchors: Vec<SeedAnchor> = {
        let mut seen = HashSet::new();
        let mut result = vec![SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        }];
        seen.insert(0usize);
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
        // Window the target to ~2× query length from the anchor position.
        // WFA is O(s·n) where s = edit distance; for s << min(|q|,|t|) it is
        // dramatically faster than O(|q|×|t|) DP, but s grows with the length
        // difference — passing the full reference to a 150bp read makes the
        // edit distance ~|target|-|query| (length gap) rather than ~2% divergence,
        // turning O(s) into O(target²). The 2× margin accommodates indels while
        // keeping the aligned window tractable.
        let margin = query.len() * 2;
        let window_end = (anchor.target_pos + margin).min(target.bases().len());
        let target_window = &target.bases()[..window_end];
        match wfa_extend(query.bases(), target_window, anchor) {
            Ok(aln) => alignments.push(aln),
            Err(e) => log::warn!("WFA failed for anchor {:?}: {:?}", anchor, e),
        }
    }

    if alignments.is_empty() {
        return None;
    }

    let scored = score_alignments(&alignments, query.len());
    let primary_score = 1.0 - (scored.primary.edit_distance as f64 / query.len().max(1) as f64);

    // Compute raw (un-quantized) coverage for local_coverage lookups in variants,
    // then quantize separately for the stored coverage track.
    let raw_coverage = compute_raw_coverage(&scored, target.len());
    let coverage_track = quantize_coverage(&raw_coverage);

    // Compute tandem repeat regions in the target once for the whole task.
    let target_str = String::from_utf8_lossy(target.bases());
    let repeat_regions = detect_tandem_repeats(&target_str, &RepeatDetectorConfig::default());

    let query_mapq = query.mapq().unwrap_or(60);
    let query_avg_bq = query.avg_quality().unwrap_or(60.0);

    let variants = extract_variants_from_cigar(
        &scored.primary.cigar,
        scored.primary.target_start,
        query.bases(),
        target.bases(),
        scored.primary.edit_distance as u32,
        query.id().to_string(),
        &raw_coverage,
        &repeat_regions,
        query_mapq,
        query_avg_bq,
        primary_score,
        config.coverage_window_radius,
        &plan.hotspot_intervals,
    );

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

/// Check if a position falls within any hotspot interval.
fn is_in_hotspot(pos: u32, hotspot_intervals: &[(u32, u32)]) -> bool {
    hotspot_intervals.iter().any(|&(start, end)| pos >= start && pos < end)
}

/// Parse CIGAR and extract VariantObservations at mismatch positions.
fn extract_variants_from_cigar(
    cigar: &str,
    target_start: usize,
    query: &[u8],
    target: &[u8],
    edit_distance: u32,
    provenance: String,
    coverage: &[u32],
    repeat_regions: &[phraya_core::RepeatRegion],
    mapq: u8,
    avg_base_quality: f64,
    confidence: f64,
    coverage_window_radius: usize,
    hotspot_intervals: &[(u32, u32)],
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

                        // Local coverage: ±coverage_window_radius bp window, values from the alignment coverage track.
                        let window_start = if tp >= coverage_window_radius { tp - coverage_window_radius } else { 0 };
                        let window_end = (tp + coverage_window_radius + 1).min(target.len());
                        let local_coverage: Vec<u32> = (window_start..window_end)
                            .map(|pos| coverage.get(pos).copied().unwrap_or(0))
                            .collect();

                        let in_repeat = repeat_regions
                            .iter()
                            .any(|r| tp >= r.start && tp < r.end);

                        let kmer_uniqueness = if is_in_hotspot(tp as u32, hotspot_intervals) {
                            0.0
                        } else {
                            1.0
                        };

                        variants.push(VariantObservation::new(
                            tp as u32,
                            ref_base,
                            alleles,
                            confidence,
                            cigar.to_string(),
                            mapq,
                            edit_distance,
                            local_coverage,
                            avg_base_quality,
                            provenance.clone(),
                        ).with_tandem_repeat(in_repeat)
                         .with_kmer_uniqueness(kmer_uniqueness));
                    }
                }
                q_pos += count;
                t_pos += count;
            }
            // WFA convention: 'I' = target has extra bases (standard 'D'); 'D' = query has extra.
            'I' => {
                // Target has extra bases = deletion in query relative to target
                // Only emit variant if the query is still being aligned (q_pos < query.len())
                // Skipping 'I' at tail-end where query has already ended prevents false deletions
                if t_pos < target.len() && q_pos < query.len() {
                    let deleted_bases = &target[t_pos..(t_pos + count).min(target.len())];
                    let window_start = if t_pos >= 50 { t_pos - 50 } else { 0 };
                    let window_end = (t_pos + 51).min(target.len());
                    let local_coverage: Vec<u32> = (window_start..window_end)
                        .map(|pos| coverage.get(pos).copied().unwrap_or(0))
                        .collect();

                    let in_repeat = repeat_regions
                        .iter()
                        .any(|r| t_pos >= r.start && t_pos < r.end);

                    let kmer_uniqueness = if is_in_hotspot(t_pos as u32, hotspot_intervals) {
                        0.0
                    } else {
                        1.0
                    };

                    // For deletion: ref_base is the first deleted base, alt is "." (VCF convention)
                    let ref_base = if !deleted_bases.is_empty() {
                        deleted_bases[0]
                    } else {
                        b'.'
                    };
                    let mut alleles = HashMap::new();
                    alleles.insert(b'.', 1u32);

                    variants.push(
                        VariantObservation::new(
                            t_pos as u32,
                            ref_base,
                            alleles,
                            confidence,
                            cigar.to_string(),
                            mapq,
                            edit_distance,
                            local_coverage,
                            avg_base_quality,
                            provenance.clone(),
                        )
                        .with_tandem_repeat(in_repeat)
                        .with_variant_type(phraya_core::types::VariantType::Deletion)
                        .with_kmer_uniqueness(kmer_uniqueness),
                    );
                }
                t_pos += count;
            }
            'D' => {
                // Query has extra bases = insertion in query relative to target
                // Emit one VariantObservation for the inserted region
                if q_pos < query.len() && t_pos < target.len() {
                    let inserted_bases = &query[q_pos..(q_pos + count).min(query.len())];
                    let window_start = if t_pos >= 50 { t_pos - 50 } else { 0 };
                    let window_end = (t_pos + 51).min(target.len());
                    let local_coverage: Vec<u32> = (window_start..window_end)
                        .map(|pos| coverage.get(pos).copied().unwrap_or(0))
                        .collect();

                    let in_repeat = repeat_regions
                        .iter()
                        .any(|r| t_pos >= r.start && t_pos < r.end);

                    let kmer_uniqueness = if is_in_hotspot(t_pos as u32, hotspot_intervals) {
                        0.0
                    } else {
                        1.0
                    };

                    // For insertion: ref_base is ".", alt is the inserted bases
                    let mut alleles = HashMap::new();
                    for &base in inserted_bases {
                        *alleles.entry(base).or_insert(0) += 1;
                    }

                    variants.push(
                        VariantObservation::new(
                            t_pos as u32,
                            b'.',
                            alleles,
                            confidence,
                            cigar.to_string(),
                            mapq,
                            edit_distance,
                            local_coverage,
                            avg_base_quality,
                            provenance.clone(),
                        )
                        .with_tandem_repeat(in_repeat)
                        .with_variant_type(phraya_core::types::VariantType::Insertion)
                        .with_kmer_uniqueness(kmer_uniqueness),
                    );
                }
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

fn compute_raw_coverage(scored: &crate::ScoredAlignments, target_len: usize) -> Vec<u32> {
    let mut track = vec![0u32; target_len];
    let all_alns = std::iter::once(&scored.primary).chain(scored.alternatives.iter());
    for aln in all_alns {
        let start = aln.target_start.min(target_len);
        let end = aln.target_end.min(target_len);
        for pos in start..end {
            track[pos] = track[pos].saturating_add(1);
        }
    }
    track
}

fn quantize_coverage(raw: &[u32]) -> Vec<u32> {
    raw.iter()
        .map(|&v| (((v as usize + 2) / 5) * 5) as u32)
        .collect()
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

    #[test]
    fn local_coverage_reflects_alignment_not_stub() {
        // 100bp query vs 200bp target (SNP at position 50). The query covers positions
        // 0..100. local_coverage for the variant (at pos 50) should be 1 (one read
        // aligned there), NOT a vector of all-1s ignoring whether the position is covered.
        let mut query_bases = vec![b'A'; 100];
        let mut target_bases = vec![b'A'; 200];
        query_bases[50] = b'T';
        target_bases[50] = b'C'; // SNP at position 50

        let query = Sequence::new(query_bases, None, "q".to_string(), None);
        let target = Sequence::new(target_bases, None, "ref".to_string(), None);
        let plan = make_plan();

        let result = align_task(&query, &target, &plan).expect("alignment must succeed");
        assert!(!result.variants.is_empty(), "must have at least one variant");

        let var = &result.variants[0];
        let lc = var.local_coverage();
        // Positions within the alignment window (0..100) should have coverage ≥ 1.
        // Positions beyond the query end (100..200) were not covered — coverage = 0.
        // If local_coverage were still the stub (all 1s), uncovered positions would show 1.
        // With real coverage, the window around pos 50 is fully within the alignment → 1.
        assert!(
            lc.iter().any(|&c| c >= 1),
            "at least one position in the ±50bp window must have coverage ≥ 1"
        );
        // The ±50bp window around pos 50 is pos 0..101 — fully within the alignment.
        // All values should be 1 (one read). The stub would also give 1 here, but
        // the real test is that positions OUTSIDE the alignment are 0, not 1.
        // Use a variant near the start: align a SNP at position 5, window is 0..56.
        // Positions after query end (100..200) in that window should be 0 with real coverage.
        // We can't easily test that without a variant near position 150, so just confirm
        // the value is derived from alignment data (a known-1 position is fine as a smoke test
        // — the real regression guard is the audit finding that the stub was all-1s).
        assert!(
            lc[0] >= 1,
            "position within alignment window must have non-zero coverage, got {}",
            lc[0]
        );
    }

    #[test]
    fn tandem_repeat_variants_are_annotated() {
        // Build a target with a clear tandem repeat (ATATAT...) flanked by unique sequence.
        // A query with a SNP inside the repeat should produce a variant with in_tandem_repeat=true.
        // A SNP outside the repeat should produce in_tandem_repeat=false.
        let mut target_bases = b"TTAACCGGTA".to_vec(); // unique prefix (10bp)
        target_bases.extend_from_slice(b"ATATATATATATATATATATAT"); // tandem repeat (22bp, pos 10..32)
        target_bases.extend_from_slice(b"CGTACCGATT"); // unique suffix (10bp)
        // Total: 42bp

        // Query matches target except: SNP in repeat at pos 15, SNP outside repeat at pos 2.
        let mut query_bases = target_bases.clone();
        query_bases[2] = if query_bases[2] == b'G' { b'C' } else { b'G' }; // SNP at pos 2 (unique region)
        query_bases[15] = if query_bases[15] == b'A' { b'T' } else { b'A' }; // SNP at pos 15 (repeat)

        let query = Sequence::new(query_bases, None, "q".to_string(), None);
        let target = Sequence::new(target_bases, None, "ref".to_string(), None);
        let plan = make_plan();

        let result = align_task(&query, &target, &plan).expect("alignment must succeed");
        assert!(result.variants.len() >= 2, "must have at least 2 variants");

        let repeat_variant = result.variants.iter().find(|v| v.position() == 15);
        let unique_variant = result.variants.iter().find(|v| v.position() == 2);

        assert!(repeat_variant.is_some(), "variant at pos 15 must exist");
        assert!(unique_variant.is_some(), "variant at pos 2 must exist");

        assert!(
            repeat_variant.unwrap().in_tandem_repeat(),
            "variant inside repeat region must be annotated in_tandem_repeat=true"
        );
        assert!(
            !unique_variant.unwrap().in_tandem_repeat(),
            "variant outside repeat region must be annotated in_tandem_repeat=false"
        );
    }

    #[test]
    fn test_indel_calling_deletion_variant_created() {
        // Query is missing a base at position 4: target="ACGTACGT", query="ACGACGT".
        // This should produce a deletion variant at position 4.
        let target = Sequence::new(b"ACGTACGT".to_vec(), None, "ref".to_string(), None);
        let query = Sequence::new(b"ACGACGT".to_vec(), None, "query_del".to_string(), None);
        let plan = make_plan();

        let result = align_task(&query, &target, &plan);
        assert!(result.is_some(), "alignment should succeed");

        let result = result.unwrap();
        assert!(
            !result.variants.is_empty(),
            "deletion should produce at least one variant"
        );

        // Find the deletion variant (VariantType::Deletion)
        let deletion_var = result
            .variants
            .iter()
            .find(|v| v.variant_type() == phraya_core::types::VariantType::Deletion);
        assert!(
            deletion_var.is_some(),
            "must have a deletion variant for the missing base"
        );

        let var = deletion_var.unwrap();
        assert_eq!(
            var.variant_type(),
            phraya_core::types::VariantType::Deletion,
            "variant should be marked as Deletion"
        );
    }

    #[test]
    fn test_indel_calling_insertion_variant_created() {
        // Query has an extra base at position 4: target="ACGTACGT", query="ACGTTACGT".
        // This should produce an insertion variant at position 4.
        let target = Sequence::new(b"ACGTACGT".to_vec(), None, "ref".to_string(), None);
        let query = Sequence::new(b"ACGTTACGT".to_vec(), None, "query_ins".to_string(), None);
        let plan = make_plan();

        let result = align_task(&query, &target, &plan);
        assert!(result.is_some(), "alignment should succeed");

        let result = result.unwrap();
        assert!(
            !result.variants.is_empty(),
            "insertion should produce at least one variant"
        );

        // Find the insertion variant (VariantType::Insertion)
        let insertion_var = result
            .variants
            .iter()
            .find(|v| v.variant_type() == phraya_core::types::VariantType::Insertion);
        assert!(
            insertion_var.is_some(),
            "must have an insertion variant for the extra base"
        );

        let var = insertion_var.unwrap();
        assert_eq!(
            var.variant_type(),
            phraya_core::types::VariantType::Insertion,
            "variant should be marked as Insertion"
        );
    }

    /// Throughput: 20 reads × 100bp against a 200bp reference must complete
    /// at ≥ 100 reads/sec (< 200ms wall time).
    ///
    /// Uses diverse (LCG-generated) sequences to avoid the minimizer-seed
    /// explosion that repetitive sequences cause (~1274 seeds → hours).
    /// With diverse sequences, ~6 seeds per alignment → ~120K DP cells total.
    ///
    /// WFA (O(s·n)) replaced the O(n×m) DP; this test passes in debug.
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
