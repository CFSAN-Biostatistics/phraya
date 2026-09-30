//! Core SNP alignment and distance matrix formatters.
//!
//! These are multi-sample output formats: they require observations from
//! **all** samples to agree on which positions are "core" (present in every
//! sample). Use with multiple input `.phraya` files:
//!
//! ```sh
//! phraya filter s1.phraya s2.phraya s3.phraya --format core-snp --output core.fasta
//! phraya filter s1.phraya s2.phraya s3.phraya --format distance-matrix --output dist.tsv
//! ```

use phraya_core::types::{VariantObservation, VariantType};
use phraya_io::phraya::PhrayaFile;
use std::collections::HashMap;

/// A core SNP position: present in all samples, with the called base per sample.
#[derive(Debug, Clone)]
pub struct CoreSnp {
    /// 0-indexed reference position.
    pub position: u32,
    /// Reference base at this position.
    pub ref_base: u8,
    /// Called base per sample, in the same order as `sample_ids`.
    pub sample_bases: Vec<u8>,
    /// Whether each sample had an observation at this position.
    pub sample_present: Vec<bool>,
}

/// Compute core SNPs from multiple filtered `.phraya` files.
///
/// A position is "core" if:
/// 1. Every sample has at least one SNP observation at that position
/// 2. All observations at that position are SNPs (indels excluded)
///
/// The called base per sample is the allele with the highest aggregate count
/// across all observations at that position. Missing → `b'N'`.
///
/// # Arguments
/// * `files` — one `PhrayaFile` per sample (already filtered)
/// * `sample_ids` — sample labels, same order as `files`
pub fn compute_core_snps(
    files: &[PhrayaFile],
    sample_ids: &[String],
) -> Vec<CoreSnp> {
    if files.is_empty() || sample_ids.len() != files.len() {
        return Vec::new();
    }

    // For each file: position → Vec<&VariantObservation> (SNP obs only)
    let mut per_file_positions: Vec<HashMap<u32, Vec<&VariantObservation>>> = Vec::new();

    for file in files {
        let mut pos_map: HashMap<u32, Vec<&VariantObservation>> = HashMap::new();
        for obs in &file.observations {
            if obs.variant_type() == VariantType::Snp {
                pos_map.entry(obs.position()).or_default().push(obs);
            }
        }
        per_file_positions.push(pos_map);
    }

    // Find positions present in ALL files
    let all_positions: Vec<u32> = per_file_positions[0]
        .keys()
        .filter(|pos| per_file_positions.iter().all(|m| m.contains_key(pos)))
        .copied()
        .collect();

    let mut core_positions = all_positions;
    core_positions.sort();

    let n_samples = files.len();
    let mut core_snps = Vec::new();

    for pos in core_positions {
        let mut sample_bases = Vec::with_capacity(n_samples);
        let mut sample_present = vec![true; n_samples];
        let mut ref_base: u8 = b'N';

        for file_idx in 0..n_samples {
            let obs_list = &per_file_positions[file_idx];
            if let Some(observations) = obs_list.get(&pos) {
                // Aggregate allele counts across all observations at this position
                let mut allele_counts: HashMap<u8, u32> = HashMap::new();
                for obs in observations {
                    for (&allele, &count) in obs.all_alleles() {
                        *allele_counts.entry(allele).or_insert(0) += count;
                    }
                    if ref_base == b'N' {
                        ref_base = obs.ref_base();
                    }
                }

                // Called base = allele with highest count
                // Tie-break: prefer the alphabetically smaller base
                let best_base = allele_counts
                    .iter()
                    .max_by_key(|(allele, &count)| (count, std::u8::MAX - *allele))
                    .map(|(&allele, _)| allele)
                    .unwrap_or(b'N');

                sample_bases.push(best_base);
            } else {
                sample_bases.push(b'N');
                sample_present[file_idx] = false;
            }
        }

        core_snps.push(CoreSnp {
            position: pos,
            ref_base,
            sample_bases,
            sample_present,
        });
    }

    core_snps
}

/// Format core SNPs as a multi-FASTA file.
///
/// Each sample gets one FASTA record. The sequence is the concatenation
/// of called bases at each core SNP position (A/C/G/T/N for missing).
///
/// # Arguments
/// * `files` — one `PhrayaFile` per sample (already filtered)
/// * `sample_ids` — FASTA headers, same order as `files`
pub fn format_core_snp(
    files: &[PhrayaFile],
    sample_ids: &[String],
) -> String {
    let core_snps = compute_core_snps(files, sample_ids);

    let mut output = String::new();

    for (i, sample_id) in sample_ids.iter().enumerate() {
        output.push_str(&format!(">{}\n", sample_id));

        let seq: String = core_snps
            .iter()
            .map(|snp| snp.sample_bases.get(i).copied().unwrap_or(b'N') as char)
            .collect();

        output.push_str(&seq);
        output.push('\n');
    }

    output
}

/// Format pairwise distances as a square TSV matrix.
///
/// Distance = Hamming distance on core SNP bases: count of positions where
/// two samples have differing called bases, **excluding** positions where
/// either sample is missing (N).
///
/// # Arguments
/// * `files` — one `PhrayaFile` per sample (already filtered)
/// * `sample_ids` — matrix row/column headers, same order as `files`
pub fn format_distance_matrix(
    files: &[PhrayaFile],
    sample_ids: &[String],
) -> String {
    let core_snps = compute_core_snps(files, sample_ids);
    let n_samples = sample_ids.len();

    // Build per-sample base vectors for efficient pairwise comparison
    let sample_bases: Vec<Vec<u8>> = (0..n_samples)
        .map(|i| {
            core_snps
                .iter()
                .map(|snp| snp.sample_bases.get(i).copied().unwrap_or(b'N'))
                .collect::<Vec<u8>>()
        })
        .collect();

    let mut output = String::new();

    // Header row
    output.push('\t');
    for id in sample_ids {
        output.push_str(id);
        output.push('\t');
    }
    output.push('\n');

    // Data rows
    for i in 0..n_samples {
        output.push_str(&sample_ids[i]);

        for j in 0..n_samples {
            let dist = compute_hamming_distance(&sample_bases[i], &sample_bases[j]);
            output.push('\t');
            output.push_str(&dist.to_string());
        }

        output.push('\n');
    }

    output
}

/// Compute Hamming distance on SNP base vectors.
///
/// Counts positions where `a` and `b` differ. Positions where either is `N`
/// (missing) are excluded from both the numerator and denominator.
fn compute_hamming_distance(a: &[u8], b: &[u8]) -> u32 {
    a.iter()
        .zip(b.iter())
        .filter(|(&x, &y)| x != b'N' && y != b'N' && x != y)
        .count() as u32
}

#[cfg(test)]
mod tests {
    use super::*;
    use phraya_core::types::Strand;
    use std::collections::HashMap;

    /// Helper: create a minimal VariantObservation for testing.
    fn make_snp_obs(
        position: u32,
        ref_base: u8,
        alt_base: u8,
        alt_count: u32,
        ref_count: u32,
    ) -> VariantObservation {
        let mut alleles = HashMap::new();
        alleles.insert(ref_base, ref_count);
        alleles.insert(alt_base, alt_count);
        VariantObservation::new(
            position,
            ref_base,
            alleles,
            0.99,
            "50M".to_string(),
            60,
            1,
            vec![10],
            35.0,
            "sample:read1".to_string(),
        )
        .with_variant_type(VariantType::Snp)
        .with_strand(Strand::Forward)
        .with_query_position(25)
        .with_coverage_window_offset(25)
    }

    #[test]
    fn test_compute_core_snps_single_sample() {
        // Single sample: all SNP positions are "core" (present in all 1 file)
        let file = PhrayaFile::new(
            100,
            "sample1".to_string(),
            "2024-01-01T00:00:00Z".to_string(),
            vec![
                make_snp_obs(10, b'A', b'T', 8, 2),
                make_snp_obs(20, b'C', b'G', 5, 5),
                make_snp_obs(30, b'G', b'A', 10, 0),
                make_snp_obs(40, b'T', b'C', 3, 7),
            ],
            phraya_core::types::CoverageTrack::uncovered(100),
        );

        let snps = compute_core_snps(&[file], &["sample1".to_string()]);

        assert_eq!(snps.len(), 4);
        assert_eq!(snps[0].position, 10);
        assert_eq!(snps[0].ref_base, b'A');
        // Called base = highest count: T(8) vs A(2) → T
        assert_eq!(snps[0].sample_bases, vec![b'T']);
    }

    #[test]
    fn test_compute_core_snps_multi_sample() {
        // Two samples, positions 10 and 20 in common
        let file1 = PhrayaFile::new(
            100,
            "sample1".to_string(),
            "2024-01-01T00:00:00Z".to_string(),
            vec![
                make_snp_obs(10, b'A', b'T', 8, 2),  // T called
                make_snp_obs(20, b'C', b'G', 5, 5),  // C called (tie: C < G alphabetically)
                make_snp_obs(30, b'G', b'A', 10, 0), // A called (sample1 only)
            ],
            phraya_core::types::CoverageTrack::uncovered(100),
        );

        let file2 = PhrayaFile::new(
            100,
            "sample2".to_string(),
            "2024-01-01T00:00:00Z".to_string(),
            vec![
                make_snp_obs(10, b'A', b'T', 1, 9),  // A called (ref, high count)
                make_snp_obs(20, b'C', b'G', 7, 3),  // G called (alt_count=7 > ref_count=3)
                make_snp_obs(40, b'T', b'C', 10, 0), // C called (sample2 only)
            ],
            phraya_core::types::CoverageTrack::uncovered(100),
        );

        let snps = compute_core_snps(&[file1, file2], &["sample1".to_string(), "sample2".to_string()]);

        // Only positions 10 and 20 are present in BOTH files → core SNPs
        assert_eq!(snps.len(), 2);
        assert_eq!(snps[0].position, 10);
        // sample1: T (8 > 2), sample2: A (9 > 1)
        assert_eq!(snps[0].sample_bases, vec![b'T', b'A']);
        assert_eq!(snps[1].position, 20);
        // sample1: C (tie 5/5, C < G alphabetically), sample2: G (7 > 3)
        assert_eq!(snps[1].sample_bases, vec![b'C', b'G']);
    }

    #[test]
    fn test_compute_core_snps_excludes_indels() {
        let file = PhrayaFile::new(
            100,
            "sample1".to_string(),
            "2024-01-01T00:00:00Z".to_string(),
            vec![
                make_snp_obs(10, b'A', b'T', 5, 5),
                VariantObservation::new(
                    15,
                    b'A',
                    {
                        let mut m = HashMap::new();
                        m.insert(b'C', 5);
                        m
                    },
                    0.99,
                    "5M1D3M".to_string(),
                    60,
                    0,
                    vec![10],
                    35.0,
                    "sample:read1".to_string(),
                )
                .with_variant_type(VariantType::Insertion)
                .with_strand(Strand::Forward)
                .with_query_position(5)
                .with_coverage_window_offset(5),
            ],
            phraya_core::types::CoverageTrack::uncovered(100),
        );

        let snps = compute_core_snps(&[file], &["sample1".to_string()]);

        // Only the SNP is a core SNP; the insertion is excluded
        assert_eq!(snps.len(), 1);
        assert_eq!(snps[0].position, 10);
        assert_eq!(snps[0].ref_base, b'A');
    }

    #[test]
    fn test_compute_core_snps_excludes_missing_positions() {
        // position 30 only in file1, not file2 → not core
        let file1 = PhrayaFile::new(
            100,
            "s1".to_string(),
            "2024-01-01T00:00:00Z".to_string(),
            vec![
                make_snp_obs(10, b'A', b'T', 5, 5),
                make_snp_obs(30, b'C', b'G', 5, 5),
            ],
            phraya_core::types::CoverageTrack::uncovered(100),
        );
        let file2 = PhrayaFile::new(
            100,
            "s2".to_string(),
            "2024-01-01T00:00:00Z".to_string(),
            vec![
                make_snp_obs(10, b'A', b'T', 5, 5),
                make_snp_obs(20, b'G', b'C', 5, 5),
            ],
            phraya_core::types::CoverageTrack::uncovered(100),
        );

        let snps = compute_core_snps(&[file1, file2], &["s1".to_string(), "s2".to_string()]);

        // Only position 10 is in both files
        assert_eq!(snps.len(), 1);
        assert_eq!(snps[0].position, 10);
    }

    #[test]
    fn test_compute_core_snps_empty() {
        let snps = compute_core_snps(&[], &[]);
        assert!(snps.is_empty());
    }

    #[test]
    fn test_format_core_snp_single_sample() {
        let file = PhrayaFile::new(
            100,
            "sample1".to_string(),
            "2024-01-01T00:00:00Z".to_string(),
            vec![
                make_snp_obs(10, b'A', b'T', 8, 2), // T
                make_snp_obs(20, b'C', b'G', 5, 5), // C (tie)
            ],
            phraya_core::types::CoverageTrack::uncovered(100),
        );

        let output = format_core_snp(&[file], &["sample1".to_string()]);
        let lines: Vec<&str> = output.lines().collect();
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0], ">sample1");
        assert_eq!(lines[1], "TC");
    }

    #[test]
    fn test_format_core_snp_multi_sample() {
        let file1 = PhrayaFile::new(
            100,
            "sample1".to_string(),
            "2024-01-01T00:00:00Z".to_string(),
            vec![
                make_snp_obs(10, b'A', b'T', 8, 2), // T
                make_snp_obs(20, b'C', b'G', 5, 5), // C
            ],
            phraya_core::types::CoverageTrack::uncovered(100),
        );
        let file2 = PhrayaFile::new(
            100,
            "sample2".to_string(),
            "2024-01-01T00:00:00Z".to_string(),
            vec![
                make_snp_obs(10, b'A', b'T', 1, 9), // A
                make_snp_obs(20, b'C', b'G', 7, 3), // G (alt_count=7 > ref_count=3)
            ],
            phraya_core::types::CoverageTrack::uncovered(100),
        );

        let output = format_core_snp(
            &[file1, file2],
            &["sample1".to_string(), "sample2".to_string()],
        );
        let lines: Vec<&str> = output.lines().collect();
        assert_eq!(lines.len(), 4);
        assert_eq!(lines[0], ">sample1");
        assert_eq!(lines[1], "TC");
        assert_eq!(lines[2], ">sample2");
        assert_eq!(lines[3], "AG");
    }

    #[test]
    fn test_format_core_snp_missing_data() {
        // Position 30 is only in s1, not in s2 → not a core SNP, excluded from output
        let file1 = PhrayaFile::new(
            100,
            "s1".to_string(),
            "2024-01-01T00:00:00Z".to_string(),
            vec![
                make_snp_obs(10, b'A', b'T', 8, 2), // T (alt_count=8 > ref_count=2)
                make_snp_obs(20, b'C', b'G', 5, 5), // C (tie, C < G alphabetically)
                make_snp_obs(30, b'G', b'A', 10, 0), // A (pos 30 only in s1 → not core)
            ],
            phraya_core::types::CoverageTrack::uncovered(100),
        );
        let file2 = PhrayaFile::new(
            100,
            "s2".to_string(),
            "2024-01-01T00:00:00Z".to_string(),
            vec![
                make_snp_obs(10, b'A', b'T', 1, 9), // A (ref_count=9 > alt_count=1)
                make_snp_obs(20, b'C', b'G', 7, 3), // G (alt_count=7 > ref_count=3)
            ],
            phraya_core::types::CoverageTrack::uncovered(100),
        );

        let output = format_core_snp(
            &[file1, file2],
            &["s1".to_string(), "s2".to_string()],
        );
        let lines: Vec<&str> = output.lines().collect();
        assert_eq!(lines.len(), 4);
        // Core SNPs = positions 10 and 20 only (30 only in s1)
        // s1: T, C → "TC"; s2: A, G → "AG"
        assert_eq!(lines[1], "TC");
        assert_eq!(lines[3], "AG");
    }

    #[test]
    fn test_format_distance_matrix_two_samples() {
        let file1 = PhrayaFile::new(
            100,
            "s1".to_string(),
            "2024-01-01T00:00:00Z".to_string(),
            vec![
                make_snp_obs(10, b'A', b'T', 8, 2), // s1 calls T
                make_snp_obs(20, b'C', b'G', 5, 5), // s1 calls C
            ],
            phraya_core::types::CoverageTrack::uncovered(100),
        );
        let file2 = PhrayaFile::new(
            100,
            "s2".to_string(),
            "2024-01-01T00:00:00Z".to_string(),
            vec![
                make_snp_obs(10, b'A', b'T', 1, 9), // s2 calls A (different from T)
                make_snp_obs(20, b'C', b'G', 7, 3), // s2 calls G (different from C)
            ],
            phraya_core::types::CoverageTrack::uncovered(100),
        );

        let output = format_distance_matrix(
            &[file1, file2],
            &["s1".to_string(), "s2".to_string()],
        );
        let lines: Vec<&str> = output.lines().collect();
        assert_eq!(lines.len(), 3); // header + 2 rows
        assert!(lines[0].starts_with('\t'));
        assert!(lines[0].contains("s1"));
        assert!(lines[0].contains("s2"));

        let row0: Vec<&str> = lines[1].split('\t').collect();
        assert_eq!(row0[0], "s1");
        assert_eq!(row0[1], "0");   // s1 vs s1
        assert_eq!(row0[2], "2");   // s1 vs s2: 2 differences

        let row1: Vec<&str> = lines[2].split('\t').collect();
        assert_eq!(row1[0], "s2");
        assert_eq!(row1[1], "2");   // s2 vs s1: 2 differences
        assert_eq!(row1[2], "0");   // s2 vs s2
    }

    #[test]
    fn test_format_distance_matrix_symmetric() {
        let file1 = PhrayaFile::new(
            100,
            "a".to_string(),
            "2024-01-01T00:00:00Z".to_string(),
            vec![
                make_snp_obs(10, b'A', b'T', 8, 2), // T
                make_snp_obs(20, b'C', b'G', 5, 5), // C
                make_snp_obs(30, b'G', b'A', 10, 0), // A
            ],
            phraya_core::types::CoverageTrack::uncovered(100),
        );
        let file2 = PhrayaFile::new(
            100,
            "b".to_string(),
            "2024-01-01T00:00:00Z".to_string(),
            vec![
                make_snp_obs(10, b'A', b'T', 1, 9), // A
                make_snp_obs(20, b'C', b'G', 7, 3), // G (alt_count=7 > ref_count=3)
                make_snp_obs(30, b'G', b'A', 10, 0), // A (same as a: alt_count=10 > ref_count=0)
            ],
            phraya_core::types::CoverageTrack::uncovered(100),
        );
        let file3 = PhrayaFile::new(
            100,
            "c".to_string(),
            "2024-01-01T00:00:00Z".to_string(),
            vec![
                make_snp_obs(10, b'A', b'T', 8, 2), // T (same as a)
                make_snp_obs(20, b'C', b'G', 5, 5), // C (same as a)
                make_snp_obs(30, b'G', b'A', 10, 0), // A (same as a)
            ],
            phraya_core::types::CoverageTrack::uncovered(100),
        );

        let output = format_distance_matrix(
            &[file1, file2, file3],
            &["a".to_string(), "b".to_string(), "c".to_string()],
        );
        let lines: Vec<&str> = output.lines().collect();
        assert_eq!(lines.len(), 4); // header + 3 rows

        // a vs b: positions 10,20 differ, pos 30 same → 2
        // a vs c: all same → 0
        // b vs c: positions 10,20 differ → 2
        let a_row: Vec<&str> = lines[1].split('\t').collect();
        assert_eq!(a_row[0], "a");
        assert_eq!(a_row[1], "0");  // a vs a
        assert_eq!(a_row[2], "2");  // a vs b
        assert_eq!(a_row[3], "0");  // a vs c

        let b_row: Vec<&str> = lines[2].split('\t').collect();
        assert_eq!(b_row[1], "2");  // b vs a
        assert_eq!(b_row[2], "0");  // b vs b
        assert_eq!(b_row[3], "2");  // b vs c

        let c_row: Vec<&str> = lines[3].split('\t').collect();
        assert_eq!(c_row[1], "0");  // c vs a
        assert_eq!(c_row[2], "2");  // c vs b
        assert_eq!(c_row[3], "0");  // c vs c
    }

    #[test]
    fn test_format_distance_matrix_excludes_missing() {
        // Position 30 only in s1 → excluded from distance (N in s2)
        let file1 = PhrayaFile::new(
            100,
            "s1".to_string(),
            "2024-01-01T00:00:00Z".to_string(),
            vec![
                make_snp_obs(10, b'A', b'T', 8, 2), // T
                make_snp_obs(30, b'G', b'A', 10, 0), // A
            ],
            phraya_core::types::CoverageTrack::uncovered(100),
        );
        let file2 = PhrayaFile::new(
            100,
            "s2".to_string(),
            "2024-01-01T00:00:00Z".to_string(),
            vec![
                make_snp_obs(10, b'A', b'T', 1, 9), // A
            ],
            phraya_core::types::CoverageTrack::uncovered(100),
        );

        // core SNPs = only position 10 (present in both)
        // Position 30 is NOT a core SNP (only in s1)
        let output = format_distance_matrix(
            &[file1, file2],
            &["s1".to_string(), "s2".to_string()],
        );

        let row1: Vec<&str> = output.lines().nth(1).unwrap().split('\t').collect();
        // Distance = 1 (T vs A at position 10)
        assert_eq!(row1[2], "1");
    }

    #[test]
    fn test_compute_hamming_distance_basic() {
        assert_eq!(compute_hamming_distance(b"ACGT", b"ACGT"), 0);
        assert_eq!(compute_hamming_distance(b"ACGT", b"AGGT"), 1);
        assert_eq!(compute_hamming_distance(b"ACGT", b"TGCA"), 4);
    }

    #[test]
    fn test_compute_hamming_distance_excludes_n() {
        // N in either position → not counted
        assert_eq!(compute_hamming_distance(b"ACGT", b"ACNT"), 0); // 4th is N in b
        assert_eq!(compute_hamming_distance(b"ACNT", b"ACGT"), 0); // 4th is N in a
        assert_eq!(compute_hamming_distance(b"ACGT", b"ACGN"), 0); // 4th is N in b
        assert_eq!(compute_hamming_distance(b"ACGT", b"TGCA"), 4); // no N, 4 differences
    }

    #[test]
    fn test_called_base_aggregation_across_observations() {
        // Two observations at the same position, same sample:
        // obs1: T(8), A(2) → T
        // obs2: G(3), A(7) → A
        // Aggregated: A(9), T(8), G(3) → A
        let obs1 = VariantObservation::new(
            50, b'A',
            { let mut m = HashMap::new(); m.insert(b'T', 8); m.insert(b'A', 2); m },
            0.95, "5M".to_string(), 60, 1, vec![10], 35.0,
            "sample:read1".to_string(),
        )
        .with_variant_type(VariantType::Snp)
        .with_strand(Strand::Forward)
        .with_query_position(2)
        .with_coverage_window_offset(5);

        let obs2 = VariantObservation::new(
            50, b'A',
            { let mut m = HashMap::new(); m.insert(b'G', 3); m.insert(b'A', 7); m },
            0.95, "5M".to_string(), 60, 1, vec![10], 35.0,
            "sample:read2".to_string(),
        )
        .with_variant_type(VariantType::Snp)
        .with_strand(Strand::Forward)
        .with_query_position(2)
        .with_coverage_window_offset(5);

        let file = PhrayaFile::new(
            100, "s1".to_string(), "2024-01-01T00:00:00Z".to_string(),
            vec![obs1, obs2],
            phraya_core::types::CoverageTrack::uncovered(100),
        );

        let snps = compute_core_snps(&[file], &["s1".to_string()]);
        assert_eq!(snps.len(), 1);
        // Aggregated: A(9) > T(8) → called base = A
        assert_eq!(snps[0].sample_bases, vec![b'A']);
    }
}
