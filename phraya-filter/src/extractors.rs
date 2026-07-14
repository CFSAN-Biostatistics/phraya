/// Feature extractors for computing derived metrics not directly stored in VariantObservation.
/// These enable filtering on computed properties like CIGAR complexity, allele frequency, and multi-mapping.
use crate::QueryIndex;
use std::collections::HashMap;

/// Extract the count of CIGAR operations (M/I/D) from a CIGAR string.
///
/// # Arguments
/// * `cigar` - CIGAR string (e.g., "50M", "10M5I10M5D25M")
///
/// # Returns
/// Number of M/I/D operations (e.g., "50M" → 1, "10M5I10M5D25M" → 5)
pub fn extract_cigar_ops(cigar: &str) -> usize {
    cigar
        .chars()
        .filter(|c| matches!(c, 'M' | 'I' | 'D'))
        .count()
}

/// Extract the allele frequency for a specific allele from the all_alleles map.
///
/// # Arguments
/// * `all_alleles` - HashMap of (base, count) for all alleles at position
/// * `allele` - the specific allele (base) to compute frequency for
///
/// # Returns
/// Frequency as f64 (0.0-1.0). Returns 0.0 if total is 0 (empty all_alleles).
pub fn extract_allele_frequency(all_alleles: &HashMap<u8, usize>, allele: u8) -> f64 {
    let total: usize = all_alleles.values().sum();
    if total == 0 {
        return 0.0;
    }
    let count = all_alleles.get(&allele).copied().unwrap_or(0);
    count as f64 / total as f64
}

/// Extract the multi-mapping fraction for a query at a specific reference position.
///
/// # Arguments
/// * `position` - the reference position (0-indexed)
/// * `query_index` - the QueryIndex (HashMap<query_id, Vec<(position, score)>>)
///
/// # Returns
/// Fraction of reads that multi-map (0.0-1.0) at this position.
/// Returns 0.0 if position has no multi-mapping reads (query_index doesn't contain it, or is empty).
pub fn extract_multi_map_fraction(position: u32, query_index: &QueryIndex) -> f64 {
    // Count total reads that have at least one alignment at this position
    let total_reads_at_position: usize = query_index
        .values()
        .filter(|alignments| alignments.iter().any(|&(p, _)| p == position))
        .count();

    if total_reads_at_position == 0 {
        return 0.0;
    }

    // Count reads that have multiple alignments overall (multi-map reads)
    let multi_map_reads: usize = query_index
        .values()
        .filter(|alignments| alignments.iter().any(|&(p, _)| p == position) && alignments.len() > 1)
        .count();

    multi_map_reads as f64 / total_reads_at_position as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    // ===== CIGAR Operation Counter Tests =====

    #[test]
    fn issue_82_extract_cigar_ops_single_match() {
        let cigar = "50M";
        assert_eq!(extract_cigar_ops(cigar), 1);
    }

    #[test]
    fn issue_82_extract_cigar_ops_multiple_operations() {
        let cigar = "10M5I10M5D25M";
        assert_eq!(extract_cigar_ops(cigar), 5);
    }

    #[test]
    fn issue_82_extract_cigar_ops_complex_cigar() {
        let cigar = "5M2D8M1I3M2D10M";
        assert_eq!(extract_cigar_ops(cigar), 7);
    }

    #[test]
    fn issue_82_extract_cigar_ops_only_matches() {
        let cigar = "100M";
        assert_eq!(extract_cigar_ops(cigar), 1);
    }

    #[test]
    fn issue_82_extract_cigar_ops_only_indels() {
        let cigar = "5I3D2I";
        assert_eq!(extract_cigar_ops(cigar), 3);
    }

    #[test]
    fn issue_82_extract_cigar_ops_empty_string() {
        let cigar = "";
        assert_eq!(extract_cigar_ops(cigar), 0);
    }

    #[test]
    fn issue_82_extract_cigar_ops_no_mdis() {
        let cigar = "10N5S";
        assert_eq!(extract_cigar_ops(cigar), 0);
    }

    #[test]
    fn issue_82_extract_cigar_ops_mixed_operations() {
        let cigar = "2M3S5I7D10N2M";
        assert_eq!(extract_cigar_ops(cigar), 4); // 2M, 5I, 7D, 2M
    }

    // ===== Allele Frequency Tests =====

    #[test]
    fn issue_82_extract_allele_frequency_simple() {
        let mut alleles = HashMap::new();
        alleles.insert(b'A', 90);
        alleles.insert(b'C', 10);

        let freq_a = extract_allele_frequency(&alleles, b'A');
        let freq_c = extract_allele_frequency(&alleles, b'C');

        assert!((freq_a - 0.9).abs() < 1e-10);
        assert!((freq_c - 0.1).abs() < 1e-10);
    }

    #[test]
    fn issue_82_extract_allele_frequency_four_bases() {
        let mut alleles = HashMap::new();
        alleles.insert(b'A', 50);
        alleles.insert(b'C', 30);
        alleles.insert(b'G', 15);
        alleles.insert(b'T', 5);

        assert!((extract_allele_frequency(&alleles, b'A') - 0.5).abs() < 1e-10);
        assert!((extract_allele_frequency(&alleles, b'C') - 0.3).abs() < 1e-10);
        assert!((extract_allele_frequency(&alleles, b'G') - 0.15).abs() < 1e-10);
        assert!((extract_allele_frequency(&alleles, b'T') - 0.05).abs() < 1e-10);
    }

    #[test]
    fn issue_82_extract_allele_frequency_missing_allele() {
        let mut alleles = HashMap::new();
        alleles.insert(b'A', 80);
        alleles.insert(b'T', 20);

        let freq_g = extract_allele_frequency(&alleles, b'G');
        assert_eq!(freq_g, 0.0);
    }

    #[test]
    fn issue_82_extract_allele_frequency_empty_alleles() {
        let alleles = HashMap::new();

        let freq = extract_allele_frequency(&alleles, b'A');
        assert_eq!(freq, 0.0);
    }

    #[test]
    fn issue_82_extract_allele_frequency_single_base() {
        let mut alleles = HashMap::new();
        alleles.insert(b'A', 100);

        let freq_a = extract_allele_frequency(&alleles, b'A');
        assert_eq!(freq_a, 1.0);
    }

    #[test]
    fn issue_82_extract_allele_frequency_single_read() {
        let mut alleles = HashMap::new();
        alleles.insert(b'C', 1);

        let freq_c = extract_allele_frequency(&alleles, b'C');
        assert_eq!(freq_c, 1.0);
    }

    #[test]
    fn issue_82_extract_allele_frequency_equal_distribution() {
        let mut alleles = HashMap::new();
        alleles.insert(b'A', 25);
        alleles.insert(b'C', 25);
        alleles.insert(b'G', 25);
        alleles.insert(b'T', 25);

        let freq = extract_allele_frequency(&alleles, b'A');
        assert!((freq - 0.25).abs() < 1e-10);
    }

    // ===== Multi-Mapping Fraction Tests =====

    #[test]
    fn issue_82_extract_multi_map_fraction_no_multimap() {
        let mut query_index = HashMap::new();
        query_index.insert("read1".to_string(), vec![(100u32, 0.98)]);
        query_index.insert("read2".to_string(), vec![(100u32, 0.97)]);
        query_index.insert("read3".to_string(), vec![(100u32, 0.96)]);

        let fraction = extract_multi_map_fraction(100, &query_index);
        assert_eq!(fraction, 0.0);
    }

    #[test]
    fn issue_82_extract_multi_map_fraction_partial_multimap() {
        let mut query_index = HashMap::new();
        query_index.insert("read1".to_string(), vec![(100u32, 0.98), (200u32, 0.95)]);
        query_index.insert("read2".to_string(), vec![(100u32, 0.97)]);
        query_index.insert("read3".to_string(), vec![(100u32, 0.96)]);

        let fraction = extract_multi_map_fraction(100, &query_index);
        assert!((fraction - 1.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn issue_82_extract_multi_map_fraction_all_multimap() {
        let mut query_index = HashMap::new();
        query_index.insert("read1".to_string(), vec![(100u32, 0.98), (200u32, 0.95)]);
        query_index.insert("read2".to_string(), vec![(100u32, 0.97), (150u32, 0.94)]);
        query_index.insert("read3".to_string(), vec![(100u32, 0.96), (250u32, 0.92)]);

        let fraction = extract_multi_map_fraction(100, &query_index);
        assert_eq!(fraction, 1.0);
    }

    #[test]
    fn issue_82_extract_multi_map_fraction_empty_query_index() {
        let query_index = HashMap::new();

        let fraction = extract_multi_map_fraction(100, &query_index);
        assert_eq!(fraction, 0.0);
    }

    #[test]
    fn issue_82_extract_multi_map_fraction_position_not_covered() {
        let mut query_index = HashMap::new();
        query_index.insert("read1".to_string(), vec![(100u32, 0.98)]);
        query_index.insert("read2".to_string(), vec![(100u32, 0.97)]);

        let fraction = extract_multi_map_fraction(500, &query_index);
        assert_eq!(fraction, 0.0);
    }

    #[test]
    fn issue_82_extract_multi_map_fraction_mixed_positions() {
        let mut query_index = HashMap::new();
        query_index.insert("read1".to_string(), vec![(100u32, 0.98), (200u32, 0.95)]);
        query_index.insert("read2".to_string(), vec![(100u32, 0.97)]);
        query_index.insert("read3".to_string(), vec![(200u32, 0.96), (300u32, 0.92)]);

        // At position 100: read1 and read2, only read1 has multimap → 1/2
        let fraction_100 = extract_multi_map_fraction(100, &query_index);
        assert!((fraction_100 - 0.5).abs() < 1e-10);

        // At position 200: read1 and read3, both have multimap → 2/2
        let fraction_200 = extract_multi_map_fraction(200, &query_index);
        assert_eq!(fraction_200, 1.0);
    }

    #[test]
    fn issue_82_extract_multi_map_fraction_single_read() {
        let mut query_index = HashMap::new();
        query_index.insert("read1".to_string(), vec![(100u32, 0.98), (200u32, 0.95)]);

        let fraction = extract_multi_map_fraction(100, &query_index);
        assert_eq!(fraction, 1.0);
    }

    #[test]
    fn issue_82_extract_multi_map_fraction_three_mappings_one_query() {
        let mut query_index = HashMap::new();
        query_index.insert(
            "read1".to_string(),
            vec![(100u32, 0.98), (200u32, 0.95), (300u32, 0.92)],
        );
        query_index.insert("read2".to_string(), vec![(100u32, 0.97)]);

        let fraction = extract_multi_map_fraction(100, &query_index);
        assert!((fraction - 0.5).abs() < 1e-10);
    }

    #[test]
    fn issue_82_extract_multi_map_fraction_large_dataset() {
        let mut query_index = HashMap::new();

        // 100 reads at position 100
        for i in 0..100 {
            let query_id = format!("read_{}", i);
            if i < 70 {
                // 70 reads map only to position 100
                query_index.insert(query_id, vec![(100u32, 0.98)]);
            } else {
                // 30 reads map to position 100 and elsewhere
                query_index.insert(query_id, vec![(100u32, 0.97), (200u32, 0.95)]);
            }
        }

        let fraction = extract_multi_map_fraction(100, &query_index);
        assert!((fraction - 0.3).abs() < 1e-10);
    }
}
