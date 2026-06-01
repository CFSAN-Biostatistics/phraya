use crate::VariantObservation;

/// TSV output column specification
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Column {
    Position,
    RefBase,
    AltBases,
    Coverage,
    Mapq,
    Cigar,
    EditDistance,
    Confidence,
    Provenance,
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    fn create_observation(
        position: u32,
        ref_base: u8,
        alleles: HashMap<u8, u32>,
        mapq: u8,
        provenance: &str,
    ) -> VariantObservation {
        // Extract coverage from the alleles map (sum of all counts)
        let coverage: u32 = alleles.values().sum();
        VariantObservation::new(
            position,
            ref_base,
            alleles,
            0.95,
            "10M".to_string(),
            mapq,
            0,
            vec![coverage],
            35.0,
            provenance.to_string(),
        )
    }

    #[test]
    fn issue_84_single_observation_all_columns() {
        let mut alleles = HashMap::new();
        alleles.insert(b'A', 10);
        alleles.insert(b'T', 5);

        let obs = create_observation(100, b'A', alleles, 60, "sample1:read1");

        let columns = [
            Column::Position,
            Column::RefBase,
            Column::AltBases,
            Column::Coverage,
            Column::Mapq,
            Column::Cigar,
            Column::EditDistance,
            Column::Confidence,
            Column::Provenance,
        ];

        let tsv = format_tsv(std::iter::once(obs), &columns);

        // Check header
        let lines: Vec<&str> = tsv.lines().collect();
        assert_eq!(lines.len(), 2, "Should have header + 1 data row");
        assert_eq!(
            lines[0],
            "Position\tRefBase\tAltBases\tCoverage\tMapq\tCigar\tEditDistance\tConfidence\tProvenance"
        );

        // Check data row - coverage is 10+5=15 from alleles
        let data_line = lines[1];
        assert!(data_line.contains("100\tA\tT\t15\t60"));
    }

    #[test]
    fn issue_84_column_selection_position_coverage_only() {
        let mut alleles = HashMap::new();
        alleles.insert(b'A', 15);

        let obs = create_observation(50, b'A', alleles, 45, "sample:read");

        let columns = [Column::Position, Column::Coverage];
        let tsv = format_tsv(std::iter::once(obs), &columns);

        let lines: Vec<&str> = tsv.lines().collect();
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0], "Position\tCoverage");
        assert_eq!(lines[1], "50\t15");
    }

    #[test]
    fn issue_84_default_columns() {
        let mut alleles = HashMap::new();
        alleles.insert(b'A', 12);
        alleles.insert(b'G', 3);

        let obs = create_observation(200, b'A', alleles, 55, "sample:read");

        let tsv = format_tsv_with_defaults(std::iter::once(obs));

        let lines: Vec<&str> = tsv.lines().collect();
        assert_eq!(lines.len(), 2);
        // Default: Position, RefBase, AltBases, Coverage, Mapq
        // Coverage is 12+3=15 from alleles
        assert_eq!(lines[0], "Position\tRefBase\tAltBases\tCoverage\tMapq");
        assert_eq!(lines[1], "200\tA\tG\t15\t55");
    }

    #[test]
    fn issue_84_special_characters_escaped_tab() {
        let mut alleles = HashMap::new();
        alleles.insert(b'A', 10);

        let provenance_with_tab = "sample1\tread42";
        let obs = create_observation(100, b'A', alleles, 60, provenance_with_tab);

        let columns = [Column::Provenance];
        let tsv = format_tsv(std::iter::once(obs), &columns);

        let lines: Vec<&str> = tsv.lines().collect();
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0], "Provenance");
        // Tab should be escaped as \\t
        assert_eq!(lines[1], "sample1\\tread42");
    }

    #[test]
    fn issue_84_special_characters_escaped_newline() {
        let mut alleles = HashMap::new();
        alleles.insert(b'A', 10);

        let provenance_with_newline = "sample1\nread42";
        let obs = create_observation(100, b'A', alleles, 60, provenance_with_newline);

        let columns = [Column::Provenance];
        let tsv = format_tsv(std::iter::once(obs), &columns);

        let lines: Vec<&str> = tsv.lines().collect();
        // Note: the actual newline is escaped, so we should get one header + one data line
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0], "Provenance");
        assert_eq!(lines[1], "sample1\\nread42");
    }

    #[test]
    fn issue_84_multi_allelic_site_sorted_alts() {
        let mut alleles = HashMap::new();
        alleles.insert(b'A', 10);
        alleles.insert(b'G', 5);
        alleles.insert(b'C', 3);

        let obs = create_observation(75, b'A', alleles, 50, "sample:read");

        let columns = [Column::AltBases];
        let tsv = format_tsv(std::iter::once(obs), &columns);

        let lines: Vec<&str> = tsv.lines().collect();
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0], "AltBases");
        // Alleles should be sorted: C, G
        assert_eq!(lines[1], "CG");
    }

    #[test]
    fn issue_84_no_alt_alleles_shows_dot() {
        let mut alleles = HashMap::new();
        alleles.insert(b'A', 20);

        let obs = create_observation(100, b'A', alleles, 60, "sample:read");

        let columns = [Column::AltBases];
        let tsv = format_tsv(std::iter::once(obs), &columns);

        let lines: Vec<&str> = tsv.lines().collect();
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0], "AltBases");
        assert_eq!(lines[1], ".");
    }

    #[test]
    fn issue_84_zero_coverage_when_empty() {
        let mut alleles = HashMap::new();
        alleles.insert(b'A', 5);

        let mut obs = create_observation(100, b'A', alleles, 60, "sample:read");
        // Manually create observation with empty local_coverage
        obs = VariantObservation::new(
            100,
            b'A',
            [(b'A', 5)].into_iter().collect(),
            0.95,
            "10M".to_string(),
            60,
            0,
            vec![], // Empty local coverage
            35.0,
            "sample:read".to_string(),
        );

        let columns = [Column::Coverage];
        let tsv = format_tsv(std::iter::once(obs), &columns);

        let lines: Vec<&str> = tsv.lines().collect();
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0], "Coverage");
        assert_eq!(lines[1], "0");
    }

    #[test]
    fn issue_84_multiple_observations() {
        let mut alleles1 = HashMap::new();
        alleles1.insert(b'A', 10);
        alleles1.insert(b'T', 5);
        let obs1 = create_observation(100, b'A', alleles1, 60, "sample1:read1");

        let mut alleles2 = HashMap::new();
        alleles2.insert(b'C', 8);
        let obs2 = create_observation(150, b'C', alleles2, 50, "sample2:read2");

        let columns = [Column::Position, Column::RefBase, Column::Coverage];
        let tsv = format_tsv(vec![obs1, obs2].into_iter(), &columns);

        let lines: Vec<&str> = tsv.lines().collect();
        assert_eq!(lines.len(), 3, "Should have header + 2 data rows");
        assert_eq!(lines[0], "Position\tRefBase\tCoverage");
        assert_eq!(lines[1], "100\tA\t15"); // Coverage is 10+5=15
        assert_eq!(lines[2], "150\tC\t8");
    }

    #[test]
    fn issue_84_empty_observations() {
        let columns = [Column::Position, Column::RefBase];
        let tsv = format_tsv(std::iter::empty(), &columns);

        let lines: Vec<&str> = tsv.lines().collect();
        // Should have header only
        assert_eq!(lines.len(), 1);
        assert_eq!(lines[0], "Position\tRefBase");
    }

    #[test]
    fn issue_84_confidence_field() {
        let mut alleles = HashMap::new();
        alleles.insert(b'A', 10);

        let obs = VariantObservation::new(
            100,
            b'A',
            alleles,
            0.85,
            "10M".to_string(),
            60,
            2,
            vec![10],
            35.0,
            "sample:read".to_string(),
        );

        let columns = [Column::Confidence];
        let tsv = format_tsv(std::iter::once(obs), &columns);

        let lines: Vec<&str> = tsv.lines().collect();
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0], "Confidence");
        assert_eq!(lines[1], "0.85");
    }

    #[test]
    fn issue_84_edit_distance_field() {
        let mut alleles = HashMap::new();
        alleles.insert(b'A', 10);

        let obs = VariantObservation::new(
            100,
            b'A',
            alleles,
            0.95,
            "8M1D1M".to_string(),
            60,
            1,
            vec![10],
            35.0,
            "sample:read".to_string(),
        );

        let columns = [Column::EditDistance];
        let tsv = format_tsv(std::iter::once(obs), &columns);

        let lines: Vec<&str> = tsv.lines().collect();
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0], "EditDistance");
        assert_eq!(lines[1], "1");
    }

    #[test]
    fn issue_84_cigar_field() {
        let mut alleles = HashMap::new();
        alleles.insert(b'A', 10);

        let obs = VariantObservation::new(
            100,
            b'A',
            alleles,
            0.95,
            "5M2I3M".to_string(),
            60,
            0,
            vec![10],
            35.0,
            "sample:read".to_string(),
        );

        let columns = [Column::Cigar];
        let tsv = format_tsv(std::iter::once(obs), &columns);

        let lines: Vec<&str> = tsv.lines().collect();
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0], "Cigar");
        assert_eq!(lines[1], "5M2I3M");
    }

    #[test]
    fn issue_84_position_is_0_indexed() {
        let mut alleles = HashMap::new();
        alleles.insert(b'A', 10);

        let obs = create_observation(0, b'A', alleles, 60, "sample:read");

        let columns = [Column::Position];
        let tsv = format_tsv(std::iter::once(obs), &columns);

        let lines: Vec<&str> = tsv.lines().collect();
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[1], "0");
    }

    #[test]
    fn issue_84_ref_base_as_char() {
        let mut alleles = HashMap::new();
        alleles.insert(b'G', 15);

        let obs = create_observation(100, b'G', alleles, 60, "sample:read");

        let columns = [Column::RefBase];
        let tsv = format_tsv(std::iter::once(obs), &columns);

        let lines: Vec<&str> = tsv.lines().collect();
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[1], "G");
    }
}
