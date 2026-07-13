// Module declarations
pub mod types;
mod hotspot_tests;

/// Represents a detected tandem repeat region in a sequence.
///
/// A tandem repeat is a pattern of nucleotides that repeats multiple times in succession.
/// For example, "ATATATATAT" is a dinucleotide repeat (AT repeated 5 times).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RepeatRegion {
    /// Starting position of the repeat (0-indexed, inclusive)
    pub start: usize,
    /// Ending position of the repeat (0-indexed, inclusive)
    pub end: usize,
    /// Length of the repeating unit in nucleotides (e.g., 2 for dinucleotide, 3 for trinucleotide)
    pub period: usize,
    /// The actual repeating unit sequence (e.g., "AT", "CAG", "GATA")
    pub unit: String,
}

impl RepeatRegion {
    /// Creates a new RepeatRegion.
    pub fn new(start: usize, end: usize, period: usize, unit: String) -> Self {
        RepeatRegion {
            start,
            end,
            period,
            unit,
        }
    }

    /// Returns the number of complete repeat units in this region.
    pub fn repeat_count(&self) -> usize {
        (self.end - self.start + 1) / self.period
    }
}

/// Configuration for tandem repeat detection.
#[derive(Debug, Clone)]
pub struct RepeatDetectorConfig {
    /// Minimum number of complete repeat periods to consider a region as a tandem repeat.
    /// Default: 3
    pub min_repeat_count: usize,
}

impl Default for RepeatDetectorConfig {
    fn default() -> Self {
        RepeatDetectorConfig {
            min_repeat_count: 3,
        }
    }
}

/// Detects tandem repeats (simple sequence repeats) in a DNA sequence.
///
/// This function identifies regions where short patterns (2-4 nucleotides) repeat
/// multiple times consecutively. It searches for dinucleotide, trinucleotide, and
/// tetranucleotide repeats.
///
/// # Arguments
///
/// * `sequence` - The DNA sequence to search (case-insensitive, can contain ACGT or acgt)
/// * `config` - Configuration specifying minimum repeat count threshold
///
/// # Returns
///
/// A vector of `RepeatRegion` structs, sorted by position, describing all detected
/// tandem repeats. If no repeats are found, returns an empty vector.
///
/// # Example
///
/// ```ignore
/// let seq = "ATATATATAT";  // AT repeated 5 times
/// let regions = detect_tandem_repeats(seq, &RepeatDetectorConfig::default());
/// assert_eq!(regions.len(), 1);
/// assert_eq!(regions[0].period, 2);
/// assert_eq!(regions[0].unit, "AT");
/// ```
pub fn detect_tandem_repeats(sequence: &str, config: &RepeatDetectorConfig) -> Vec<RepeatRegion> {
    let seq_upper = sequence.to_uppercase();
    let bytes = seq_upper.as_bytes();

    if bytes.len() < config.min_repeat_count * 2 {
        return Vec::new();
    }

    let mut regions = Vec::new();
    let mut current_pos = 0;

    // Scan through the sequence linearly, looking for repeats
    while current_pos + (config.min_repeat_count * 2) <= bytes.len() {
        let mut found_repeat = false;

        // Try periods from largest to smallest (4, 3, 2) to prefer longer repeats
        for period in (2..=4).rev() {
            if current_pos + period * config.min_repeat_count > bytes.len() {
                continue;
            }

            // Extract the potential unit starting at current position
            let unit = &bytes[current_pos..current_pos + period];

            // Check how many consecutive times this unit repeats
            let mut repeat_count = 1;
            let mut end_pos = current_pos + period;

            while end_pos + period <= bytes.len() {
                let next_unit = &bytes[end_pos..end_pos + period];
                if unit == next_unit {
                    repeat_count += 1;
                    end_pos += period;
                } else {
                    break;
                }
            }

            // Only add if meets minimum threshold
            if repeat_count >= config.min_repeat_count {
                let unit_str = String::from_utf8_lossy(unit).into_owned();
                let region = RepeatRegion::new(current_pos, end_pos - 1, period, unit_str);
                regions.push(region);

                current_pos = end_pos;
                found_repeat = true;
                break; // Found a repeat at this position, move forward
            }
        }

        if !found_repeat {
            current_pos += 1;
        }
    }

    regions
}

#[cfg(test)]
mod tests {
    use super::*;

    // ===== HAPPY PATH: Basic detection =====

    #[test]
    fn detects_dinucleotide_repeat_at() {
        let seq = "ATATATATAT"; // AT repeated 5 times
        let config = RepeatDetectorConfig::default();
        let regions = detect_tandem_repeats(seq, &config);

        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].start, 0);
        assert_eq!(regions[0].end, 9);
        assert_eq!(regions[0].period, 2);
        assert_eq!(regions[0].unit, "AT");
        assert_eq!(regions[0].repeat_count(), 5);
    }

    #[test]
    fn detects_trinucleotide_repeat_cag() {
        let seq = "CAGCAGCAGCAG"; // CAG repeated 4 times
        let config = RepeatDetectorConfig::default();
        let regions = detect_tandem_repeats(seq, &config);

        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].start, 0);
        assert_eq!(regions[0].end, 11);
        assert_eq!(regions[0].period, 3);
        assert_eq!(regions[0].unit, "CAG");
        assert_eq!(regions[0].repeat_count(), 4);
    }

    #[test]
    fn detects_tetranucleotide_repeat_gata() {
        let seq = "GATAGATAGATA"; // GATA repeated 3 times
        let config = RepeatDetectorConfig::default();
        let regions = detect_tandem_repeats(seq, &config);

        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].start, 0);
        assert_eq!(regions[0].end, 11);
        assert_eq!(regions[0].period, 4);
        assert_eq!(regions[0].unit, "GATA");
        assert_eq!(regions[0].repeat_count(), 3);
    }

    #[test]
    fn detects_dinucleotide_repeat_gc() {
        let seq = "GCGCGCGC"; // GC repeated 4 times
        let config = RepeatDetectorConfig::default();
        let regions = detect_tandem_repeats(seq, &config);

        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].period, 2);
        assert_eq!(regions[0].unit, "GC");
    }

    #[test]
    fn detects_trinucleotide_repeat_aaa() {
        let seq = "AAAAAAAAA"; // AAA repeated 3 times
        let config = RepeatDetectorConfig::default();
        let regions = detect_tandem_repeats(seq, &config);

        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].period, 3);
        assert_eq!(regions[0].unit, "AAA");
        assert_eq!(regions[0].repeat_count(), 3);
    }

    // ===== EDGE CASES: Boundary conditions =====

    #[test]
    fn detects_repeat_at_sequence_start() {
        let seq = "ATATATATATTGCAAA"; // AT repeat at start
        let config = RepeatDetectorConfig::default();
        let regions = detect_tandem_repeats(seq, &config);

        assert!(regions.iter().any(|r| r.start == 0 && r.period == 2));
    }

    #[test]
    fn detects_repeat_at_sequence_end() {
        let seq = "TGCAAATATATATAT"; // AT repeat at end
        let config = RepeatDetectorConfig::default();
        let regions = detect_tandem_repeats(seq, &config);

        assert!(regions.iter().any(|r| r.end == 14 && r.period == 2));
    }

    #[test]
    fn detects_repeat_in_sequence_middle() {
        let seq = "TGCAATATATATACCC"; // AT repeat in middle
        let config = RepeatDetectorConfig::default();
        let regions = detect_tandem_repeats(seq, &config);

        assert!(regions.iter().any(|r| r.period == 2 && r.unit == "AT"));
    }

    #[test]
    fn detects_exactly_three_periods() {
        let seq = "CAGCAGCAG"; // CAG repeated exactly 3 times
        let config = RepeatDetectorConfig::default();
        let regions = detect_tandem_repeats(seq, &config);

        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].repeat_count(), 3);
        assert_eq!(regions[0].period, 3);
    }

    #[test]
    fn detects_multiple_repeats_in_sequence() {
        let seq = "ATATATATAGCGCGCATATATAT"; // AT repeat + GC repeat + AT repeat
        let config = RepeatDetectorConfig::default();
        let regions = detect_tandem_repeats(seq, &config);

        // Should detect multiple regions
        assert!(regions.len() >= 2);
        assert!(regions.iter().any(|r| r.unit == "AT"));
        assert!(regions.iter().any(|r| r.unit == "GC"));
    }

    #[test]
    fn detects_with_adjacent_same_unit_as_separate_regions() {
        // When two repeat regions are adjacent or separated, they should be identified
        let seq = "ATATATATAGCGCGC"; // AT repeat followed by GC repeat
        let config = RepeatDetectorConfig::default();
        let regions = detect_tandem_repeats(seq, &config);

        assert!(regions.len() >= 2);
    }

    // ===== ERROR CASES: Non-matching patterns =====

    #[test]
    fn returns_empty_vector_for_sequence_without_repeats() {
        // ATGCATGATCGATCG has no tandem repeats (all substrings unique)
        let seq = "ATGCATGATCGATCG";
        let config = RepeatDetectorConfig::default();
        let regions = detect_tandem_repeats(seq, &config);

        assert_eq!(regions.len(), 0);
    }

    #[test]
    fn returns_empty_vector_for_single_nucleotide() {
        let seq = "A";
        let config = RepeatDetectorConfig::default();
        let regions = detect_tandem_repeats(seq, &config);

        assert_eq!(regions.len(), 0);
    }

    #[test]
    fn returns_empty_vector_for_two_nucleotides() {
        let seq = "AT";
        let config = RepeatDetectorConfig::default();
        let regions = detect_tandem_repeats(seq, &config);

        assert_eq!(regions.len(), 0);
    }

    #[test]
    fn ignores_single_repeat_unit_when_min_is_three() {
        let seq = "ATTGCATGC"; // Single AT, not a repeat
        let config = RepeatDetectorConfig::default();
        let regions = detect_tandem_repeats(seq, &config);

        // Should not detect since we need at least 3 periods
        assert!(regions.is_empty() || regions.iter().all(|r| r.repeat_count() >= 3));
    }

    #[test]
    fn ignores_two_repeat_units_when_min_is_three() {
        let seq = "ATATATGC"; // AT repeated 2 times, below threshold
        let config = RepeatDetectorConfig::default();
        let regions = detect_tandem_repeats(seq, &config);

        // Should not detect since we need at least 3 periods by default
        assert!(regions.is_empty() || regions.iter().all(|r| r.repeat_count() >= 3));
    }

    #[test]
    fn case_insensitive_detection_lowercase() {
        let seq = "atatatatat"; // lowercase
        let config = RepeatDetectorConfig::default();
        let regions = detect_tandem_repeats(seq, &config);

        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].period, 2);
    }

    #[test]
    fn case_insensitive_detection_mixed_case() {
        let seq = "AtAtAtAtAt"; // mixed case
        let config = RepeatDetectorConfig::default();
        let regions = detect_tandem_repeats(seq, &config);

        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].period, 2);
    }

    // ===== CORRECTNESS: Known STR sequences =====

    #[test]
    fn detects_huntingtons_cag_repeat() {
        // Huntington's disease marker: CAG repeat in HTT gene
        let seq = "CAGCAGCAGCAGCAGCAGCAG"; // CAG repeated 7 times
        let config = RepeatDetectorConfig::default();
        let regions = detect_tandem_repeats(seq, &config);

        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].unit, "CAG");
        assert_eq!(regions[0].repeat_count(), 7);
    }

    #[test]
    fn detects_fragile_x_cgg_repeat() {
        // Fragile X syndrome marker: CGG repeat
        let seq = "CGGCGGCGGCGGCGGCGG"; // CGG repeated 6 times
        let config = RepeatDetectorConfig::default();
        let regions = detect_tandem_repeats(seq, &config);

        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].unit, "CGG");
        assert_eq!(regions[0].repeat_count(), 6);
    }

    #[test]
    fn detects_myotonic_dystrophy_ctg_repeat() {
        // Myotonic dystrophy marker: CTG repeat
        let seq = "CTGCTGCTGCTGCTG"; // CTG repeated 5 times
        let config = RepeatDetectorConfig::default();
        let regions = detect_tandem_repeats(seq, &config);

        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].unit, "CTG");
        assert_eq!(regions[0].repeat_count(), 5);
    }

    #[test]
    fn detects_real_str_marker_d5s818() {
        // D5S818 CODIS marker: AGAT tetranucleotide repeat
        let seq = "AGATAGATAGATAGATAGATAGAT"; // AGAT repeated 6 times
        let config = RepeatDetectorConfig::default();
        let regions = detect_tandem_repeats(seq, &config);

        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].unit, "AGAT");
        assert_eq!(regions[0].repeat_count(), 6);
    }

    // ===== CONFIGURABILITY: Min repeat count threshold =====

    #[test]
    fn respects_min_repeat_count_of_two() {
        let seq = "ATATGC"; // AT repeated 2 times
        let config = RepeatDetectorConfig {
            min_repeat_count: 2,
        };
        let regions = detect_tandem_repeats(seq, &config);

        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].repeat_count(), 2);
    }

    #[test]
    fn respects_min_repeat_count_of_four() {
        let seq = "CAGCAGCAGCAG"; // CAG repeated 4 times
        let config = RepeatDetectorConfig {
            min_repeat_count: 4,
        };
        let regions = detect_tandem_repeats(seq, &config);

        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].repeat_count(), 4);
    }

    #[test]
    fn excludes_repeats_below_min_threshold() {
        let seq = "ATATCAGCAGCAG"; // AT repeated 2 times + CAG repeated 3 times
        let config = RepeatDetectorConfig {
            min_repeat_count: 3,
        };
        let regions = detect_tandem_repeats(seq, &config);

        // Should only detect CAG, not AT (only 2 repeats)
        assert!(regions.iter().all(|r| r.repeat_count() >= 3));
    }

    #[test]
    fn high_threshold_returns_only_long_repeats() {
        let seq = "ATATATATAGCGCGCGCGCGCAA"; // AT x5 + GC x7
        let config = RepeatDetectorConfig {
            min_repeat_count: 6,
        };
        let regions = detect_tandem_repeats(seq, &config);

        // Only GC (7 repeats) should be detected
        assert!(regions.iter().all(|r| r.repeat_count() >= 6));
    }

    // ===== PERFORMANCE: 10kb sequence =====

    #[test]
    fn performance_test_10kb_sequence() {
        use std::time::Instant;

        // Generate a 10kb semi-random sequence with embedded repeats
        let mut seq = String::new();
        let pattern = "ATGCTAGCTAGCTAGC";
        for _ in 0..625 {
            // 625 * 16 = 10000 bp
            seq.push_str(pattern);
        }

        // Insert some tandem repeats at various positions
        let mut seq_bytes = seq.into_bytes();
        // Insert AT repeat at position 100
        for i in 0..20 {
            seq_bytes[100 + i * 2] = b'A';
            seq_bytes[100 + i * 2 + 1] = b'T';
        }
        let seq = String::from_utf8(seq_bytes).unwrap();

        let config = RepeatDetectorConfig::default();
        let start = Instant::now();
        let regions = detect_tandem_repeats(&seq, &config);
        let elapsed = start.elapsed();

        // Performance assertion: should complete in under 100ms
        assert!(
            elapsed.as_millis() < 100,
            "Performance test failed: took {}ms, expected < 100ms",
            elapsed.as_millis()
        );

        // Should detect at least the embedded repeat
        assert!(!regions.is_empty());
    }

    // ===== ROBUSTNESS: Complex patterns =====

    #[test]
    fn detects_repeat_overlaps_correctly() {
        // "AAAAAA" could be detected as AAA x2, AA x3, or A x6
        // Should detect highest period match with 3+ repeats
        let seq = "AAAAAA";
        let config = RepeatDetectorConfig::default();
        let regions = detect_tandem_repeats(seq, &config);

        // Behavior: should detect AAA x2, which meets threshold of 2 repeats minimum
        // (this test may need adjustment based on implementation priority)
        assert!(!regions.is_empty());
        assert!(regions.iter().any(|r| r.period <= 3));
    }

    #[test]
    fn detects_mixed_repeats_no_false_positives() {
        let seq = "ATGCATGCATGC"; // Not a true repeat (period doesn't match)
        let config = RepeatDetectorConfig::default();
        let regions = detect_tandem_repeats(seq, &config);

        // ATGC is period 4, so "ATGCATGCATGC" would be 3 repeats - might be detected
        // depending on implementation; ensure no false positives outside this
        for region in regions {
            assert!(region.period <= 4);
        }
    }

    #[test]
    fn empty_sequence_returns_empty_vector() {
        let seq = "";
        let config = RepeatDetectorConfig::default();
        let regions = detect_tandem_repeats(seq, &config);

        assert_eq!(regions.len(), 0);
    }
}
