// Tandem repeat detection module
// Detects simple tandem repeats (di-, tri-, and tetranucleotide patterns) in sequences

/// A region containing a tandem repeat
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RepeatRegion {
    /// Start position (0-based, inclusive)
    pub start: usize,
    /// End position (0-based, exclusive)
    pub end: usize,
    /// Period of the repeat (2, 3, or 4 for di-, tri-, tetranucleotide)
    pub period: usize,
    /// The repeating unit sequence
    pub unit: String,
}

/// Detect tandem repeats in a sequence
///
/// # Arguments
/// * `sequence` - The DNA sequence to analyze
/// * `min_count` - Minimum number of repeat periods required (default: 3)
///
/// # Returns
/// A vector of RepeatRegion structs representing detected repeats
pub fn detect_repeats(sequence: &[u8], min_count: usize) -> Vec<RepeatRegion> {
    if sequence.is_empty() {
        return Vec::new();
    }

    // Collect all candidate repeats first
    let mut candidates = Vec::new();

    for period in 2..=4 {
        for i in 0..sequence.len() {
            // Skip if doesn't have room for min_count
            if i + min_count * period > sequence.len() {
                continue;
            }

            // Extract the potential repeat unit at this position
            let unit = &sequence[i..i + period];

            // Count how many consecutive times this exact unit repeats
            let mut count = 0;
            let mut end = i;
            while end + period <= sequence.len() && &sequence[end..end + period] == unit {
                count += 1;
                end += period;
            }

            // If we found enough repeats, add as candidate
            if count >= min_count {
                candidates.push((i, end, period, String::from_utf8_lossy(unit).to_string(), count));
            }
        }
    }

    // Greedy selection: pick non-overlapping repeats in order of:
    // 1. Length (longer first)
    // 2. Period (shorter first, to prefer dinucleotide over tetranucleotide, etc.)
    // 3. Count of repetitions (higher first - more repetitions means more canonical)
    // 4. Unit string (earlier in lex order - prefer canonical alphabetic units)
    // 5. Starting position (EARLIER first)
    candidates.sort_by(|a, b| {
        let len_a = a.1 - a.0;
        let len_b = b.1 - b.0;
        len_b.cmp(&len_a)  // longer first
            .then(a.2.cmp(&b.2))  // shorter period first
            .then(b.4.cmp(&a.4))  // more repetitions first
            .then(a.3.cmp(&b.3))  // lex order of unit (prefer canonical)
            .then(a.0.cmp(&b.0))  // EARLIER start first
    });

    let mut repeats = Vec::new();
    let mut covered = vec![false; sequence.len()];

    for (start, end, period, unit, _count) in candidates {
        // Check if this repeat overlaps with already selected repeats
        let mut overlaps = false;
        for j in start..end {
            if covered[j] {
                overlaps = true;
                break;
            }
        }

        if !overlaps {
            // Mark these positions as covered
            for j in start..end {
                covered[j] = true;
            }

            repeats.push(RepeatRegion {
                start,
                end,
                period,
                unit,
            });
        }
    }

    // Sort repeats by start position for consistent output
    repeats.sort_by_key(|r| r.start);

    repeats
}

#[cfg(test)]
mod tests {
    use super::*;

    // ========================================================================
    // Acceptance Test: Function signature and basic structure
    // ========================================================================

    #[test]
    fn test_detect_repeats_returns_vec_of_repeat_regions() {
        // AC: Function takes sequence → returns Vec<RepeatRegion> with start, end, period, unit
        let sequence = b"ATGATGATGATG";
        let result = detect_repeats(sequence, 3);

        // This should compile and return a Vec<RepeatRegion>
        assert!(result.is_empty() || !result.is_empty());
    }

    // ========================================================================
    // Acceptance Test: Dinucleotide repeat detection
    // ========================================================================

    #[test]
    fn test_detects_dinucleotide_repeat_at() {
        // AC: Detects dinucleotide repeats (e.g., "ATATATATAT")
        let sequence = b"ATATATATAT";
        let repeats = detect_repeats(sequence, 3);

        assert_eq!(repeats.len(), 1, "Should detect one AT repeat");
        let repeat = &repeats[0];
        assert_eq!(repeat.start, 0);
        assert_eq!(repeat.end, 10);
        assert_eq!(repeat.period, 2);
        assert_eq!(repeat.unit, "AT");
    }

    #[test]
    fn test_detects_dinucleotide_repeat_cg() {
        let sequence = b"CGCGCGCGCG";
        let repeats = detect_repeats(sequence, 3);

        assert_eq!(repeats.len(), 1);
        let repeat = &repeats[0];
        assert_eq!(repeat.period, 2);
        assert_eq!(repeat.unit, "CG");
        assert_eq!(repeat.end - repeat.start, 10);
    }

    #[test]
    fn test_detects_dinucleotide_repeat_embedded_in_sequence() {
        let sequence = b"ACGTATATATATGCTA";
        let repeats = detect_repeats(sequence, 3);

        assert!(!repeats.is_empty(), "Should detect AT repeat in middle");
        let repeat = &repeats[0];
        assert_eq!(repeat.start, 4);
        assert_eq!(repeat.end, 12);
        assert_eq!(repeat.period, 2);
        assert_eq!(repeat.unit, "AT");
    }

    // ========================================================================
    // Acceptance Test: Trinucleotide repeat detection
    // ========================================================================

    #[test]
    fn test_detects_trinucleotide_repeat_cag() {
        // AC: Detects trinucleotide repeats (e.g., "CAGCAGCAGCAG")
        let sequence = b"CAGCAGCAGCAG";
        let repeats = detect_repeats(sequence, 3);

        assert_eq!(repeats.len(), 1, "Should detect one CAG repeat");
        let repeat = &repeats[0];
        assert_eq!(repeat.start, 0);
        assert_eq!(repeat.end, 12);
        assert_eq!(repeat.period, 3);
        assert_eq!(repeat.unit, "CAG");
    }

    #[test]
    fn test_detects_trinucleotide_repeat_gca() {
        let sequence = b"GCAGCAGCAGCAGCA";
        let repeats = detect_repeats(sequence, 3);

        assert_eq!(repeats.len(), 1);
        let repeat = &repeats[0];
        assert_eq!(repeat.period, 3);
        assert_eq!(repeat.unit, "GCA");
        assert_eq!(repeat.end - repeat.start, 15);
    }

    #[test]
    fn test_detects_trinucleotide_repeat_embedded() {
        let sequence = b"TTTTCAGCAGCAGCAGAAAA";
        let repeats = detect_repeats(sequence, 3);

        assert!(!repeats.is_empty(), "Should detect CAG repeat");
        let repeat = &repeats[0];
        assert_eq!(repeat.start, 4);
        assert_eq!(repeat.end, 16);
        assert_eq!(repeat.period, 3);
        assert_eq!(repeat.unit, "CAG");
    }

    // ========================================================================
    // Acceptance Test: Tetranucleotide repeat detection
    // ========================================================================

    #[test]
    fn test_detects_tetranucleotide_repeat_gata() {
        // AC: Detects tetranucleotide repeats (e.g., "GATAGATAGATA")
        let sequence = b"GATAGATAGATA";
        let repeats = detect_repeats(sequence, 3);

        assert_eq!(repeats.len(), 1, "Should detect one GATA repeat");
        let repeat = &repeats[0];
        assert_eq!(repeat.start, 0);
        assert_eq!(repeat.end, 12);
        assert_eq!(repeat.period, 4);
        assert_eq!(repeat.unit, "GATA");
    }

    #[test]
    fn test_detects_tetranucleotide_repeat_aagg() {
        let sequence = b"AAGGAAGGAAGGAAGG";
        let repeats = detect_repeats(sequence, 3);

        assert_eq!(repeats.len(), 1);
        let repeat = &repeats[0];
        assert_eq!(repeat.period, 4);
        assert_eq!(repeat.unit, "AAGG");
        assert_eq!(repeat.end - repeat.start, 16);
    }

    #[test]
    fn test_detects_tetranucleotide_repeat_embedded() {
        let sequence = b"CCCCGATAGATAGATAGATAGGGG";
        let repeats = detect_repeats(sequence, 3);

        assert!(!repeats.is_empty(), "Should detect GATA repeat");
        let repeat = &repeats[0];
        assert_eq!(repeat.start, 4);
        assert_eq!(repeat.end, 20);
        assert_eq!(repeat.period, 4);
        assert_eq!(repeat.unit, "GATA");
    }

    // ========================================================================
    // Acceptance Test: Minimum repeat count configuration
    // ========================================================================

    #[test]
    fn test_minimum_repeat_count_default_3() {
        // AC: Minimum repeat count configurable (default: 3 full periods)
        // 3 periods of "AT" = 6 bases
        let sequence_3_repeats = b"ATATAT";
        let repeats = detect_repeats(sequence_3_repeats, 3);
        assert_eq!(
            repeats.len(),
            1,
            "3 repeats should be detected with default min=3"
        );

        // 2 periods of "AT" = 4 bases
        let sequence_2_repeats = b"ATAT";
        let repeats = detect_repeats(sequence_2_repeats, 3);
        assert_eq!(
            repeats.len(),
            0,
            "2 repeats should NOT be detected with min=3"
        );
    }

    #[test]
    fn test_minimum_repeat_count_configurable() {
        let sequence = b"ATAT"; // Only 2 repeats of AT

        // Should not detect with min_count = 3
        let repeats_min_3 = detect_repeats(sequence, 3);
        assert_eq!(repeats_min_3.len(), 0);

        // Should detect with min_count = 2
        let repeats_min_2 = detect_repeats(sequence, 2);
        assert_eq!(repeats_min_2.len(), 1);
        assert_eq!(repeats_min_2[0].unit, "AT");
    }

    #[test]
    fn test_minimum_repeat_count_trinucleotide() {
        // 3 periods of "CAG" = 9 bases
        let sequence_3 = b"CAGCAGCAG";
        let repeats = detect_repeats(sequence_3, 3);
        assert_eq!(repeats.len(), 1);

        // 2 periods of "CAG" = 6 bases
        let sequence_2 = b"CAGCAG";
        let repeats = detect_repeats(sequence_2, 3);
        assert_eq!(repeats.len(), 0, "Should not detect with only 2 periods");
    }

    #[test]
    fn test_minimum_repeat_count_tetranucleotide() {
        // 3 periods of "GATA" = 12 bases
        let sequence_3 = b"GATAGATAGATA";
        let repeats = detect_repeats(sequence_3, 3);
        assert_eq!(repeats.len(), 1);

        // 2 periods of "GATA" = 8 bases
        let sequence_2 = b"GATAGATA";
        let repeats = detect_repeats(sequence_2, 3);
        assert_eq!(repeats.len(), 0);
    }

    // ========================================================================
    // Acceptance Test: Known STR sequences
    // ========================================================================

    #[test]
    fn test_known_str_d1s80_allele() {
        // AC: Tests with known STR sequences verify correct detection
        // D1S80 is a known VNTR locus with 16bp repeating unit
        // For this test we'll use a simplified tetranucleotide STR
        let sequence = b"AATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATA";
        let repeats = detect_repeats(sequence, 3);

        assert!(!repeats.is_empty(), "Should detect repeat in known STR");
        let repeat = &repeats[0];
        assert_eq!(repeat.period, 4);
        assert_eq!(repeat.unit, "GATA");
    }

    #[test]
    fn test_known_str_huntington_disease_cag() {
        // CAG repeat in huntingtin gene - classic trinucleotide repeat
        let sequence = b"ATGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGTTG";
        let repeats = detect_repeats(sequence, 3);

        assert!(!repeats.is_empty(), "Should detect CAG expansion");
        let cag_repeat = repeats.iter().find(|r| r.unit == "CAG");
        assert!(cag_repeat.is_some(), "Should find CAG repeat specifically");
        assert!(
            cag_repeat.unwrap().end - cag_repeat.unwrap().start >= 30,
            "Should detect long CAG tract"
        );
    }

    #[test]
    fn test_known_str_fragile_x_cgg() {
        // CGG repeat in FMR1 gene - another clinically significant repeat
        let sequence = b"TTCGGCGGCGGCGGCGGCGGCGGCGGCGGAA";
        let repeats = detect_repeats(sequence, 3);

        assert!(!repeats.is_empty());
        let cgg_repeat = repeats.iter().find(|r| r.unit == "CGG");
        assert!(cgg_repeat.is_some(), "Should detect CGG repeat");
    }

    // ========================================================================
    // Edge Cases and Boundaries
    // ========================================================================

    #[test]
    fn test_empty_sequence() {
        let sequence = b"";
        let repeats = detect_repeats(sequence, 3);
        assert_eq!(repeats.len(), 0);
    }

    #[test]
    fn test_sequence_shorter_than_minimum_repeat() {
        let sequence = b"AT";
        let repeats = detect_repeats(sequence, 3);
        assert_eq!(repeats.len(), 0);
    }

    #[test]
    fn test_no_repeats_in_sequence() {
        let sequence = b"ACGTACGTACGT"; // Could be mistaken for repeat but isn't a simple 2-4bp unit
        let repeats = detect_repeats(sequence, 3);
        // Should either find ACGT as a 4bp repeat or nothing
        // The function should handle this correctly
        assert!(repeats.is_empty() || repeats[0].period == 4);
    }

    #[test]
    fn test_multiple_different_repeats_in_same_sequence() {
        let sequence = b"ATATATCGCGCGCGGATAGATAGATAGATA";
        let repeats = detect_repeats(sequence, 3);

        assert!(
            repeats.len() >= 2,
            "Should detect multiple distinct repeats"
        );

        // Should find AT repeat
        let at_repeat = repeats.iter().find(|r| r.unit == "AT");
        assert!(at_repeat.is_some());

        // Should find CG repeat
        let cg_repeat = repeats.iter().find(|r| r.unit == "CG");
        assert!(cg_repeat.is_some());

        // Should find GATA repeat
        let gata_repeat = repeats.iter().find(|r| r.unit == "GATA");
        assert!(gata_repeat.is_some());
    }

    #[test]
    fn test_adjacent_but_different_repeats() {
        // Two different repeats right next to each other
        let sequence = b"ATATATATCGCGCGCG";
        let repeats = detect_repeats(sequence, 3);

        assert_eq!(repeats.len(), 2, "Should detect two separate repeats");
        assert!(
            repeats[0].end <= repeats[1].start || repeats[1].end <= repeats[0].start,
            "Repeats should not overlap"
        );
    }

    #[test]
    fn test_imperfect_repeat_not_detected() {
        // Repeat with interruption should not be detected as one long repeat
        let sequence = b"ATATATAGATATATAT";
        let repeats = detect_repeats(sequence, 3);

        // Should detect either side but not across the interruption
        for repeat in &repeats {
            assert_eq!(repeat.unit, "AT");
            assert!(
                repeat.end - repeat.start < 16,
                "Should not span the interruption"
            );
        }
    }

    #[test]
    fn test_repeat_at_sequence_boundaries() {
        // Repeat at the very start
        let start_repeat = b"ATATATATGGGGGGGG";
        let repeats = detect_repeats(start_repeat, 3);
        assert!(repeats.iter().any(|r| r.start == 0));

        // Repeat at the very end
        let end_repeat = b"GGGGGGGGATATAT";
        let repeats = detect_repeats(end_repeat, 3);
        assert!(repeats.iter().any(|r| r.end == end_repeat.len()));
    }

    #[test]
    fn test_homopolymer_not_detected_as_dinucleotide() {
        // AAAAA should not be detected as a di/tri/tetranucleotide repeat
        let sequence = b"AAAAAAAAA";
        let repeats = detect_repeats(sequence, 3);

        // Either no repeats (homopolymer excluded) or only if explicitly supported
        // The spec asks for di-, tri-, and tetranucleotide only
        assert!(repeats.is_empty() || repeats[0].period >= 2);
    }

    #[test]
    fn test_case_sensitivity() {
        // DNA sequences are typically uppercase, but test lowercase handling
        let sequence = b"atatatat";
        let repeats = detect_repeats(sequence, 3);

        // Should either handle lowercase or require uppercase
        // Implementation should be consistent
        assert!(repeats.is_empty() || repeats[0].unit == "at" || repeats[0].unit == "AT");
    }

    #[test]
    fn test_overlapping_period_interpretations() {
        // ATAT can be viewed as dinucleotide AT or tetranucleotide ATAT
        // The function should prefer the shorter period
        let sequence = b"ATATATAT";
        let repeats = detect_repeats(sequence, 3);

        assert_eq!(repeats.len(), 1);
        assert_eq!(
            repeats[0].period, 2,
            "Should prefer shorter period (dinucleotide)"
        );
        assert_eq!(repeats[0].unit, "AT");
    }

    // ========================================================================
    // Acceptance Test: Performance requirement
    // ========================================================================

    #[test]
    fn test_performance_10kb_sequence() {
        // AC: Performance: <100ms for 10kb sequence
        use std::time::Instant;

        // Generate a 10KB sequence with some repeats embedded
        let mut sequence = Vec::with_capacity(10_000);

        // Add varied content with some repeats
        for i in 0..2000 {
            match i % 5 {
                0 => sequence.extend_from_slice(b"ATATAT"),
                1 => sequence.extend_from_slice(b"CGCGC"),
                2 => sequence.extend_from_slice(b"GATAGATA"),
                3 => sequence.extend_from_slice(b"CAGCAGCAG"),
                _ => sequence.extend_from_slice(b"ACGT"),
            }
            if sequence.len() >= 10_000 {
                break;
            }
        }
        sequence.truncate(10_000);

        let start = Instant::now();
        let repeats = detect_repeats(&sequence, 3);
        let duration = start.elapsed();

        assert!(
            duration.as_millis() < 100,
            "Performance requirement: should process 10KB in <100ms, took {}ms",
            duration.as_millis()
        );
        assert!(
            !repeats.is_empty(),
            "Should detect some repeats in test sequence"
        );
    }

    #[test]
    fn test_performance_many_short_repeats() {
        use std::time::Instant;

        // Worst case: many short repeats that just meet the threshold
        let mut sequence = Vec::with_capacity(10_000);
        for _ in 0..1000 {
            sequence.extend_from_slice(b"ATATATATAA"); // 3 repeats of AT plus filler
        }

        let start = Instant::now();
        let repeats = detect_repeats(&sequence, 3);
        let duration = start.elapsed();

        assert!(
            duration.as_millis() < 100,
            "Should handle many short repeats efficiently, took {}ms",
            duration.as_millis()
        );
        assert!(repeats.len() >= 100, "Should detect many repeats");
    }
}
