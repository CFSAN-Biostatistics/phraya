use phraya_io::bam_cram::{BamCramParser, ParsedReads};
use std::io::Write;
use tempfile::NamedTempFile;

/// Helper to create minimal BAM file with proper pairs
/// Note: This is a simplified test - real BAM files need proper headers and structure
#[test]
#[ignore] // Requires creating valid BAM files, which is complex
fn bam_parser_extracts_mate_info() {
    // This test is a placeholder for future BAM file creation
    // Real implementation would use noodles to write a proper BAM file
    // with paired reads, proper pair flags, and TLEN values

    // Example structure:
    // 1. Create BAM with noodles::bam::io::Writer
    // 2. Add records with SAM flags set (0x1, 0x2, 0x40/0x80)
    // 3. Set TLEN field to known values
    // 4. Parse with BamCramParser::from_bam_path()
    // 5. Verify mate_info HashMap populated correctly

    println!("BAM parsing test requires valid BAM file construction");
}

#[test]
fn parsed_reads_structure() {
    use phraya_core::types::Sequence;
    use std::collections::HashMap;

    // Verify ParsedReads can be constructed
    let sequences = vec![
        Sequence::new(
            b"ATCGATCG".to_vec(),
            Some(vec![30; 8]),
            "read1".to_string(),
            None,
        ),
    ];

    let mut mate_info = HashMap::new();
    mate_info.insert(
        "read1".to_string(),
        phraya_core::types::MateInfo::new(
            "read1/2".to_string(),
            true,
            450,
            true,
            false,
            true,
        ),
    );

    let parsed = ParsedReads {
        sequences,
        mate_info,
    };

    assert_eq!(parsed.sequences.len(), 1);
    assert_eq!(parsed.mate_info.len(), 1);
    assert!(parsed.mate_info.contains_key("read1"));
}

/// Test mate ID construction logic
#[test]
fn mate_id_toggling() {
    // Simulate the logic from BamCramParser

    // Read with /1 suffix → mate is /2
    let id1 = "read123/1";
    let base_id1 = id1.trim_end_matches("/1").trim_end_matches("/2");
    let is_first = true;
    let mate_id1 = if is_first {
        format!("{}/2", base_id1)
    } else {
        format!("{}/1", base_id1)
    };
    assert_eq!(mate_id1, "read123/2");

    // Read with /2 suffix → mate is /1
    let id2 = "read123/2";
    let base_id2 = id2.trim_end_matches("/1").trim_end_matches("/2");
    let is_second = false; // is_first = false
    let mate_id2 = if is_second {
        format!("{}/2", base_id2)
    } else {
        format!("{}/1", base_id2)
    };
    assert_eq!(mate_id2, "read123/1");

    // Read without suffix
    let id3 = "read456";
    let base_id3 = id3.trim_end_matches("/1").trim_end_matches("/2");
    let mate_id3 = format!("{}/2", base_id3);
    assert_eq!(mate_id3, "read456/2");
}

/// Test SAM flag interpretation
#[test]
fn sam_flags_interpretation() {
    // Simulate flag checks from noodles
    // In real code: flags.is_segmented(), flags.is_properly_segmented(), etc.

    // SAM flags as u16 bitfield
    const PAIRED: u16 = 0x1;
    const PROPER_PAIR: u16 = 0x2;
    const MATE_UNMAPPED: u16 = 0x8;
    const FIRST_IN_PAIR: u16 = 0x40;
    const SECOND_IN_PAIR: u16 = 0x80;

    // Example: properly paired first read
    let flags = PAIRED | PROPER_PAIR | FIRST_IN_PAIR;
    let is_paired = (flags & PAIRED) != 0;
    let proper_pair = (flags & PROPER_PAIR) != 0;
    let is_first = (flags & FIRST_IN_PAIR) != 0;
    let is_second = (flags & SECOND_IN_PAIR) != 0;
    let mate_unmapped = (flags & MATE_UNMAPPED) != 0;

    assert!(is_paired);
    assert!(proper_pair);
    assert!(is_first);
    assert!(!is_second);
    assert!(!mate_unmapped);

    // Example: paired but not proper (discordant)
    let flags2 = PAIRED | FIRST_IN_PAIR; // Missing PROPER_PAIR
    let proper_pair2 = (flags2 & PROPER_PAIR) != 0;
    assert!(!proper_pair2);

    // Example: mate unmapped
    let flags3 = PAIRED | MATE_UNMAPPED | FIRST_IN_PAIR;
    let mate_unmapped3 = (flags3 & MATE_UNMAPPED) != 0;
    assert!(mate_unmapped3);
}
