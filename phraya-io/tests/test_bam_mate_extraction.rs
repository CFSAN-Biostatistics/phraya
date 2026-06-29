use phraya_core::types::MateInfo;
use phraya_io::bam_cram::BamCramParser;
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

/// Verify that BamCramParser rejects files that are plainly not BAM.
/// This calls the real parser entry point; a non-BAM file must return an error,
/// not silently produce empty output.
#[test]
fn bam_parser_returns_err_on_non_bam_input() {
    let mut f = NamedTempFile::new().unwrap();
    write!(f, ">read1\nACGTACGT\n").unwrap();
    f.flush().unwrap();
    let result = BamCramParser::from_bam_path(f.path());
    assert!(result.is_err(), "FASTA content must not be accepted as BAM");
}

/// MateInfo::is_discordant uses the absolute value of insert_size, treats
/// zero as "unmapped / unpaired" (never discordant), and applies the sigma
/// multiplier to std_dev. All four behaviours tested here.
#[test]
fn mate_info_is_discordant_semantics() {
    // mean=400, std_dev=50, sigma=3.0 → threshold deviation = 150 → window [250, 550]
    let ok = MateInfo::new("r/2".to_string(), true, 450, true, false, true);
    assert!(!ok.is_discordant(400, 50, 3.0), "450 within 3σ must be concordant");

    let out = MateInfo::new("r/2".to_string(), true, 700, true, false, true);
    assert!(out.is_discordant(400, 50, 3.0), "700 beyond 3σ must be discordant");

    // Negative insert_size: mate is upstream; absolute value must be used.
    let neg_ok = MateInfo::new("r/2".to_string(), true, -450, true, false, true);
    assert!(!neg_ok.is_discordant(400, 50, 3.0), "|-450| within 3σ must be concordant");

    let neg_out = MateInfo::new("r/2".to_string(), true, -700, true, false, true);
    assert!(neg_out.is_discordant(400, 50, 3.0), "|-700| beyond 3σ must be discordant");

    // Zero insert_size means unmapped mate: must never be considered discordant.
    let zero = MateInfo::new("r/2".to_string(), false, 0, true, false, false);
    assert!(!zero.is_discordant(400, 50, 3.0), "insert_size=0 must not be discordant");

    // Wider sigma accepts the same value that narrow sigma rejects.
    let borderline = MateInfo::new("r/2".to_string(), true, 600, true, false, true);
    assert!(borderline.is_discordant(400, 50, 3.0), "600 rejected at 3σ (window [250,550])");
    assert!(!borderline.is_discordant(400, 50, 5.0), "600 accepted at 5σ (window [150,650])");
}

/// All six fields passed to MateInfo::new must round-trip correctly.
/// This guards against silent field-ordering bugs in the constructor.
#[test]
fn mate_info_fields_roundtrip_through_constructor() {
    let first = MateInfo::new(
        "read123/2".to_string(),
        true,   // proper_pair
        -450,   // insert_size (negative = mate upstream)
        true,   // is_first_in_pair
        false,  // is_second_in_pair
        true,   // mate_mapped
    );
    assert_eq!(first.mate_id, "read123/2");
    assert!(first.proper_pair);
    assert_eq!(first.insert_size, -450);
    assert!(first.is_first_in_pair);
    assert!(!first.is_second_in_pair);
    assert!(first.mate_mapped);

    // Second-in-pair with unmapped mate.
    let second = MateInfo::new("read123/1".to_string(), false, 0, false, true, false);
    assert_eq!(second.mate_id, "read123/1");
    assert!(!second.proper_pair);
    assert_eq!(second.insert_size, 0);
    assert!(!second.is_first_in_pair);
    assert!(second.is_second_in_pair);
    assert!(!second.mate_mapped);
}
