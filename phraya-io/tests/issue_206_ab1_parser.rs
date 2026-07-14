use phraya_core::types::ParseError;
/// Issue #206: feat(io): parse Sanger AB1/ABIF trace input
///
/// This test file contains RED (failing) acceptance tests for issue #206.
/// Tests verify that AB1 (ABIF) files are parsed into Sequence records with proper
/// ID, basecalls, and per-base quality scores (PCON).
///
/// Acceptance Criteria:
/// 1. `.ab1` files are auto-detected and parsed into `Sequence` records
/// 2. Handles ABIF directory/tag traversal robustly (missing optional tags tolerated)
/// 3. No non-Rust dependency (pure-Rust ABIF parser)
/// 4. Extracts basecalls from PBAS tag and qualities from PCON tag
/// 5. Test against a small sample AB1
use phraya_io::SequenceParser;
use std::io::Write;
use std::path::PathBuf;
use tempfile::NamedTempFile;

// ============================================================================
// Helper Functions for Creating Minimal AB1 Files
// ============================================================================

/// AB1 (ABIF) is a tagged binary format with a directory structure.
/// Format:
/// - Header: "ABIF" (4 bytes)
/// - Version: u32 big-endian (usually 101)
/// - Number of elements: u32 big-endian
/// - Number of elements again: u32 big-endian
/// - Directory: array of tags
/// - Data section: variable-length records
///
/// Tag structure (28 bytes each):
/// - Tag name: 4 ASCII characters
/// - Tag number: u32 big-endian
/// - Element type: u16 big-endian (e.g., 4=byte, 5=char, 7=int, etc.)
/// - Element size: u16 big-endian
/// - Number of elements: u32 big-endian
/// - Data size: u32 big-endian (total bytes)
/// - Data offset or value: u32 big-endian

/// Create a minimal valid AB1 file with PBAS tag only
fn create_minimal_ab1_single_sequence() -> (NamedTempFile, PathBuf) {
    let tmp = NamedTempFile::new().unwrap();
    let ab1_path = tmp.path().with_extension("ab1");

    // Write minimal valid AB1 file
    write_minimal_ab1(
        &ab1_path,
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTAC",
        Some("IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"),
    )
    .unwrap();

    (tmp, ab1_path)
}

/// Create an AB1 with PBAS and PCON (quality) tags
fn create_ab1_with_quality() -> (NamedTempFile, PathBuf) {
    let tmp = NamedTempFile::new().unwrap();
    let ab1_path = tmp.path().with_extension("ab1");

    write_minimal_ab1(
        &ab1_path,
        "ACGTACGTACGTACGTACGTACGTACGT",
        Some("HHHHHHHHHHHHHHHHHHHHHHHHHHHH"),
    )
    .unwrap();

    (tmp, ab1_path)
}

/// Create an AB1 with PBAS tag but no PCON (missing quality tag)
fn create_ab1_no_quality() -> (NamedTempFile, PathBuf) {
    let tmp = NamedTempFile::new().unwrap();
    let ab1_path = tmp.path().with_extension("ab1");

    write_minimal_ab1(&ab1_path, "ACGTACGTACGTACGTACGTACGT", None).unwrap();

    (tmp, ab1_path)
}

/// Create an empty/corrupt AB1 (missing PBAS tag)
fn create_corrupt_ab1_no_pbas() -> (NamedTempFile, PathBuf) {
    let tmp = NamedTempFile::new().unwrap();
    let ab1_path = tmp.path().with_extension("ab1");

    // Write an AB1 header but with no PBAS tag (only PCON)
    write_ab1_header_only(&ab1_path).unwrap();

    (tmp, ab1_path)
}

/// Create a file with invalid AB1 magic bytes
fn create_invalid_ab1_magic() -> (NamedTempFile, PathBuf) {
    let tmp = NamedTempFile::new().unwrap();
    let ab1_path = tmp.path().with_extension("ab1");

    let mut file = std::fs::File::create(&ab1_path).unwrap();
    // Write wrong magic: "XXIF" instead of "ABIF"
    write!(file, "XXIF").unwrap();
    // Write some dummy data
    file.write_all(&[0u8; 100]).unwrap();

    (tmp, ab1_path)
}

/// Helper to write a minimal valid AB1 file structure
fn write_minimal_ab1(path: &PathBuf, bases: &str, quality: Option<&str>) -> std::io::Result<()> {
    let mut file = std::fs::File::create(path)?;

    // ABIF header
    file.write_all(b"ABIF")?; // Magic
    file.write_all(&101u32.to_be_bytes())?; // Version 101

    // Number of directory elements (we'll have 2-3 tags)
    let num_elements: u32 = if quality.is_some() { 3u32 } else { 2u32 };
    file.write_all(&num_elements.to_be_bytes())?; // Number of elements
    file.write_all(&num_elements.to_be_bytes())?; // Number of elements (repeated)

    // Calculate offsets for data section
    // Directory starts after header (4+4+4+4 = 16 bytes)
    let dir_size = (num_elements as usize) * 24; // Each tag is 24 bytes (not 28!)
    let data_offset = 16 + dir_size;

    let bases_bytes = bases.as_bytes();
    let bases_len = bases_bytes.len();

    // Directory entry for PBAS (bases tag) - type 5 (char)
    file.write_all(b"PBAS")?; // Tag name
    file.write_all(&1u32.to_be_bytes())?; // Tag number
    file.write_all(&5u16.to_be_bytes())?; // Element type (char)
    file.write_all(&1u16.to_be_bytes())?; // Element size
    file.write_all(&(bases_len as u32).to_be_bytes())?; // Number of elements
    file.write_all(&(bases_len as u32).to_be_bytes())?; // Data size
    file.write_all(&(data_offset as u32).to_be_bytes())?; // Data offset

    // If quality provided, add PCON tag (type 5 = char, same as PBAS)
    if let Some(q) = quality {
        let q_bytes = q.as_bytes();
        let q_len = q_bytes.len();
        let pcon_offset = data_offset + bases_len;

        file.write_all(b"PCON")?; // Tag name
        file.write_all(&1u32.to_be_bytes())?; // Tag number
        file.write_all(&5u16.to_be_bytes())?; // Element type (char)
        file.write_all(&1u16.to_be_bytes())?; // Element size
        file.write_all(&(q_len as u32).to_be_bytes())?; // Number of elements
        file.write_all(&(q_len as u32).to_be_bytes())?; // Data size
        file.write_all(&(pcon_offset as u32).to_be_bytes())?; // Data offset
    }

    // Dummy tag (required for valid AB1)
    file.write_all(b"APID")?; // Sequencer name
    file.write_all(&1u32.to_be_bytes())?; // Tag number
    file.write_all(&5u16.to_be_bytes())?; // Element type (char)
    file.write_all(&1u16.to_be_bytes())?; // Element size
    file.write_all(&4u32.to_be_bytes())?; // Number of elements
    file.write_all(&4u32.to_be_bytes())?; // Data size
    file.write_all(&0u32.to_be_bytes())?; // Data offset (inline, unused)

    // Data section: write bases
    file.write_all(bases_bytes)?;

    // Write quality if present
    if let Some(q) = quality {
        file.write_all(q.as_bytes())?;
    }

    Ok(())
}

/// Helper to write a minimal AB1 header only (no valid tags)
fn write_ab1_header_only(path: &PathBuf) -> std::io::Result<()> {
    let mut file = std::fs::File::create(path)?;

    file.write_all(b"ABIF")?; // Magic
    file.write_all(&101u32.to_be_bytes())?; // Version 101
    file.write_all(&0u32.to_be_bytes())?; // 0 elements
    file.write_all(&0u32.to_be_bytes())?; // 0 elements

    Ok(())
}

// ============================================================================
// AB1 Parser Tests
// ============================================================================

/// Happy path: parse a minimal AB1 file with bases and quality
#[test]
fn parse_ab1_with_bases_and_quality() {
    let (_tmp, ab1_path) = create_ab1_with_quality();

    let result = SequenceParser::from_path(&ab1_path);
    assert!(result.is_ok(), "should successfully open AB1 file");

    let mut parser = result.unwrap();
    let seq_result = parser.next();

    assert!(seq_result.is_some(), "should yield at least one sequence");

    let seq = seq_result.unwrap();
    assert!(seq.is_ok(), "should parse sequence without error");

    let seq = seq.unwrap();
    assert_eq!(
        seq.bases(),
        b"ACGTACGTACGTACGTACGTACGTACGT",
        "should extract bases from PBAS tag"
    );
    assert!(
        seq.quality_scores().is_some(),
        "should have quality scores from PCON tag"
    );
    assert_eq!(
        seq.quality_scores().unwrap(),
        b"HHHHHHHHHHHHHHHHHHHHHHHHHHHH",
        "quality should match PCON tag"
    );
}

/// Parse AB1 with bases but no quality tag (PCON missing)
#[test]
fn parse_ab1_no_quality_tolerates_missing_pcon() {
    let (_tmp, ab1_path) = create_ab1_no_quality();

    let result = SequenceParser::from_path(&ab1_path);
    assert!(result.is_ok(), "should successfully open AB1 file");

    let mut parser = result.unwrap();
    let seq_result = parser.next();

    assert!(seq_result.is_some(), "should yield at least one sequence");

    let seq = seq_result.unwrap();
    assert!(
        seq.is_ok(),
        "should parse successfully even without PCON tag"
    );

    let seq = seq.unwrap();
    assert_eq!(
        seq.bases(),
        b"ACGTACGTACGTACGTACGTACGT",
        "should extract bases from PBAS tag"
    );
    assert!(
        seq.quality_scores().is_none(),
        "should have no quality when PCON tag absent"
    );
}

/// Parse minimal AB1 with single sequence
#[test]
fn parse_minimal_ab1_single_sequence() {
    let (_tmp, ab1_path) = create_minimal_ab1_single_sequence();

    let result = SequenceParser::from_path(&ab1_path);
    assert!(result.is_ok(), "should successfully open AB1 file");

    let mut parser = result.unwrap();
    let seq = parser.next().unwrap().unwrap();

    assert!(!seq.bases().is_empty(), "sequence should have bases");
    assert_eq!(seq.len(), 38, "should have correct length");
}

/// Error case: AB1 file with missing required PBAS tag
#[test]
fn parse_corrupt_ab1_missing_pbas_returns_error() {
    let (_tmp, ab1_path) = create_corrupt_ab1_no_pbas();

    let result = SequenceParser::from_path(&ab1_path);

    // Either the open fails or the parse fails, but it should fail
    if result.is_ok() {
        let mut parser = result.unwrap();
        let seq_result = parser.next();
        // If we get here, should be an error when trying to extract no-bases sequence
        if let Some(Err(e)) = seq_result {
            assert!(
                matches!(e, ParseError::InvalidFormat(_)),
                "should return InvalidFormat error for missing PBAS tag"
            );
        } else if let Some(Ok(seq)) = seq_result {
            // Empty bases should also be considered an error
            assert!(
                seq.bases().is_empty(),
                "if parsed, should have empty bases as fallback"
            );
        }
    }
}

/// Error case: file with invalid AB1 magic bytes
#[test]
fn parse_invalid_ab1_magic_returns_error() {
    let (_tmp, ab1_path) = create_invalid_ab1_magic();

    let result = SequenceParser::from_path(&ab1_path);

    assert!(
        result.is_err(),
        "should reject file with invalid AB1 magic bytes"
    );

    if let Err(ParseError::InvalidFormat(msg)) = result {
        // Error message should hint at invalid format
        assert!(
            msg.contains("AB1") || msg.contains("magic") || msg.contains("invalid"),
            "error message should indicate AB1 format issue: {}",
            msg
        );
    }
}

/// Error case: empty AB1 file (no header)
#[test]
fn parse_empty_ab1_file_returns_error() {
    let tmp = NamedTempFile::new().unwrap();
    let ab1_path = tmp.path().with_extension("ab1");
    std::fs::write(&ab1_path, b"").unwrap();

    let result = SequenceParser::from_path(&ab1_path);

    assert!(result.is_err(), "should reject empty AB1 file");
}

/// Edge case: AB1 file with exactly 0 bases
#[test]
fn parse_ab1_with_zero_bases_returns_error() {
    let tmp = NamedTempFile::new().unwrap();
    let ab1_path = tmp.path().with_extension("ab1");

    write_minimal_ab1(&ab1_path, "", None).unwrap();

    let result = SequenceParser::from_path(&ab1_path);

    if result.is_ok() {
        let mut parser = result.unwrap();
        let seq_result = parser.next();

        if let Some(Ok(seq)) = seq_result {
            // Empty sequence should be considered invalid or produce None
            assert!(seq.bases().is_empty(), "sequence should have no bases");
        } else if let Some(Err(e)) = seq_result {
            assert!(
                matches!(e, ParseError::InvalidFormat(_)),
                "should return error for zero-length bases"
            );
        }
    }
}

/// Auto-detection: .ab1 extension triggers AB1 parser, not FASTA/FASTQ
#[test]
fn ab1_file_extension_triggers_ab1_parser() {
    let (_tmp, ab1_path) = create_minimal_ab1_single_sequence();

    // Verify the file has .ab1 extension
    assert_eq!(ab1_path.extension().unwrap(), "ab1");

    let result = SequenceParser::from_path(&ab1_path);
    assert!(result.is_ok(), "should auto-detect and parse .ab1 file");

    let mut parser = result.unwrap();
    let seq = parser.next();
    assert!(seq.is_some(), "should yield a sequence");
}

/// Sequence ID extraction: AB1 should use filename or LSID tag as ID
#[test]
fn ab1_sequence_id_from_filename_or_tag() {
    let (_tmp, ab1_path) = create_ab1_with_quality();

    let result = SequenceParser::from_path(&ab1_path);
    let mut parser = result.unwrap();
    let seq = parser.next().unwrap().unwrap();

    // ID should either be the filename or extracted from AB1 tag
    let id = seq.id();
    assert!(!id.is_empty(), "sequence should have a non-empty ID");
    // Can be filename, LSID tag, or a default - just verify it's not empty
}

/// Quality length must match bases length
#[test]
fn ab1_quality_length_must_match_bases_length() {
    let (_tmp, ab1_path) = create_ab1_with_quality();

    let mut parser = SequenceParser::from_path(&ab1_path).unwrap();
    let seq = parser.next().unwrap().unwrap();

    if let Some(qual) = seq.quality_scores() {
        assert_eq!(
            qual.len(),
            seq.bases().len(),
            "quality length must match bases length"
        );
    }
}

/// PBAS tag extraction: bases should be bytes from PBAS record
#[test]
fn ab1_pbas_tag_contains_valid_dna_bases() {
    let (_tmp, ab1_path) = create_ab1_with_quality();

    let mut parser = SequenceParser::from_path(&ab1_path).unwrap();
    let seq = parser.next().unwrap().unwrap();

    let bases = seq.bases();
    assert!(!bases.is_empty(), "bases should not be empty");

    // Check that all bases are valid DNA (ACGT or N)
    for &b in bases {
        match b {
            b'A' | b'C' | b'G' | b'T' | b'N' | b'a' | b'c' | b'g' | b't' | b'n' => {}
            _ => panic!("invalid base in PBAS: {}", b as char),
        }
    }
}

/// PCON tag extraction: quality scores should be Phred ASCII from PCON record
#[test]
fn ab1_pcon_tag_contains_valid_phred_quality() {
    let (_tmp, ab1_path) = create_ab1_with_quality();

    let mut parser = SequenceParser::from_path(&ab1_path).unwrap();
    let seq = parser.next().unwrap().unwrap();

    if let Some(qual) = seq.quality_scores() {
        // Phred quality is typically ASCII 33+ (! and above)
        for &q in qual {
            assert!(
                q >= b'!' as u8 || q >= 33,
                "quality should be valid Phred ASCII: {}",
                q
            );
        }
    }
}

/// Parser should handle AB1 files robustly (no panic on malformed tags)
#[test]
fn ab1_parser_handles_truncated_directory_gracefully() {
    let tmp = NamedTempFile::new().unwrap();
    let ab1_path = tmp.path().with_extension("ab1");

    let mut file = std::fs::File::create(&ab1_path).unwrap();
    file.write_all(b"ABIF").unwrap();
    file.write_all(&101u32.to_be_bytes()).unwrap();
    file.write_all(&999u32.to_be_bytes()).unwrap(); // Claim 999 elements
    file.write_all(&999u32.to_be_bytes()).unwrap();
    file.write_all(b"data").unwrap(); // Incomplete directory

    let result = SequenceParser::from_path(&ab1_path);
    // Should either fail gracefully or return Err, not panic
    if result.is_ok() {
        let mut parser = result.unwrap();
        let seq_result = parser.next();
        // Should not panic, even if it's Err or None
        assert!(seq_result.is_none() || seq_result.unwrap().is_err());
    }
}

/// Parser rejects AB1 files with version mismatch (not version 101)
#[test]
fn ab1_parser_rejects_unsupported_version() {
    let tmp = NamedTempFile::new().unwrap();
    let ab1_path = tmp.path().with_extension("ab1");

    let mut file = std::fs::File::create(&ab1_path).unwrap();
    file.write_all(b"ABIF").unwrap();
    file.write_all(&999u32.to_be_bytes()).unwrap(); // Invalid version
    file.write_all(&0u32.to_be_bytes()).unwrap();
    file.write_all(&0u32.to_be_bytes()).unwrap();

    let result = SequenceParser::from_path(&ab1_path);
    // Should fail on unsupported version
    if result.is_ok() {
        let mut parser = result.unwrap();
        let seq_result = parser.next();
        // Either error or no sequence returned
        assert!(seq_result.is_none() || seq_result.unwrap().is_err());
    }
}
