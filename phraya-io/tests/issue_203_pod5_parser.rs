/// Issue #203: feat(io): parse Oxford Nanopore POD5 input (basecalled reads)
///
/// This test file contains RED (failing) acceptance tests for issue #203.
/// Tests verify that POD5 files are parsed into Sequence records with proper
/// ID, basecalls, and optional quality scores.
///
/// Acceptance Criteria:
/// 1. `.pod5` files are auto-detected and parsed into `Sequence` records
/// 2. Multi-read POD5 files iterate all reads
/// 3. Implemented without non-Rust runtime dependency (arrow-rs)
/// 4. Round-trip test against a small real/sample POD5

use phraya_io::SequenceParser;
use phraya_core::types::ParseError;
use std::io::Write;
use std::path::{Path, PathBuf};
use tempfile::NamedTempFile;

// ============================================================================
// Helper Functions for Creating Minimal POD5 Files
// ============================================================================

/// Create a minimal valid POD5 file with a single read
/// POD5 is an Apache Arrow IPC file format with Parquet serialization
fn create_minimal_pod5_single_read() -> (NamedTempFile, PathBuf) {
    let tmp = NamedTempFile::new().unwrap();
    let pod5_path = tmp.path().with_extension("pod5");

    // Write minimal POD5 Arrow file with one read
    // This is a placeholder that will be replaced with real POD5 when arrow-rs is available
    // For now, we create the basic structure
    write_minimal_arrow_pod5(&pod5_path, vec![
        PodRead {
            read_id: "read_00001".to_string(),
            basecall: "ACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_string(),
            quality: Some("IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII".to_string()),
        }
    ]).unwrap();

    (tmp, pod5_path)
}

/// Create a minimal POD5 with multiple reads
fn create_minimal_pod5_multiple_reads() -> (NamedTempFile, PathBuf) {
    let tmp = NamedTempFile::new().unwrap();
    let pod5_path = tmp.path().with_extension("pod5");

    write_minimal_arrow_pod5(&pod5_path, vec![
        PodRead {
            read_id: "read_00001".to_string(),
            basecall: "ACGTACGTACGTACGTACGTACGTACGT".to_string(),
            quality: Some("IIIIIIIIIIIIIIIIIIIIIIII".to_string()),
        },
        PodRead {
            read_id: "read_00002".to_string(),
            basecall: "TGCATGCATGCATGCATGCATGCA".to_string(),
            quality: Some("HHHHHHHHHHHHHHHHHHHHHHHH".to_string()),
        },
        PodRead {
            read_id: "read_00003".to_string(),
            basecall: "AAATTTGGGGCCCCAAATTTGGGG".to_string(),
            quality: Some("JJJJJJJJJJJJJJJJJJJJJJJJ".to_string()),
        }
    ]).unwrap();

    (tmp, pod5_path)
}

/// Create a POD5 with reads but no quality scores
fn create_minimal_pod5_no_quality() -> (NamedTempFile, PathBuf) {
    let tmp = NamedTempFile::new().unwrap();
    let pod5_path = tmp.path().with_extension("pod5");

    write_minimal_arrow_pod5(&pod5_path, vec![
        PodRead {
            read_id: "read_00001".to_string(),
            basecall: "ACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_string(),
            quality: None,
        }
    ]).unwrap();

    (tmp, pod5_path)
}

/// Create an empty POD5 file with no reads
fn create_minimal_pod5_empty() -> (NamedTempFile, PathBuf) {
    let tmp = NamedTempFile::new().unwrap();
    let pod5_path = tmp.path().with_extension("pod5");

    write_minimal_arrow_pod5(&pod5_path, vec![]).unwrap();

    (tmp, pod5_path)
}

/// Helper struct representing a POD5 read record
struct PodRead {
    read_id: String,
    basecall: String,
    quality: Option<String>,
}

/// Write a minimal Arrow POD5 file
/// This is a simplified representation; actual implementation will use arrow-rs
fn write_minimal_arrow_pod5(path: &Path, reads: Vec<PodRead>) -> std::io::Result<()> {
    use std::fs::File;

    let mut file = File::create(path)?;

    // Write POD5 magic number (Arrow IPC format)
    // Arrow format starts with specific bytes that identify it as Arrow
    // This is a placeholder that will be implemented with arrow-rs
    file.write_all(b"POD5")?;  // POD5 signature

    // In the real implementation, this will use arrow-rs to construct:
    // - RecordBatch with read_id (string), basecall (string), quality (optional string)
    // - Write as Apache Arrow IPC format
    // For testing purposes, we're creating a minimal placeholder

    for read in &reads {
        writeln!(file, "read_id:{}", read.read_id)?;
        writeln!(file, "basecall:{}", read.basecall)?;
        if let Some(q) = &read.quality {
            writeln!(file, "quality:{}", q)?;
        }
    }

    Ok(())
}

// ============================================================================
// POD5 Parser Tests - Happy Path
// ============================================================================

#[test]
#[ignore = "RED: awaiting real #203 POD5 implementation"]
fn issue_203_parse_single_read_pod5_file() {
    // Test: Parse a POD5 file with a single read
    // Expected: One Sequence record with correct ID and bases
    let (_tmp, path) = create_minimal_pod5_single_read();

    let result = SequenceParser::from_path(&path);

    // Should successfully parse POD5 file
    assert!(result.is_ok(), "POD5 file should parse without error");

    let mut parser = result.unwrap();
    let seq = parser.next();

    // Should yield at least one sequence
    assert!(seq.is_some(), "POD5 should yield at least one read");

    let seq_result = seq.unwrap();
    assert!(seq_result.is_ok(), "First read should parse successfully");

    let seq = seq_result.unwrap();

    // Verify read ID
    assert_eq!(seq.id(), "read_00001", "Read ID should match POD5 record");

    // Verify basecall sequence
    assert_eq!(
        seq.bases(),
        b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT",
        "Basecall sequence should match POD5 record"
    );

    // Verify sequence length
    assert_eq!(seq.len(), 38, "Sequence length should match basecall length");

    // Verify quality scores present
    assert!(seq.quality_at(0).is_some(), "Quality scores should be present");
    assert_eq!(seq.quality_at(0), Some(b'I'), "First quality score should be 'I'");
}

#[test]
#[ignore = "RED: awaiting real #203 POD5 implementation"]
fn issue_203_parse_multiple_reads_from_pod5() {
    // Test: Parse a POD5 file with multiple reads
    // Expected: All reads extracted in order
    let (_tmp, path) = create_minimal_pod5_multiple_reads();

    let result = SequenceParser::from_path(&path);
    assert!(result.is_ok(), "POD5 with multiple reads should parse");

    let parser = result.unwrap();

    // Collect all sequences
    let sequences: Vec<_> = parser
        .map(|r| r.expect("Each read should parse"))
        .collect();

    // Should have exactly 3 reads
    assert_eq!(sequences.len(), 3, "POD5 should yield 3 reads");

    // Verify first read
    assert_eq!(sequences[0].id(), "read_00001");
    assert_eq!(sequences[0].bases(), b"ACGTACGTACGTACGTACGTACGT");

    // Verify second read
    assert_eq!(sequences[1].id(), "read_00002");
    assert_eq!(sequences[1].bases(), b"TGCATGCATGCATGCATGCATGCA");
    assert_eq!(sequences[1].quality_at(0), Some(b'H'));

    // Verify third read
    assert_eq!(sequences[2].id(), "read_00003");
    assert_eq!(sequences[2].bases(), b"AAATTTGGGGCCCCAAATTTGGGG");
    assert_eq!(sequences[2].len(), 25);
}

#[test]
#[ignore = "RED: awaiting real #203 POD5 implementation"]
fn issue_203_parse_pod5_with_quality_scores() {
    // Test: Quality scores are properly extracted from POD5
    // Expected: quality_at() returns correct Phred scores
    let (_tmp, path) = create_minimal_pod5_single_read();

    let mut parser = SequenceParser::from_path(&path).unwrap();
    let seq = parser.next().unwrap().unwrap();

    // Verify quality scores are present
    for i in 0..seq.len() {
        let qual = seq.quality_at(i);
        assert!(qual.is_some(), "Quality score at position {} should exist", i);
        assert_eq!(qual, Some(b'I'), "Quality score at position {} should be 'I'", i);
    }
}

#[test]
#[ignore = "RED: awaiting real #203 POD5 implementation"]
fn issue_203_pod5_file_auto_detection_by_extension() {
    // Test: .pod5 extension triggers POD5 parser
    // Expected: File is recognized and parsed as POD5
    let (_tmp, path) = create_minimal_pod5_single_read();

    // Verify file has .pod5 extension
    assert!(
        path.extension().map(|e| e == "pod5").unwrap_or(false),
        "Test file should have .pod5 extension"
    );

    // Should parse correctly via SequenceParser
    let result = SequenceParser::from_path(&path);
    assert!(
        result.is_ok(),
        ".pod5 file should be auto-detected and parsed"
    );
}

#[test]
#[ignore = "RED: awaiting real #203 POD5 implementation"]
fn issue_203_parse_pod5_without_quality_scores() {
    // Test: POD5 reads without quality scores are handled correctly
    // Expected: Sequence created with quality=None
    let (_tmp, path) = create_minimal_pod5_no_quality();

    let mut parser = SequenceParser::from_path(&path).unwrap();
    let seq = parser.next().unwrap().unwrap();

    // Verify read ID and sequence
    assert_eq!(seq.id(), "read_00001");
    assert_eq!(seq.bases(), b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT");

    // Quality should be absent (not None for individual positions)
    // Depending on implementation, either quality_at() returns None consistently,
    // or Sequence::new was called with quality=None
    let qual_status = seq.quality_at(0);
    assert!(
        qual_status.is_none(),
        "Quality should be absent when not provided in POD5"
    );
}

// ============================================================================
// POD5 Parser Tests - Edge Cases
// ============================================================================

#[test]
#[ignore = "RED: awaiting real #203 POD5 implementation"]
fn issue_203_empty_pod5_file_returns_empty_iterator() {
    // Test: Empty POD5 with no reads
    // Expected: Parser returns successfully but yields no sequences
    let (_tmp, path) = create_minimal_pod5_empty();

    let result = SequenceParser::from_path(&path);
    assert!(result.is_ok(), "Empty POD5 should parse without error");

    let mut parser = result.unwrap();
    assert_eq!(
        parser.next(),
        None,
        "Empty POD5 should yield no sequences"
    );
}

#[test]
#[ignore = "RED: awaiting real #203 POD5 implementation"]
fn issue_203_pod5_read_id_extraction() {
    // Test: Read IDs are correctly extracted from POD5
    // Expected: seq.id() matches the POD5 read_id field
    let (_tmp, path) = create_minimal_pod5_multiple_reads();

    let parser = SequenceParser::from_path(&path).unwrap();

    let ids: Vec<_> = parser
        .map(|r| r.map(|seq| seq.id().to_string()))
        .collect::<Result<Vec<_>, _>>()
        .unwrap();

    assert_eq!(ids, vec!["read_00001", "read_00002", "read_00003"]);
}

#[test]
#[ignore = "RED: awaiting real #203 POD5 implementation"]
fn issue_203_pod5_sequence_length_matches_basecall_length() {
    // Test: seq.len() matches the basecall string length
    // Expected: All bases count correctly
    let (_tmp, path) = create_minimal_pod5_multiple_reads();

    let mut parser = SequenceParser::from_path(&path).unwrap();

    // Expected lengths from the test data
    let expected_lengths = vec![25, 25, 25];

    for (i, expected_len) in expected_lengths.iter().enumerate() {
        let seq = parser.next().unwrap().unwrap();
        assert_eq!(
            seq.len(),
            *expected_len,
            "Sequence {} should have length {}",
            i,
            expected_len
        );
    }
}

#[test]
#[ignore = "RED: awaiting real #203 POD5 implementation"]
fn issue_203_pod5_quality_length_matches_basecall_length() {
    // Test: Quality score count matches basecall length
    // Expected: Each base has exactly one quality score
    let (_tmp, path) = create_minimal_pod5_single_read();

    let mut parser = SequenceParser::from_path(&path).unwrap();
    let seq = parser.next().unwrap().unwrap();

    // Count quality scores by checking quality_at() for each position
    let mut quality_count = 0;
    for i in 0..seq.len() {
        if seq.quality_at(i).is_some() {
            quality_count += 1;
        }
    }

    assert_eq!(
        quality_count, seq.len(),
        "Quality score count should match sequence length"
    );
}

// ============================================================================
// POD5 Parser Tests - Error Handling
// ============================================================================

#[test]
fn issue_203_invalid_pod5_file_returns_error() {
    // Test: Malformed POD5 file triggers ParseError
    // Expected: ParseError with descriptive message
    let mut tmp = NamedTempFile::new().unwrap();
    let pod5_path = tmp.path().with_extension("pod5");

    // Write invalid POD5 data
    std::io::Write::write_all(&mut tmp, b"INVALID_POD5_DATA").unwrap();
    tmp.flush().unwrap();

    let result = SequenceParser::from_path(&pod5_path);

    assert!(
        result.is_err(),
        "Invalid POD5 should return ParseError"
    );

    if let Err(ParseError::InvalidFormat(msg)) = result {
        // Error message should indicate format problem
        assert!(
            msg.len() > 0,
            "Error message should describe the format issue"
        );
    } else {
        panic!("Expected ParseError::InvalidFormat");
    }
}

#[test]
fn issue_203_nonexistent_pod5_file_returns_error() {
    // Test: Non-existent file triggers appropriate error
    // Expected: ParseError indicating file not found
    use std::path::Path;

    let nonexistent = Path::new("/nonexistent/path/file.pod5");
    let result = SequenceParser::from_path(nonexistent);

    assert!(result.is_err(), "Non-existent file should return error");

    if let Err(ParseError::InvalidFormat(msg)) = result {
        assert!(
            msg.contains("failed to open") || msg.contains("not found") || msg.len() > 0,
            "Error should indicate file open failure"
        );
    } else {
        panic!("Expected ParseError::InvalidFormat");
    }
}

#[test]
fn issue_203_pod5_zero_length_basecall_handled() {
    // Test: POD5 with empty basecall string
    // Expected: Either skipped or created as empty Sequence
    // (Implementation detail; test that it doesn't panic)
    let tmp = NamedTempFile::new().unwrap();
    let pod5_path = tmp.path().with_extension("pod5");

    write_minimal_arrow_pod5(&pod5_path, vec![
        PodRead {
            read_id: "empty_read".to_string(),
            basecall: "".to_string(),  // Empty basecall
            quality: None,
        }
    ]).unwrap();

    let result = SequenceParser::from_path(&pod5_path);

    // Should either:
    // 1. Parse successfully and yield empty sequence (which then gets filtered)
    // 2. Parse successfully and skip the empty read
    // 3. Parse but yield empty iterator (empty Sequences are filtered)
    // The test just verifies no panic occurs
    if let Ok(parser) = result {
        let sequences: Vec<_> = parser.collect::<Result<Vec<_>, _>>().unwrap_or_default();
        // Empty sequences might be filtered; this is OK
        assert!(sequences.is_empty() || sequences.iter().all(|s| s.len() == 0));
    }
}

// ============================================================================
// POD5 Format Detection Tests
// ============================================================================

#[test]
fn issue_203_pod5_format_detected_before_fasta_fastq() {
    // Test: When parsing a .pod5 file, POD5 parser is tried first
    // Expected: POD5 parser is invoked, not FASTA/FASTQ parser
    let (_tmp, path) = create_minimal_pod5_single_read();

    // The parser should attempt POD5 first due to extension
    let result = SequenceParser::from_path(&path);

    // If it errors, the error should be from POD5 parsing, not from
    // trying to interpret it as FASTA/FASTQ
    if let Err(ParseError::InvalidFormat(msg)) = result {
        // Should be a POD5-specific error or a generic parsing error
        // Not an error about invalid FASTA/FASTQ format characters
        assert!(
            !msg.contains("must start with '>' (FASTA) or '@' (FASTQ)"),
            "Should not fail with FASTA/FASTQ error message"
        );
    }
}

// ============================================================================
// POD5 Round-Trip and Integration Tests
// ============================================================================

#[test]
#[ignore = "RED: awaiting real #203 POD5 implementation"]
fn issue_203_pod5_round_trip_preserves_all_metadata() {
    // Test: Reading a POD5 and examining all read properties
    // Expected: All metadata preserved through parsing
    let (_tmp, path) = create_minimal_pod5_single_read();

    let mut parser = SequenceParser::from_path(&path).unwrap();
    let seq = parser.next().unwrap().unwrap();

    // Verify all key metadata
    assert_eq!(seq.id(), "read_00001");
    assert_eq!(seq.bases(), b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT");
    assert_eq!(seq.len(), 38);

    // Quality should be consistent for all bases
    for i in 0..seq.len() {
        assert_eq!(seq.quality_at(i), Some(b'I'));
    }
}

#[test]
#[ignore = "RED: awaiting real #203 POD5 implementation"]
fn issue_203_pod5_large_read_handling() {
    // Test: POD5 with longer reads (closer to real ONT reads)
    // Expected: Long sequences parsed correctly without truncation
    let tmp = NamedTempFile::new().unwrap();
    let pod5_path = tmp.path().with_extension("pod5");

    // Create a read ~5000bp (realistic for ONT)
    let long_basecall = "ACGT".repeat(1250);  // 5000bp
    let long_quality = "I".repeat(5000);

    write_minimal_arrow_pod5(&pod5_path, vec![
        PodRead {
            read_id: "long_read".to_string(),
            basecall: long_basecall.clone(),
            quality: Some(long_quality),
        }
    ]).unwrap();

    let mut parser = SequenceParser::from_path(&pod5_path).unwrap();
    let seq = parser.next().unwrap().unwrap();

    assert_eq!(seq.len(), 5000, "Long read should be fully parsed");
    assert_eq!(seq.id(), "long_read");
    assert_eq!(seq.bases(), long_basecall.as_bytes());
}

#[test]
#[ignore = "RED: awaiting real #203 POD5 implementation"]
fn issue_203_pod5_many_reads_iteration() {
    // Test: POD5 with many reads yields all of them
    // Expected: All reads iterated in order
    let tmp = NamedTempFile::new().unwrap();
    let pod5_path = tmp.path().with_extension("pod5");

    // Create a POD5 with 100 reads
    let reads: Vec<_> = (0..100)
        .map(|i| PodRead {
            read_id: format!("read_{:05}", i),
            basecall: format!("ACGT"),
            quality: Some("IIII".to_string()),
        })
        .collect();

    write_minimal_arrow_pod5(&pod5_path, reads).unwrap();

    let parser = SequenceParser::from_path(&pod5_path).unwrap();
    let count = parser.count();

    assert_eq!(count, 100, "Should iterate exactly 100 reads");
}

// ============================================================================
// POD5 Integration with Phraya Pipeline
// ============================================================================

#[test]
#[ignore = "RED: awaiting real #203 POD5 implementation"]
fn issue_203_pod5_sequence_compatible_with_alignment_pipeline() {
    // Test: POD5-parsed Sequences work with downstream code
    // Expected: Sequences have all required fields for alignment
    let (_tmp, path) = create_minimal_pod5_single_read();

    let mut parser = SequenceParser::from_path(&path).unwrap();
    let seq = parser.next().unwrap().unwrap();

    // Sequences must support:
    // 1. .bases() accessor
    assert!(!seq.bases().is_empty(), "Must have bases");

    // 2. .id() accessor
    assert!(!seq.id().is_empty(), "Must have ID");

    // 3. .len() method
    assert!(seq.len() > 0, "Must have non-zero length");

    // 4. Optional quality scores
    assert!(seq.quality_at(0).is_some(), "Should have quality if provided");

    // 5. Should be cloneable/movable for use in alignment
    let seq_clone = seq.clone();
    assert_eq!(seq.id(), seq_clone.id());
}

#[test]
#[ignore = "RED: awaiting real #203 POD5 implementation"]
fn issue_203_pod5_dna_characters_valid() {
    // Test: POD5 basecalls contain only valid DNA bases
    // Expected: Characters are ACGT (case may vary, but test expects uppercase)
    let (_tmp, path) = create_minimal_pod5_single_read();

    let mut parser = SequenceParser::from_path(&path).unwrap();
    let seq = parser.next().unwrap().unwrap();

    for &base in seq.bases() {
        match base {
            b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't' | b'N' | b'n' => {
                // Valid DNA base
            }
            _ => panic!("Invalid DNA base: {}", base as char),
        }
    }
}
