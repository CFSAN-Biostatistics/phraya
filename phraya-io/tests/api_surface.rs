//! API surface tests for phraya-io
//!
//! These tests document and verify the expected public API for FASTA/FASTQ I/O.
//! All tests should FAIL initially (RED phase) until implementation is complete.

use phraya_core::Sequence;
use phraya_io::{ParseError, SequenceFormat, SequenceReader, parse_sequences};
use std::path::Path;

/// Test that the main parsing function exists and has the correct signature
#[test]
fn test_parse_sequences_function_exists() {
    // parse_sequences should accept a Path-like argument
    let path = Path::new("dummy.fasta");

    // Should return Result<SequenceReader, ParseError>
    let _result: Result<SequenceReader, ParseError> = parse_sequences(path);
}

/// Test that SequenceReader is an iterator
#[test]
fn test_sequence_reader_is_iterator() {
    let path = Path::new("dummy.fasta");
    let reader = parse_sequences(path).unwrap();

    // Should implement Iterator<Item = Result<Sequence, ParseError>>
    for result in reader {
        match result {
            Ok(seq) => {
                // Should have access to Sequence fields
                let _id = seq.id;
                let _description = seq.description;
                let _data = seq.data;
                let _quality = seq.quality;
            }
            Err(e) => {
                // Should have structured error types
                match e {
                    ParseError::IoError { path: _, source: _ } => {}
                    ParseError::MalformedEntry { line: _, reason: _ } => {}
                    ParseError::QualityLengthMismatch {
                        seq_len: _,
                        qual_len: _,
                    } => {}
                    ParseError::UnsupportedFormat { path: _, reason: _ } => {}
                }
            }
        }
    }
}

/// Test that SequenceFormat::detect exists and returns format information
#[test]
fn test_sequence_format_detect() {
    let path = Path::new("test.fasta");

    // Should return Result<SequenceFormat, ParseError>
    let _result: Result<SequenceFormat, ParseError> = SequenceFormat::detect(path);
}

/// Test that SequenceFormat provides format query methods
#[test]
fn test_sequence_format_query_methods() {
    let path = Path::new("test.fasta");
    let format = SequenceFormat::detect(path).unwrap();

    // Should have query methods
    let _is_fasta: bool = format.is_fasta();
    let _is_fastq: bool = format.is_fastq();
    let _is_compressed: bool = format.is_compressed();
}

/// Test Sequence type structure from phraya-core
#[test]
fn test_sequence_type_structure() {
    // Sequence should have these fields
    let seq = Sequence {
        id: Some("seq1".to_string()),
        description: Some("test sequence".to_string()),
        data: vec![b'A', b'C', b'T', b'G'],
        quality: Some(vec![b'I', b'I', b'I', b'I']),
        pairing_info: None,
    };

    assert_eq!(seq.id, Some("seq1".to_string()));
    assert_eq!(seq.data.len(), 4);
    assert!(seq.quality.is_some());
}

/// Test ParseError variants are comprehensive
#[test]
fn test_parse_error_variants() {
    // IoError for file access problems
    let _io_err = ParseError::IoError {
        path: "test.fasta".into(),
        source: std::io::Error::new(std::io::ErrorKind::NotFound, "file not found"),
    };

    // MalformedEntry for parse errors
    let _malformed_err = ParseError::MalformedEntry {
        line: 42,
        reason: "Invalid header".to_string(),
    };

    // QualityLengthMismatch for FASTQ-specific errors
    let _quality_err = ParseError::QualityLengthMismatch {
        seq_len: 100,
        qual_len: 98,
    };

    // UnsupportedFormat for format detection failures
    let _format_err = ParseError::UnsupportedFormat {
        path: "test.txt".into(),
        reason: "Not a FASTA or FASTQ file".to_string(),
    };
}

/// Test that ParseError implements std::error::Error
#[test]
fn test_parse_error_is_std_error() {
    fn assert_is_error<T: std::error::Error>() {}
    assert_is_error::<ParseError>();
}

/// Test that ParseError has Display implementation
#[test]
fn test_parse_error_display() {
    let err = ParseError::MalformedEntry {
        line: 10,
        reason: "Missing header".to_string(),
    };

    let display_string = format!("{}", err);
    assert!(display_string.contains("line"));
    assert!(display_string.contains("10"));
}

/// Test that Sequence supports Debug
#[test]
fn test_sequence_debug() {
    let seq = Sequence {
        id: Some("test".to_string()),
        description: None,
        data: vec![b'A', b'C', b'T', b'G'],
        quality: None,
        pairing_info: None,
    };

    let debug_string = format!("{:?}", seq);
    assert!(debug_string.contains("test"));
}

/// Test that Sequence supports Clone
#[test]
fn test_sequence_clone() {
    let seq = Sequence {
        id: Some("test".to_string()),
        description: None,
        data: vec![b'A', b'C', b'T', b'G'],
        quality: None,
        pairing_info: None,
    };

    let seq2 = seq.clone();
    assert_eq!(seq.id, seq2.id);
    assert_eq!(seq.data, seq2.data);
}
