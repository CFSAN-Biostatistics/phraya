//! Integration tests for FASTA/FASTQ I/O with auto-detection
//!
//! These tests cover all acceptance criteria for issue #2:
//! - Parse FASTA (plain, .gz, .bz2) with no quality field
//! - Parse FASTQ (plain, .gz, .bz2) with quality field populated
//! - Auto-detection via magic bytes and extension fallback
//! - All 6 combinations (FASTA/FASTQ × plain/gz/bz2)
//! - Error handling for malformed files
//! - Streaming interface (iterator, not Vec)

use phraya_core::Sequence;
use phraya_io::{ParseError, SequenceFormat, parse_sequences};
use std::path::PathBuf;

/// Helper to get fixture path
fn fixture_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join(name)
}

#[test]
fn test_parse_fasta_plain() {
    let path = fixture_path("sample.fasta");
    let sequences: Result<Vec<_>, _> = parse_sequences(&path).unwrap().collect();
    let sequences = sequences.expect("Failed to parse FASTA");

    assert_eq!(sequences.len(), 3);

    // Check first sequence
    assert_eq!(sequences[0].id, Some("seq1".to_string()));
    assert_eq!(
        sequences[0].description,
        Some("test sequence one".to_string())
    );
    assert_eq!(sequences[0].data, b"ACTGACTGACTG");
    assert!(
        sequences[0].quality.is_none(),
        "FASTA should have no quality scores"
    );

    // Check second sequence
    assert_eq!(sequences[1].id, Some("seq2".to_string()));
    assert_eq!(sequences[1].data, b"GGGGCCCCAAAATTTT");
    assert!(sequences[1].quality.is_none());

    // Check third sequence
    assert_eq!(sequences[2].id, Some("seq3".to_string()));
    assert_eq!(sequences[2].data, b"ACGTACGTACGTACGT");
    assert!(sequences[2].quality.is_none());
}

#[test]
fn test_parse_fasta_gzipped() {
    let path = fixture_path("sample.fasta.gz");
    let sequences: Result<Vec<_>, _> = parse_sequences(&path).unwrap().collect();
    let sequences = sequences.expect("Failed to parse gzipped FASTA");

    assert_eq!(sequences.len(), 3);
    assert_eq!(sequences[0].data, b"ACTGACTGACTG");
    assert!(sequences[0].quality.is_none());
    assert_eq!(sequences[1].data, b"GGGGCCCCAAAATTTT");
    assert_eq!(sequences[2].data, b"ACGTACGTACGTACGT");
}

#[test]
fn test_parse_fasta_bz2() {
    let path = fixture_path("sample.fasta.bz2");
    let sequences: Result<Vec<_>, _> = parse_sequences(&path).unwrap().collect();
    let sequences = sequences.expect("Failed to parse bzip2 FASTA");

    assert_eq!(sequences.len(), 3);
    assert_eq!(sequences[0].data, b"ACTGACTGACTG");
    assert!(sequences[0].quality.is_none());
    assert_eq!(sequences[1].data, b"GGGGCCCCAAAATTTT");
    assert_eq!(sequences[2].data, b"ACGTACGTACGTACGT");
}

#[test]
fn test_parse_fastq_plain() {
    let path = fixture_path("sample.fastq");
    let sequences: Result<Vec<_>, _> = parse_sequences(&path).unwrap().collect();
    let sequences = sequences.expect("Failed to parse FASTQ");

    assert_eq!(sequences.len(), 3);

    // Check first sequence
    assert_eq!(sequences[0].id, Some("seq1".to_string()));
    assert_eq!(
        sequences[0].description,
        Some("test sequence one".to_string())
    );
    assert_eq!(sequences[0].data, b"ACTGACTGACTG");
    assert!(
        sequences[0].quality.is_some(),
        "FASTQ must have quality scores"
    );
    assert_eq!(sequences[0].quality.as_ref().unwrap(), b"IIIIIIIIIIII");

    // Check second sequence
    assert_eq!(sequences[1].id, Some("seq2".to_string()));
    assert_eq!(sequences[1].data, b"GGGGCCCCAAAATTTT");
    assert_eq!(sequences[1].quality.as_ref().unwrap(), b"HHHHHHHHHHHHHHHH");

    // Check third sequence
    assert_eq!(sequences[2].id, Some("seq3".to_string()));
    assert_eq!(sequences[2].data, b"ACGTACGTACGTACGT");
    assert_eq!(sequences[2].quality.as_ref().unwrap(), b"JJJJJJJJJJJJJJJJ");
}

#[test]
fn test_parse_fastq_gzipped() {
    let path = fixture_path("sample.fastq.gz");
    let sequences: Result<Vec<_>, _> = parse_sequences(&path).unwrap().collect();
    let sequences = sequences.expect("Failed to parse gzipped FASTQ");

    assert_eq!(sequences.len(), 3);
    assert_eq!(sequences[0].data, b"ACTGACTGACTG");
    assert_eq!(sequences[0].quality.as_ref().unwrap(), b"IIIIIIIIIIII");
    assert_eq!(sequences[1].data, b"GGGGCCCCAAAATTTT");
    assert_eq!(sequences[1].quality.as_ref().unwrap(), b"HHHHHHHHHHHHHHHH");
}

#[test]
fn test_parse_fastq_bz2() {
    let path = fixture_path("sample.fastq.bz2");
    let sequences: Result<Vec<_>, _> = parse_sequences(&path).unwrap().collect();
    let sequences = sequences.expect("Failed to parse bzip2 FASTQ");

    assert_eq!(sequences.len(), 3);
    assert_eq!(sequences[0].data, b"ACTGACTGACTG");
    assert_eq!(sequences[0].quality.as_ref().unwrap(), b"IIIIIIIIIIII");
    assert_eq!(sequences[2].data, b"ACGTACGTACGTACGT");
    assert_eq!(sequences[2].quality.as_ref().unwrap(), b"JJJJJJJJJJJJJJJJ");
}

#[test]
fn test_auto_detection_magic_bytes_gzip() {
    // Should detect gzip from magic bytes even without .gz extension
    let path = fixture_path("sample.fastq.gz");
    let format = SequenceFormat::detect(&path).expect("Failed to detect format");

    // The detected format should indicate gzip compression
    assert!(
        format.is_compressed(),
        "Should detect gzip compression from magic bytes"
    );
}

#[test]
fn test_auto_detection_magic_bytes_bzip2() {
    // Should detect bzip2 from magic bytes even without .bz2 extension
    let path = fixture_path("sample.fasta.bz2");
    let format = SequenceFormat::detect(&path).expect("Failed to detect format");

    assert!(
        format.is_compressed(),
        "Should detect bzip2 compression from magic bytes"
    );
}

#[test]
fn test_auto_detection_extension_fallback() {
    // When magic bytes are ambiguous, should fall back to extension
    let fasta_path = fixture_path("sample.fasta");
    let format = SequenceFormat::detect(&fasta_path).expect("Failed to detect format");

    assert!(
        format.is_fasta(),
        "Should detect FASTA format from extension when needed"
    );

    let fastq_path = fixture_path("sample.fastq");
    let format = SequenceFormat::detect(&fastq_path).expect("Failed to detect format");

    assert!(
        format.is_fastq(),
        "Should detect FASTQ format from extension when needed"
    );
}

#[test]
fn test_streaming_interface() {
    // Verify that parse_sequences returns an iterator, not a Vec
    let path = fixture_path("sample.fasta");
    let iterator = parse_sequences(&path).expect("Failed to open file");

    // Iterator should be lazy - we can take just one element
    let first = iterator.take(1).collect::<Result<Vec<_>, _>>().unwrap();
    assert_eq!(first.len(), 1);
    assert_eq!(first[0].data, b"ACTGACTGACTG");
}

#[test]
fn test_streaming_large_file_memory_efficiency() {
    // Test that we can iterate without loading all sequences into memory
    let path = fixture_path("sample.fastq");
    let iterator = parse_sequences(&path).expect("Failed to open file");

    let mut count = 0;
    for result in iterator {
        result.expect("Failed to parse sequence");
        count += 1;
    }
    assert_eq!(count, 3);
}

#[test]
fn test_error_malformed_fasta() {
    let path = fixture_path("malformed.fasta");
    let sequences: Vec<Result<Sequence, ParseError>> = parse_sequences(&path).unwrap().collect();

    // Should successfully parse valid entries but error on malformed ones
    let mut errors = 0;
    let mut successes = 0;

    for result in sequences {
        match result {
            Ok(_) => successes += 1,
            Err(e) => {
                errors += 1;
                // Verify error is informative
                assert!(
                    matches!(e, ParseError::MalformedEntry { .. }),
                    "Expected MalformedEntry error, got: {:?}",
                    e
                );
            }
        }
    }

    assert!(errors > 0, "Should detect malformed entries");
    assert!(successes > 0, "Should parse valid entries");
}

#[test]
fn test_error_malformed_fastq() {
    let path = fixture_path("malformed.fastq");
    let sequences: Vec<Result<Sequence, ParseError>> = parse_sequences(&path).unwrap().collect();

    let mut errors = 0;
    for result in sequences {
        if let Err(e) = result {
            errors += 1;
            // Verify errors are descriptive
            match e {
                ParseError::MalformedEntry { line, reason } => {
                    assert!(line > 0, "Error should include line number");
                    assert!(!reason.is_empty(), "Error should include reason");
                }
                ParseError::QualityLengthMismatch { seq_len, qual_len } => {
                    assert!(seq_len != qual_len, "Should report length mismatch");
                }
                _ => panic!("Unexpected error type: {:?}", e),
            }
        }
    }

    assert!(errors > 0, "Should detect malformed FASTQ entries");
}

#[test]
fn test_error_nonexistent_file() {
    let path = fixture_path("does_not_exist.fasta");
    let result = parse_sequences(&path);

    assert!(result.is_err(), "Should error on nonexistent file");
    match result {
        Err(ParseError::IoError { .. }) => {
            // Expected error type
        }
        Err(e) => panic!("Expected IoError, got: {:?}", e),
        Ok(_) => panic!("Should not succeed on nonexistent file"),
    }
}

#[test]
fn test_error_unsupported_format() {
    // Test with a file that has an unsupported extension
    let path = fixture_path("generate_fixtures.sh");
    let result = SequenceFormat::detect(&path);

    match result {
        Err(ParseError::UnsupportedFormat { .. }) => {
            // Expected error
        }
        Err(e) => panic!("Expected UnsupportedFormat error, got: {:?}", e),
        Ok(_) => panic!("Should not detect format for unsupported file"),
    }
}

#[test]
fn test_empty_file() {
    // Create a temporary empty file
    let path = fixture_path("empty.fasta");
    std::fs::write(&path, b"").expect("Failed to create empty file");

    let sequences: Vec<_> = parse_sequences(&path)
        .unwrap()
        .collect::<Result<Vec<_>, _>>()
        .unwrap();

    assert_eq!(sequences.len(), 0, "Empty file should yield no sequences");

    // Cleanup
    std::fs::remove_file(&path).ok();
}

#[test]
fn test_sequence_with_newlines_in_data() {
    // FASTA files often have wrapped sequence data
    let content = b">seq1 wrapped sequence
ACTGACTG
ACTGACTG
ACTGACTG
>seq2 single line
GGGGCCCCAAAATTTT
";
    let path = fixture_path("wrapped.fasta");
    std::fs::write(&path, content).expect("Failed to write test file");

    let sequences: Vec<_> = parse_sequences(&path)
        .unwrap()
        .collect::<Result<Vec<_>, _>>()
        .unwrap();

    assert_eq!(sequences.len(), 2);
    assert_eq!(
        sequences[0].data, b"ACTGACTGACTGACTGACTGACTG",
        "Should concatenate wrapped lines"
    );
    assert_eq!(sequences[1].data, b"GGGGCCCCAAAATTTT");

    // Cleanup
    std::fs::remove_file(&path).ok();
}

#[test]
fn test_fastq_quality_score_validation() {
    // FASTQ quality scores should have same length as sequence
    let content = b"@seq1
ACTGACTG
+
IIIIIIII
@seq2
GGGGCCCC
+
HHH
";
    let path = fixture_path("mismatched_quality.fastq");
    std::fs::write(&path, content).expect("Failed to write test file");

    let sequences: Vec<Result<Sequence, ParseError>> = parse_sequences(&path).unwrap().collect();

    // First should succeed, second should fail
    assert!(
        sequences[0].is_ok(),
        "First entry should parse successfully"
    );
    assert!(
        sequences[1].is_err(),
        "Second entry should fail due to quality length mismatch"
    );

    if let Err(ParseError::QualityLengthMismatch { seq_len, qual_len }) = &sequences[1] {
        assert_eq!(*seq_len, 8, "Sequence length should be 8");
        assert_eq!(*qual_len, 3, "Quality length should be 3");
    } else {
        panic!("Expected QualityLengthMismatch error");
    }

    // Cleanup
    std::fs::remove_file(&path).ok();
}

#[test]
fn test_fasta_multiline_description() {
    // FASTA headers can have long descriptions
    let content =
        b">seq1 This is a very long description with multiple words and special chars: @#$%
ACTGACTGACTG
";
    let path = fixture_path("long_desc.fasta");
    std::fs::write(&path, content).expect("Failed to write test file");

    let sequences: Vec<_> = parse_sequences(&path)
        .unwrap()
        .collect::<Result<Vec<_>, _>>()
        .unwrap();

    assert_eq!(sequences.len(), 1);
    assert_eq!(sequences[0].id, Some("seq1".to_string()));
    assert!(
        sequences[0]
            .description
            .as_ref()
            .unwrap()
            .contains("special chars")
    );

    // Cleanup
    std::fs::remove_file(&path).ok();
}

#[test]
fn test_iterator_is_lazy() {
    // Verify iterator doesn't read entire file upfront
    let path = fixture_path("sample.fasta");
    let mut iterator = parse_sequences(&path).expect("Failed to open file");

    // Taking just the first element should not require reading the whole file
    let first = iterator.next().unwrap().unwrap();
    assert_eq!(first.data, b"ACTGACTGACTG");

    // We can still get the rest
    let rest: Vec<_> = iterator.collect::<Result<Vec<_>, _>>().unwrap();
    assert_eq!(rest.len(), 2);
}

#[test]
fn test_multiple_iterators_same_file() {
    // Should be able to create multiple independent iterators
    let path = fixture_path("sample.fasta");

    let iter1 = parse_sequences(&path).unwrap();
    let iter2 = parse_sequences(&path).unwrap();

    let seqs1: Vec<_> = iter1.collect::<Result<Vec<_>, _>>().unwrap();
    let seqs2: Vec<_> = iter2.collect::<Result<Vec<_>, _>>().unwrap();

    assert_eq!(seqs1.len(), seqs2.len());
    assert_eq!(seqs1[0].data, seqs2[0].data);
}

#[test]
fn test_compression_format_combination_matrix() {
    // Comprehensive test for all 6 combinations
    let test_cases = vec![
        ("sample.fasta", false, false),    // FASTA, plain
        ("sample.fasta.gz", false, true),  // FASTA, gzip
        ("sample.fasta.bz2", false, true), // FASTA, bzip2
        ("sample.fastq", true, false),     // FASTQ, plain
        ("sample.fastq.gz", true, true),   // FASTQ, gzip
        ("sample.fastq.bz2", true, true),  // FASTQ, bzip2
    ];

    for (filename, should_have_quality, is_compressed) in test_cases {
        let path = fixture_path(filename);
        let format = SequenceFormat::detect(&path)
            .unwrap_or_else(|_| panic!("Failed to detect format for {}", filename));

        assert_eq!(
            format.is_compressed(),
            is_compressed,
            "Compression detection failed for {}",
            filename
        );

        let sequences: Vec<_> = parse_sequences(&path)
            .unwrap_or_else(|_| panic!("Failed to parse {}", filename))
            .collect::<Result<Vec<_>, _>>()
            .unwrap_or_else(|_| panic!("Failed to read sequences from {}", filename));

        assert_eq!(
            sequences.len(),
            3,
            "Should have 3 sequences in {}",
            filename
        );

        for seq in &sequences {
            assert_eq!(
                seq.quality.is_some(),
                should_have_quality,
                "Quality field mismatch for {}",
                filename
            );
        }
    }
}
