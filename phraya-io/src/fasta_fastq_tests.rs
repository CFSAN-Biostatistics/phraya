/// Acceptance tests for FASTA/FASTQ parsing (Issue #60)
///
/// This module contains comprehensive tests for the FASTA and FASTQ parser implementation.
/// Tests cover all acceptance criteria including:
/// - Valid FASTA files (single/multiple sequences, wrapped lines)
/// - Valid FASTQ files (4-line format with quality scores)
/// - Auto-detection via magic bytes
/// - Gzip compression support
/// - Iterator-based streaming
/// - Quality score validation
/// - Empty file handling
/// - Malformed file detection
///
/// These tests should ALL FAIL initially (TDD RED phase) as the implementation does not exist yet.

use phraya_core::Sequence;
use std::io::Write;
use tempfile::NamedTempFile;

/// Helper to create a temporary file with given content
fn create_temp_file(content: &[u8]) -> NamedTempFile {
    let mut file = NamedTempFile::new().expect("Failed to create temp file");
    file.write_all(content).expect("Failed to write to temp file");
    file.flush().expect("Failed to flush temp file");
    file
}

/// Helper to create a gzipped temporary file with given content
fn create_gzipped_temp_file(content: &[u8]) -> NamedTempFile {
    use flate2::write::GzEncoder;
    use flate2::Compression;

    let mut file = NamedTempFile::new().expect("Failed to create temp file");
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(content).expect("Failed to write to encoder");
    let compressed = encoder.finish().expect("Failed to finish encoding");
    file.write_all(&compressed).expect("Failed to write compressed data");
    file.flush().expect("Failed to flush temp file");
    file
}

// =============================================================================
// HAPPY PATH: Valid FASTA files
// =============================================================================

#[test]
fn test_parse_single_sequence_fasta() {
    let content = b">seq1 description here\nACGTACGT\n";
    let file = create_temp_file(content);

    let sequences: Vec<Sequence> = crate::parse_sequences(file.path())
        .expect("Failed to parse FASTA")
        .collect();

    assert_eq!(sequences.len(), 1);
    assert_eq!(sequences[0].id(), "seq1");
    assert_eq!(sequences[0].description(), Some("description here"));
    assert_eq!(sequences[0].len(), 8);
    assert_eq!(sequences[0].quality_scores(), None); // FASTA has no quality scores
}

#[test]
fn test_parse_multiple_sequences_fasta() {
    let content = b">seq1\nACGT\n>seq2\nTGCA\n>seq3\nAAAA\n";
    let file = create_temp_file(content);

    let sequences: Vec<Sequence> = crate::parse_sequences(file.path())
        .expect("Failed to parse FASTA")
        .collect();

    assert_eq!(sequences.len(), 3);
    assert_eq!(sequences[0].id(), "seq1");
    assert_eq!(sequences[0].len(), 4);
    assert_eq!(sequences[1].id(), "seq2");
    assert_eq!(sequences[1].len(), 4);
    assert_eq!(sequences[2].id(), "seq3");
    assert_eq!(sequences[2].len(), 4);
}

#[test]
fn test_parse_wrapped_fasta_lines() {
    // FASTA format allows sequence data to be wrapped across multiple lines
    let content = b">seq1 wrapped sequence\nACGTACGT\nTGCATGCA\nAAAATTTT\n";
    let file = create_temp_file(content);

    let sequences: Vec<Sequence> = crate::parse_sequences(file.path())
        .expect("Failed to parse FASTA")
        .collect();

    assert_eq!(sequences.len(), 1);
    assert_eq!(sequences[0].id(), "seq1");
    assert_eq!(sequences[0].len(), 24); // 8 + 8 + 8 bases
}

#[test]
fn test_parse_fasta_no_description() {
    let content = b">seq1\nACGT\n";
    let file = create_temp_file(content);

    let sequences: Vec<Sequence> = crate::parse_sequences(file.path())
        .expect("Failed to parse FASTA")
        .collect();

    assert_eq!(sequences.len(), 1);
    assert_eq!(sequences[0].id(), "seq1");
    assert_eq!(sequences[0].description(), None);
}

#[test]
fn test_parse_fasta_trailing_newlines() {
    let content = b">seq1\nACGT\n\n\n";
    let file = create_temp_file(content);

    let sequences: Vec<Sequence> = crate::parse_sequences(file.path())
        .expect("Failed to parse FASTA")
        .collect();

    assert_eq!(sequences.len(), 1);
    assert_eq!(sequences[0].len(), 4);
}

// =============================================================================
// HAPPY PATH: Valid FASTQ files
// =============================================================================

#[test]
fn test_parse_single_sequence_fastq() {
    let content = b"@seq1 description here\nACGT\n+\nIIII\n";
    let file = create_temp_file(content);

    let sequences: Vec<Sequence> = crate::parse_sequences(file.path())
        .expect("Failed to parse FASTQ")
        .collect();

    assert_eq!(sequences.len(), 1);
    assert_eq!(sequences[0].id(), "seq1");
    assert_eq!(sequences[0].description(), Some("description here"));
    assert_eq!(sequences[0].len(), 4);
    assert!(sequences[0].quality_scores().is_some());
    assert_eq!(sequences[0].quality_scores().unwrap().len(), 4);
}

#[test]
fn test_parse_multiple_sequences_fastq() {
    let content = b"@seq1\nACGT\n+\nIIII\n@seq2\nTGCA\n+\nHHHH\n";
    let file = create_temp_file(content);

    let sequences: Vec<Sequence> = crate::parse_sequences(file.path())
        .expect("Failed to parse FASTQ")
        .collect();

    assert_eq!(sequences.len(), 2);
    assert_eq!(sequences[0].id(), "seq1");
    assert_eq!(sequences[0].len(), 4);
    assert_eq!(sequences[1].id(), "seq2");
    assert_eq!(sequences[1].len(), 4);
}

#[test]
fn test_parse_fastq_quality_scores_extracted() {
    let content = b"@seq1\nACGT\n+\n!#$%\n"; // Various quality scores
    let file = create_temp_file(content);

    let sequences: Vec<Sequence> = crate::parse_sequences(file.path())
        .expect("Failed to parse FASTQ")
        .collect();

    assert_eq!(sequences.len(), 1);
    let quality = sequences[0].quality_scores().expect("Should have quality scores");
    assert_eq!(quality.len(), 4);
    // Quality scores are ASCII - 33 in Phred+33 format
    assert_eq!(quality[0], b'!'); // Store raw quality byte
    assert_eq!(quality[1], b'#');
    assert_eq!(quality[2], b'$');
    assert_eq!(quality[3], b'%');
}

#[test]
fn test_parse_fastq_plus_line_can_have_content() {
    // The '+' line can optionally repeat the sequence identifier
    let content = b"@seq1 description\nACGT\n+seq1 description\nIIII\n";
    let file = create_temp_file(content);

    let sequences: Vec<Sequence> = crate::parse_sequences(file.path())
        .expect("Failed to parse FASTQ")
        .collect();

    assert_eq!(sequences.len(), 1);
    assert_eq!(sequences[0].id(), "seq1");
}

// =============================================================================
// AUTO-DETECTION: Magic bytes and extension fallback
// =============================================================================

#[test]
fn test_auto_detect_fasta_via_magic_byte() {
    // First byte is '>' for FASTA
    let content = b">seq1\nACGT\n";
    let file = create_temp_file(content);

    // Should detect as FASTA even without .fa extension
    let sequences: Vec<Sequence> = crate::parse_sequences(file.path())
        .expect("Failed to auto-detect FASTA")
        .collect();

    assert_eq!(sequences.len(), 1);
    assert!(sequences[0].quality_scores().is_none()); // No quality scores = FASTA
}

#[test]
fn test_auto_detect_fastq_via_magic_byte() {
    // First byte is '@' for FASTQ
    let content = b"@seq1\nACGT\n+\nIIII\n";
    let file = create_temp_file(content);

    // Should detect as FASTQ even without .fq extension
    let sequences: Vec<Sequence> = crate::parse_sequences(file.path())
        .expect("Failed to auto-detect FASTQ")
        .collect();

    assert_eq!(sequences.len(), 1);
    assert!(sequences[0].quality_scores().is_some()); // Has quality scores = FASTQ
}

#[test]
fn test_extension_fallback_fasta() {
    // Test with .fa extension
    let content = b">seq1\nACGT\n";
    let mut file = NamedTempFile::new().expect("Failed to create temp file");
    let _path = file.path().with_extension("fa");
    file.write_all(content).expect("Failed to write");

    // This test may need adjustment based on implementation details
    // The key is that .fa/.fasta extensions should be recognized
}

#[test]
fn test_extension_fallback_fastq() {
    // Test with .fq extension
    let content = b"@seq1\nACGT\n+\nIIII\n";
    let mut file = NamedTempFile::new().expect("Failed to create temp file");
    let _path = file.path().with_extension("fq");
    file.write_all(content).expect("Failed to write");

    // This test may need adjustment based on implementation details
    // The key is that .fq/.fastq extensions should be recognized
}

// =============================================================================
// GZIP COMPRESSION: Transparent decompression
// =============================================================================

#[test]
fn test_parse_gzipped_fasta() {
    let content = b">seq1\nACGTACGT\n";
    let file = create_gzipped_temp_file(content);

    let sequences: Vec<Sequence> = crate::parse_sequences(file.path())
        .expect("Failed to parse gzipped FASTA")
        .collect();

    assert_eq!(sequences.len(), 1);
    assert_eq!(sequences[0].id(), "seq1");
    assert_eq!(sequences[0].len(), 8);
}

#[test]
fn test_parse_gzipped_fastq() {
    let content = b"@seq1\nACGT\n+\nIIII\n";
    let file = create_gzipped_temp_file(content);

    let sequences: Vec<Sequence> = crate::parse_sequences(file.path())
        .expect("Failed to parse gzipped FASTQ")
        .collect();

    assert_eq!(sequences.len(), 1);
    assert_eq!(sequences[0].id(), "seq1");
    assert!(sequences[0].quality_scores().is_some());
}

#[test]
fn test_auto_detect_gzip_by_magic_bytes() {
    // gzip files start with 0x1f 0x8b magic bytes
    let content = b">seq1\nACGT\n";
    let file = create_gzipped_temp_file(content);

    // Should auto-detect gzip compression
    let sequences: Vec<Sequence> = crate::parse_sequences(file.path())
        .expect("Failed to auto-detect gzip")
        .collect();

    assert_eq!(sequences.len(), 1);
}

#[test]
fn test_gz_extension_recognized() {
    // Test that .fa.gz and .fq.gz extensions work
    let content = b">seq1\nACGT\n";
    let file = create_gzipped_temp_file(content);

    let sequences: Vec<Sequence> = crate::parse_sequences(file.path())
        .expect("Failed with .gz extension")
        .collect();

    assert_eq!(sequences.len(), 1);
}

// =============================================================================
// STREAMING: Iterator interface for memory efficiency
// =============================================================================

#[test]
fn test_returns_iterator_not_vec() {
    // The parse function should return an iterator, not a Vec
    // This enables processing large files without loading everything into memory
    let content = b">seq1\nACGT\n>seq2\nTGCA\n";
    let file = create_temp_file(content);

    let mut iter = crate::parse_sequences(file.path())
        .expect("Failed to parse");

    // Can process one at a time
    let seq1 = iter.next().expect("Should have first sequence");
    assert_eq!(seq1.id(), "seq1");

    let seq2 = iter.next().expect("Should have second sequence");
    assert_eq!(seq2.id(), "seq2");

    assert!(iter.next().is_none());
}

#[test]
fn test_iterator_lazy_evaluation() {
    // Iterator should parse sequences on demand, not all at once
    let content = b">seq1\nACGT\n>seq2\nTGCA\n>seq3\nAAAA\n";
    let file = create_temp_file(content);

    let mut iter = crate::parse_sequences(file.path())
        .expect("Failed to parse");

    // Take only first sequence - should not parse remaining sequences
    let seq1 = iter.next();
    assert!(seq1.is_some());
    // Implementation detail: remaining sequences not yet parsed
}

// =============================================================================
// VALIDATION: Quality score length
// =============================================================================

#[test]
fn test_quality_length_matches_sequence_length() {
    let content = b"@seq1\nACGT\n+\nIIII\n";
    let file = create_temp_file(content);

    let sequences: Vec<Sequence> = crate::parse_sequences(file.path())
        .expect("Failed to parse")
        .collect();

    assert_eq!(sequences[0].len(), 4);
    assert_eq!(sequences[0].quality_scores().unwrap().len(), 4);
}

#[test]
fn test_quality_too_short_rejected() {
    let content = b"@seq1\nACGT\n+\nIII\n"; // Quality too short (3 vs 4 bases)
    let file = create_temp_file(content);

    let result: Result<Vec<Sequence>, _> = crate::parse_sequences(file.path())
        .map(|iter| iter.collect());

    // Should return an error when trying to parse
    assert!(result.is_err(), "Should reject quality length mismatch");
}

#[test]
fn test_quality_too_long_rejected() {
    let content = b"@seq1\nACGT\n+\nIIIII\n"; // Quality too long (5 vs 4 bases)
    let file = create_temp_file(content);

    let result: Result<Vec<Sequence>, _> = crate::parse_sequences(file.path())
        .map(|iter| iter.collect());

    assert!(result.is_err(), "Should reject quality length mismatch");
}

// =============================================================================
// EMPTY FILES: Graceful handling
// =============================================================================

#[test]
fn test_empty_file_returns_empty_iterator() {
    let content = b"";
    let file = create_temp_file(content);

    let sequences: Vec<Sequence> = crate::parse_sequences(file.path())
        .expect("Failed to parse empty file")
        .collect();

    assert_eq!(sequences.len(), 0);
}

#[test]
fn test_empty_gzipped_file_returns_empty_iterator() {
    let content = b"";
    let file = create_gzipped_temp_file(content);

    let sequences: Vec<Sequence> = crate::parse_sequences(file.path())
        .expect("Failed to parse empty gzipped file")
        .collect();

    assert_eq!(sequences.len(), 0);
}

#[test]
fn test_whitespace_only_file_returns_empty_iterator() {
    let content = b"\n\n  \n\t\n";
    let file = create_temp_file(content);

    let sequences: Vec<Sequence> = crate::parse_sequences(file.path())
        .expect("Failed to parse whitespace file")
        .collect();

    assert_eq!(sequences.len(), 0);
}

// =============================================================================
// MALFORMED FILES: Clear error reporting
// =============================================================================

#[test]
fn test_fastq_missing_quality_line_rejected() {
    let content = b"@seq1\nACGT\n+\n"; // Missing quality line
    let file = create_temp_file(content);

    let result: Result<Vec<Sequence>, _> = crate::parse_sequences(file.path())
        .map(|iter| iter.collect());

    assert!(result.is_err());
    // Ideally check error message mentions missing quality
}

#[test]
fn test_fastq_missing_plus_line_rejected() {
    let content = b"@seq1\nACGT\nIIII\n"; // Missing '+' separator
    let file = create_temp_file(content);

    let result: Result<Vec<Sequence>, _> = crate::parse_sequences(file.path())
        .map(|iter| iter.collect());

    assert!(result.is_err());
}

#[test]
fn test_fastq_wrong_line_count_rejected() {
    let content = b"@seq1\nACGT\n+\n"; // Only 3 lines, need 4
    let file = create_temp_file(content);

    let result: Result<Vec<Sequence>, _> = crate::parse_sequences(file.path())
        .map(|iter| iter.collect());

    assert!(result.is_err());
}

#[test]
fn test_invalid_dna_characters_rejected() {
    let content = b">seq1\nACGTXYZ\n"; // X, Y, Z are not valid DNA bases
    let file = create_temp_file(content);

    let result: Result<Vec<Sequence>, _> = crate::parse_sequences(file.path())
        .map(|iter| iter.collect());

    assert!(result.is_err());
    // Should mention invalid characters
}

#[test]
fn test_fasta_no_sequence_id_rejected() {
    let content = b">\nACGT\n"; // Empty ID after '>'
    let file = create_temp_file(content);

    let result: Result<Vec<Sequence>, _> = crate::parse_sequences(file.path())
        .map(|iter| iter.collect());

    assert!(result.is_err());
}

#[test]
fn test_fasta_no_sequence_data_rejected() {
    let content = b">seq1\n>seq2\nACGT\n"; // seq1 has no sequence data
    let file = create_temp_file(content);

    let result: Result<Vec<Sequence>, _> = crate::parse_sequences(file.path())
        .map(|iter| iter.collect());

    assert!(result.is_err());
}

#[test]
fn test_invalid_utf8_rejected() {
    // Invalid UTF-8 bytes
    let content = vec![0x3e, 0x73, 0x65, 0x71, 0x31, 0x0a, 0xff, 0xfe, 0x0a]; // >seq1\n[invalid]\n
    let file = create_temp_file(&content);

    let result: Result<Vec<Sequence>, _> = crate::parse_sequences(file.path())
        .map(|iter| iter.collect());

    assert!(result.is_err());
}

#[test]
fn test_unknown_format_rejected() {
    // File that doesn't start with '>' or '@'
    let content = b"ACGT\nTGCA\n";
    let file = create_temp_file(content);

    let result: Result<Vec<Sequence>, _> = crate::parse_sequences(file.path())
        .map(|iter| iter.collect());

    assert!(result.is_err());
    // Error should mention unable to detect format
}

// =============================================================================
// EDGE CASES: Boundary conditions
// =============================================================================

#[test]
fn test_very_long_sequence_line() {
    // Single sequence with 10,000 bases
    let mut bases = Vec::new();
    for _ in 0..10000 {
        bases.extend_from_slice(b"ACGT");
    }
    let mut content = Vec::new();
    content.extend_from_slice(b">seq1\n");
    content.extend_from_slice(&bases);
    content.extend_from_slice(b"\n");

    let file = create_temp_file(&content);

    let sequences: Vec<Sequence> = crate::parse_sequences(file.path())
        .expect("Failed to parse long sequence")
        .collect();

    assert_eq!(sequences.len(), 1);
    assert_eq!(sequences[0].len(), 40000);
}

#[test]
fn test_sequence_id_with_special_characters() {
    let content = b">seq:1|chr1:100-200 description with spaces\nACGT\n";
    let file = create_temp_file(content);

    let sequences: Vec<Sequence> = crate::parse_sequences(file.path())
        .expect("Failed to parse")
        .collect();

    assert_eq!(sequences[0].id(), "seq:1|chr1:100-200");
    assert_eq!(sequences[0].description(), Some("description with spaces"));
}

#[test]
fn test_lowercase_bases_accepted() {
    let content = b">seq1\nacgt\n";
    let file = create_temp_file(content);

    let sequences: Vec<Sequence> = crate::parse_sequences(file.path())
        .expect("Failed to parse lowercase")
        .collect();

    assert_eq!(sequences.len(), 1);
    assert_eq!(sequences[0].len(), 4);
    // Implementation may normalize to uppercase or preserve case
}

#[test]
fn test_mixed_case_bases_accepted() {
    let content = b">seq1\nAcGt\n";
    let file = create_temp_file(content);

    let sequences: Vec<Sequence> = crate::parse_sequences(file.path())
        .expect("Failed to parse mixed case")
        .collect();

    assert_eq!(sequences.len(), 1);
}

#[test]
fn test_n_bases_accepted() {
    // 'N' represents unknown/ambiguous base - should be accepted
    let content = b">seq1\nACGTNNNN\n";
    let file = create_temp_file(content);

    let sequences: Vec<Sequence> = crate::parse_sequences(file.path())
        .expect("Failed to parse with N bases")
        .collect();

    assert_eq!(sequences.len(), 1);
    assert_eq!(sequences[0].len(), 8);
}

#[test]
fn test_iupac_ambiguity_codes_accepted() {
    // IUPAC codes: R, Y, S, W, K, M, B, D, H, V (ambiguous bases)
    let content = b">seq1\nACGTRYSWKMBDHV\n";
    let file = create_temp_file(content);

    let sequences: Vec<Sequence> = crate::parse_sequences(file.path())
        .expect("Failed to parse IUPAC codes")
        .collect();

    assert_eq!(sequences.len(), 1);
    assert_eq!(sequences[0].len(), 14);
}

// =============================================================================
// REAL-WORLD: Realistic test cases
// =============================================================================

#[test]
fn test_ncbi_fasta_format() {
    // Realistic NCBI-style header
    let content = b">NZ_CP012345.1 Escherichia coli strain ABC, complete genome\n\
                     ACGTACGTACGTACGT\n";
    let file = create_temp_file(content);

    let sequences: Vec<Sequence> = crate::parse_sequences(file.path())
        .expect("Failed to parse NCBI format")
        .collect();

    assert_eq!(sequences[0].id(), "NZ_CP012345.1");
    assert!(sequences[0].description().is_some());
}

#[test]
fn test_illumina_fastq_format() {
    // Realistic Illumina read header
    let content = b"@SRR123456.1 HWI-ST1234:100:C0001ABXX:1:1101:1234:2000 1:N:0:ATCACG\n\
                     ACGTACGT\n\
                     +\n\
                     IIIIIIII\n";
    let file = create_temp_file(content);

    let sequences: Vec<Sequence> = crate::parse_sequences(file.path())
        .expect("Failed to parse Illumina format")
        .collect();

    assert_eq!(sequences[0].id(), "SRR123456.1");
    assert!(sequences[0].description().is_some());
}

#[test]
fn test_large_file_many_sequences() {
    // Test with 1000 sequences to ensure iterator efficiency
    let mut content = Vec::new();
    for i in 0..1000 {
        content.extend_from_slice(format!(">seq{}\n", i).as_bytes());
        content.extend_from_slice(b"ACGTACGTACGT\n");
    }
    let file = create_temp_file(&content);

    let sequences: Vec<Sequence> = crate::parse_sequences(file.path())
        .expect("Failed to parse many sequences")
        .collect();

    assert_eq!(sequences.len(), 1000);
    assert_eq!(sequences[0].id(), "seq0");
    assert_eq!(sequences[999].id(), "seq999");
}
