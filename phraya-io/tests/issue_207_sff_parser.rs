/// Issue #207: feat(io): parse SFF (Standard Flowgram Format) input
///
/// This test file contains RED (failing) acceptance tests for issue #207.
/// Tests verify that SFF files are parsed into Sequence records with proper
/// ID, basecall sequences, quality scores, and clip points applied.
///
/// Acceptance Criteria:
/// 1. `.sff` files are auto-detected and parsed into `Sequence` records (bases + qualities)
/// 2. Read clip points (qual/adapter left/right) are applied to the emitted sequence
/// 3. Multi-read SFF iterates all reads
/// 4. No non-Rust dependency (zero binary deps principle)
/// 5. Test against a small sample SFF

use phraya_io::SequenceParser;
use phraya_core::types::ParseError;
use std::io::Write;
use std::path::Path;
use tempfile::NamedTempFile;

// ============================================================================
// Helper Functions for Creating Minimal SFF Files
// ============================================================================

/// Reference: https://trace.ncbi.nlm.nih.gov/Traces/trace.cgi?cmd=show&f=formats&m=doc&s=formats

/// Helper: Write SFF file with minimal valid structure
/// Creates a valid SFF header + read records
fn write_minimal_sff(path: &Path, reads: Vec<(&str, &str, &str, u16, u16)>) -> std::io::Result<()> {
    // Format: (read_id, basecalls, quality_scores, clip_qual_right, clip_adapter_right)
    let mut file = std::fs::File::create(path)?;

    // SFF header (simplified: write magic, version, num_reads, basic structure)
    file.write_all(b".sff")?;           // magic (4 bytes)
    file.write_all(&1_u32.to_be_bytes())?;  // version (4 bytes)
    file.write_all(&0_u64.to_be_bytes())?;  // index_offset (8 bytes, placeholder)
    file.write_all(&0_u32.to_be_bytes())?;  // index_length (4 bytes)
    file.write_all(&(reads.len() as u32).to_be_bytes())?;  // num_reads (4 bytes)
    file.write_all(&31_u16.to_be_bytes())?; // header_length (16 bytes after this field, standard SFF is 31+4)
    file.write_all(&4_u16.to_be_bytes())?;  // key_length (4 bytes for "TCAG")
    file.write_all(&100_u16.to_be_bytes())?; // num_flows (100 flows typical for 454)
    file.write_all(&1_u8.to_be_bytes())?;   // flowgram_format (1 = uint16)
    file.write_all(b"TACG")?;           // flow_chars (TACG is standard for 454)
    file.write_all(&[0_u8; 252])?;      // pad flow_chars to 256
    file.write_all(b"TCAG")?;           // key_sequence (TCAG is standard)
    file.write_all(&[0_u8; 252])?;      // pad key_sequence to 256

    // SFF read records
    for (read_id, bases, qualities, clip_qual_right, clip_adapter_right) in reads {
        // Read header
        let id_bytes = read_id.as_bytes();
        let name_length = id_bytes.len() as u16;
        let num_bases = bases.len() as u32;

        file.write_all(&16_u16.to_be_bytes())?; // read_header_length (fixed at 16 for simple case)
        file.write_all(&name_length.to_be_bytes())?; // name_length
        file.write_all(&num_bases.to_be_bytes())?;   // num_bases
        file.write_all(&0_u16.to_be_bytes())?;  // clip_qual_left (no clipping on left)
        file.write_all(&clip_qual_right.to_be_bytes())?;  // clip_qual_right
        file.write_all(&0_u16.to_be_bytes())?;  // clip_adapter_left
        file.write_all(&clip_adapter_right.to_be_bytes())?; // clip_adapter_right

        // Pad read header to 8-byte boundary if needed
        let header_with_name_len = 16 + name_length as usize;
        let padding = (8 - (header_with_name_len % 8)) % 8;
        if padding > 0 {
            file.write_all(&vec![0_u8; padding])?;
        }

        // Read name
        file.write_all(id_bytes)?;

        // Pad name to 8-byte boundary
        let padding = (8 - (name_length as usize % 8)) % 8;
        if padding > 0 {
            file.write_all(&vec![0_u8; padding])?;
        }

        // Flowgram values (placeholder: write as u16 values, typical 100 flows)
        for _ in 0..100 {
            file.write_all(&50_u16.to_be_bytes())?;
        }

        // Flow indices (one per base, 1 byte each)
        for _ in 0..num_bases {
            file.write_all(&1_u8.to_be_bytes())?;
        }

        // Bases (DNA characters as ASCII)
        file.write_all(bases.as_bytes())?;

        // Quality scores (ASCII, typically Phred+33)
        file.write_all(qualities.as_bytes())?;

        // Pad entire read to 8-byte boundary
        let read_size = 16 + name_length as usize + padding + (100 * 2) + num_bases as usize + num_bases as usize;
        let padding = (8 - (read_size % 8)) % 8;
        if padding > 0 {
            file.write_all(&vec![0_u8; padding])?;
        }
    }

    Ok(())
}

/// Create a minimal valid SFF file with a single read, no clipping
fn create_minimal_sff_single_read() -> NamedTempFile {
    let tmp = NamedTempFile::new().unwrap();
    write_minimal_sff(tmp.path(), vec![
        ("read_00001", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTAC", "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII", 38, 38),
    ]).unwrap();
    tmp
}

/// Create SFF file with multiple reads
fn create_minimal_sff_multiple_reads() -> NamedTempFile {
    let tmp = NamedTempFile::new().unwrap();
    write_minimal_sff(tmp.path(), vec![
        ("read_00001", "ACGTACGTACGTACGTACGT", "IIIIIIIIIIIIIIIIIIII", 20, 20),
        ("read_00002", "TGCATGCATGCATGCATGCA", "HHHHHHHHHHHHHHHHHHHH", 20, 20),
        ("read_00003", "AAATTTGGGGCCCCAAATTT", "JJJJJJJJJJJJJJJJJJJJ", 20, 20),
    ]).unwrap();
    tmp
}

/// Create SFF file with quality clipping applied (clip_qual_right < num_bases)
fn create_minimal_sff_with_quality_clipping() -> NamedTempFile {
    let tmp = NamedTempFile::new().unwrap();
    write_minimal_sff(tmp.path(), vec![
        // Full sequence: "ACGTACGTACGTACGTACGTACGTACGTACGTACGT" (38bp)
        // After qual clip right=30: "ACGTACGTACGTACGTACGTACGTACGTAC" (30bp)
        ("read_clipped", "ACGTACGTACGTACGTACGTACGTACGTACGTACGT", "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII", 30, 38),
    ]).unwrap();
    tmp
}

/// Create SFF file with adapter clipping applied
fn create_minimal_sff_with_adapter_clipping() -> NamedTempFile {
    let tmp = NamedTempFile::new().unwrap();
    write_minimal_sff(tmp.path(), vec![
        // Full sequence: "ACGTACGTACGTACGTACGTACGTACGTACGTACGTAC" (38bp)
        // After adapter clip right=35: "ACGTACGTACGTACGTACGTACGTACGTACGTAC" (35bp)
        ("read_adapter", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTAC", "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII", 38, 35),
    ]).unwrap();
    tmp
}

/// Create SFF file with combined clipping (qual_right and adapter_right both constrain)
fn create_minimal_sff_with_combined_clipping() -> NamedTempFile {
    let tmp = NamedTempFile::new().unwrap();
    write_minimal_sff(tmp.path(), vec![
        // Full sequence: "ACGTACGTACGTACGTACGTACGTACGTACGTACGT" (38bp)
        // clip_qual_right=35, clip_adapter_right=32 -> final clip = min(35, 32) = 32bp
        ("read_both", "ACGTACGTACGTACGTACGTACGTACGTACGTACGT", "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII", 35, 32),
    ]).unwrap();
    tmp
}

/// Create SFF file with invalid magic bytes (should fail parsing)
fn create_invalid_sff_bad_magic() -> NamedTempFile {
    let tmp = NamedTempFile::new().unwrap();
    let mut file = std::fs::File::create(tmp.path()).unwrap();
    file.write_all(b"XXXX").unwrap();  // Invalid magic
    file.write_all(&[0_u8; 100]).unwrap();
    tmp
}

/// Create SFF file with truncated header
fn create_invalid_sff_truncated_header() -> NamedTempFile {
    let tmp = NamedTempFile::new().unwrap();
    let mut file = std::fs::File::create(tmp.path()).unwrap();
    file.write_all(b".sff").unwrap();  // Valid magic
    file.write_all(&[0_u8; 10]).unwrap();  // Only 10 bytes instead of full header
    tmp
}

// ============================================================================
// Acceptance Tests
// ============================================================================

#[test]
fn sff_auto_detect_single_read() {
    // Test: SFF files are auto-detected by extension or magic and parsed
    let sff_file = create_minimal_sff_single_read();

    let mut parser = SequenceParser::from_path(sff_file.path()).unwrap();
    let seq = parser.next().unwrap().unwrap();

    assert_eq!(seq.id(), "read_00001");
    assert_eq!(seq.len(), 38);
    assert_eq!(seq.bases(), b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTAC");
    assert_eq!(seq.quality_scores(), Some(&vec![b'I'; 38]));
}

#[test]
fn sff_multi_read_iteration() {
    // Test: Multi-read SFF files iterate all reads
    let sff_file = create_minimal_sff_multiple_reads();

    let mut parser = SequenceParser::from_path(sff_file.path()).unwrap();

    let seq1 = parser.next().unwrap().unwrap();
    assert_eq!(seq1.id(), "read_00001");
    assert_eq!(seq1.len(), 20);
    assert_eq!(seq1.bases(), b"ACGTACGTACGTACGTACGT");

    let seq2 = parser.next().unwrap().unwrap();
    assert_eq!(seq2.id(), "read_00002");
    assert_eq!(seq2.len(), 20);
    assert_eq!(seq2.bases(), b"TGCATGCATGCATGCATGCA");

    let seq3 = parser.next().unwrap().unwrap();
    assert_eq!(seq3.id(), "read_00003");
    assert_eq!(seq3.len(), 20);
    assert_eq!(seq3.bases(), b"AAATTTGGGGCCCCAAATTT");

    assert_eq!(parser.next(), None);
}

#[test]
fn sff_quality_clip_right_boundary() {
    // Test: Quality clipping right boundary is applied (clip_qual_right)
    let sff_file = create_minimal_sff_with_quality_clipping();

    let mut parser = SequenceParser::from_path(sff_file.path()).unwrap();
    let seq = parser.next().unwrap().unwrap();

    // Full bases: 38bp, clipped to 30bp
    assert_eq!(seq.len(), 30);
    assert_eq!(seq.bases(), b"ACGTACGTACGTACGTACGTACGTACGTAC");
    assert_eq!(seq.quality_scores(), Some(&vec![b'I'; 30]));
}

#[test]
fn sff_adapter_clip_right_boundary() {
    // Test: Adapter clipping right boundary is applied (clip_adapter_right)
    let sff_file = create_minimal_sff_with_adapter_clipping();

    let mut parser = SequenceParser::from_path(sff_file.path()).unwrap();
    let seq = parser.next().unwrap().unwrap();

    // Full bases: 38bp, clipped to 35bp
    assert_eq!(seq.len(), 35);
    assert_eq!(seq.bases(), b"ACGTACGTACGTACGTACGTACGTACGTACGTACG");
    assert_eq!(seq.quality_scores(), Some(&vec![b'I'; 35]));
}

#[test]
fn sff_combined_quality_and_adapter_clipping() {
    // Test: Both quality and adapter clipping are applied (min of both)
    let sff_file = create_minimal_sff_with_combined_clipping();

    let mut parser = SequenceParser::from_path(sff_file.path()).unwrap();
    let seq = parser.next().unwrap().unwrap();

    // clip_qual_right=35, clip_adapter_right=32 -> min = 32
    assert_eq!(seq.len(), 32);
    assert_eq!(seq.bases(), b"ACGTACGTACGTACGTACGTACGTACGTACGT");
    assert_eq!(seq.quality_scores(), Some(&vec![b'I'; 32]));
}

#[test]
fn sff_invalid_magic_bytes_returns_error() {
    // Test: Invalid SFF magic bytes are rejected with ParseError
    let sff_file = create_invalid_sff_bad_magic();

    let result = SequenceParser::from_path(sff_file.path());

    assert!(result.is_err());
    if let Err(ParseError::InvalidFormat(msg)) = result {
        assert!(msg.contains("SFF") || msg.contains("magic"));
    }
}

#[test]
fn sff_truncated_header_returns_error() {
    // Test: Truncated SFF header is rejected
    let sff_file = create_invalid_sff_truncated_header();

    let result = SequenceParser::from_path(sff_file.path());

    assert!(result.is_err());
    if let Err(ParseError::InvalidFormat(msg)) = result {
        assert!(msg.contains("truncated") || msg.contains("unexpected EOF"));
    }
}

#[test]
fn sff_preserves_sequence_id() {
    // Test: SFF read names are correctly extracted as sequence IDs
    let sff_file = create_minimal_sff_single_read();

    let mut parser = SequenceParser::from_path(sff_file.path()).unwrap();
    let seq = parser.next().unwrap().unwrap();

    assert_eq!(seq.id(), "read_00001");
    // ID should not contain extra padding or artifacts
    assert!(!seq.id().contains("\0"));
}

#[test]
fn sff_preserves_base_sequence() {
    // Test: DNA bases in SFF are correctly extracted (bases field)
    let sff_file = create_minimal_sff_single_read();

    let mut parser = SequenceParser::from_path(sff_file.path()).unwrap();
    let seq = parser.next().unwrap().unwrap();

    assert_eq!(seq.bases(), b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTAC");
    // Verify only ACGT characters (after clipping applied)
    for &base in seq.bases() {
        assert!(matches!(base, b'A' | b'C' | b'G' | b'T' | b'N'));
    }
}

#[test]
fn sff_preserves_quality_scores() {
    // Test: Quality scores in SFF are correctly extracted and stored
    let sff_file = create_minimal_sff_single_read();

    let mut parser = SequenceParser::from_path(sff_file.path()).unwrap();
    let seq = parser.next().unwrap().unwrap();

    assert!(seq.quality_scores().is_some());
    let quality = seq.quality_scores().unwrap();
    assert_eq!(quality.len(), 38);
    // All quality scores should be printable ASCII (Phred+33 encoding)
    for &q in quality.iter() {
        assert!(q >= 33 && q <= 126, "quality score {} is out of Phred+33 range", q);
    }
}

#[test]
fn sff_quality_and_bases_length_match() {
    // Test: After clipping, bases and quality scores have matching length
    let sff_file = create_minimal_sff_with_quality_clipping();

    let mut parser = SequenceParser::from_path(sff_file.path()).unwrap();
    let seq = parser.next().unwrap().unwrap();

    let bases_len = seq.bases().len();
    let quality_len = seq.quality_scores().map(|q| q.len()).unwrap_or(0);

    assert_eq!(bases_len, quality_len);
}

#[test]
fn sff_extension_auto_detection() {
    // Test: Files with .sff extension are recognized and parsed
    // (not treated as FASTA/FASTQ)
    let tmp = NamedTempFile::new().unwrap();
    let sff_path = tmp.path().with_extension("sff");
    write_minimal_sff(
        &sff_path,
        vec![(
            "read_00001",
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTAC",
            "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
            38,
            38,
        )],
    )
    .unwrap();
    let path_str = sff_path.to_string_lossy();

    assert!(path_str.contains(".sff"));
    // Parser should succeed without format errors
    let result = SequenceParser::from_path(&sff_path);
    assert!(result.is_ok());
}

#[test]
fn sff_empty_file_returns_empty_iterator() {
    // Test: Empty or minimal SFF file returns empty iterator
    let tmp = NamedTempFile::new().unwrap();
    let mut file = std::fs::File::create(tmp.path()).unwrap();
    // Write only magic + minimal header, no reads
    file.write_all(b".sff").unwrap();
    file.write_all(&0_u32.to_be_bytes()).unwrap();  // version
    file.write_all(&0_u64.to_be_bytes()).unwrap();  // index_offset
    file.write_all(&0_u32.to_be_bytes()).unwrap();  // index_length
    file.write_all(&0_u32.to_be_bytes()).unwrap();  // num_reads = 0
    drop(file);

    let mut parser = SequenceParser::from_path(tmp.path()).unwrap();
    assert_eq!(parser.next(), None);
}

#[test]
fn sff_zero_clip_returns_zero_length_sequence() {
    // Test: Clipping to zero returns valid but empty Sequence
    let tmp = NamedTempFile::new().unwrap();
    write_minimal_sff(tmp.path(), vec![
        ("read_empty", "ACGTACGTACGTACGTACGTACGTACGTACGTACGT", "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII", 0, 0),
    ]).unwrap();

    let mut parser = SequenceParser::from_path(tmp.path()).unwrap();
    let seq = parser.next().unwrap().unwrap();

    assert_eq!(seq.len(), 0);
    assert_eq!(seq.id(), "read_empty");
}

#[test]
fn sff_dna_bases_validation() {
    // Test: Only valid DNA bases (ACGTN) appear in output
    let sff_file = create_minimal_sff_single_read();

    let mut parser = SequenceParser::from_path(sff_file.path()).unwrap();
    let seq = parser.next().unwrap().unwrap();

    for &base in seq.bases() {
        assert!(matches!(base, b'A' | b'C' | b'G' | b'T' | b'N'),
                "Invalid base character: {} ({})", base as char, base);
    }
}

#[test]
fn sff_quality_encoding_validation() {
    // Test: Quality scores are valid Phred+33 values
    let sff_file = create_minimal_sff_single_read();

    let mut parser = SequenceParser::from_path(sff_file.path()).unwrap();
    let seq = parser.next().unwrap().unwrap();

    if let Some(quality) = seq.quality_scores() {
        for &q in quality.iter() {
            // Phred+33: ASCII 33-126 correspond to quality 0-93
            assert!(q >= 33, "Quality score {} is below Phred+33 minimum", q);
            assert!(q <= 126, "Quality score {} is above Phred+33 maximum", q);
        }
    }
}
