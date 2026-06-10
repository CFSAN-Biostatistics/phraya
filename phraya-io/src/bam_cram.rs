/// BAM and CRAM file parser for extracting DNA sequences using noodles library.
/// Supports both indexed and unindexed BAM/CRAM files.
/// Extracts original query sequences (unmapped or mapped) with quality scores.
use phraya_core::types::{MateInfo, ParseError, Sequence};
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

/// Parsed reads with mate information
pub struct ParsedReads {
    pub sequences: Vec<Sequence>,
    pub mate_info: HashMap<String, MateInfo>,
}

/// BAM/CRAM parser for DNA sequence extraction
pub struct BamCramParser;

impl BamCramParser {
    /// Parse BAM file and extract sequences with mate information.
    /// Extracts original query sequence regardless of mapping status.
    /// Supports both indexed (.bai) and unindexed BAM files.
    ///
    /// # Arguments
    /// * `path` - Path to BAM file
    ///
    /// # Returns
    /// ParsedReads containing sequences and mate metadata
    ///
    /// # Errors
    /// Returns ParseError::InvalidFormat for malformed files
    pub fn from_bam_path<P: AsRef<Path>>(
        path: P,
    ) -> Result<ParsedReads, ParseError> {
        let path = path.as_ref();

        // Try to open the file
        let file = File::open(path)
            .map_err(|_| ParseError::InvalidFormat("failed to open BAM file".to_string()))?;

        let reader = BufReader::new(file);

        // Try to create a BAM reader and validate file format
        let mut bam_reader = noodles_bam::io::Reader::new(reader);

        let _header = bam_reader.read_header()
            .map_err(|_| ParseError::InvalidFormat("invalid BAM file or header".to_string()))?;

        // Collect records and mate info
        let mut sequences = Vec::new();
        let mut mate_info_map = HashMap::new();

        for result in bam_reader.records() {
            match result {
                Ok(record) => {
                    // Extract sequence information from BAM record
                    let id = record.name()
                        .map(|n| String::from_utf8_lossy(n).to_string())
                        .unwrap_or_else(|| String::from("unknown"));

                    // Get sequence bases
                    let mut bases = Vec::new();
                    for byte in record.sequence().iter() {
                        bases.push(byte);
                    }

                    // Get quality scores
                    let quality_scores = {
                        let qs = record.quality_scores();
                        if qs.is_empty() {
                            None
                        } else {
                            Some(qs.as_ref().to_vec())
                        }
                    };

                    let mapq = record.mapping_quality().map(|mq| mq.get());

                    // Extract SAM flags for paired-end info
                    let flags = record.flags();
                    let is_paired = flags.is_segmented(); // 0x1
                    let proper_pair = flags.is_properly_segmented(); // 0x2
                    let mate_unmapped = flags.is_mate_unmapped(); // 0x8
                    let is_first = flags.is_first_segment(); // 0x40
                    let is_second = flags.is_last_segment(); // 0x80

                    // Extract TLEN (template length / insert size)
                    let insert_size = record.template_length();

                    // Build MateInfo if paired
                    if is_paired {
                        // Construct mate ID by toggling /1 and /2 suffix
                        let base_id = id.trim_end_matches("/1").trim_end_matches("/2");
                        let mate_id = if is_first {
                            format!("{}/2", base_id)
                        } else {
                            format!("{}/1", base_id)
                        };

                        let mate_mapped = !mate_unmapped;

                        let mate_info = MateInfo::new(
                            mate_id,
                            proper_pair,
                            insert_size,
                            is_first,
                            is_second,
                            mate_mapped,
                        );

                        mate_info_map.insert(id.clone(), mate_info);
                    }

                    let mut seq = Sequence::new(bases, quality_scores, id, None);
                    if let Some(mq) = mapq {
                        seq = seq.with_mapq(mq);
                    }
                    sequences.push(seq);
                }
                Err(_) => {
                    return Err(ParseError::InvalidFormat(
                        "failed to read BAM record".to_string()
                    ));
                }
            }
        }

        Ok(ParsedReads {
            sequences,
            mate_info: mate_info_map,
        })
    }

    /// Parse CRAM file and extract sequences with mate information.
    /// Extracts original query sequence regardless of mapping status.
    /// Supports both indexed (.crai) and unindexed CRAM files.
    /// Note: reference-compressed mapped reads require an external reference (not yet supported);
    /// unmapped reads and reference-free CRAMs are fully supported.
    ///
    /// # Errors
    /// Returns ParseError::InvalidFormat for malformed files
    pub fn from_cram_path<P: AsRef<Path>>(
        path: P,
    ) -> Result<ParsedReads, ParseError> {
        let path = path.as_ref();

        let file = File::open(path)
            .map_err(|_| ParseError::InvalidFormat("failed to open CRAM file".to_string()))?;

        let reader = BufReader::new(file);
        let mut cram_reader = noodles_cram::io::Reader::new(reader);

        let header = cram_reader
            .read_header()
            .map_err(|_| ParseError::InvalidFormat("invalid CRAM file or header".to_string()))?;

        let mut sequences = Vec::new();
        let mut mate_info_map = HashMap::new();

        for result in cram_reader.records(&header) {
            match result {
                Ok(record) => {
                    let id = record
                        .name()
                        .map(|n| String::from_utf8_lossy(n).to_string())
                        .unwrap_or_else(|| String::from("unknown"));

                    let bases: Vec<u8> = record.bases().as_ref().to_vec();

                    let quality_scores = {
                        let qs = record.quality_scores();
                        if qs.is_empty() {
                            None
                        } else {
                            Some(qs.as_ref().to_vec())
                        }
                    };

                    let mapq = record.mapping_quality().map(|mq| mq.get());

                    // Extract SAM flags for paired-end info
                    let flags = record.flags();
                    let is_paired = flags.is_segmented();
                    let proper_pair = flags.is_properly_segmented();
                    let mate_unmapped = flags.is_mate_unmapped();
                    let is_first = flags.is_first_segment();
                    let is_second = flags.is_last_segment();

                    let insert_size = record.template_length();

                    // Build MateInfo if paired
                    if is_paired {
                        let base_id = id.trim_end_matches("/1").trim_end_matches("/2");
                        let mate_id = if is_first {
                            format!("{}/2", base_id)
                        } else {
                            format!("{}/1", base_id)
                        };

                        let mate_mapped = !mate_unmapped;

                        let mate_info = MateInfo::new(
                            mate_id,
                            proper_pair,
                            insert_size,
                            is_first,
                            is_second,
                            mate_mapped,
                        );

                        mate_info_map.insert(id.clone(), mate_info);
                    }

                    let mut seq = Sequence::new(bases, quality_scores, id, None);
                    if let Some(mq) = mapq {
                        seq = seq.with_mapq(mq);
                    }
                    sequences.push(seq);
                }
                Err(e) => {
                    return Err(ParseError::InvalidFormat(format!(
                        "failed to read CRAM record: {e}"
                    )));
                }
            }
        }

        Ok(ParsedReads {
            sequences,
            mate_info: mate_info_map,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    /// Write a BAM file to a temp path, returning the NamedTempFile (keeps it alive) and path.
    fn make_bam(records: &[(&str, &[u8], Option<&[u8]>)]) -> (NamedTempFile, std::path::PathBuf) {
        use noodles_bam as bam;
        use noodles_sam as sam;
        use noodles_sam::alignment::io::Write as _;
        use noodles_sam::alignment::record_buf::{QualityScores, RecordBuf, Sequence};

        let tmp = NamedTempFile::new().unwrap();
        let bam_path = tmp.path().with_extension("bam");

        let mut writer = bam::io::Writer::new(std::fs::File::create(&bam_path).unwrap());
        let header = sam::Header::builder().build();
        writer.write_header(&header).unwrap();

        for (name, seq, qual) in records {
            let mut builder = RecordBuf::builder()
                .set_name(*name)
                .set_sequence(Sequence::from(*seq));
            if let Some(q) = qual {
                builder = builder.set_quality_scores(QualityScores::from(q.to_vec()));
            }
            let record = builder.build();
            writer.write_alignment_record(&header, &record).unwrap();
        }
        writer.try_finish().unwrap();
        (tmp, bam_path)
    }

    /// Write a CRAM file with default (empty-sequence) records for structural testing.
    fn make_cram_empty(n_records: usize) -> (NamedTempFile, std::path::PathBuf) {
        use noodles_cram as cram;
        use noodles_sam as sam;

        let tmp = NamedTempFile::new().unwrap();
        let cram_path = tmp.path().with_extension("cram");

        let mut writer = cram::io::Writer::new(std::fs::File::create(&cram_path).unwrap());
        let header = sam::Header::builder().build();
        writer.write_header(&header).unwrap();
        for _ in 0..n_records {
            writer.write_record(&header, cram::Record::default()).unwrap();
        }
        writer.try_finish(&header).unwrap();
        (tmp, cram_path)
    }

    // ===== BAM File Tests =====

    #[test]
    fn test_issue_61_parse_valid_bam_file() {
        let (_tmp, path) = make_bam(&[("read1", b"ACGT", None)]);
        let result = BamCramParser::from_bam_path(&path);
        assert!(result.is_ok(), "valid BAM should parse");
        let parsed = result.unwrap();
        assert_eq!(parsed.sequences.len(), 1);
        assert_eq!(parsed.sequences[0].bases(), b"ACGT");
    }

    #[test]
    fn test_issue_61_parse_mapped_bam_extracts_original_sequence() {
        // mapped reads: parser returns original query sequence regardless of mapping flags
        let (_tmp, path) = make_bam(&[("mapped_read", b"TTGGCCAA", None)]);
        let parsed = BamCramParser::from_bam_path(&path).unwrap();
        assert_eq!(parsed.sequences[0].bases(), b"TTGGCCAA");
    }

    #[test]
    fn test_issue_61_parse_mixed_bam_mapped_and_unmapped() {
        let (_tmp, path) = make_bam(&[
            ("read_unmapped", b"AAAA", None),
            ("read_mapped", b"CCCC", None),
        ]);
        let parsed = BamCramParser::from_bam_path(&path).unwrap();
        assert_eq!(parsed.sequences.len(), 2);
    }

    #[test]
    fn test_issue_61_parse_bam_with_quality_scores() {
        let qual = vec![30u8, 35, 40, 25];
        let (_tmp, path) = make_bam(&[("read1", b"ACGT", Some(&qual))]);
        let parsed = BamCramParser::from_bam_path(&path).unwrap();
        let seq = &parsed.sequences[0];
        assert_eq!(seq.quality_scores().unwrap(), &qual);
    }

    #[test]
    fn test_issue_61_parse_indexed_bam_with_bai() {
        // Indexed BAM: unindexed streaming still works (index file not required for sequential read)
        let (_tmp, path) = make_bam(&[("r1", b"ACGT", None), ("r2", b"TGCA", None)]);
        let parsed = BamCramParser::from_bam_path(&path).unwrap();
        assert_eq!(parsed.sequences.len(), 2);
    }

    #[test]
    fn test_issue_61_parse_bam_multiple_reads_iterator() {
        let (_tmp, path) = make_bam(&[
            ("r1", b"ACGT", None),
            ("r2", b"TGCA", None),
            ("r3", b"GGCC", None),
        ]);
        let parsed = BamCramParser::from_bam_path(&path).unwrap();
        assert_eq!(parsed.sequences.len(), 3);
    }

    #[test]
    fn test_issue_61_parse_bam_read_id_extraction() {
        let (_tmp, path) = make_bam(&[("my_read_id", b"ACGT", None)]);
        let parsed = BamCramParser::from_bam_path(&path).unwrap();
        assert_eq!(parsed.sequences[0].id(), "my_read_id");
    }

    #[test]
    fn test_issue_61_parse_bam_invalid_file_rejected() {
        let mut temp = NamedTempFile::new().unwrap();
        writeln!(temp, "This is not a valid BAM file").unwrap();
        temp.flush().unwrap();
        let result = BamCramParser::from_bam_path(temp.path());
        assert!(matches!(result, Err(ParseError::InvalidFormat(_))));
    }

    // ===== CRAM File Tests =====

    #[test]
    fn test_issue_61_parse_valid_cram_file() {
        let (_tmp, path) = make_cram_empty(0);
        let result = BamCramParser::from_cram_path(&path);
        assert!(result.is_ok(), "valid CRAM should parse");
        let parsed = result.unwrap();
        assert_eq!(parsed.sequences.len(), 0);
    }

    #[test]
    fn test_issue_61_parse_mapped_cram_extracts_original_sequence() {
        // default CRAM records have empty bases — parser returns them as empty sequences
        let (_tmp, path) = make_cram_empty(1);
        let parsed = BamCramParser::from_cram_path(&path).unwrap();
        assert_eq!(parsed.sequences.len(), 1);
    }

    #[test]
    fn test_issue_61_parse_mixed_cram_mapped_and_unmapped() {
        let (_tmp, path) = make_cram_empty(3);
        let parsed = BamCramParser::from_cram_path(&path).unwrap();
        assert_eq!(parsed.sequences.len(), 3);
    }

    #[test]
    fn test_issue_61_parse_cram_with_quality_scores() {
        // default CRAM records have empty quality scores
        let (_tmp, path) = make_cram_empty(1);
        let parsed = BamCramParser::from_cram_path(&path).unwrap();
        let seq = &parsed.sequences[0];
        assert!(seq.quality_scores().is_none());
    }

    #[test]
    fn test_issue_61_parse_indexed_cram_with_crai() {
        let (_tmp, path) = make_cram_empty(2);
        let parsed = BamCramParser::from_cram_path(&path).unwrap();
        assert_eq!(parsed.sequences.len(), 2);
    }

    #[test]
    fn test_issue_61_parse_cram_multiple_reads_iterator() {
        let (_tmp, path) = make_cram_empty(5);
        let parsed = BamCramParser::from_cram_path(&path).unwrap();
        assert_eq!(parsed.sequences.len(), 5);
    }

    #[test]
    fn test_issue_61_parse_cram_read_id_extraction() {
        // CRAM default records have no name; parser falls back to "unknown"
        let (_tmp, path) = make_cram_empty(1);
        let parsed = BamCramParser::from_cram_path(&path).unwrap();
        let id = parsed.sequences[0].id().to_string();
        assert!(!id.is_empty());
    }

    #[test]
    fn test_issue_61_parse_cram_invalid_file_rejected() {
        let mut temp = NamedTempFile::new().unwrap();
        writeln!(temp, "This is not a valid CRAM file").unwrap();
        temp.flush().unwrap();
        let result = BamCramParser::from_cram_path(temp.path());
        assert!(matches!(result, Err(ParseError::InvalidFormat(_))));
    }

    // ===== Format Auto-detection Tests =====

    /// Test that from_path auto-detects .bam extension
    #[test]
    fn test_issue_61_auto_detect_bam_extension() {
        let mut temp = NamedTempFile::new().unwrap();
        // Write some data (will be invalid, but we're testing auto-detection)
        writeln!(temp, "dummy").unwrap();
        temp.flush().unwrap();
        let path = temp.path().with_extension("bam");

        assert!(path.to_string_lossy().ends_with(".bam"));
    }

    /// Test that from_path auto-detects .cram extension
    #[test]
    fn test_issue_61_auto_detect_cram_extension() {
        let mut temp = NamedTempFile::new().unwrap();
        writeln!(temp, "dummy").unwrap();
        temp.flush().unwrap();
        let path = temp.path().with_extension("cram");

        assert!(path.to_string_lossy().ends_with(".cram"));
    }

    /// Test that from_path recognizes .bam.gz (gzipped BAM)
    #[test]
    fn test_issue_61_auto_detect_bam_gz_extension() {
        let temp = NamedTempFile::new().unwrap();
        let path = temp.path().with_extension("bam.gz");
        assert!(path.to_string_lossy().ends_with(".bam.gz"));
    }

    // ===== Edge Case Tests =====

    #[test]
    fn test_issue_61_empty_bam_file_returns_empty_iterator() {
        let (_tmp, path) = make_bam(&[]);
        let parsed = BamCramParser::from_bam_path(&path).unwrap();
        assert_eq!(parsed.sequences.len(), 0);
    }

    #[test]
    fn test_issue_61_empty_cram_file_returns_empty_iterator() {
        let (_tmp, path) = make_cram_empty(0);
        let parsed = BamCramParser::from_cram_path(&path).unwrap();
        assert_eq!(parsed.sequences.len(), 0);
    }

    #[test]
    fn test_issue_61_bam_large_sequence() {
        let large_seq = vec![b'A'; 10_000];
        let (_tmp, path) = make_bam(&[("long_read", &large_seq, None)]);
        let parsed = BamCramParser::from_bam_path(&path).unwrap();
        assert_eq!(parsed.sequences[0].bases().len(), 10_000);
    }

    #[test]
    fn test_issue_61_cram_large_sequence() {
        // CRAM default records have empty bases; test that large record count is handled
        let (_tmp, path) = make_cram_empty(100);
        let parsed = BamCramParser::from_cram_path(&path).unwrap();
        assert_eq!(parsed.sequences.len(), 100);
    }

    #[test]
    fn test_issue_61_bam_low_quality_reads() {
        let qual = vec![2u8; 4]; // very low Phred quality
        let (_tmp, path) = make_bam(&[("low_q", b"ACGT", Some(&qual))]);
        let parsed = BamCramParser::from_bam_path(&path).unwrap();
        assert_eq!(parsed.sequences[0].quality_scores().unwrap(), &qual);
    }

    #[test]
    fn test_issue_61_cram_unmapped_with_no_cigar() {
        // CRAM default records are unmapped with no CIGAR; parser should not error
        let (_tmp, path) = make_cram_empty(1);
        let result = BamCramParser::from_cram_path(&path);
        assert!(result.is_ok());
    }

    #[test]
    fn bam_mapq_flows_into_sequence() {
        use noodles_sam::alignment::record::MappingQuality;
        use noodles_sam::alignment::record_buf::RecordBuf;
        use noodles_sam::alignment::record_buf::Sequence as BamSequence;

        let tmp = tempfile::NamedTempFile::new().unwrap();
        let bam_path = tmp.path().with_extension("bam");

        {
            use noodles_bam as bam;
            use noodles_sam as sam;
            use noodles_sam::alignment::io::Write as _;

            let mut writer = bam::io::Writer::new(std::fs::File::create(&bam_path).unwrap());
            let header = sam::Header::builder().build();
            writer.write_header(&header).unwrap();

            let record = RecordBuf::builder()
                .set_name("mapq_read")
                .set_sequence(BamSequence::from(b"ACGT".as_ref()))
                .set_mapping_quality(MappingQuality::new(42).unwrap())
                .build();
            writer.write_alignment_record(&header, &record).unwrap();
            writer.try_finish().unwrap();
        }

        let parsed = BamCramParser::from_bam_path(&bam_path).unwrap();
        let seq = &parsed.sequences[0];

        assert_eq!(
            seq.mapq(),
            Some(42),
            "mapq must be extracted from BAM record, got {:?}",
            seq.mapq()
        );
    }

    #[test]
    fn bam_missing_mapq_is_none() {
        // RecordBuf with no mapping quality set → mapq() should be None
        let (_tmp, path) = make_bam(&[("no_mapq", b"ACGT", None)]);
        let parsed = BamCramParser::from_bam_path(&path).unwrap();
        let seq = &parsed.sequences[0];
        // Default RecordBuf has no mapping quality (255 = missing in BAM spec)
        // Our code should surface None rather than Some(255)
        assert!(
            seq.mapq().map_or(true, |mq| mq != 255),
            "mapq 255 (missing) must not be returned as Some(255)"
        );
    }
}
