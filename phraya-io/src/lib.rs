pub mod ab1;
pub mod bam_cram;
pub mod phraya;
pub mod plan;
pub mod queries;
pub mod sff;
pub mod use_case;

pub use use_case::{classify_input, InputType, UseCaseError};

use flate2::read::GzDecoder;
use phraya_core::types::{ParseError, Sequence};
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

/// Parser for FASTA, FASTQ, and SFF files with auto-detection.
pub struct SequenceParser;

impl SequenceParser {
    /// Parse sequences from a file (FASTA, FASTQ, AB1, or SFF auto-detected, with optional gzip).
    /// Supports .fa/.fasta/.fq/.fastq/.ab1/.sff and .gz variants.
    /// Returns an iterator of Sequence objects.
    pub fn from_path<P: AsRef<Path>>(
        path: P,
    ) -> Result<Box<dyn Iterator<Item = Result<Sequence, ParseError>>>, ParseError> {
        let path = path.as_ref();
        let path_str = path.as_os_str().to_string_lossy();

        // Check for AB1 extension
        if path_str.ends_with(".ab1") {
            let seq = ab1::parse_ab1_file(path)?;
            return Ok(Box::new(std::iter::once(Ok(seq))));
        }

        // Check for SFF format (.sff extension)
        if path_str.ends_with(".sff") {
            return sff::parse_sff_file(path);
        }

        // Check for SFF by magic bytes
        let mut file = File::open(path)
            .map_err(|e| ParseError::InvalidFormat(format!("failed to open file: {}", e)))?;

        let mut magic = [0u8; 4];
        match file.read_exact(&mut magic) {
            Ok(()) => {
                if &magic == b".sff" {
                    // Reopen file for SFF parsing
                    let file = File::open(path)
                        .map_err(|e| ParseError::InvalidFormat(format!("failed to open file: {}", e)))?;
                    return sff::parse_sff_file_from_reader(Box::new(file));
                }
            }
            Err(_) => {
                // File too short to contain magic bytes, not SFF
            }
        }

        // Reopen file for text parsing
        let file = File::open(path)
            .map_err(|e| ParseError::InvalidFormat(format!("failed to open file: {}", e)))?;
        let is_gzipped = path_str.ends_with(".gz");

        if is_gzipped {
            let reader = GzDecoder::new(file);
            Self::parse_reader(reader)
        } else {
            Self::parse_reader(file)
        }
    }

    fn parse_reader<R: Read + 'static>(
        reader: R,
    ) -> Result<Box<dyn Iterator<Item = Result<Sequence, ParseError>>>, ParseError> {
        let buf_reader = BufReader::new(reader);
        let mut lines = buf_reader.lines();

        let first_line = match lines.next() {
            Some(Ok(line)) => line,
            Some(Err(e)) => return Err(ParseError::InvalidFormat(format!("io error: {}", e))),
            None => return Ok(Box::new(std::iter::empty())),
        };

        let format = if first_line.starts_with('>') {
            Format::Fasta
        } else if first_line.starts_with('@') {
            Format::Fastq
        } else if first_line.starts_with(".sff") {
            // This shouldn't happen (SFF files should have been caught by magic detection)
            // But provide a better error message
            return Err(ParseError::InvalidFormat(
                "file appears to be SFF format (starts with .sff) but could not be parsed as valid SFF".to_string(),
            ));
        } else {
            return Err(ParseError::InvalidFormat(
                format!("invalid file format; magic bytes should be '>' (FASTA), '@' (FASTQ), or '.sff' (SFF), got: '{}'",
                        first_line.chars().take(4).collect::<String>())
            ));
        };

        Ok(Box::new(SequenceIterator::new(format, first_line, lines)))
    }
}

fn parse_header(line: &str) -> (String, Option<String>) {
    let parts: Vec<&str> = line[1..].split_whitespace().collect();
    let id = parts.get(0).map(|s| s.to_string()).unwrap_or_default();
    let description = if parts.len() > 1 {
        Some(parts[1..].join(" "))
    } else {
        None
    };
    (id, description)
}

enum Format {
    Fasta,
    Fastq,
}

struct SequenceIterator {
    format: Format,
    lines: Box<dyn Iterator<Item = std::io::Result<String>>>,
    next_header: Option<String>,
}

impl SequenceIterator {
    fn new<R: BufRead + 'static>(
        format: Format,
        first_line: String,
        lines: std::io::Lines<R>,
    ) -> Self {
        SequenceIterator {
            format,
            lines: Box::new(lines),
            next_header: Some(first_line),
        }
    }
}

impl Iterator for SequenceIterator {
    type Item = Result<Sequence, ParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.format {
            Format::Fasta => self.next_fasta(),
            Format::Fastq => self.next_fastq(),
        }
    }
}

impl SequenceIterator {
    fn next_fasta(&mut self) -> Option<Result<Sequence, ParseError>> {
        let header = self.next_header.take()?;
        let (id, description) = parse_header(&header);

        let mut bases = Vec::new();

        for line_result in self.lines.by_ref() {
            match line_result {
                Ok(line) => {
                    if line.starts_with('>') {
                        self.next_header = Some(line);
                        break;
                    }
                    bases.extend_from_slice(line.as_bytes());
                }
                Err(e) => return Some(Err(ParseError::InvalidFormat(format!("io error: {}", e)))),
            }
        }

        if bases.is_empty() && self.next_header.is_none() {
            return None;
        }

        Some(Ok(Sequence::new(bases, None, id, description)))
    }

    fn next_fastq(&mut self) -> Option<Result<Sequence, ParseError>> {
        let header = self.next_header.take()?;
        let (id, description) = parse_header(&header);

        let seq_line = match self.lines.next() {
            Some(Ok(line)) => line,
            Some(Err(e)) => {
                return Some(Err(ParseError::InvalidFormat(format!("io error: {}", e))))
            }
            None => {
                return Some(Err(ParseError::InvalidFormat(
                    "unexpected EOF: missing sequence".to_string(),
                )))
            }
        };

        let bases = seq_line.as_bytes().to_vec();

        match self.lines.next() {
            Some(Ok(_)) => {}
            Some(Err(e)) => {
                return Some(Err(ParseError::InvalidFormat(format!("io error: {}", e))))
            }
            None => {
                return Some(Err(ParseError::InvalidFormat(
                    "unexpected EOF: missing plus".to_string(),
                )))
            }
        };

        let qual_line = match self.lines.next() {
            Some(Ok(line)) => line,
            Some(Err(e)) => {
                return Some(Err(ParseError::InvalidFormat(format!("io error: {}", e))))
            }
            None => {
                return Some(Err(ParseError::InvalidFormat(
                    "unexpected EOF: missing quality".to_string(),
                )))
            }
        };

        let quality = qual_line.as_bytes().to_vec();

        if quality.len() != bases.len() {
            return Some(Err(ParseError::InvalidFormat(format!(
                "quality score length ({}) != sequence length ({})",
                quality.len(),
                bases.len()
            ))));
        }

        if let Some(Ok(line)) = self.lines.next() {
            if line.starts_with('@') {
                self.next_header = Some(line);
            }
        }

        Some(Ok(Sequence::new(bases, Some(quality), id, description)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn parse_single_sequence_fasta() {
        let mut temp = NamedTempFile::new().unwrap();
        writeln!(temp, ">seq1 description").unwrap();
        writeln!(temp, "ACGT").unwrap();
        temp.flush().unwrap();

        let mut parser = SequenceParser::from_path(temp.path()).unwrap();
        let seq = parser.next().unwrap().unwrap();

        assert_eq!(seq.id(), "seq1");
        assert_eq!(seq.description(), Some("description"));
        assert_eq!(seq.len(), 4);
    }

    #[test]
    fn parse_multiple_sequences_fasta() {
        let mut temp = NamedTempFile::new().unwrap();
        writeln!(temp, ">seq1").unwrap();
        writeln!(temp, "ACGT").unwrap();
        writeln!(temp, ">seq2").unwrap();
        writeln!(temp, "TGCA").unwrap();
        temp.flush().unwrap();

        let mut parser = SequenceParser::from_path(temp.path()).unwrap();
        let seq1 = parser.next().unwrap().unwrap();
        let seq2 = parser.next().unwrap().unwrap();

        assert_eq!(seq1.id(), "seq1");
        assert_eq!(seq2.id(), "seq2");
        assert_eq!(parser.next(), None);
    }

    #[test]
    fn parse_wrapped_fasta_sequence() {
        let mut temp = NamedTempFile::new().unwrap();
        writeln!(temp, ">seq1").unwrap();
        writeln!(temp, "ACGT").unwrap();
        writeln!(temp, "TGCA").unwrap();
        writeln!(temp, "AAAA").unwrap();
        temp.flush().unwrap();

        let mut parser = SequenceParser::from_path(temp.path()).unwrap();
        let seq = parser.next().unwrap().unwrap();

        assert_eq!(seq.len(), 12);
    }

    #[test]
    fn parse_single_sequence_fastq() {
        let mut temp = NamedTempFile::new().unwrap();
        writeln!(temp, "@seq1 description").unwrap();
        writeln!(temp, "ACGT").unwrap();
        writeln!(temp, "+").unwrap();
        writeln!(temp, "IIII").unwrap();
        temp.flush().unwrap();

        let mut parser = SequenceParser::from_path(temp.path()).unwrap();
        let seq = parser.next().unwrap().unwrap();

        assert_eq!(seq.id(), "seq1");
        assert_eq!(seq.description(), Some("description"));
        assert_eq!(seq.len(), 4);
        assert_eq!(seq.quality_at(0), Some(b'I'));
    }

    #[test]
    fn parse_multiple_sequences_fastq() {
        let mut temp = NamedTempFile::new().unwrap();
        writeln!(temp, "@seq1").unwrap();
        writeln!(temp, "ACGT").unwrap();
        writeln!(temp, "+").unwrap();
        writeln!(temp, "IIII").unwrap();
        writeln!(temp, "@seq2").unwrap();
        writeln!(temp, "TGCA").unwrap();
        writeln!(temp, "+").unwrap();
        writeln!(temp, "HHHH").unwrap();
        temp.flush().unwrap();

        let mut parser = SequenceParser::from_path(temp.path()).unwrap();
        let seq1 = parser.next().unwrap().unwrap();
        let seq2 = parser.next().unwrap().unwrap();

        assert_eq!(seq1.id(), "seq1");
        assert_eq!(seq2.id(), "seq2");
    }

    #[test]
    fn fastq_quality_length_validation() {
        let mut temp = NamedTempFile::new().unwrap();
        writeln!(temp, "@seq1").unwrap();
        writeln!(temp, "ACGT").unwrap();
        writeln!(temp, "+").unwrap();
        writeln!(temp, "II").unwrap(); // Too short
        temp.flush().unwrap();

        let mut parser = SequenceParser::from_path(temp.path()).unwrap();
        let result = parser.next().unwrap();

        assert!(result.is_err());
        if let Err(ParseError::InvalidFormat(msg)) = result {
            assert!(msg.contains("quality score length"));
        }
    }

    #[test]
    fn empty_file_returns_empty_iterator() {
        let temp = NamedTempFile::new().unwrap();
        let mut parser = SequenceParser::from_path(temp.path()).unwrap();
        assert_eq!(parser.next(), None);
    }

    #[test]
    fn invalid_format_magic_byte() {
        let mut temp = NamedTempFile::new().unwrap();
        writeln!(temp, "ACGT").unwrap(); // No > or @
        temp.flush().unwrap();

        let result = SequenceParser::from_path(temp.path());
        assert!(result.is_err());
    }

    #[test]
    fn fastq_missing_quality_line() {
        let mut temp = NamedTempFile::new().unwrap();
        writeln!(temp, "@seq1").unwrap();
        writeln!(temp, "ACGT").unwrap();
        writeln!(temp, "+").unwrap();
        temp.flush().unwrap();

        let mut parser = SequenceParser::from_path(temp.path()).unwrap();
        let result = parser.next().unwrap();

        assert!(matches!(result, Err(ParseError::InvalidFormat(ref msg)) if msg.contains("missing quality")));
    }

    #[test]
    fn parse_gzipped_fasta() {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::Write;

        let temp = NamedTempFile::new().unwrap();
        let mut encoder = GzEncoder::new(temp.as_file(), Compression::default());
        writeln!(encoder, ">seq1 description").unwrap();
        writeln!(encoder, "ACGT").unwrap();
        let _file = encoder.finish().unwrap();

        let path = std::path::PathBuf::from(temp.path()).with_extension("fa.gz");
        std::fs::rename(temp.path(), &path).unwrap();

        let mut parser = SequenceParser::from_path(&path).unwrap();
        let seq = parser.next().unwrap().unwrap();

        assert_eq!(seq.id(), "seq1");
        assert_eq!(seq.description(), Some("description"));
        assert_eq!(seq.len(), 4);

        std::fs::remove_file(&path).unwrap();
    }

    #[test]
    fn parse_gzipped_fastq() {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::Write;

        let temp = NamedTempFile::new().unwrap();
        let mut encoder = GzEncoder::new(temp.as_file(), Compression::default());
        writeln!(encoder, "@seq1 description").unwrap();
        writeln!(encoder, "ACGT").unwrap();
        writeln!(encoder, "+").unwrap();
        writeln!(encoder, "IIII").unwrap();
        let _file = encoder.finish().unwrap();

        let path = std::path::PathBuf::from(temp.path()).with_extension("fq.gz");
        std::fs::rename(temp.path(), &path).unwrap();

        let mut parser = SequenceParser::from_path(&path).unwrap();
        let seq = parser.next().unwrap().unwrap();

        assert_eq!(seq.id(), "seq1");
        assert_eq!(seq.len(), 4);
        assert_eq!(seq.quality_at(0), Some(b'I'));

        std::fs::remove_file(&path).unwrap();
    }

    #[test]
    fn first_line_invalid_utf8_returns_io_error() {
        let mut temp = NamedTempFile::new().unwrap();
        temp.write_all(&[0xFF, 0xFE, b'\n']).unwrap();
        temp.flush().unwrap();

        let result = SequenceParser::from_path(temp.path());
        assert!(matches!(result, Err(ParseError::InvalidFormat(_))));
    }

    #[test]
    fn fasta_trailing_header_with_no_body_returns_none() {
        // A FASTA file ending in a header line with nothing after it: that
        // dangling header has no bases and no further header follows it, so it
        // must be dropped (None) rather than emitted as an empty-bases Sequence.
        let mut temp = NamedTempFile::new().unwrap();
        writeln!(temp, ">seq1").unwrap();
        writeln!(temp, "ACGT").unwrap();
        writeln!(temp, ">seq2").unwrap();
        temp.flush().unwrap();

        let mut parser = SequenceParser::from_path(temp.path()).unwrap();
        let seq1 = parser.next().unwrap().unwrap();
        assert_eq!(seq1.id(), "seq1");

        assert!(parser.next().is_none());
    }

    #[test]
    fn fastq_invalid_utf8_sequence_line_returns_io_error() {
        let mut temp = NamedTempFile::new().unwrap();
        writeln!(temp, "@seq1").unwrap();
        temp.write_all(&[0xFF, 0xFE, b'\n']).unwrap();
        temp.flush().unwrap();

        let mut parser = SequenceParser::from_path(temp.path()).unwrap();
        let result = parser.next().unwrap();
        assert!(matches!(result, Err(ParseError::InvalidFormat(_))));
    }

    #[test]
    fn fastq_missing_sequence_line_returns_eof_error() {
        let mut temp = NamedTempFile::new().unwrap();
        writeln!(temp, "@seq1").unwrap();
        temp.flush().unwrap();

        let mut parser = SequenceParser::from_path(temp.path()).unwrap();
        let result = parser.next().unwrap();
        assert!(matches!(result, Err(ParseError::InvalidFormat(msg)) if msg.contains("missing sequence")));
    }

    #[test]
    fn fastq_invalid_utf8_plus_line_returns_io_error() {
        let mut temp = NamedTempFile::new().unwrap();
        writeln!(temp, "@seq1").unwrap();
        writeln!(temp, "ACGT").unwrap();
        temp.write_all(&[0xFF, 0xFE, b'\n']).unwrap();
        temp.flush().unwrap();

        let mut parser = SequenceParser::from_path(temp.path()).unwrap();
        let result = parser.next().unwrap();
        assert!(matches!(result, Err(ParseError::InvalidFormat(_))));
    }

    #[test]
    fn fastq_missing_plus_line_returns_eof_error() {
        let mut temp = NamedTempFile::new().unwrap();
        writeln!(temp, "@seq1").unwrap();
        writeln!(temp, "ACGT").unwrap();
        temp.flush().unwrap();

        let mut parser = SequenceParser::from_path(temp.path()).unwrap();
        let result = parser.next().unwrap();
        assert!(matches!(result, Err(ParseError::InvalidFormat(msg)) if msg.contains("missing plus")));
    }

    #[test]
    fn fastq_invalid_utf8_quality_line_returns_io_error() {
        let mut temp = NamedTempFile::new().unwrap();
        writeln!(temp, "@seq1").unwrap();
        writeln!(temp, "ACGT").unwrap();
        writeln!(temp, "+").unwrap();
        temp.write_all(&[0xFF, 0xFE, b'\n']).unwrap();
        temp.flush().unwrap();

        let mut parser = SequenceParser::from_path(temp.path()).unwrap();
        let result = parser.next().unwrap();
        assert!(matches!(result, Err(ParseError::InvalidFormat(_))));
    }

}
