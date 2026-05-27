use crate::{ParseError, SequenceFormat};
use phraya_core::Sequence;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::path::Path;

/// Iterator over sequences in a FASTA/FASTQ file
pub struct SequenceReader {
    reader: Box<dyn BufRead>,
    format: SequenceFormat,
    line_number: usize,
    /// Buffer for peeked line (for handling multi-record parsing)
    peeked_line: Option<String>,
}

impl SequenceReader {
    /// Create a new SequenceReader
    pub fn new(path: &Path, format: SequenceFormat) -> Result<Self, ParseError> {
        let path_str = path.to_string_lossy().to_string();
        let file = File::open(path).map_err(|e| ParseError::IoError {
            path: path_str.clone(),
            source: e,
        })?;

        let boxed_reader: Box<dyn BufRead> = match format {
            SequenceFormat::FastaPlain | SequenceFormat::FastqPlain => {
                Box::new(BufReader::new(file))
            }
            SequenceFormat::FastaGzip | SequenceFormat::FastqGzip => {
                let decoder = flate2::read::GzDecoder::new(file);
                Box::new(BufReader::new(decoder))
            }
            SequenceFormat::FastaBzip2 | SequenceFormat::FastqBzip2 => {
                let decoder = bzip2::read::BzDecoder::new(file);
                Box::new(BufReader::new(decoder))
            }
        };

        Ok(SequenceReader {
            reader: boxed_reader,
            format,
            line_number: 0,
            peeked_line: None,
        })
    }
}

impl Iterator for SequenceReader {
    type Item = Result<Sequence, ParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.format.is_fasta() {
            self.next_fasta()
        } else {
            self.next_fastq()
        }
    }
}

impl SequenceReader {
    fn next_line(&mut self) -> std::io::Result<(String, bool)> {
        // Check if we have a peeked line
        if let Some(line) = self.peeked_line.take() {
            return Ok((line, true)); // true = was peeked
        }

        let mut line = String::new();
        let bytes_read = self.reader.read_line(&mut line)?;

        if bytes_read == 0 {
            // EOF
            Ok((String::new(), false))
        } else {
            self.line_number += 1;
            Ok((line, false))
        }
    }

    fn next_fasta(&mut self) -> Option<Result<Sequence, ParseError>> {
        let mut line;

        // Skip until we find a header line
        loop {
            match self.next_line() {
                Ok((l, _peeked)) => {
                    if l.is_empty() {
                        return None; // EOF
                    }
                    line = l;
                    let trimmed = line.trim();

                    // Skip empty lines
                    if trimmed.is_empty() {
                        continue;
                    }

                    // Check for header
                    if trimmed.starts_with('>') {
                        line = trimmed.to_string();
                        break;
                    } else {
                        return Some(Err(ParseError::MalformedEntry {
                            line: self.line_number,
                            reason: "Expected '>' for FASTA header".to_string(),
                        }));
                    }
                }
                Err(e) => {
                    return Some(Err(ParseError::IoError {
                        path: "unknown".to_string(),
                        source: e,
                    }));
                }
            }
        }

        // Parse header
        let header = line.strip_prefix('>').unwrap_or("");
        let (id, description) = parse_header(header);

        // Read sequence data
        let mut data = Vec::new();
        loop {
            match self.next_line() {
                Ok((l, peeked)) => {
                    if l.is_empty() {
                        // EOF - return the sequence
                        return Some(Ok(Sequence {
                            id,
                            description,
                            data,
                            quality: None,
                            pairing_info: None,
                        }));
                    }

                    let trimmed = l.trim();

                    // Empty lines can appear in FASTA files
                    if trimmed.is_empty() {
                        continue;
                    }

                    // Check if this is a new header
                    if trimmed.starts_with('>') {
                        // If this was peeked, we don't need to save it
                        if !peeked {
                            // Save this line for the next iteration
                            self.peeked_line = Some(trimmed.to_string());
                        }
                        return Some(Ok(Sequence {
                            id,
                            description,
                            data,
                            quality: None,
                            pairing_info: None,
                        }));
                    }

                    // Add sequence data
                    data.extend_from_slice(trimmed.as_bytes());
                }
                Err(e) => {
                    return Some(Err(ParseError::IoError {
                        path: "unknown".to_string(),
                        source: e,
                    }));
                }
            }
        }
    }

    fn next_fastq(&mut self) -> Option<Result<Sequence, ParseError>> {
        // Skip until we find a header line
        loop {
            match self.next_line() {
                Ok((l, _peeked)) => {
                    if l.is_empty() {
                        return None; // EOF
                    }
                    let trimmed = l.trim();

                    // Skip empty lines
                    if trimmed.is_empty() {
                        continue;
                    }

                    // Check for header
                    if trimmed.starts_with('@') {
                        // Parse header
                        let header = trimmed.strip_prefix('@').unwrap_or("");
                        let (id, description) = parse_header(header);

                        // Read sequence data
                        let data = match self.next_line() {
                            Ok((seq_line, _)) => {
                                if seq_line.is_empty() {
                                    return Some(Err(ParseError::MalformedEntry {
                                        line: self.line_number,
                                        reason: "Unexpected EOF while reading sequence".to_string(),
                                    }));
                                }
                                seq_line.trim().as_bytes().to_vec()
                            }
                            Err(e) => {
                                return Some(Err(ParseError::IoError {
                                    path: "unknown".to_string(),
                                    source: e,
                                }));
                            }
                        };

                        // Read separator line
                        match self.next_line() {
                            Ok((sep_line, _)) => {
                                if sep_line.is_empty() {
                                    return Some(Err(ParseError::MalformedEntry {
                                        line: self.line_number,
                                        reason: "Unexpected EOF while reading separator"
                                            .to_string(),
                                    }));
                                }
                                let trimmed = sep_line.trim();
                                if !trimmed.starts_with('+') {
                                    return Some(Err(ParseError::MalformedEntry {
                                        line: self.line_number,
                                        reason: "Expected '+' for FASTQ separator".to_string(),
                                    }));
                                }
                            }
                            Err(e) => {
                                return Some(Err(ParseError::IoError {
                                    path: "unknown".to_string(),
                                    source: e,
                                }));
                            }
                        }

                        // Read quality line
                        let quality = match self.next_line() {
                            Ok((qual_line, _)) => {
                                if qual_line.is_empty() {
                                    return Some(Err(ParseError::MalformedEntry {
                                        line: self.line_number,
                                        reason: "Unexpected EOF while reading quality".to_string(),
                                    }));
                                }
                                qual_line.trim().as_bytes().to_vec()
                            }
                            Err(e) => {
                                return Some(Err(ParseError::IoError {
                                    path: "unknown".to_string(),
                                    source: e,
                                }));
                            }
                        };

                        // Validate quality length
                        if data.len() != quality.len() {
                            return Some(Err(ParseError::QualityLengthMismatch {
                                seq_len: data.len(),
                                qual_len: quality.len(),
                            }));
                        }

                        return Some(Ok(Sequence {
                            id,
                            description,
                            data,
                            quality: Some(quality),
                            pairing_info: None,
                        }));
                    } else {
                        return Some(Err(ParseError::MalformedEntry {
                            line: self.line_number,
                            reason: "Expected '@' for FASTQ header".to_string(),
                        }));
                    }
                }
                Err(e) => {
                    return Some(Err(ParseError::IoError {
                        path: "unknown".to_string(),
                        source: e,
                    }));
                }
            }
        }
    }
}

/// Parse header into ID and optional description
fn parse_header(header: &str) -> (Option<String>, Option<String>) {
    let parts: Vec<&str> = header.splitn(2, ' ').collect();

    let id = if parts[0].is_empty() {
        None
    } else {
        Some(parts[0].to_string())
    };

    let description = if parts.len() > 1 && !parts[1].is_empty() {
        Some(parts[1].to_string())
    } else {
        None
    };

    (id, description)
}
