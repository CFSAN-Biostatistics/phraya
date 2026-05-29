use phraya_core::Sequence;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

/// Parse sequences from a FASTA or FASTQ file (auto-detecting format).
/// Returns an iterator of Sequence objects for memory-efficient streaming.
///
/// Validates the entire file format before returning. If any errors are detected,
/// returns Err. Otherwise, returns Ok(iterator) where the iterator yields all valid sequences.
pub fn parse_sequences<P: AsRef<Path>>(
    path: P,
) -> Result<Box<dyn Iterator<Item = Sequence>>, phraya_core::ParseError> {
    let path = path.as_ref();

    // Open file
    let file = File::open(path).map_err(|e| {
        phraya_core::ParseError::InvalidFormat(format!("Failed to open file: {}", e))
    })?;

    // Detect if gzipped via magic bytes
    let reader: Box<dyn Read> = if is_gzip(&path)? {
        Box::new(flate2::read::GzDecoder::new(file))
    } else {
        Box::new(file)
    };

    let buf_reader = BufReader::new(reader);
    let mut lines_iter = buf_reader.lines();

    // Read all lines to detect format and create iterator
    let mut lines = Vec::new();
    let mut first_non_empty_char: Option<char> = None;

    loop {
        match lines_iter.next() {
            Some(Ok(line)) => {
                if first_non_empty_char.is_none() && !line.trim().is_empty() {
                    first_non_empty_char = line.trim().chars().next();
                }
                lines.push(line);
            }
            Some(Err(e)) => {
                return Err(phraya_core::ParseError::InvalidUtf8(format!(
                    "IO error reading file: {}",
                    e
                )))
            }
            None => break,
        }
    }

    // Detect format from first non-empty character
    match first_non_empty_char {
        Some('>') => {
            // Parse all sequences upfront and validate
            let mut iter = FastaIterator::new(lines);
            let mut sequences = Vec::new();
            while let Some(result) = iter.next() {
                sequences.push(result?);
            }
            Ok(Box::new(sequences.into_iter()))
        }
        Some('@') => {
            // Parse all sequences upfront and validate
            let mut iter = FastqIterator::new(lines);
            let mut sequences = Vec::new();
            while let Some(result) = iter.next() {
                sequences.push(result?);
            }
            Ok(Box::new(sequences.into_iter()))
        }
        Some(_) => Err(phraya_core::ParseError::InvalidFormat(
            "Cannot detect format: file does not start with '>' (FASTA) or '@' (FASTQ)"
                .to_string(),
        )),
        None => {
            // Empty file - return an iterator that yields nothing
            Ok(Box::new(std::iter::empty::<Sequence>()))
        }
    }
}

struct FastaIterator {
    lines: Vec<String>,
    index: usize,
    current_line: Option<String>,
    error: Option<phraya_core::ParseError>,
}

impl FastaIterator {
    fn new(lines: Vec<String>) -> Self {
        FastaIterator {
            lines,
            index: 0,
            current_line: None,
            error: None,
        }
    }

    fn skip_empty_lines(&mut self) -> Option<String> {
        while self.index < self.lines.len() {
            let line = self.lines[self.index].clone();
            self.index += 1;
            if !line.trim().is_empty() {
                return Some(line);
            }
        }
        None
    }
}

impl Iterator for FastaIterator {
    type Item = Result<Sequence, phraya_core::ParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        // If we previously encountered an error, stop iteration
        if self.error.is_some() {
            return None;
        }

        let header = self.current_line.take().or_else(|| self.skip_empty_lines())?;

        if !header.starts_with('>') {
            let err = phraya_core::ParseError::InvalidFormat(
                "Invalid FASTA header".to_string(),
            );
            self.error = None; // Consumed the error
            return Some(Err(err));
        }

        // Parse header: ">id description"
        let header_content = &header[1..]; // Skip '>'
        let (id, description) = if let Some(space_pos) =
            header_content.find(|c: char| c.is_whitespace())
        {
            let id = header_content[..space_pos].to_string();
            let desc = header_content[space_pos..].trim().to_string();
            (
                id,
                if desc.is_empty() { None } else { Some(desc) },
            )
        } else {
            (header_content.to_string(), None)
        };

        if id.is_empty() {
            let err = phraya_core::ParseError::InvalidFormat(
                "Empty sequence ID".to_string(),
            );
            return Some(Err(err));
        }

        // Read sequence lines until next header or EOF
        let mut bases = Vec::new();
        loop {
            match self.skip_empty_lines() {
                Some(line) => {
                    if line.starts_with('>') {
                        // Next header found, save it for next iteration
                        self.current_line = Some(line);
                        break;
                    } else {
                        // Sequence data
                        let trimmed = line.trim();
                        if validate_dna_bases(trimmed).is_err() {
                            // Invalid characters found
                            let err = phraya_core::ParseError::InvalidFormat(
                                format!("Invalid DNA base in sequence: {}", trimmed),
                            );
                            return Some(Err(err));
                        }
                        bases.extend_from_slice(trimmed.as_bytes());
                    }
                }
                None => break,
            }
        }

        if bases.is_empty() {
            let err = phraya_core::ParseError::InvalidFormat(
                "No sequence data found".to_string(),
            );
            return Some(Err(err));
        }

        Some(Ok(Sequence::new(bases, None, id, description)))
    }
}

struct FastqIterator {
    lines: Vec<String>,
    index: usize,
    error: Option<phraya_core::ParseError>,
}

impl FastqIterator {
    fn new(lines: Vec<String>) -> Self {
        FastqIterator {
            lines,
            index: 0,
            error: None,
        }
    }

    fn skip_empty_lines(&mut self) -> Option<String> {
        while self.index < self.lines.len() {
            let line = self.lines[self.index].clone();
            self.index += 1;
            if !line.trim().is_empty() {
                return Some(line);
            }
        }
        None
    }

    fn next_line(&mut self) -> Option<String> {
        if self.index < self.lines.len() {
            let line = self.lines[self.index].clone();
            self.index += 1;
            Some(line)
        } else {
            None
        }
    }
}

impl Iterator for FastqIterator {
    type Item = Result<Sequence, phraya_core::ParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        // If we previously encountered an error, stop iteration
        if self.error.is_some() {
            return None;
        }

        let header = self.skip_empty_lines()?;

        if !header.starts_with('@') {
            let err = phraya_core::ParseError::InvalidFormat(
                "Invalid FASTQ header".to_string(),
            );
            return Some(Err(err));
        }

        // Parse header: "@id description"
        let header_content = &header[1..]; // Skip '@'
        let (id, description) = if let Some(space_pos) =
            header_content.find(|c: char| c.is_whitespace())
        {
            let id = header_content[..space_pos].to_string();
            let desc = header_content[space_pos..].trim().to_string();
            (
                id,
                if desc.is_empty() { None } else { Some(desc) },
            )
        } else {
            (header_content.to_string(), None)
        };

        if id.is_empty() {
            let err = phraya_core::ParseError::InvalidFormat(
                "Empty sequence ID".to_string(),
            );
            return Some(Err(err));
        }

        // Line 2: sequence
        let sequence = match self.next_line() {
            Some(line) => line,
            None => {
                let err = phraya_core::ParseError::InvalidFormat(
                    "Missing sequence line in FASTQ".to_string(),
                );
                return Some(Err(err));
            }
        };

        let sequence = sequence.trim();
        if validate_dna_bases(sequence).is_err() {
            let err = phraya_core::ParseError::InvalidFormat(
                format!("Invalid DNA base in sequence: {}", sequence),
            );
            return Some(Err(err));
        }

        // Line 3: "+"
        let sep = match self.next_line() {
            Some(line) => line,
            None => {
                let err = phraya_core::ParseError::InvalidFormat(
                    "Missing separator line in FASTQ".to_string(),
                );
                return Some(Err(err));
            }
        };

        if !sep.trim().starts_with('+') {
            let err = phraya_core::ParseError::InvalidFormat(
                "Invalid separator line in FASTQ (expected '+')".to_string(),
            );
            return Some(Err(err));
        }

        // Line 4: quality scores
        let quality = match self.next_line() {
            Some(line) => line,
            None => {
                let err = phraya_core::ParseError::InvalidFormat(
                    "Missing quality line in FASTQ".to_string(),
                );
                return Some(Err(err));
            }
        };

        let quality = quality.trim();

        // Validate quality length matches sequence length
        if quality.len() != sequence.len() {
            let err = phraya_core::ParseError::InvalidFormat(
                format!(
                    "Quality score length ({}) does not match sequence length ({})",
                    quality.len(),
                    sequence.len()
                ),
            );
            return Some(Err(err));
        }

        Some(Ok(Sequence::new(
            sequence.as_bytes().to_vec(),
            Some(quality.as_bytes().to_vec()),
            id,
            description,
        )))
    }
}

fn validate_dna_bases(seq: &str) -> Result<(), ()> {
    for ch in seq.chars() {
        match ch {
            'A' | 'a' | 'C' | 'c' | 'G' | 'g' | 'T' | 't' | 'U' | 'u' | 'N' | 'n' | 'R'
            | 'r' | 'Y' | 'y' | 'S' | 's' | 'W' | 'w' | 'K' | 'k' | 'M' | 'm' | 'B' | 'b'
            | 'D' | 'd' | 'H' | 'h' | 'V' | 'v' => {}
            _ => return Err(()),
        }
    }
    Ok(())
}

fn is_gzip(path: &Path) -> Result<bool, phraya_core::ParseError> {
    // Check file extension first
    if let Some(ext) = path.extension() {
        if let Some(ext_str) = ext.to_str() {
            if ext_str == "gz" {
                return Ok(true);
            }
        }
    }

    // Check magic bytes: 0x1f 0x8b
    let file = File::open(path).map_err(|e| {
        phraya_core::ParseError::InvalidFormat(format!("Failed to check gzip: {}", e))
    })?;

    let mut reader = BufReader::new(file);
    let mut magic_bytes = [0u8; 2];

    match reader.read_exact(&mut magic_bytes) {
        Ok(()) => Ok(magic_bytes == [0x1f, 0x8b]),
        Err(_) => Ok(false), // Not enough bytes = not gzip
    }
}

#[cfg(test)]
mod fasta_fastq_tests;
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn placeholder_now_works() {
        // Just verify the function signature works
        let result = parse_sequences("nonexistent.fa");
        assert!(result.is_err());
    }

    #[test]
    fn test_iterator_yields_sequences() {
        // Verify the iterator yields Sequence items
        use tempfile::NamedTempFile;
        use std::io::Write;

        let mut file = NamedTempFile::new().unwrap();
        file.write_all(b">seq1\nACGT\n").unwrap();
        file.flush().unwrap();

        let result = parse_sequences(file.path());
        assert!(result.is_ok());

        let iter = result.unwrap();
        let items: Vec<_> = iter.collect();

        // Should have one sequence
        assert_eq!(items.len(), 1);
        assert_eq!(items[0].id(), "seq1");
    }
}
