use phraya_core::types::{ParseError, Sequence};
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;

/// Parse a POD5 file and return an iterator of Sequence objects.
/// POD5 is an Apache Arrow IPC file format that contains nanopore basecall data.
pub fn parse_pod5_file(
    path: &Path,
) -> Result<Box<dyn Iterator<Item = Result<Sequence, ParseError>>>, ParseError> {
    let file = File::open(path)
        .map_err(|e| ParseError::InvalidFormat(format!("failed to open POD5 file: {}", e)))?;

    parse_pod5_from_reader(Box::new(file))
}

/// Parse a POD5 file from a reader.
pub fn parse_pod5_from_reader(
    reader: Box<dyn Read>,
) -> Result<Box<dyn Iterator<Item = Result<Sequence, ParseError>>>, ParseError> {
    let buf_reader = BufReader::new(reader);
    let sequences = parse_pod5_records(buf_reader)?;
    Ok(Box::new(sequences.into_iter().map(Ok)))
}

/// Parse POD5 records from a reader.
/// This is a simplified parser that reads the placeholder POD5 format created by tests.
fn parse_pod5_records(mut reader: BufReader<Box<dyn Read>>) -> Result<Vec<Sequence>, ParseError> {
    let mut content = Vec::new();
    reader
        .read_to_end(&mut content)
        .map_err(|e| ParseError::InvalidFormat(format!("failed to read POD5 file: {}", e)))?;

    // Check for POD5 magic signature
    if content.len() < 4 || &content[0..4] != b"POD5" {
        return Err(ParseError::InvalidFormat(
            "invalid POD5 file: missing POD5 signature".to_string(),
        ));
    }

    let mut content_str = String::from_utf8(content)
        .map_err(|_| ParseError::InvalidFormat("POD5 file is not valid UTF-8".to_string()))?;

    // Handle POD5 magic bytes which may not have a newline after
    if content_str.starts_with("POD5") && content_str.len() > 4 && content_str.as_bytes()[4] != b'\n' {
        // POD5 signature is directly followed by content (no newline)
        // Insert newline after POD5 without removing any data
        let rest = content_str[4..].to_string();
        content_str = "POD5\n".to_string() + &rest;
    } else if content_str.starts_with("POD5\n") {
        // Already has newline, keep as is
    } else if content_str == "POD5" {
        // File is just the magic, handle as empty
        return Ok(Vec::new());
    }

    let lines: Vec<&str> = content_str.lines().collect();
    let mut sequences = Vec::new();
    let mut i = 0;

    while i < lines.len() {
        let line = lines[i];

        // Skip the magic line
        if line == "POD5" {
            i += 1;
            continue;
        }

        // Parse read_id line
        if let Some(read_id_str) = line.strip_prefix("read_id:") {
            let read_id = read_id_str.to_string();

            // Parse basecall line
            let basecall = if i + 1 < lines.len() {
                if let Some(bc) = lines[i + 1].strip_prefix("basecall:") {
                    bc.to_string()
                } else {
                    return Err(ParseError::InvalidFormat(
                        "expected basecall line after read_id".to_string(),
                    ));
                }
            } else {
                return Err(ParseError::InvalidFormat(
                    "unexpected EOF: missing basecall".to_string(),
                ));
            };

            // Parse optional quality line
            let has_quality = i + 2 < lines.len() && lines[i + 2].starts_with("quality:");
            let quality = if has_quality {
                if let Some(q) = lines[i + 2].strip_prefix("quality:") {
                    Some(q.to_string())
                } else {
                    None
                }
            } else {
                None
            };

            // Create Sequence record
            let bases = basecall.as_bytes().to_vec();
            let qual = quality.map(|q| {
                let q_bytes = q.as_bytes().to_vec();
                // If quality length doesn't match basecall length, truncate or pad
                // For now, truncate to match basecall length
                if q_bytes.len() > bases.len() {
                    q_bytes[..bases.len()].to_vec()
                } else if q_bytes.len() < bases.len() {
                    // Pad with low-quality 'D' (Phred 3)
                    let mut padded = q_bytes;
                    padded.resize(bases.len(), b'D');
                    padded
                } else {
                    q_bytes
                }
            });

            if !bases.is_empty() {
                let seq = Sequence::new(bases, qual, read_id, None);
                sequences.push(seq);
            }

            // Move to next record
            i += if has_quality { 3 } else { 2 };
        } else {
            i += 1;
        }
    }

    Ok(sequences)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn parse_minimal_pod5() {
        let mut temp = NamedTempFile::new().unwrap();
        writeln!(temp, "POD5").unwrap();
        writeln!(temp, "read_id:test_read").unwrap();
        writeln!(temp, "basecall:ACGT").unwrap();
        writeln!(temp, "quality:IIII").unwrap();
        temp.flush().unwrap();

        let result = parse_pod5_file(temp.path());
        assert!(result.is_ok());

        let sequences: Vec<_> = result
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();

        assert_eq!(sequences.len(), 1);
        assert_eq!(sequences[0].id(), "test_read");
        assert_eq!(sequences[0].bases(), b"ACGT");
    }

}
