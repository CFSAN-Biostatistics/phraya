use crate::SequenceParser;
use std::path::Path;
use thiserror::Error;

/// Detected input classification
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputType {
    /// Short sequence (<5kb, likely from sequencing reads)
    Read,
    /// Long sequence (≥5kb, likely an assembly contig)
    Contig,
}

/// Errors during input classification
#[derive(Debug, Clone, Error)]
pub enum UseCaseError {
    /// I/O errors during file parsing
    #[error("io error: {0}")]
    IoError(String),

    /// No sequences found in input
    #[error("no sequences found in input")]
    NoSequences,
}

/// Detect input type (read or contig) from the first sequence in a file.
///
/// Returns InputType::Read if first sequence is <5kb, InputType::Contig if ≥5kb.
pub fn classify_input(path: &Path) -> Result<InputType, UseCaseError> {
    let mut parser =
        SequenceParser::from_path(path).map_err(|e| UseCaseError::IoError(e.to_string()))?;

    match parser.next() {
        Some(Ok(seq)) => {
            if seq.bases().len() >= 5000 {
                Ok(InputType::Contig)
            } else {
                Ok(InputType::Read)
            }
        }
        Some(Err(e)) => Err(UseCaseError::IoError(e.to_string())),
        None => Err(UseCaseError::NoSequences),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_fasta(content: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        write!(f, "{}", content).unwrap();
        f.flush().unwrap();
        f
    }

    #[test]
    fn classify_short_sequence_as_read() {
        let f = write_fasta(">read1\nACGT\n");
        assert_eq!(classify_input(f.path()).unwrap(), InputType::Read);
    }

    #[test]
    fn classify_long_sequence_as_contig() {
        let fasta = format!(">contig1\n{}\n", "A".repeat(5000));
        let f = write_fasta(&fasta);
        assert_eq!(classify_input(f.path()).unwrap(), InputType::Contig);
    }

    #[test]
    fn classify_exactly_5kb_as_contig() {
        let fasta = format!(">boundary\n{}\n", "A".repeat(5000));
        let f = write_fasta(&fasta);
        assert_eq!(classify_input(f.path()).unwrap(), InputType::Contig);
    }

    #[test]
    fn classify_4999bp_as_read() {
        let fasta = format!(">just_under\n{}\n", "A".repeat(4999));
        let f = write_fasta(&fasta);
        assert_eq!(classify_input(f.path()).unwrap(), InputType::Read);
    }

    #[test]
    fn classify_empty_file_returns_no_sequences() {
        let f = NamedTempFile::new().unwrap();
        assert!(matches!(
            classify_input(f.path()),
            Err(UseCaseError::NoSequences)
        ));
    }

    #[test]
    fn classify_nonexistent_file_returns_io_error() {
        use std::path::Path;
        assert!(matches!(
            classify_input(Path::new("/nonexistent/path/file.fa")),
            Err(UseCaseError::IoError(_))
        ));
    }

    #[test]
    fn classify_invalid_utf8_body_returns_io_error() {
        // Header line is valid so the format detects as FASTA, but the sequence
        // body line contains invalid UTF-8, which BufRead::lines() surfaces as
        // an io::Error when the iterator advances past the header.
        let mut f = NamedTempFile::new().unwrap();
        f.write_all(b">seq1\n").unwrap();
        f.write_all(&[b'A', 0xFF, 0xFE, b'\n']).unwrap();
        f.flush().unwrap();

        assert!(matches!(
            classify_input(f.path()),
            Err(UseCaseError::IoError(_))
        ));
    }
}
