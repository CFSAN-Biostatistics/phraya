use std::io;
use std::path::Path;
use thiserror::Error;

mod format;
mod reader;

pub use format::SequenceFormat;
pub use reader::SequenceReader;

/// Error type for sequence parsing operations
#[derive(Error, Debug)]
pub enum ParseError {
    #[error("I/O error reading {path}: {source}")]
    IoError { path: String, source: io::Error },

    #[error("Malformed entry at line {line}: {reason}")]
    MalformedEntry { line: usize, reason: String },

    #[error("Quality score length mismatch: sequence length {seq_len}, quality length {qual_len}")]
    QualityLengthMismatch { seq_len: usize, qual_len: usize },

    #[error("Unsupported format for {path}: {reason}")]
    UnsupportedFormat { path: String, reason: String },
}

/// Parse sequences from a file with auto-detection of format and compression
///
/// # Arguments
/// * `path` - Path to FASTA or FASTQ file (compressed or plain)
///
/// # Returns
/// * `Result<SequenceReader, ParseError>` - Iterator over sequences or error
pub fn parse_sequences<P: AsRef<Path>>(path: P) -> Result<SequenceReader, ParseError> {
    let path_ref = path.as_ref();

    // Detect format
    let format = SequenceFormat::detect(path_ref)?;

    // Create reader
    SequenceReader::new(path_ref, format)
}
