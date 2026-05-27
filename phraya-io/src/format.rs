use std::fs::File;
use std::io::Read;
use std::path::Path;

use crate::ParseError;

/// Represents the format and compression of a sequence file
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SequenceFormat {
    /// FASTA, uncompressed
    FastaPlain,
    /// FASTA, gzip-compressed
    FastaGzip,
    /// FASTA, bzip2-compressed
    FastaBzip2,
    /// FASTQ, uncompressed
    FastqPlain,
    /// FASTQ, gzip-compressed
    FastqGzip,
    /// FASTQ, bzip2-compressed
    FastqBzip2,
}

impl SequenceFormat {
    /// Detect format from magic bytes and file extension
    pub fn detect(path: &Path) -> Result<Self, ParseError> {
        let path_str = path.to_string_lossy().to_string();

        // Try to read magic bytes
        let mut file = File::open(path).map_err(|e| ParseError::IoError {
            path: path_str.clone(),
            source: e,
        })?;

        let mut magic = [0u8; 3];
        let bytes_read = file.read(&mut magic).map_err(|e| ParseError::IoError {
            path: path_str.clone(),
            source: e,
        })?;

        // Check for gzip magic bytes (1f 8b)
        if bytes_read >= 2 && magic[0] == 0x1f && magic[1] == 0x8b {
            let format_type = Self::detect_format_type(path)?;
            return Ok(match format_type {
                FormatType::Fasta => SequenceFormat::FastaGzip,
                FormatType::Fastq => SequenceFormat::FastqGzip,
            });
        }

        // Check for bzip2 magic bytes (42 5a)
        if bytes_read >= 2 && magic[0] == 0x42 && magic[1] == 0x5a {
            let format_type = Self::detect_format_type(path)?;
            return Ok(match format_type {
                FormatType::Fasta => SequenceFormat::FastaBzip2,
                FormatType::Fastq => SequenceFormat::FastqBzip2,
            });
        }

        // Fall back to extension-based detection
        let format_type = Self::detect_format_type(path)?;
        Ok(match format_type {
            FormatType::Fasta => SequenceFormat::FastaPlain,
            FormatType::Fastq => SequenceFormat::FastqPlain,
        })
    }

    /// Detect whether format is FASTA or FASTQ from file extension
    fn detect_format_type(path: &Path) -> Result<FormatType, ParseError> {
        let path_str = path.to_string_lossy().to_string();
        let path_lower = path_str.to_lowercase();

        if path_lower.ends_with(".fasta")
            || path_lower.ends_with(".fasta.gz")
            || path_lower.ends_with(".fasta.bz2")
            || path_lower.ends_with(".fa")
            || path_lower.ends_with(".fa.gz")
            || path_lower.ends_with(".fa.bz2")
        {
            Ok(FormatType::Fasta)
        } else if path_lower.ends_with(".fastq")
            || path_lower.ends_with(".fastq.gz")
            || path_lower.ends_with(".fastq.bz2")
            || path_lower.ends_with(".fq")
            || path_lower.ends_with(".fq.gz")
            || path_lower.ends_with(".fq.bz2")
        {
            Ok(FormatType::Fastq)
        } else {
            Err(ParseError::UnsupportedFormat {
                path: path_str,
                reason: "Unknown file extension. Expected .fasta/.fa or .fastq/.fq with optional .gz/.bz2 compression".to_string(),
            })
        }
    }

    /// Check if this format is FASTA
    pub fn is_fasta(self) -> bool {
        matches!(
            self,
            SequenceFormat::FastaPlain | SequenceFormat::FastaGzip | SequenceFormat::FastaBzip2
        )
    }

    /// Check if this format is FASTQ
    pub fn is_fastq(self) -> bool {
        matches!(
            self,
            SequenceFormat::FastqPlain | SequenceFormat::FastqGzip | SequenceFormat::FastqBzip2
        )
    }

    /// Check if this format is compressed
    pub fn is_compressed(self) -> bool {
        matches!(
            self,
            SequenceFormat::FastaGzip
                | SequenceFormat::FastaBzip2
                | SequenceFormat::FastqGzip
                | SequenceFormat::FastqBzip2
        )
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FormatType {
    Fasta,
    Fastq,
}
