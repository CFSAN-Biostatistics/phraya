//! I/O operations for Phraya
//!
//! This crate handles reading and writing various genomic file formats:
//! - FASTQ/FASTA input (planned)
//! - VCF output (planned)
//! - BAM/CRAM output (planned)
//! - `.phraya` binary format (MessagePack-based native format)

mod phraya_format;

pub use phraya_format::{PhrayaMetadata, read_phraya, read_phraya_metadata, write_phraya};

use thiserror::Error;

#[derive(Debug, Error)]
pub enum IoError {
    #[error("Failed to write file: {0}")]
    WriteError(String),

    #[error("Failed to read file: {0}")]
    ReadError(String),

    #[error("Serialization error: {0}")]
    SerializationError(String),

    #[error("Deserialization error: {0}")]
    DeserializationError(String),

    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
}

pub type Result<T> = std::result::Result<T, IoError>;
