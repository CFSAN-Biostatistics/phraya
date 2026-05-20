//! Phraya Core Types
//!
//! Defines foundational types for variant evidence, sequence representation, and confidence scoring.

use serde::{Deserialize, Serialize};
use thiserror::Error;

/// Represents a DNA/RNA sequence with optional quality scores.
///
/// A sequence is the fundamental unit of sequence data in Phraya.
/// It can represent assembly contigs, short reads, or long reads,
/// with optional per-base quality scores (typically Phred scores for reads).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct Sequence {
    /// The sequence identifier (e.g., contig name, read ID)
    pub id: String,
    /// The sequence string (ACGT for DNA, ACGU for RNA, N for unknown)
    pub sequence: String,
    /// Optional per-base quality scores (Phred scale, typically 0-60)
    pub quality: Option<Vec<u8>>,
}

impl Sequence {
    /// Create a new sequence without quality scores.
    pub fn new(id: String, sequence: String) -> Self {
        todo!("Sequence::new")
    }

    /// Create a new sequence with quality scores.
    pub fn with_quality(id: String, sequence: String, quality: Vec<u8>) -> Self {
        todo!("Sequence::with_quality")
    }

    /// Return the length of the sequence.
    pub fn len(&self) -> usize {
        todo!("Sequence::len")
    }

    /// Return true if the sequence is empty.
    pub fn is_empty(&self) -> bool {
        todo!("Sequence::is_empty")
    }
}

/// Represents a per-base confidence score (0.0 to 1.0).
///
/// Confidence scores integrate multiple evidence sources (alignment quality,
/// population support, repeat context, etc.) into a single reliability metric.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub struct BaseConfidence(f64);

impl BaseConfidence {
    /// Create a new confidence score with validation (0.0 <= value <= 1.0).
    pub fn new(value: f64) -> Result<Self, String> {
        todo!("BaseConfidence::new")
    }

    /// Create a confidence score of 0.0 (no confidence).
    pub fn zero() -> Self {
        todo!("BaseConfidence::zero")
    }

    /// Create a confidence score of 1.0 (full confidence).
    pub fn one() -> Self {
        todo!("BaseConfidence::one")
    }

    /// Retrieve the confidence value (0.0 to 1.0).
    pub fn value(&self) -> f64 {
        todo!("BaseConfidence::value")
    }
}

/// Represents a variant observation at a specific genomic position.
///
/// Captures evidence from a single alignment: the position, reference/alternate bases,
/// depth of coverage, and confidence score informed by alignment context and population support.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct VariantObservation {
    /// Position in the reference sequence (0-based, inclusive)
    pub position: usize,
    /// Reference base at this position
    pub reference_base: char,
    /// Alternate base (SNP), or '-' for indel anchor
    pub alternate_base: char,
    /// Total depth of coverage at this position
    pub depth: u32,
    /// Confidence in this observation
    pub confidence: BaseConfidence,
    /// Source identifier (e.g., "align_000" or sample name), optional
    pub source: Option<String>,
}

impl VariantObservation {
    /// Create a new variant observation.
    pub fn new(
        position: usize,
        reference_base: char,
        alternate_base: char,
        depth: u32,
        confidence: BaseConfidence,
    ) -> Self {
        todo!("VariantObservation::new")
    }

    /// Create a variant observation with a source identifier.
    pub fn with_source(
        position: usize,
        reference_base: char,
        alternate_base: char,
        depth: u32,
        confidence: BaseConfidence,
        source: String,
    ) -> Self {
        todo!("VariantObservation::with_source")
    }

    /// Calculate the allele frequency (depth) associated with this variant.
    pub fn allele_frequency(&self) -> f64 {
        todo!("VariantObservation::allele_frequency")
    }
}

/// Errors that occur during alignment operations.
#[derive(Error, Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum AlignError {
    /// Alignment computation failed
    #[error("alignment failed: {0}")]
    AlignmentFailed(String),

    /// Invalid input sequence (e.g., invalid bases, mismatched quality length)
    #[error("invalid sequence: {0}")]
    InvalidSequence(String),

    /// Insufficient coverage to produce confident results
    #[error("insufficient coverage: {0}")]
    InsufficientCoverage(String),

    /// Invalid parameters (e.g., negative gap penalty)
    #[error("invalid parameters: {0}")]
    InvalidParameters(String),
}

/// Errors that occur during index operations.
#[derive(Error, Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum IndexError {
    /// Index construction failed
    #[error("index building failed: {0}")]
    BuildFailed(String),

    /// Index not found or could not be loaded
    #[error("index not found: {0}")]
    NotFound(String),

    /// Index file format is invalid
    #[error("invalid index format: {0}")]
    InvalidFormat(String),

    /// Query not found in index
    #[error("query not found in index: {0}")]
    QueryNotFound(String),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[should_panic]
    fn test_sequence_new() {
        let _ = Sequence::new("id".to_string(), "ACGT".to_string());
    }

    #[test]
    #[should_panic]
    fn test_sequence_with_quality() {
        let _ = Sequence::with_quality("id".to_string(), "ACGT".to_string(), vec![30, 31, 32, 33]);
    }

    #[test]
    #[should_panic]
    fn test_sequence_len() {
        let seq = Sequence {
            id: "id".to_string(),
            sequence: "ACGT".to_string(),
            quality: None,
        };
        let _ = seq.len();
    }

    #[test]
    #[should_panic]
    fn test_sequence_is_empty() {
        let seq = Sequence {
            id: "id".to_string(),
            sequence: "ACGT".to_string(),
            quality: None,
        };
        let _ = seq.is_empty();
    }

    #[test]
    #[should_panic]
    fn test_base_confidence_new() {
        let _ = BaseConfidence::new(0.5);
    }

    #[test]
    #[should_panic]
    fn test_base_confidence_zero() {
        let _ = BaseConfidence::zero();
    }

    #[test]
    #[should_panic]
    fn test_base_confidence_one() {
        let _ = BaseConfidence::one();
    }

    #[test]
    #[should_panic]
    fn test_base_confidence_value() {
        let conf = BaseConfidence(0.5);
        let _ = conf.value();
    }

    #[test]
    #[should_panic]
    fn test_variant_observation_new() {
        let _ = VariantObservation::new(100, 'A', 'T', 50, BaseConfidence(0.9));
    }

    #[test]
    #[should_panic]
    fn test_variant_observation_with_source() {
        let _ = VariantObservation::with_source(100, 'A', 'T', 50, BaseConfidence(0.9), "src".to_string());
    }

    #[test]
    #[should_panic]
    fn test_variant_observation_allele_frequency() {
        let obs = VariantObservation {
            position: 100,
            reference_base: 'A',
            alternate_base: 'T',
            depth: 50,
            confidence: BaseConfidence(0.9),
            source: None,
        };
        let _ = obs.allele_frequency();
    }
}
