//! Phraya Core Types
//!
//! Defines foundational types for variant evidence, sequence representation, and confidence scoring.

use serde::{Deserialize, Serialize};
use thiserror::Error;

/// Errors that occur during alignment operations.
#[derive(Error, Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum AlignError {
    #[error("alignment failed: {0}")]
    AlignmentFailed(String),

    #[error("invalid sequence: {0}")]
    InvalidSequence(String),

    #[error("insufficient coverage: {0}")]
    InsufficientCoverage(String),

    #[error("invalid parameters: {0}")]
    InvalidParameters(String),
}

/// Errors that occur during index operations.
#[derive(Error, Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum IndexError {
    #[error("index building failed: {0}")]
    BuildFailed(String),

    #[error("index not found: {0}")]
    NotFound(String),

    #[error("invalid index format: {0}")]
    InvalidFormat(String),

    #[error("query not found in index: {0}")]
    QueryNotFound(String),
}

/// Represents a DNA/RNA sequence with optional quality scores.
///
/// A sequence is the fundamental unit of sequence data in Phraya.
/// It can represent assembly contigs, short reads, or long reads,
/// with optional per-base quality scores (typically Phred scores for reads).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct Sequence {
    /// The sequence string (ACGT for DNA, ACGU for RNA, N for unknown)
    pub sequence: String,

    /// Optional per-base quality scores (Phred scale)
    /// If present, length must match sequence length
    pub quality: Option<Vec<u8>>,

    /// Sequence identifier (contig name, read name, etc.)
    pub id: String,
}

impl Sequence {
    /// Creates a new sequence without quality scores.
    pub fn new(id: String, sequence: String) -> Self {
        Sequence {
            sequence,
            quality: None,
            id,
        }
    }

    /// Creates a new sequence with quality scores.
    pub fn with_quality(id: String, sequence: String, quality: Vec<u8>) -> Self {
        Sequence {
            sequence,
            quality: Some(quality),
            id,
        }
    }

    /// Returns the length of the sequence.
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// Returns true if the sequence is empty.
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }
}

/// Represents a confidence score for a base call.
///
/// Confidence scores quantify how certain we are about a base call
/// or variant observation. Scores range from 0.0 (no confidence)
/// to 1.0 (complete confidence).
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, PartialOrd)]
pub struct BaseConfidence(f64);

impl BaseConfidence {
    /// Creates a new confidence score from a value between 0.0 and 1.0.
    /// Returns None if the value is outside the valid range.
    pub fn new(value: f64) -> Option<Self> {
        if value >= 0.0 && value <= 1.0 && !value.is_nan() {
            Some(BaseConfidence(value))
        } else {
            None
        }
    }

    /// Returns the confidence value as an f64.
    pub fn value(&self) -> f64 {
        self.0
    }

    /// Creates a confidence of 0.0 (no confidence).
    pub fn zero() -> Self {
        BaseConfidence(0.0)
    }

    /// Creates a confidence of 1.0 (complete confidence).
    pub fn one() -> Self {
        BaseConfidence(1.0)
    }
}

/// Represents an observation of a variant at a specific position.
///
/// A variant observation captures evidence that a specific alternate base
/// occurs at a particular position in a sequence. Multiple observations
/// at the same position may come from different reads, assemblies, or regions.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct VariantObservation {
    /// Reference position (0-based, inclusive)
    pub position: u64,

    /// Reference base at this position
    pub reference_base: char,

    /// Alternate (observed) base
    pub alternate_base: char,

    /// Count of observations supporting this variant
    pub observation_count: u32,

    /// Count of observations supporting the reference
    pub reference_count: u32,

    /// Optional source identifier (e.g., sample name, read name, contig name)
    pub source: Option<String>,
}

impl VariantObservation {
    /// Creates a new variant observation.
    pub fn new(
        position: u64,
        reference_base: char,
        alternate_base: char,
        observation_count: u32,
        reference_count: u32,
    ) -> Self {
        VariantObservation {
            position,
            reference_base,
            alternate_base,
            observation_count,
            reference_count,
            source: None,
        }
    }

    /// Creates a new variant observation with a source identifier.
    pub fn with_source(
        position: u64,
        reference_base: char,
        alternate_base: char,
        observation_count: u32,
        reference_count: u32,
        source: String,
    ) -> Self {
        VariantObservation {
            position,
            reference_base,
            alternate_base,
            observation_count,
            reference_count,
            source: Some(source),
        }
    }

    /// Returns the total coverage (observations + reference_count).
    pub fn total_depth(&self) -> u64 {
        (self.observation_count as u64) + (self.reference_count as u64)
    }

    /// Returns the allele frequency as a fraction (0.0 to 1.0).
    pub fn allele_frequency(&self) -> f64 {
        let total = self.total_depth();
        if total == 0 {
            0.0
        } else {
            self.observation_count as f64 / total as f64
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ===== Sequence Tests =====

    #[test]
    fn test_sequence_new_without_quality() {
        let seq = Sequence::new("read001".to_string(), "ACGTACGT".to_string());
        assert_eq!(seq.id, "read001");
        assert_eq!(seq.sequence, "ACGTACGT");
        assert_eq!(seq.quality, None);
        assert_eq!(seq.len(), 8);
        assert!(!seq.is_empty());
    }

    #[test]
    fn test_sequence_with_quality() {
        let quality = vec![30, 31, 32, 33, 30, 31, 32, 33];
        let seq = Sequence::with_quality(
            "read001".to_string(),
            "ACGTACGT".to_string(),
            quality.clone(),
        );
        assert_eq!(seq.id, "read001");
        assert_eq!(seq.sequence, "ACGTACGT");
        assert_eq!(seq.quality, Some(quality));
        assert_eq!(seq.len(), 8);
    }

    #[test]
    fn test_sequence_empty() {
        let seq = Sequence::new("empty".to_string(), String::new());
        assert!(seq.is_empty());
        assert_eq!(seq.len(), 0);
    }

    #[test]
    fn test_sequence_serialize_deserialize() {
        let seq = Sequence::with_quality(
            "test_seq".to_string(),
            "ACGTNN".to_string(),
            vec![30, 31, 32, 33, 0, 0],
        );

        let json = serde_json::to_string(&seq).expect("serialization failed");
        let deserialized: Sequence = serde_json::from_str(&json).expect("deserialization failed");

        assert_eq!(seq, deserialized);
    }

    #[test]
    fn test_sequence_serialize_without_quality() {
        let seq = Sequence::new("test".to_string(), "ACGT".to_string());
        let json = serde_json::to_string(&seq).expect("serialization failed");
        let deserialized: Sequence = serde_json::from_str(&json).expect("deserialization failed");

        assert_eq!(seq, deserialized);
        assert_eq!(deserialized.quality, None);
    }

    // ===== BaseConfidence Tests =====

    #[test]
    fn test_base_confidence_valid_range() {
        assert_eq!(BaseConfidence::new(0.0).unwrap().value(), 0.0);
        assert_eq!(BaseConfidence::new(0.5).unwrap().value(), 0.5);
        assert_eq!(BaseConfidence::new(1.0).unwrap().value(), 1.0);
    }

    #[test]
    fn test_base_confidence_zero() {
        let conf = BaseConfidence::zero();
        assert_eq!(conf.value(), 0.0);
    }

    #[test]
    fn test_base_confidence_one() {
        let conf = BaseConfidence::one();
        assert_eq!(conf.value(), 1.0);
    }

    #[test]
    fn test_base_confidence_invalid_negative() {
        assert!(BaseConfidence::new(-0.1).is_none());
    }

    #[test]
    fn test_base_confidence_invalid_above_one() {
        assert!(BaseConfidence::new(1.1).is_none());
    }

    #[test]
    fn test_base_confidence_invalid_nan() {
        assert!(BaseConfidence::new(f64::NAN).is_none());
    }

    #[test]
    fn test_base_confidence_serialize_deserialize() {
        let conf = BaseConfidence::new(0.75).unwrap();
        let json = serde_json::to_string(&conf).expect("serialization failed");
        let deserialized: BaseConfidence = serde_json::from_str(&json).expect("deserialization failed");

        assert_eq!(conf, deserialized);
        assert_eq!(deserialized.value(), 0.75);
    }

    #[test]
    fn test_base_confidence_min_value() {
        let conf = BaseConfidence::new(0.0).unwrap();
        assert_eq!(conf.value(), 0.0);
    }

    #[test]
    fn test_base_confidence_max_value() {
        let conf = BaseConfidence::new(1.0).unwrap();
        assert_eq!(conf.value(), 1.0);
    }

    #[test]
    fn test_base_confidence_comparison() {
        let low = BaseConfidence::new(0.3).unwrap();
        let high = BaseConfidence::new(0.7).unwrap();

        assert!(low < high);
        assert!(high > low);
        assert!(low <= high);
        assert!(high >= low);
    }

    // ===== VariantObservation Tests =====

    #[test]
    fn test_variant_observation_new() {
        let obs = VariantObservation::new(100, 'A', 'T', 5, 10);
        assert_eq!(obs.position, 100);
        assert_eq!(obs.reference_base, 'A');
        assert_eq!(obs.alternate_base, 'T');
        assert_eq!(obs.observation_count, 5);
        assert_eq!(obs.reference_count, 10);
        assert_eq!(obs.source, None);
    }

    #[test]
    fn test_variant_observation_with_source() {
        let obs = VariantObservation::with_source(
            100,
            'A',
            'T',
            5,
            10,
            "sample001".to_string(),
        );
        assert_eq!(obs.source, Some("sample001".to_string()));
    }

    #[test]
    fn test_variant_observation_total_depth() {
        let obs = VariantObservation::new(100, 'A', 'T', 5, 10);
        assert_eq!(obs.total_depth(), 15);
    }

    #[test]
    fn test_variant_observation_total_depth_zero() {
        let obs = VariantObservation::new(100, 'A', 'T', 0, 0);
        assert_eq!(obs.total_depth(), 0);
    }

    #[test]
    fn test_variant_observation_allele_frequency() {
        let obs = VariantObservation::new(100, 'A', 'T', 5, 10);
        assert!((obs.allele_frequency() - 0.333_333_33).abs() < 0.001);
    }

    #[test]
    fn test_variant_observation_allele_frequency_fixed_ratio() {
        let obs = VariantObservation::new(100, 'A', 'T', 1, 1);
        assert_eq!(obs.allele_frequency(), 0.5);
    }

    #[test]
    fn test_variant_observation_allele_frequency_all_alternate() {
        let obs = VariantObservation::new(100, 'A', 'T', 10, 0);
        assert_eq!(obs.allele_frequency(), 1.0);
    }

    #[test]
    fn test_variant_observation_allele_frequency_all_reference() {
        let obs = VariantObservation::new(100, 'A', 'T', 0, 10);
        assert_eq!(obs.allele_frequency(), 0.0);
    }

    #[test]
    fn test_variant_observation_allele_frequency_zero_depth() {
        let obs = VariantObservation::new(100, 'A', 'T', 0, 0);
        assert_eq!(obs.allele_frequency(), 0.0);
    }

    #[test]
    fn test_variant_observation_serialize_deserialize() {
        let obs = VariantObservation::with_source(
            100,
            'A',
            'T',
            5,
            10,
            "sample001".to_string(),
        );

        let json = serde_json::to_string(&obs).expect("serialization failed");
        let deserialized: VariantObservation = serde_json::from_str(&json).expect("deserialization failed");

        assert_eq!(obs, deserialized);
    }

    #[test]
    fn test_variant_observation_serialize_without_source() {
        let obs = VariantObservation::new(100, 'A', 'T', 5, 10);
        let json = serde_json::to_string(&obs).expect("serialization failed");
        let deserialized: VariantObservation = serde_json::from_str(&json).expect("deserialization failed");

        assert_eq!(obs, deserialized);
        assert_eq!(deserialized.source, None);
    }

    // ===== AlignError Tests =====

    #[test]
    fn test_align_error_display() {
        let err = AlignError::AlignmentFailed("test error".to_string());
        assert!(err.to_string().contains("alignment failed"));
        assert!(err.to_string().contains("test error"));
    }

    #[test]
    fn test_align_error_invalid_sequence() {
        let err = AlignError::InvalidSequence("contains invalid chars".to_string());
        assert!(err.to_string().contains("invalid sequence"));
    }

    #[test]
    fn test_align_error_insufficient_coverage() {
        let err = AlignError::InsufficientCoverage("depth < 5".to_string());
        assert!(err.to_string().contains("insufficient coverage"));
    }

    #[test]
    fn test_align_error_invalid_parameters() {
        let err = AlignError::InvalidParameters("bandwidth negative".to_string());
        assert!(err.to_string().contains("invalid parameters"));
    }

    #[test]
    fn test_align_error_serialize_deserialize() {
        let err = AlignError::AlignmentFailed("test error".to_string());
        let json = serde_json::to_string(&err).expect("serialization failed");
        let deserialized: AlignError = serde_json::from_str(&json).expect("deserialization failed");

        assert_eq!(err, deserialized);
    }

    #[test]
    fn test_align_error_clone() {
        let err = AlignError::InvalidSequence("original".to_string());
        let cloned = err.clone();
        assert_eq!(err, cloned);
    }

    // ===== IndexError Tests =====

    #[test]
    fn test_index_error_display() {
        let err = IndexError::BuildFailed("out of memory".to_string());
        assert!(err.to_string().contains("index building failed"));
    }

    #[test]
    fn test_index_error_not_found() {
        let err = IndexError::NotFound("no such index".to_string());
        assert!(err.to_string().contains("index not found"));
    }

    #[test]
    fn test_index_error_invalid_format() {
        let err = IndexError::InvalidFormat("wrong magic bytes".to_string());
        assert!(err.to_string().contains("invalid index format"));
    }

    #[test]
    fn test_index_error_query_not_found() {
        let err = IndexError::QueryNotFound("ACGTACGT".to_string());
        assert!(err.to_string().contains("query not found"));
    }

    #[test]
    fn test_index_error_serialize_deserialize() {
        let err = IndexError::BuildFailed("out of memory".to_string());
        let json = serde_json::to_string(&err).expect("serialization failed");
        let deserialized: IndexError = serde_json::from_str(&json).expect("deserialization failed");

        assert_eq!(err, deserialized);
    }

    #[test]
    fn test_index_error_clone() {
        let err = IndexError::NotFound("original".to_string());
        let cloned = err.clone();
        assert_eq!(err, cloned);
    }

    // ===== Edge Cases and Integration Tests =====

    #[test]
    fn test_sequence_round_trip_all_bases() {
        let seq = Sequence::new("all_bases".to_string(), "ACGTNNACGT".to_string());
        let json = serde_json::to_string(&seq).expect("serialization failed");
        let deserialized: Sequence = serde_json::from_str(&json).expect("deserialization failed");
        assert_eq!(seq, deserialized);
    }

    #[test]
    fn test_sequence_round_trip_max_quality() {
        // Typical Phred max is 93 (93rd percentile confidence)
        let quality = vec![60; 100];
        let seq = Sequence::with_quality(
            "high_quality".to_string(),
            "A".repeat(100),
            quality,
        );
        let json = serde_json::to_string(&seq).expect("serialization failed");
        let deserialized: Sequence = serde_json::from_str(&json).expect("deserialization failed");
        assert_eq!(seq, deserialized);
    }

    #[test]
    fn test_base_confidence_serialization_precision() {
        let original = BaseConfidence::new(0.999_999_9).unwrap();
        let json = serde_json::to_string(&original).expect("serialization failed");
        let deserialized: BaseConfidence = serde_json::from_str(&json).expect("deserialization failed");

        // Allow for small floating point differences
        assert!((original.value() - deserialized.value()).abs() < 1e-6);
    }

    #[test]
    fn test_variant_observation_large_depth() {
        let obs = VariantObservation::new(100, 'A', 'T', 1_000_000, 9_000_000);
        assert_eq!(obs.total_depth(), 10_000_000);
        assert!((obs.allele_frequency() - 0.1).abs() < 0.000_001);
    }

    #[test]
    fn test_variant_observation_multiple_sources_roundtrip() {
        let sources = vec!["sample_1", "sample_2", "contig_A"];
        for source in sources {
            let obs = VariantObservation::with_source(
                42,
                'C',
                'G',
                3,
                7,
                source.to_string(),
            );
            let json = serde_json::to_string(&obs).expect("serialization failed");
            let deserialized: VariantObservation = serde_json::from_str(&json).expect("deserialization failed");
            assert_eq!(obs, deserialized);
            assert_eq!(deserialized.source.as_deref(), Some(source));
        }
    }
}
