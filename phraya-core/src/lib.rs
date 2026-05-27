//! Core types for Phraya
//!
//! This crate defines foundational types used throughout Phraya:
//! - `VariantObservation`: Variant call with position and evidence
//! - `BaseConfidence`: Confidence scoring metadata
//! - `Sequence`: DNA/RNA sequence with optional quality scores
//!
//! These types are serializable for the .phraya binary format.

use serde::{Serialize, Deserialize};

/// Temporary stub for VariantObservation
///
/// This is a minimal stub to allow tests in phraya-io to compile.
/// Full implementation will be provided by issue #1.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct VariantObservation {
    pub position: u64,
    pub ref_base: u8,
    pub alt_base: u8,
}

impl VariantObservation {
    /// Temporary constructor stub
    pub fn new(position: u64, ref_base: u8, alt_base: u8) -> Self {
        Self {
            position,
            ref_base,
            alt_base,
        }
    }
}
