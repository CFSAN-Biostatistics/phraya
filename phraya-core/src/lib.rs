use serde::{Deserialize, Serialize};
use std::collections::HashMap;

mod evidence;

/// Core sequence type for evidence extraction
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Sequence {
    pub id: String,
    pub bases: Vec<u8>,
    pub quality: Option<Vec<u8>>,
}

/// Evidence layer summarizing cross-sample patterns
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct EvidenceLayer {
    /// K-mer uniqueness scores per position (from FM-index)
    pub kmer_uniqueness: HashMap<usize, f64>,

    /// Positions where all samples match reference
    pub invariant_positions: Vec<usize>,

    /// Positions with variation across samples
    pub polymorphic_sites: HashMap<usize, PolymorphicSite>,

    /// Number of samples analyzed
    pub sample_count: usize,

    /// Length of reference sequence
    pub reference_length: usize,
}

/// Information about a polymorphic site
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct PolymorphicSite {
    /// Reference base at this position (if reference mode)
    pub reference_base: u8,

    /// Count of each observed allele
    pub allele_counts: HashMap<u8, usize>,
}

/// Extract evidence layer from sequences
///
/// Scans all input sequences and builds a compact summary of cross-sample patterns.
/// For Phase 1, computes:
/// - Reference k-mer uniqueness map (from FM-index)
/// - Invariant positions (all samples match reference)
/// - Basic polymorphic site catalog (positions that vary, with allele frequencies)
///
/// # Arguments
/// * `sequences` - Input sequences to analyze
/// * `reference` - Optional reference sequence (None for MSA mode)
///
/// # Returns
/// Compact evidence layer suitable for parallel alignment
pub fn extract_evidence(_sequences: &[Sequence], _reference: Option<&Sequence>) -> EvidenceLayer {
    // Placeholder that will be implemented
    unimplemented!("extract_evidence not yet implemented")
}
