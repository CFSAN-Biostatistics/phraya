// This module will contain core types for Phraya.
// Tests are written first (TDD RED phase).

use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use thiserror::Error;

/// Variant type classification for variant observations.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
pub enum VariantType {
    /// Single nucleotide polymorphism (SNP)
    #[serde(rename = "snp")]
    Snp,
    /// Insertion (query has bases, reference does not)
    #[serde(rename = "insertion")]
    Insertion,
    /// Deletion (reference has bases, query does not)
    #[serde(rename = "deletion")]
    Deletion,
}

impl Default for VariantType {
    fn default() -> Self {
        VariantType::Snp
    }
}

/// Sequence type with DNA bytes, optional per-base quality scores (Phred), metadata (id, description).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Sequence {
    /// Raw DNA bases as bytes (e.g., b"ACGT")
    bases: Vec<u8>,
    /// Optional per-base quality scores (Phred format)
    quality_scores: Option<Vec<u8>>,
    /// Sequence identifier (e.g., "seq1")
    id: String,
    /// Optional description
    description: Option<String>,
    /// Mapping quality from BAM/CRAM record (0–254); None for FASTA/FASTQ input
    #[serde(default)]
    mapq: Option<u8>,
}

impl Sequence {
    /// Create a new sequence. Panics if quality_scores.len() != bases.len() when Some.
    pub fn new(
        bases: Vec<u8>,
        quality_scores: Option<Vec<u8>>,
        id: String,
        description: Option<String>,
    ) -> Self {
        if let Some(ref scores) = quality_scores {
            assert_eq!(
                scores.len(),
                bases.len(),
                "quality scores length must match sequence length"
            );
        }
        Sequence {
            bases,
            quality_scores,
            id,
            description,
            mapq: None,
        }
    }

    /// Set the mapping quality (from BAM/CRAM records).
    pub fn with_mapq(mut self, mapq: u8) -> Self {
        self.mapq = Some(mapq);
        self
    }

    /// Return the mapping quality if available.
    pub fn mapq(&self) -> Option<u8> {
        self.mapq
    }

    /// Return the length of the sequence
    pub fn len(&self) -> usize {
        self.bases.len()
    }

    /// Check if the sequence is empty
    pub fn is_empty(&self) -> bool {
        self.bases.is_empty()
    }

    /// Get the sequence ID
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Get the sequence description if present
    pub fn description(&self) -> Option<&str> {
        self.description.as_deref()
    }

    /// Get the raw DNA bases
    pub fn bases(&self) -> &[u8] {
        &self.bases
    }

    /// Get the quality score at a specific position (0-indexed), or None if out of bounds
    pub fn quality_at(&self, pos: usize) -> Option<u8> {
        self.quality_scores
            .as_ref()
            .and_then(|q| q.get(pos).copied())
    }

    /// Get all quality scores if present
    pub fn quality_scores(&self) -> Option<&Vec<u8>> {
        self.quality_scores.as_ref()
    }

    /// Calculate the average quality score, or None if no quality scores
    pub fn avg_quality(&self) -> Option<f64> {
        self.quality_scores.as_ref().map(|scores| {
            if scores.is_empty() {
                0.0
            } else {
                let sum: u32 = scores.iter().map(|&q| q as u32).sum();
                sum as f64 / scores.len() as f64
            }
        })
    }
}

/// K-mer uniqueness default value for VariantObservation
fn default_kmer_uniqueness() -> f64 {
    1.0
}

/// Variant observation at a genomic position with full alignment metadata.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct VariantObservation {
    /// Position in the reference (0-indexed)
    position: u32,
    /// Reference base
    ref_base: u8,
    /// All alleles with their counts
    all_alleles: HashMap<u8, u32>,
    /// Confidence score (0.0-1.0)
    confidence: f64,
    /// CIGAR string for this alignment
    cigar: String,
    /// Mapping quality (0-60)
    mapq: u8,
    /// Edit distance between query and reference
    edit_distance: u32,
    /// Local coverage (±50bp window) as a vector of counts
    local_coverage: Vec<u32>,
    /// Average per-base quality at this position
    avg_base_quality: f64,
    /// Provenance: identifies the sample and read (e.g., "sample1:read42")
    provenance: String,
    /// Whether this variant falls within a tandem repeat region
    #[serde(default)]
    in_tandem_repeat: bool,
    /// Variant type (SNP, insertion, or deletion)
    #[serde(default)]
    variant_type: VariantType,
    /// K-mer uniqueness score at this position (0.0-1.0)
    #[serde(default = "default_kmer_uniqueness")]
    kmer_uniqueness: f64,
    /// Number of reads at this position that are paired (SAM flag 0x1)
    #[serde(default)]
    total_paired_count: u32,
    /// Number of reads at this position that are properly paired (SAM flag 0x2)
    #[serde(default)]
    proper_pair_count: u32,
    /// Sum of absolute insert sizes from paired reads at this position (for mean computation).
    #[serde(default)]
    insert_size_sum: i64,
    /// Number of paired reads contributing an insert size at this position.
    #[serde(default)]
    insert_size_count: u32,
    /// Number of paired reads at this position whose mate was unmapped (SAM flag 0x8).
    #[serde(default)]
    unmapped_mate_count: u32,
    /// Mate relationship metadata for paired-end reads (insert size filters).
    /// Kept as the final field so that omitting it when `None` (the common case) only
    /// shortens the trailing end of the MessagePack array, leaving all earlier fields
    /// correctly positioned; `#[serde(default)]` restores it to `None` on read.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    mate_info: Option<MateInfo>,
}

impl VariantObservation {
    /// Create a new variant observation
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        position: u32,
        ref_base: u8,
        all_alleles: HashMap<u8, u32>,
        confidence: f64,
        cigar: String,
        mapq: u8,
        edit_distance: u32,
        local_coverage: Vec<u32>,
        avg_base_quality: f64,
        provenance: String,
    ) -> Self {
        VariantObservation {
            position,
            ref_base,
            all_alleles,
            confidence,
            cigar,
            mapq,
            edit_distance,
            local_coverage,
            avg_base_quality,
            provenance,
            in_tandem_repeat: false,
            variant_type: VariantType::default(),
            kmer_uniqueness: default_kmer_uniqueness(),
            mate_info: None,
            total_paired_count: 0,
            proper_pair_count: 0,
            insert_size_sum: 0,
            insert_size_count: 0,
            unmapped_mate_count: 0,
        }
    }

    /// Mark whether this variant falls in a tandem repeat region.
    pub fn with_tandem_repeat(mut self, value: bool) -> Self {
        self.in_tandem_repeat = value;
        self
    }

    /// Whether this variant is in a tandem repeat region.
    pub fn in_tandem_repeat(&self) -> bool {
        self.in_tandem_repeat
    }

    /// Set the variant type.
    pub fn with_variant_type(mut self, variant_type: VariantType) -> Self {
        self.variant_type = variant_type;
        self
    }

    /// Get the variant type.
    pub fn variant_type(&self) -> VariantType {
        self.variant_type
    }

    /// Set the k-mer uniqueness score.
    pub fn with_kmer_uniqueness(mut self, value: f64) -> Self {
        self.kmer_uniqueness = value;
        self
    }

    /// Get the k-mer uniqueness score.
    pub fn kmer_uniqueness(&self) -> f64 {
        self.kmer_uniqueness
    }

    /// Set mate information for paired-end reads.
    pub fn with_mate_info(mut self, mate_info: MateInfo) -> Self {
        self.mate_info = Some(mate_info);
        self
    }

    /// Get mate information.
    pub fn mate_info(&self) -> Option<&MateInfo> {
        self.mate_info.as_ref()
    }

    /// Set aggregate pair counts (reads at this position that are paired / properly paired).
    pub fn with_pair_counts(mut self, total_paired: u32, proper_paired: u32) -> Self {
        self.total_paired_count = total_paired;
        self.proper_pair_count = proper_paired;
        self
    }

    /// Raw pair counts: (total_paired, proper_paired).
    pub fn pair_counts(&self) -> (u32, u32) {
        (self.total_paired_count, self.proper_pair_count)
    }

    /// Add pair counts from another observation at the same position (used during merge).
    pub fn add_pair_counts(&mut self, total_paired: u32, proper_paired: u32) {
        self.total_paired_count += total_paired;
        self.proper_pair_count += proper_paired;
    }

    /// Fraction of paired reads at this position that are properly paired (0.0–1.0).
    /// Returns None if no paired reads cover this position.
    pub fn proper_pair_fraction(&self) -> Option<f64> {
        if self.total_paired_count == 0 {
            None
        } else {
            Some(self.proper_pair_count as f64 / self.total_paired_count as f64)
        }
    }

    /// Set aggregate insert-size stats from contributing reads.
    pub fn with_insert_stats(mut self, insert_size_sum: i64, insert_size_count: u32) -> Self {
        self.insert_size_sum = insert_size_sum;
        self.insert_size_count = insert_size_count;
        self
    }

    /// Add insert-size stats from another read at the same position (used during merge).
    pub fn add_insert_stats(&mut self, insert_size_sum: i64, insert_size_count: u32) {
        self.insert_size_sum += insert_size_sum;
        self.insert_size_count += insert_size_count;
    }

    /// Mean absolute insert size across paired reads at this position.
    /// Returns None if no paired reads contributed insert sizes.
    pub fn mean_insert_size(&self) -> Option<f64> {
        if self.insert_size_count == 0 {
            None
        } else {
            Some(self.insert_size_sum as f64 / self.insert_size_count as f64)
        }
    }

    /// Raw insert-size stats: (sum, count).
    pub fn insert_stats(&self) -> (i64, u32) {
        (self.insert_size_sum, self.insert_size_count)
    }

    /// Increment the count of paired reads whose mate was unmapped at this position.
    pub fn add_unmapped_mate_count(&mut self, count: u32) {
        self.unmapped_mate_count += count;
    }

    /// Set unmapped-mate count (builder pattern).
    pub fn with_unmapped_mate_count(mut self, count: u32) -> Self {
        self.unmapped_mate_count = count;
        self
    }

    /// Fraction of paired reads at this position whose mate was unmapped.
    /// Returns None if no paired reads cover this position.
    pub fn unmapped_mate_fraction(&self) -> Option<f64> {
        if self.total_paired_count == 0 {
            None
        } else {
            Some(self.unmapped_mate_count as f64 / self.total_paired_count as f64)
        }
    }

    /// Raw unmapped-mate count.
    pub fn unmapped_mate_count(&self) -> u32 {
        self.unmapped_mate_count
    }

    /// Get the position of this variant
    pub fn position(&self) -> u32 {
        self.position
    }

    /// Get the reference base
    pub fn ref_base(&self) -> u8 {
        self.ref_base
    }

    /// Get all alleles with their counts
    pub fn all_alleles(&self) -> &HashMap<u8, u32> {
        &self.all_alleles
    }

    /// Get the confidence score
    pub fn confidence(&self) -> f64 {
        self.confidence
    }

    /// Get the CIGAR string
    pub fn cigar(&self) -> &str {
        &self.cigar
    }

    /// Get the mapping quality
    pub fn mapq(&self) -> u8 {
        self.mapq
    }

    /// Get the edit distance
    pub fn edit_distance(&self) -> u32 {
        self.edit_distance
    }

    /// Get the local coverage vector
    pub fn local_coverage(&self) -> &Vec<u32> {
        &self.local_coverage
    }

    /// Get the average base quality
    pub fn avg_base_quality(&self) -> f64 {
        self.avg_base_quality
    }

    /// Get the provenance string
    pub fn provenance(&self) -> &str {
        &self.provenance
    }
}

/// Evidence layer containing k-mer and alignment-derived evidence
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct EvidenceLayer {
    /// K-mer uniqueness at each position (value 0.0-1.0)
    kmer_uniqueness: HashMap<u32, f64>,
    /// Polymorphic sites and their alleles
    polymorphic_sites: HashMap<u32, Vec<u8>>,
    /// Invariant positions (where all samples match reference)
    invariant_positions: HashSet<u32>,
    /// Multi-mapping fraction at each position (0.0-1.0)
    multi_map_fraction: HashMap<u32, f64>,
    /// Average score ratio gap to next-best alignment
    avg_score_ratio_gap: HashMap<u32, f64>,
}

impl EvidenceLayer {
    /// Create a new evidence layer
    pub fn new(
        kmer_uniqueness: HashMap<u32, f64>,
        polymorphic_sites: HashMap<u32, Vec<u8>>,
        invariant_positions: HashSet<u32>,
        multi_map_fraction: HashMap<u32, f64>,
        avg_score_ratio_gap: HashMap<u32, f64>,
    ) -> Self {
        EvidenceLayer {
            kmer_uniqueness,
            polymorphic_sites,
            invariant_positions,
            multi_map_fraction,
            avg_score_ratio_gap,
        }
    }

    /// Get the k-mer uniqueness map
    pub fn kmer_uniqueness(&self) -> &HashMap<u32, f64> {
        &self.kmer_uniqueness
    }

    /// Get the polymorphic sites map
    pub fn polymorphic_sites(&self) -> &HashMap<u32, Vec<u8>> {
        &self.polymorphic_sites
    }

    /// Get the invariant positions set
    pub fn invariant_positions(&self) -> &HashSet<u32> {
        &self.invariant_positions
    }

    /// Get the multi-mapping fraction map
    pub fn multi_map_fraction(&self) -> &HashMap<u32, f64> {
        &self.multi_map_fraction
    }

    /// Get the average score ratio gap map
    pub fn avg_score_ratio_gap(&self) -> &HashMap<u32, f64> {
        &self.avg_score_ratio_gap
    }
}

/// Coverage track with RLE compression and quantization to nearest 5.
/// Stores (value, length) pairs for efficient representation of coverage across reference.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct CoverageTrack {
    runs: Vec<(u8, u32)>, // (quantized coverage, run length)
    total_length: u32,    // total reference length
}

impl CoverageTrack {
    /// Create a CoverageTrack from an array of coverage values.
    /// Values are quantized to nearest 5 (0, 5, 10, 15, ..., 255).
    pub fn new(coverage: Vec<usize>) -> Self {
        let total_length = coverage.len() as u32;
        let quantized: Vec<u8> = coverage.iter().map(|&c| Self::quantize(c)).collect();

        let mut runs = Vec::new();
        let mut current_val = 0u8;
        let mut current_len = 0u32;

        for &val in &quantized {
            if val == current_val {
                current_len += 1;
            } else {
                if current_len > 0 {
                    runs.push((current_val, current_len));
                }
                current_val = val;
                current_len = 1;
            }
        }

        if current_len > 0 {
            runs.push((current_val, current_len));
        }

        CoverageTrack { runs, total_length }
    }

    /// Quantize a coverage value to the nearest multiple of 5.
    /// 0-2 → 0, 3-7 → 5, 8-12 → 10, etc.
    pub fn quantize(value: usize) -> u8 {
        let rounded = ((value + 2) / 5) * 5;
        (rounded.min(255)) as u8
    }

    /// Get coverage at a specific position via binary search on runs.
    pub fn coverage_at(&self, pos: u32) -> Option<u8> {
        if pos >= self.total_length {
            return None;
        }

        let mut current_pos = 0u32;
        for &(value, length) in &self.runs {
            if pos < current_pos + length {
                return Some(value);
            }
            current_pos += length;
        }
        None
    }

    /// Get total reference length.
    pub fn total_length(&self) -> u32 {
        self.total_length
    }

    /// Decompress to full coverage array (for validation or downstream processing).
    pub fn decompress(&self) -> Vec<u8> {
        let mut result = Vec::with_capacity(self.total_length as usize);
        for &(value, length) in &self.runs {
            for _ in 0..length {
                result.push(value);
            }
        }
        result
    }

    /// Iterate over (position, coverage) pairs.
    pub fn iter(&self) -> impl Iterator<Item = (u32, u8)> + '_ {
        let mut pos = 0u32;
        self.runs.iter().flat_map(move |&(value, length)| {
            let start_pos = pos;
            pos += length;
            (0..length).map(move |i| (start_pos + i, value))
        })
    }

    /// Get compression ratio as (compressed_size, original_size).
    pub fn compression_ratio(&self) -> (usize, usize) {
        let compressed = self.runs.len() * std::mem::size_of::<(u8, u32)>();
        let original = self.total_length as usize * std::mem::size_of::<u8>();
        (compressed, original)
    }
}

/// Parse errors for FASTA, FASTQ, and other input formats
#[derive(Debug, Clone, Error, Serialize, Deserialize, PartialEq)]
pub enum ParseError {
    #[error("invalid format: {0}")]
    InvalidFormat(String),
    #[error("invalid UTF-8: {0}")]
    InvalidUtf8(String),
}

/// I/O errors for file operations
#[derive(Debug, Clone, Error, Serialize, Deserialize, PartialEq)]
pub enum IoError {
    #[error("file not found: {0}")]
    FileNotFound(String),
    #[error("permission denied: {0}")]
    PermissionDenied(String),
    #[error("read error: {0}")]
    ReadError(String),
}

/// Alignment errors during seeding, extension, or scoring
#[derive(Debug, Clone, Error, Serialize, Deserialize, PartialEq)]
pub enum AlignmentError {
    #[error("no seeds found for sequence: {0}")]
    NoSeeds(String),
    #[error("alignment failed: {0}")]
    AlignmentFailed(String),
}

/// Filter errors during threshold or expression evaluation
#[derive(Debug, Clone, Error, Serialize, Deserialize, PartialEq)]
pub enum FilterError {
    #[error("invalid threshold for {0}: {1}")]
    InvalidThreshold(String, f64),
    #[error("invalid filter expression: {0}")]
    InvalidExpression(String),
}

#[cfg(test)]
mod tests {
    use super::*;

    // ===== Sequence type tests =====

    #[test]
    fn sequence_creation_without_quality() {
        let seq = Sequence::new(
            b"ACGT".to_vec(),
            None,
            "seq1".to_string(),
            Some("test sequence".to_string()),
        );
        assert_eq!(seq.len(), 4);
        assert_eq!(seq.id(), "seq1");
        assert_eq!(seq.description(), Some("test sequence"));
        assert!(seq.quality_scores().is_none());
    }

    #[test]
    fn sequence_creation_with_quality() {
        let seq = Sequence::new(
            b"ACGT".to_vec(),
            Some(vec![30, 35, 40, 38]),
            "seq2".to_string(),
            None,
        );
        assert_eq!(seq.len(), 4);
        assert_eq!(seq.quality_at(0), Some(30));
        assert_eq!(seq.quality_at(1), Some(35));
        assert_eq!(seq.quality_at(2), Some(40));
        assert_eq!(seq.quality_at(3), Some(38));
        assert_eq!(seq.quality_at(4), None);
    }

    #[test]
    fn sequence_avg_quality_calculation() {
        let seq = Sequence::new(
            b"ACGT".to_vec(),
            Some(vec![20, 30, 40, 30]),
            "seq3".to_string(),
            None,
        );
        assert_eq!(seq.avg_quality(), Some(30.0));
    }

    #[test]
    fn sequence_avg_quality_without_scores() {
        let seq = Sequence::new(b"ACGT".to_vec(), None, "seq4".to_string(), None);
        assert_eq!(seq.avg_quality(), None);
    }

    #[test]
    fn sequence_empty() {
        let seq = Sequence::new(b"".to_vec(), None, "empty".to_string(), None);
        assert_eq!(seq.len(), 0);
        assert_eq!(seq.avg_quality(), None);
    }

    #[test]
    #[should_panic(expected = "quality scores length must match sequence length")]
    fn sequence_quality_length_mismatch_panics() {
        let _seq = Sequence::new(
            b"ACGT".to_vec(),
            Some(vec![30, 35]), // Only 2 quality scores for 4 bases
            "bad".to_string(),
            None,
        );
    }

    #[test]
    fn sequence_serialization() {
        let seq = Sequence::new(
            b"ACGT".to_vec(),
            Some(vec![30, 35, 40, 38]),
            "seq5".to_string(),
            Some("description".to_string()),
        );

        let json = serde_json::to_string(&seq).expect("serialization failed");
        let deserialized: Sequence = serde_json::from_str(&json).expect("deserialization failed");

        assert_eq!(deserialized.len(), 4);
        assert_eq!(deserialized.id(), "seq5");
        assert_eq!(deserialized.quality_at(0), Some(30));
    }

    // ===== VariantObservation type tests =====

    #[test]
    fn variant_observation_creation() {
        let mut alleles = std::collections::HashMap::new();
        alleles.insert(b'A', 10);
        alleles.insert(b'T', 5);

        let obs = VariantObservation::new(
            100,
            b'A',
            alleles,
            0.95,
            "10M".to_string(),
            60,
            2,
            vec![10, 12, 15, 18, 20],
            35.5,
            "sample1:read42".to_string(),
        );

        assert_eq!(obs.position(), 100);
        assert_eq!(obs.ref_base(), b'A');
        assert_eq!(obs.confidence(), 0.95);
        assert_eq!(obs.mapq(), 60);
        assert_eq!(obs.edit_distance(), 2);
        assert_eq!(obs.avg_base_quality(), 35.5);
        assert_eq!(obs.provenance(), "sample1:read42");
    }

    #[test]
    fn variant_observation_allele_counts() {
        let mut alleles = std::collections::HashMap::new();
        alleles.insert(b'A', 10);
        alleles.insert(b'T', 5);
        alleles.insert(b'G', 2);

        let obs = VariantObservation::new(
            200,
            b'A',
            alleles.clone(),
            0.98,
            "5M1I4M".to_string(),
            50,
            1,
            vec![20, 22, 25],
            40.0,
            "sample2:read99".to_string(),
        );

        let all_alleles = obs.all_alleles();
        assert_eq!(all_alleles.get(&b'A'), Some(&10));
        assert_eq!(all_alleles.get(&b'T'), Some(&5));
        assert_eq!(all_alleles.get(&b'G'), Some(&2));

        let total: u32 = all_alleles.values().sum();
        assert_eq!(total, 17);
    }

    #[test]
    fn variant_observation_local_coverage() {
        let obs = VariantObservation::new(
            150,
            b'C',
            [(b'C', 8), (b'T', 2)].into_iter().collect(),
            0.90,
            "20M".to_string(),
            55,
            0,
            vec![8, 9, 10, 10, 10, 12, 15, 18, 20, 22],
            38.0,
            "sample3:read1".to_string(),
        );

        let coverage = obs.local_coverage();
        assert_eq!(coverage.len(), 10);
        assert_eq!(coverage[0], 8);
        assert_eq!(coverage[9], 22);
    }

    #[test]
    fn variant_observation_serialization() {
        let obs = VariantObservation::new(
            300,
            b'G',
            [(b'G', 15), (b'A', 3)].into_iter().collect(),
            0.99,
            "25M".to_string(),
            60,
            0,
            vec![15, 16, 18, 20],
            42.0,
            "sample4:read5".to_string(),
        );

        let json = serde_json::to_string(&obs).expect("serialization failed");
        let deserialized: VariantObservation =
            serde_json::from_str(&json).expect("deserialization failed");

        assert_eq!(deserialized.position(), 300);
        assert_eq!(deserialized.ref_base(), b'G');
        assert_eq!(deserialized.mapq(), 60);
    }

    // ===== EvidenceLayer type tests =====

    #[test]
    fn evidence_layer_creation() {
        let mut kmer_uniqueness = std::collections::HashMap::new();
        kmer_uniqueness.insert(100, 1.0);
        kmer_uniqueness.insert(200, 0.5);

        let mut polymorphic_sites = std::collections::HashMap::new();
        polymorphic_sites.insert(150, vec![b'A', b'T']);

        let mut invariant_positions = std::collections::HashSet::new();
        invariant_positions.insert(50);
        invariant_positions.insert(51);

        let mut multi_map_fraction = std::collections::HashMap::new();
        multi_map_fraction.insert(100, 0.2);

        let mut avg_score_ratio_gap = std::collections::HashMap::new();
        avg_score_ratio_gap.insert(100, 0.15);

        let evidence = EvidenceLayer::new(
            kmer_uniqueness,
            polymorphic_sites,
            invariant_positions,
            multi_map_fraction,
            avg_score_ratio_gap,
        );

        assert_eq!(evidence.kmer_uniqueness().get(&100), Some(&1.0));
        assert_eq!(evidence.kmer_uniqueness().get(&200), Some(&0.5));
        assert_eq!(
            evidence.polymorphic_sites().get(&150),
            Some(&vec![b'A', b'T'])
        );
        assert!(evidence.invariant_positions().contains(&50));
        assert!(evidence.invariant_positions().contains(&51));
        assert_eq!(evidence.multi_map_fraction().get(&100), Some(&0.2));
        assert_eq!(evidence.avg_score_ratio_gap().get(&100), Some(&0.15));
    }

    #[test]
    fn evidence_layer_empty() {
        let evidence = EvidenceLayer::new(
            std::collections::HashMap::new(),
            std::collections::HashMap::new(),
            std::collections::HashSet::new(),
            std::collections::HashMap::new(),
            std::collections::HashMap::new(),
        );

        assert!(evidence.kmer_uniqueness().is_empty());
        assert!(evidence.polymorphic_sites().is_empty());
        assert!(evidence.invariant_positions().is_empty());
        assert!(evidence.multi_map_fraction().is_empty());
        assert!(evidence.avg_score_ratio_gap().is_empty());
    }

    #[test]
    fn evidence_layer_serialization() {
        let mut kmer_uniqueness = std::collections::HashMap::new();
        kmer_uniqueness.insert(100, 1.0);

        let evidence = EvidenceLayer::new(
            kmer_uniqueness,
            std::collections::HashMap::new(),
            std::collections::HashSet::new(),
            std::collections::HashMap::new(),
            std::collections::HashMap::new(),
        );

        let json = serde_json::to_string(&evidence).expect("serialization failed");
        let deserialized: EvidenceLayer =
            serde_json::from_str(&json).expect("deserialization failed");

        assert_eq!(deserialized.kmer_uniqueness().get(&100), Some(&1.0));
    }

    // ===== CoverageTrack type tests =====

    #[test]
    fn coverage_track_quantization_zeros() {
        assert_eq!(CoverageTrack::quantize(0), 0);
        assert_eq!(CoverageTrack::quantize(1), 0);
        assert_eq!(CoverageTrack::quantize(2), 0);
    }

    #[test]
    fn coverage_track_quantization_fives() {
        assert_eq!(CoverageTrack::quantize(3), 5);
        assert_eq!(CoverageTrack::quantize(4), 5);
        assert_eq!(CoverageTrack::quantize(5), 5);
        assert_eq!(CoverageTrack::quantize(6), 5);
        assert_eq!(CoverageTrack::quantize(7), 5);
    }

    #[test]
    fn coverage_track_quantization_tens() {
        assert_eq!(CoverageTrack::quantize(8), 10);
        assert_eq!(CoverageTrack::quantize(12), 10);
        assert_eq!(CoverageTrack::quantize(13), 15);
    }

    #[test]
    fn coverage_track_uniform_coverage() {
        let coverage = vec![10, 10, 10, 10];
        let track = CoverageTrack::new(coverage);

        assert_eq!(track.total_length(), 4);
        assert_eq!(track.coverage_at(0), Some(10));
        assert_eq!(track.coverage_at(2), Some(10));
        assert_eq!(track.coverage_at(3), Some(10));
        assert_eq!(track.coverage_at(4), None);
    }

    #[test]
    fn coverage_track_alternating_coverage() {
        let coverage = vec![10, 5, 10, 5, 10];
        let track = CoverageTrack::new(coverage);

        assert_eq!(track.coverage_at(0), Some(10));
        assert_eq!(track.coverage_at(1), Some(5));
        assert_eq!(track.coverage_at(2), Some(10));
        assert_eq!(track.coverage_at(3), Some(5));
        assert_eq!(track.coverage_at(4), Some(10));
    }

    #[test]
    fn coverage_track_zero_coverage() {
        let coverage = vec![0, 0, 5, 5, 0, 0];
        let track = CoverageTrack::new(coverage);

        assert_eq!(track.coverage_at(0), Some(0));
        assert_eq!(track.coverage_at(1), Some(0));
        assert_eq!(track.coverage_at(2), Some(5));
        assert_eq!(track.coverage_at(3), Some(5));
        assert_eq!(track.coverage_at(4), Some(0));
    }

    #[test]
    fn coverage_track_round_trip_encoding() {
        let coverage = vec![10, 10, 20, 20, 5, 5, 15, 15];
        let track = CoverageTrack::new(coverage.clone());
        let decompressed = track.decompress();

        let quantized_expected: Vec<u8> = coverage
            .iter()
            .map(|&c| CoverageTrack::quantize(c))
            .collect();

        assert_eq!(decompressed, quantized_expected);
    }

    #[test]
    fn coverage_track_quantization_idempotence() {
        // Quantizing twice should equal quantizing once
        let val = 17usize; // Rounds to 15
        let once = CoverageTrack::quantize(val);
        let twice = CoverageTrack::quantize(once as usize);
        assert_eq!(once, twice);
    }

    #[test]
    fn coverage_track_iterator() {
        let coverage = vec![10, 10, 5, 5];
        let track = CoverageTrack::new(coverage);

        let positions_and_coverage: Vec<_> = track.iter().collect();
        assert_eq!(positions_and_coverage.len(), 4);
        assert_eq!(positions_and_coverage[0], (0, 10));
        assert_eq!(positions_and_coverage[1], (1, 10));
        assert_eq!(positions_and_coverage[2], (2, 5));
        assert_eq!(positions_and_coverage[3], (3, 5));
    }

    #[test]
    fn coverage_track_compression_ratio() {
        let coverage = vec![10; 1000]; // 1000 positions with same coverage
        let track = CoverageTrack::new(coverage);

        let (compressed, original) = track.compression_ratio();
        // Should compress to 1 run of 1000 positions
        assert!(compressed < original / 2); // At least 50% compression
    }

    #[test]
    fn coverage_track_alternating_pattern_rle_overhead() {
        // Alternating patterns cause RLE overhead, not compression
        let coverage = (0..200)
            .map(|i| if i % 2 == 0 { 10 } else { 5 })
            .collect::<Vec<_>>();
        let track = CoverageTrack::new(coverage);

        // 200 positions, alternating between 10 and 5, creates 200 runs
        // Each run is 5 bytes (u8 + u32), so 1000 bytes compressed vs 200 bytes original
        // This is expected overhead for alternating patterns
        let (compressed, original) = track.compression_ratio();
        assert!(compressed > original); // RLE expands for alternating patterns
    }

    #[test]
    fn coverage_track_single_position() {
        let coverage = vec![15];
        let track = CoverageTrack::new(coverage);

        assert_eq!(track.total_length(), 1);
        assert_eq!(track.coverage_at(0), Some(15));
        assert_eq!(track.coverage_at(1), None);
    }

    #[test]
    fn coverage_track_empty() {
        let coverage = vec![];
        let track = CoverageTrack::new(coverage);

        assert_eq!(track.total_length(), 0);
        assert_eq!(track.coverage_at(0), None);
    }

    #[test]
    fn coverage_track_serialization() {
        let coverage = vec![10, 10, 5, 5];
        let track = CoverageTrack::new(coverage);

        let json = serde_json::to_string(&track).expect("serialization failed");
        let deserialized: CoverageTrack =
            serde_json::from_str(&json).expect("deserialization failed");

        assert_eq!(track, deserialized);
    }

    // ===== Error type tests =====

    #[test]
    fn parse_error_creation() {
        let err = ParseError::InvalidFormat("bad FASTA".to_string());
        assert_eq!(format!("{}", err), "invalid format: bad FASTA");
    }

    #[test]
    fn io_error_creation() {
        let err = IoError::FileNotFound("/path/to/file".to_string());
        assert_eq!(format!("{}", err), "file not found: /path/to/file");
    }

    #[test]
    fn alignment_error_creation() {
        let err = AlignmentError::NoSeeds("seq1".to_string());
        assert_eq!(format!("{}", err), "no seeds found for sequence: seq1");
    }

    #[test]
    fn filter_error_creation() {
        let err = FilterError::InvalidThreshold("coverage".to_string(), -5.0);
        assert_eq!(format!("{}", err), "invalid threshold for coverage: -5");
    }

    #[test]
    fn error_serialization() {
        // Errors should be serializable for structured logging
        let err = ParseError::InvalidFormat("test".to_string());
        let json = serde_json::to_string(&err).expect("serialization failed");
        let deserialized: ParseError = serde_json::from_str(&json).expect("deserialization failed");
        assert_eq!(format!("{}", deserialized), format!("{}", err));
    }

    // ===== Property tests =====

    #[test]
    fn property_quality_scores_match_sequence_length() {
        // Property: if quality scores are present, they must match sequence length
        let sequences = vec![
            (b"ACGT".to_vec(), vec![30, 35, 40, 38]),
            (b"A".to_vec(), vec![40]),
            (b"ACGTACGT".to_vec(), vec![30, 30, 30, 30, 40, 40, 40, 40]),
        ];

        for (bases, quality) in sequences {
            let seq = Sequence::new(
                bases.clone(),
                Some(quality.clone()),
                "test".to_string(),
                None,
            );
            assert_eq!(bases.len(), quality.len());
            assert_eq!(seq.len(), quality.len());
        }
    }

    #[test]
    fn property_allele_counts_are_positive() {
        // Property: all allele counts must be positive
        let mut alleles = std::collections::HashMap::new();
        alleles.insert(b'A', 10);
        alleles.insert(b'T', 5);
        alleles.insert(b'G', 1);

        let obs = VariantObservation::new(
            100,
            b'A',
            alleles.clone(),
            0.95,
            "10M".to_string(),
            60,
            0,
            vec![10],
            35.0,
            "test".to_string(),
        );

        for count in obs.all_alleles().values() {
            assert!(*count > 0, "all allele counts must be positive");
        }
    }

    #[test]
    fn property_allele_counts_sum_correctly() {
        // Property: sum of allele counts should equal total coverage at position
        let mut alleles = std::collections::HashMap::new();
        alleles.insert(b'A', 10);
        alleles.insert(b'T', 5);
        alleles.insert(b'G', 2);
        alleles.insert(b'C', 3);

        let obs = VariantObservation::new(
            100,
            b'A',
            alleles.clone(),
            0.95,
            "10M".to_string(),
            60,
            0,
            vec![20],
            35.0,
            "test".to_string(),
        );

        let sum: u32 = obs.all_alleles().values().sum();
        assert_eq!(sum, 20, "sum of allele counts should match total coverage");
    }

    #[test]
    fn property_ref_base_should_be_in_alleles() {
        // Property: reference base should typically be present in alleles map
        let mut alleles = std::collections::HashMap::new();
        alleles.insert(b'A', 10);
        alleles.insert(b'T', 5);

        let obs = VariantObservation::new(
            100,
            b'A',
            alleles.clone(),
            0.95,
            "10M".to_string(),
            60,
            0,
            vec![15],
            35.0,
            "test".to_string(),
        );

        assert!(
            obs.all_alleles().contains_key(&obs.ref_base()),
            "reference base should be in alleles map"
        );
    }

    #[test]
    fn property_confidence_in_valid_range() {
        // Property: confidence should be between 0.0 and 1.0
        let obs = VariantObservation::new(
            100,
            b'A',
            [(b'A', 10)].into_iter().collect(),
            0.95,
            "10M".to_string(),
            60,
            0,
            vec![10],
            35.0,
            "test".to_string(),
        );

        let conf = obs.confidence();
        assert!(
            conf >= 0.0 && conf <= 1.0,
            "confidence must be in range [0.0, 1.0]"
        );
    }

    #[test]
    fn property_mapq_in_valid_range() {
        // Property: MAPQ should be in valid range [0, 60]
        let obs = VariantObservation::new(
            100,
            b'A',
            [(b'A', 10)].into_iter().collect(),
            0.95,
            "10M".to_string(),
            60,
            0,
            vec![10],
            35.0,
            "test".to_string(),
        );

        let mapq = obs.mapq();
        assert!(mapq <= 60, "MAPQ should be in valid range [0, 60]");
    }

    #[test]
    fn property_kmer_uniqueness_in_valid_range() {
        // Property: k-mer uniqueness values should be between 0.0 and 1.0
        let mut kmer_uniqueness = std::collections::HashMap::new();
        kmer_uniqueness.insert(100, 1.0);
        kmer_uniqueness.insert(200, 0.5);
        kmer_uniqueness.insert(300, 0.0);

        let evidence = EvidenceLayer::new(
            kmer_uniqueness.clone(),
            std::collections::HashMap::new(),
            std::collections::HashSet::new(),
            std::collections::HashMap::new(),
            std::collections::HashMap::new(),
        );

        for value in evidence.kmer_uniqueness().values() {
            assert!(
                *value >= 0.0 && *value <= 1.0,
                "k-mer uniqueness must be in range [0.0, 1.0]"
            );
        }
    }

    #[test]
    fn property_multi_map_fraction_in_valid_range() {
        // Property: multi-mapping fraction should be between 0.0 and 1.0
        let mut multi_map_fraction = std::collections::HashMap::new();
        multi_map_fraction.insert(100, 0.0);
        multi_map_fraction.insert(200, 0.5);
        multi_map_fraction.insert(300, 1.0);

        let evidence = EvidenceLayer::new(
            std::collections::HashMap::new(),
            std::collections::HashMap::new(),
            std::collections::HashSet::new(),
            multi_map_fraction.clone(),
            std::collections::HashMap::new(),
        );

        for value in evidence.multi_map_fraction().values() {
            assert!(
                *value >= 0.0 && *value <= 1.0,
                "multi-map fraction must be in range [0.0, 1.0]"
            );
        }
    }
}

// ============================================================================
// Minimizer sketching
// ============================================================================

/// Default k-mer length for minimizer sketching (standard for bacterial genomics)
pub const DEFAULT_K: usize = 21;

/// Default window length for minimizer sketching
pub const DEFAULT_W: usize = 11;

/// A minimizer sketch: sparse representation of a sequence as (hash, position) pairs.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct MinimizerSketch {
    /// (minimizer_hash, position) pairs — position is 0-indexed into the original sequence
    pub minimizers: Vec<(u64, u32)>,
    /// k-mer length
    pub k: usize,
    /// window length
    pub w: usize,
}

impl MinimizerSketch {
    pub fn len(&self) -> usize {
        self.minimizers.len()
    }

    pub fn is_empty(&self) -> bool {
        self.minimizers.is_empty()
    }

    /// Find minimizers shared between self and other.
    /// Returns Vec of (hash, self_pos, other_pos).
    pub fn find_shared_minimizers(&self, other: &MinimizerSketch) -> Vec<(u64, u32, u32)> {
        let mut other_map: HashMap<u64, Vec<u32>> = HashMap::new();
        for &(val, pos) in &other.minimizers {
            other_map.entry(val).or_default().push(pos);
        }
        let mut shared = Vec::new();
        for &(val, qpos) in &self.minimizers {
            if let Some(tposs) = other_map.get(&val) {
                for &tpos in tposs {
                    shared.push((val, qpos, tpos));
                }
            }
        }
        shared.sort_by_key(|&(_, qpos, _)| qpos);
        shared
    }
}

/// Sketch a raw byte sequence using simd-minimizers canonical minimizers.
pub fn sketch(sequence: &[u8], k: usize, w: usize) -> MinimizerSketch {
    use packed_seq::AsciiSeq;
    let mut positions = Vec::new();
    let output = simd_minimizers::canonical_minimizers(k, w)
        .run(AsciiSeq(sequence), &mut positions);
    let minimizers: Vec<(u64, u32)> = output
        .pos_and_values_u64()
        .map(|(pos, val)| (val, pos))
        .collect();
    MinimizerSketch { minimizers, k, w }
}

/// Sketch a Sequence using given k and w.
pub fn sketch_sequence(seq: &Sequence, k: usize, w: usize) -> MinimizerSketch {
    sketch(seq.bases(), k, w)
}

/// Sketch a Sequence using default parameters (k=21, w=11).
pub fn sketch_sequence_default(seq: &Sequence) -> MinimizerSketch {
    sketch_sequence(seq, DEFAULT_K, DEFAULT_W)
}

fn jaccard_similarity(a: &MinimizerSketch, b: &MinimizerSketch) -> f64 {
    let set_a: HashSet<u64> = a.minimizers.iter().map(|&(val, _)| val).collect();
    let set_b: HashSet<u64> = b.minimizers.iter().map(|&(val, _)| val).collect();
    if set_a.is_empty() && set_b.is_empty() {
        return 1.0;
    }
    let intersection = set_a.intersection(&set_b).count();
    let union = set_a.union(&set_b).count();
    if union == 0 { 0.0 } else { intersection as f64 / union as f64 }
}

/// Select the centroid sketch — the one with median Jaccard similarity to all others.
/// Returns the index of the centroid in the input slice, or None if empty.
pub fn select_centroid(sketches: &[MinimizerSketch]) -> Option<usize> {
    if sketches.is_empty() {
        return None;
    }
    if sketches.len() == 1 {
        return Some(0);
    }
    let mut avg_sims: Vec<f64> = sketches
        .iter()
        .enumerate()
        .map(|(i, sk_i)| {
            let mut sims: Vec<f64> = sketches
                .iter()
                .enumerate()
                .filter(|&(j, _)| j != i)
                .map(|(_, sk_j)| jaccard_similarity(sk_i, sk_j))
                .collect();
            sims.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let n = sims.len();
            if n % 2 == 0 { (sims[n / 2 - 1] + sims[n / 2]) / 2.0 } else { sims[n / 2] }
        })
        .collect();
    let mut indexed: Vec<(usize, f64)> = avg_sims.drain(..).enumerate().collect();
    indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    Some(indexed[indexed.len() / 2].0)
}

/// Compute k-mer uniqueness scores from multiple sketches.
/// Returns position → uniqueness score in (0.0, 1.0].
pub fn compute_kmer_uniqueness(sketches: &[MinimizerSketch]) -> HashMap<u32, f64> {
    if sketches.is_empty() {
        return HashMap::new();
    }
    let mut counts: HashMap<(u64, u32), usize> = HashMap::new();
    for sk in sketches {
        for &(val, pos) in &sk.minimizers {
            *counts.entry((val, pos)).or_insert(0) += 1;
        }
    }
    let mut uniqueness: HashMap<u32, f64> = HashMap::new();
    for ((_, pos), count) in counts {
        let score = 1.0 / count as f64;
        uniqueness
            .entry(pos)
            .and_modify(|s| *s = s.min(score))
            .or_insert(score);
    }
    uniqueness
}

/// Detect variation hotspot intervals from k-mer uniqueness scores.
///
/// Scans the uniqueness map for positions where the score is below the threshold,
/// merging contiguous/adjacent positions into intervals (start, end).
/// Returns sorted intervals by start position.
///
/// # Arguments
///
/// * `uniqueness` - HashMap mapping position (u32) to uniqueness score (f64)
/// * `threshold` - Uniqueness threshold; positions with score < threshold are considered hotspots
///
/// # Returns
///
/// Vec of (start, end) intervals, sorted by start position. Returns empty vec if map is empty
/// or no positions fall below threshold.
pub fn detect_hotspot_intervals(uniqueness: &HashMap<u32, f64>, threshold: f64) -> Vec<(u32, u32)> {
    if uniqueness.is_empty() {
        return Vec::new();
    }

    // Collect positions below threshold and sort them
    let mut positions: Vec<u32> = uniqueness
        .iter()
        .filter(|(_, score)| **score < threshold)
        .map(|(&pos, _)| pos)
        .collect();

    if positions.is_empty() {
        return Vec::new();
    }

    positions.sort_unstable();

    // Merge contiguous/adjacent positions into intervals
    let mut intervals = Vec::new();
    let mut start = positions[0];
    let mut end = positions[0];

    for &pos in &positions[1..] {
        if pos == end + 1 {
            // Adjacent position, extend current interval
            end = pos;
        } else {
            // Gap detected, save current interval and start new one
            intervals.push((start, end));
            start = pos;
            end = pos;
        }
    }

    // Don't forget the last interval
    intervals.push((start, end));

    intervals
}

/// Mate relationship metadata for paired-end reads
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct MateInfo {
    /// Read mate identifier (e.g., "read123/2" for the mate of "read123/1")
    pub mate_id: String,

    /// SAM flag 0x2: read mapped in proper pair
    pub proper_pair: bool,

    /// Template length (TLEN field) - signed, 0 for unpaired/unmapped mates
    pub insert_size: i32,

    /// SAM flag 0x40: first read in pair
    pub is_first_in_pair: bool,

    /// SAM flag 0x80: second read in pair
    pub is_second_in_pair: bool,

    /// Mate is mapped (inverse of SAM flag 0x8)
    pub mate_mapped: bool,
}

impl MateInfo {
    pub fn new(
        mate_id: String,
        proper_pair: bool,
        insert_size: i32,
        is_first_in_pair: bool,
        is_second_in_pair: bool,
        mate_mapped: bool,
    ) -> Self {
        Self {
            mate_id,
            proper_pair,
            insert_size,
            is_first_in_pair,
            is_second_in_pair,
            mate_mapped,
        }
    }

    /// Check if insert size indicates discordant pair (beyond expected distribution)
    pub fn is_discordant(&self, mean: i32, std_dev: i32, sigma: f64) -> bool {
        if self.insert_size == 0 {
            return false; // Unpaired or unmapped mate
        }
        let threshold = (std_dev as f64 * sigma) as i32;
        // Use absolute value of insert_size (negative = mate upstream)
        let deviation = (self.insert_size.abs() - mean).abs();
        deviation > threshold
    }
}
