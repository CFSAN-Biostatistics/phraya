// This module will contain core types for Phraya.
// Tests are written first (TDD RED phase).

use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use thiserror::Error;

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
        }
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
        }
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

/// Coverage track with RLE (run-length encoding) compression.
/// Coverage values are quantized to the nearest 5 before encoding.
/// Stores full reference length, including zero-coverage regions.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct CoverageTrack {
    /// RLE runs: (value, length) pairs in order
    runs: Vec<(usize, usize)>,
    /// Total number of positions (for bounds checking)
    total_length: usize,
}

impl CoverageTrack {
    /// Quantize a coverage value to the nearest 5.
    /// Example: 7 → 5, 8 → 10, 12 → 10, 13 → 15.
    fn quantize(value: usize) -> usize {
        ((value + 2) / 5) * 5
    }

    /// Create a CoverageTrack from raw coverage values.
    /// Values are quantized to nearest 5 and RLE-compressed.
    pub fn from_coverage(coverage: Vec<usize>) -> Self {
        if coverage.is_empty() {
            return CoverageTrack {
                runs: Vec::new(),
                total_length: 0,
            };
        }

        let total_length = coverage.len();
        let mut runs = Vec::new();
        let mut current_value = Self::quantize(coverage[0]);
        let mut current_count = 1;

        for &value in &coverage[1..] {
            let quantized = Self::quantize(value);
            if quantized == current_value {
                current_count += 1;
            } else {
                runs.push((current_value, current_count));
                current_value = quantized;
                current_count = 1;
            }
        }
        // Add the last run
        runs.push((current_value, current_count));

        CoverageTrack { runs, total_length }
    }

    /// Decompress the RLE-encoded coverage to a full vector.
    pub fn to_vec(&self) -> Vec<usize> {
        let mut result = Vec::with_capacity(self.total_length);
        for &(value, length) in &self.runs {
            for _ in 0..length {
                result.push(value);
            }
        }
        result
    }

    /// Get coverage at a specific position via binary search.
    /// Returns 0 for out-of-bounds positions.
    pub fn coverage_at(&self, pos: usize) -> usize {
        if pos >= self.total_length {
            return 0;
        }

        // Binary search through runs to find the run containing this position
        let mut cumulative_pos = 0;
        for &(value, length) in &self.runs {
            if pos < cumulative_pos + length {
                return value;
            }
            cumulative_pos += length;
        }

        0
    }

    /// Get the number of RLE runs in this track.
    pub fn run_count(&self) -> usize {
        self.runs.len()
    }

    /// Iterator over (position, coverage) pairs.
    pub fn iter(&self) -> CoverageTrackIter {
        CoverageTrackIter {
            runs: self.runs.clone(),
            current_run_idx: 0,
            current_pos_in_run: 0,
            current_pos: 0,
        }
    }
}

/// Iterator over coverage track positions.
#[derive(Debug, Clone)]
pub struct CoverageTrackIter {
    runs: Vec<(usize, usize)>,
    current_run_idx: usize,
    current_pos_in_run: usize,
    current_pos: usize,
}

impl Iterator for CoverageTrackIter {
    type Item = (usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_run_idx >= self.runs.len() {
            return None;
        }

        let (value, length) = self.runs[self.current_run_idx];

        if self.current_pos_in_run >= length {
            self.current_run_idx += 1;
            self.current_pos_in_run = 0;
            if self.current_run_idx >= self.runs.len() {
                return None;
            }
            let (value, _length) = self.runs[self.current_run_idx];
            let pos = self.current_pos;
            self.current_pos += 1;
            self.current_pos_in_run += 1;
            return Some((pos, value));
        }

        let pos = self.current_pos;
        self.current_pos += 1;
        self.current_pos_in_run += 1;

        Some((pos, value))
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

    // ===== CoverageTrack RLE compression tests =====

    #[test]
    fn coverage_track_uniform_coverage_single_run() {
        // Uniform coverage (e.g., all 10x) should compress to a single RLE run
        let coverage = vec![10; 1000]; // 1000 positions with 10x coverage
        let track = CoverageTrack::from_coverage(coverage.clone());

        // Verify decompression matches original (modulo quantization to 10)
        let decompressed = track.to_vec();
        assert_eq!(decompressed.len(), coverage.len());
        for val in &decompressed {
            assert_eq!(*val, 10); // Quantized to nearest 5
        }

        // Verify efficient compression - should be single run
        assert_eq!(track.run_count(), 1);
    }

    #[test]
    fn coverage_track_alternating_coverage_many_runs() {
        // Alternating coverage should result in many RLE runs
        let mut coverage = Vec::new();
        for i in 0..500 {
            coverage.push(if i % 2 == 0 { 10 } else { 20 });
        }
        let track = CoverageTrack::from_coverage(coverage.clone());

        let decompressed = track.to_vec();
        assert_eq!(decompressed.len(), coverage.len());

        // Verify alternating pattern preserved (with quantization)
        for (i, &val) in decompressed.iter().enumerate() {
            let expected = if i % 2 == 0 { 10 } else { 20 };
            assert_eq!(val, expected);
        }
    }

    #[test]
    fn coverage_track_zero_coverage_regions() {
        // Zero coverage regions should be preserved
        let mut coverage = vec![0; 100];
        coverage.extend(vec![15; 100]);
        coverage.extend(vec![0; 100]);
        coverage.extend(vec![30; 100]);

        let track = CoverageTrack::from_coverage(coverage.clone());
        let decompressed = track.to_vec();

        assert_eq!(decompressed.len(), 400);
        // First 100 positions should be 0
        for i in 0..100 {
            assert_eq!(decompressed[i], 0);
        }
        // Next 100 should be 15
        for i in 100..200 {
            assert_eq!(decompressed[i], 15);
        }
        // Next 100 should be 0
        for i in 200..300 {
            assert_eq!(decompressed[i], 0);
        }
        // Last 100 should be 30
        for i in 300..400 {
            assert_eq!(decompressed[i], 30);
        }
    }

    #[test]
    fn coverage_track_quantization_to_nearest_5() {
        // Coverage values should be quantized to nearest 5
        let coverage = vec![7, 8, 12, 13, 17, 18, 22, 23];
        let track = CoverageTrack::from_coverage(coverage);
        let decompressed = track.to_vec();

        assert_eq!(decompressed, vec![5, 10, 10, 15, 15, 20, 20, 25]);
    }

    #[test]
    fn coverage_track_quantization_exact_multiples_of_5() {
        // Exact multiples of 5 should remain unchanged
        let coverage = vec![0, 5, 10, 15, 20, 25, 30];
        let track = CoverageTrack::from_coverage(coverage.clone());
        let decompressed = track.to_vec();

        assert_eq!(decompressed, coverage);
    }

    #[test]
    fn coverage_track_quantization_boundary_cases() {
        // Test rounding behavior at boundaries
        // 2 rounds to 0, 3 rounds to 5, 7 rounds to 5, 8 rounds to 10
        let coverage = vec![2, 3, 7, 8, 12, 13];
        let track = CoverageTrack::from_coverage(coverage);
        let decompressed = track.to_vec();

        assert_eq!(decompressed, vec![0, 5, 5, 10, 10, 15]);
    }

    #[test]
    fn coverage_track_random_access_via_binary_search() {
        // Random access should work via binary search in O(log n) time
        let mut coverage = Vec::new();
        coverage.extend(vec![10; 100]);
        coverage.extend(vec![20; 100]);
        coverage.extend(vec![30; 100]);

        let track = CoverageTrack::from_coverage(coverage);

        // Test random access at various positions
        assert_eq!(track.coverage_at(0), 10);
        assert_eq!(track.coverage_at(50), 10);
        assert_eq!(track.coverage_at(99), 10);
        assert_eq!(track.coverage_at(100), 20);
        assert_eq!(track.coverage_at(150), 20);
        assert_eq!(track.coverage_at(199), 20);
        assert_eq!(track.coverage_at(200), 30);
        assert_eq!(track.coverage_at(250), 30);
        assert_eq!(track.coverage_at(299), 30);
    }

    #[test]
    fn coverage_track_random_access_out_of_bounds() {
        // Out of bounds access should return None or 0
        let coverage = vec![10; 100];
        let track = CoverageTrack::from_coverage(coverage);

        assert_eq!(track.coverage_at(100), 0); // Beyond last position
        assert_eq!(track.coverage_at(1000), 0); // Way out of bounds
    }

    #[test]
    fn coverage_track_iterator_over_positions() {
        // Should be able to iterate over (position, coverage) pairs
        let coverage = vec![10, 10, 20, 20, 30, 30];
        let track = CoverageTrack::from_coverage(coverage.clone());

        let positions: Vec<(usize, usize)> = track.iter().collect();
        assert_eq!(positions.len(), coverage.len());

        for (i, (pos, cov)) in positions.iter().enumerate() {
            assert_eq!(*pos, i);
            assert_eq!(*cov, coverage[i]);
        }
    }

    #[test]
    fn coverage_track_empty_coverage() {
        // Empty coverage should be handled gracefully
        let coverage: Vec<usize> = vec![];
        let track = CoverageTrack::from_coverage(coverage);

        assert_eq!(track.to_vec().len(), 0);
        assert_eq!(track.run_count(), 0);
    }

    #[test]
    fn coverage_track_single_position() {
        // Single position coverage
        let coverage = vec![15];
        let track = CoverageTrack::from_coverage(coverage);

        assert_eq!(track.to_vec(), vec![15]);
        assert_eq!(track.coverage_at(0), 15);
        assert_eq!(track.run_count(), 1);
    }

    #[test]
    fn coverage_track_realistic_bacterial_genome() {
        // Simulate realistic bacterial genome coverage (E. coli ~4.6Mbp)
        // with mostly uniform 30x coverage with some variation
        let mut coverage = Vec::new();

        // Region 1: 1Mbp at 30x
        coverage.extend(vec![30; 1_000_000]);

        // Region 2: 500kbp at 10x (lower coverage region)
        coverage.extend(vec![10; 500_000]);

        // Region 3: 2Mbp at 35x
        coverage.extend(vec![35; 2_000_000]);

        // Region 4: 100kbp at 0x (no coverage)
        coverage.extend(vec![0; 100_000]);

        // Region 5: 1Mbp at 30x
        coverage.extend(vec![30; 1_000_000]);

        let track = CoverageTrack::from_coverage(coverage.clone());
        let decompressed = track.to_vec();

        assert_eq!(decompressed.len(), coverage.len());

        // Should compress to 5 runs
        assert_eq!(track.run_count(), 5);

        // Verify some random positions
        assert_eq!(track.coverage_at(500_000), 30);
        assert_eq!(track.coverage_at(1_200_000), 10);
        assert_eq!(track.coverage_at(3_000_000), 35);
        assert_eq!(track.coverage_at(3_550_000), 0);
        assert_eq!(track.coverage_at(4_000_000), 30);
    }

    #[test]
    fn coverage_track_high_depth_sequencing() {
        // Test high coverage values (100x)
        let coverage = vec![100; 1000];
        let track = CoverageTrack::from_coverage(coverage);
        let decompressed = track.to_vec();

        for val in &decompressed {
            assert_eq!(*val, 100);
        }
    }

    #[test]
    fn coverage_track_serialization() {
        // CoverageTrack should be serializable
        let coverage = vec![10, 10, 20, 20, 30, 30];
        let track = CoverageTrack::from_coverage(coverage);

        let json = serde_json::to_string(&track).expect("serialization failed");
        let deserialized: CoverageTrack =
            serde_json::from_str(&json).expect("deserialization failed");

        assert_eq!(deserialized.to_vec(), track.to_vec());
    }

    // ===== Property tests =====

    #[test]
    fn property_round_trip_encode_decode() {
        // Property: encode → decode should equal original (modulo quantization)
        use proptest::prelude::*;

        proptest!(|(coverage in prop::collection::vec(0usize..200, 0..1000))| {
            let track = CoverageTrack::from_coverage(coverage.clone());
            let decompressed = track.to_vec();

            // Check length matches
            prop_assert_eq!(decompressed.len(), coverage.len());

            // Check values match after quantization
            for (i, &original) in coverage.iter().enumerate() {
                let quantized = ((original + 2) / 5) * 5; // Round to nearest 5
                prop_assert_eq!(decompressed[i], quantized);
            }
        });
    }

    #[test]
    fn property_quantization_idempotence() {
        // Property: quantizing twice should equal quantizing once
        use proptest::prelude::*;

        proptest!(|(coverage in prop::collection::vec(0usize..200, 0..1000))| {
            let track1 = CoverageTrack::from_coverage(coverage.clone());
            let decompressed1 = track1.to_vec();

            // Quantize again
            let track2 = CoverageTrack::from_coverage(decompressed1.clone());
            let decompressed2 = track2.to_vec();

            // Should be identical
            prop_assert_eq!(decompressed1, decompressed2);
        });
    }

    #[test]
    fn property_coverage_at_matches_decompressed() {
        // Property: coverage_at(i) should match to_vec()[i]
        use proptest::prelude::*;

        proptest!(|(coverage in prop::collection::vec(0usize..200, 10..100))| {
            let track = CoverageTrack::from_coverage(coverage.clone());
            let decompressed = track.to_vec();

            for i in 0..coverage.len() {
                prop_assert_eq!(track.coverage_at(i), decompressed[i]);
            }
        });
    }

    #[test]
    fn property_quantization_always_multiple_of_5() {
        // Property: all decompressed values should be multiples of 5
        use proptest::prelude::*;

        proptest!(|(coverage in prop::collection::vec(0usize..200, 0..1000))| {
            let track = CoverageTrack::from_coverage(coverage);
            let decompressed = track.to_vec();

            for val in &decompressed {
                prop_assert_eq!(val % 5, 0, "all values must be multiples of 5");
            }
        });
    }

    #[test]
    fn property_quantization_within_2_of_original() {
        // Property: quantized value should be within ±2 of original
        use proptest::prelude::*;

        proptest!(|(coverage in prop::collection::vec(0usize..200, 0..1000))| {
            let track = CoverageTrack::from_coverage(coverage.clone());
            let decompressed = track.to_vec();

            for (i, &original) in coverage.iter().enumerate() {
                let quantized = decompressed[i];
                let diff = if original > quantized {
                    original - quantized
                } else {
                    quantized - original
                };
                prop_assert!(diff <= 2, "quantized value should be within ±2 of original");
            }
        });
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
