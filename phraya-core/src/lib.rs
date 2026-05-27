use phraya_index::FmIndex;

/// Sequence represents a DNA sequence with optional quality scores and alignment information.
#[derive(Clone, Debug)]
pub struct Sequence {
    pub bases: Vec<u8>,
    pub quality: Option<Vec<u8>>,
    pub metadata: Option<String>,
}

impl Sequence {
    /// Create a new Sequence
    pub fn new(bases: Vec<u8>, quality: Option<Vec<u8>>, metadata: Option<String>) -> Self {
        Sequence {
            bases,
            quality,
            metadata,
        }
    }

    /// Get the average quality score
    pub fn average_quality(&self) -> f64 {
        match &self.quality {
            None => 30.0, // Default if no quality scores
            Some(quals) => {
                if quals.is_empty() {
                    30.0
                } else {
                    let sum: f64 = quals.iter().map(|&q| q as f64).sum();
                    sum / quals.len() as f64
                }
            }
        }
    }
}

/// Feature weights for AlignmentContextScorer
#[derive(Clone, Debug)]
pub struct FeatureWeights {
    pub edge_distance_weight: f64,
    pub gc_content_weight: f64,
    pub kmer_uniqueness_weight: f64,
    pub repeat_penalty_weight: f64,
    pub snp_density_weight: f64,
    pub alignment_identity_weight: f64,
}

impl FeatureWeights {
    /// Create default normalized weights
    pub fn default_normalized() -> Self {
        // Weights are normalized to sum to 1.0
        // Alignment identity is the strongest signal for confidence
        FeatureWeights {
            edge_distance_weight: 0.1,
            gc_content_weight: 0.08,
            kmer_uniqueness_weight: 0.12,
            repeat_penalty_weight: 0.1,
            snp_density_weight: 0.1,
            alignment_identity_weight: 0.5,
        }
    }
}

/// BaseConfidence represents confidence metrics for a position
#[derive(Clone, Debug)]
pub struct BaseConfidence {
    combined_confidence: f64,
    edge_distance_score: f64,
    local_gc_content: f64,
    local_gc_content_score: f64,
    kmer_uniqueness_score: f64,
    in_repeat_region: bool,
    in_homopolymer: bool,
    repeat_period: Option<u32>,
    snp_density_15bp: f64,
    snp_density_125bp: f64,
    snp_density_1000bp: f64,
    alignment_identity_score: f64,
}

impl BaseConfidence {
    /// Combined confidence score [0.0, 1.0]
    pub fn combined_confidence(&self) -> f64 {
        self.combined_confidence
    }

    /// Edge distance score [0.0, 1.0]
    pub fn edge_distance_score(&self) -> f64 {
        self.edge_distance_score
    }

    /// Local GC content as fraction [0.0, 1.0]
    pub fn local_gc_content(&self) -> f64 {
        self.local_gc_content
    }

    /// Local GC content score [0.0, 1.0]
    pub fn local_gc_content_score(&self) -> f64 {
        self.local_gc_content_score
    }

    /// K-mer uniqueness score [0.0, 1.0]
    pub fn kmer_uniqueness_score(&self) -> f64 {
        self.kmer_uniqueness_score
    }

    /// Whether position is in a repeat region
    pub fn in_repeat_region(&self) -> bool {
        self.in_repeat_region
    }

    /// Whether position is in a homopolymer
    pub fn in_homopolymer(&self) -> bool {
        self.in_homopolymer
    }

    /// Repeat period if in repeat region
    pub fn repeat_period(&self) -> Option<u32> {
        self.repeat_period
    }

    /// SNP density in 15bp window
    pub fn snp_density_15bp(&self) -> f64 {
        self.snp_density_15bp
    }

    /// SNP density in 125bp window
    pub fn snp_density_125bp(&self) -> f64 {
        self.snp_density_125bp
    }

    /// SNP density in 1000bp window
    pub fn snp_density_1000bp(&self) -> f64 {
        self.snp_density_1000bp
    }

    /// Alignment identity score [0.0, 1.0]
    pub fn alignment_identity_score(&self) -> f64 {
        self.alignment_identity_score
    }
}

/// AlignmentContextScorer integrates multiple confidence features
pub struct AlignmentContextScorer {
    reference: Vec<u8>,
    _fm_index: FmIndex,
    weights: FeatureWeights,
}

impl AlignmentContextScorer {
    /// Create a new scorer
    pub fn new(reference: &[u8], fm_index: &FmIndex) -> Self {
        AlignmentContextScorer {
            reference: reference.to_vec(),
            _fm_index: fm_index.clone(),
            weights: FeatureWeights::default_normalized(),
        }
    }

    /// Score a position with an aligned sequence and identity
    pub fn score(
        &self,
        position: usize,
        aligned_sequence: &Sequence,
        alignment_identity: f64,
    ) -> BaseConfidence {
        // Calculate individual features
        let edge_distance_score = self.compute_edge_distance_score(position);
        let local_gc_content = self.compute_local_gc_content_from_sequence(aligned_sequence);
        let local_gc_content_score = self.compute_gc_content_score(local_gc_content);
        let kmer_uniqueness_score = self.compute_kmer_uniqueness_from_sequence(aligned_sequence);
        let (in_repeat_region, repeat_period) = self.detect_repeat_region(position);
        let in_homopolymer = self.detect_homopolymer(position);

        // Compute SNP density (based on alignment identity)
        let snp_density_15bp = self.estimate_snp_density(position, 15, alignment_identity);
        let snp_density_125bp = self.estimate_snp_density(position, 125, alignment_identity);
        let snp_density_1000bp = self.estimate_snp_density(position, 1000, alignment_identity);

        let alignment_identity_score = alignment_identity; // Direct mapping for now

        // Compute repeat penalty
        let repeat_penalty = if in_repeat_region { 0.7 } else { 1.0 };

        // Apply base quality boost if available
        let base_quality_boost = self.compute_base_quality_boost(aligned_sequence);

        // Combine features with weights
        let snp_density_score = 1.0 - (snp_density_125bp.min(1.0));
        let combined = self.weights.edge_distance_weight * edge_distance_score
            + self.weights.gc_content_weight * local_gc_content_score
            + self.weights.kmer_uniqueness_weight * kmer_uniqueness_score
            + self.weights.snp_density_weight * snp_density_score
            + self.weights.alignment_identity_weight * alignment_identity_score;

        // Apply repeat penalty
        let combined_with_penalty = combined * repeat_penalty * base_quality_boost;

        let combined_confidence = combined_with_penalty.max(0.0).min(1.0);

        BaseConfidence {
            combined_confidence,
            edge_distance_score,
            local_gc_content,
            local_gc_content_score,
            kmer_uniqueness_score,
            in_repeat_region,
            in_homopolymer,
            repeat_period,
            snp_density_15bp,
            snp_density_125bp,
            snp_density_1000bp,
            alignment_identity_score,
        }
    }

    /// Get feature weights for transparency
    pub fn feature_weights(&self) -> FeatureWeights {
        self.weights.clone()
    }

    // Helper methods for computing features

    fn compute_edge_distance_score(&self, position: usize) -> f64 {
        let ref_len = self.reference.len();
        let window = 50; // Consider 50bp margin as edge

        let dist_from_start = position.min(window) as f64;
        let dist_from_end = if ref_len > position {
            (ref_len - position - 1).min(window) as f64
        } else {
            0.0
        };

        let min_dist = dist_from_start.min(dist_from_end);
        let score = (min_dist / window as f64).min(1.0);

        score
    }

    fn compute_local_gc_content_from_sequence(&self, sequence: &Sequence) -> f64 {
        if sequence.bases.is_empty() {
            return 0.5;
        }

        let gc_count = sequence.bases.iter().filter(|&&b| b == b'G' || b == b'C').count();
        gc_count as f64 / sequence.bases.len() as f64
    }

    fn compute_gc_content_score(&self, gc_content: f64) -> f64 {
        // Score is higher near 0.5 (optimal for PCR and sequencing)
        // Lower at extremes
        let optimal = 0.5;
        let distance = (gc_content - optimal).abs();
        1.0 - (distance * 2.0).min(1.0)
    }

    fn compute_kmer_uniqueness_from_sequence(&self, sequence: &Sequence) -> f64 {
        // Compute k-mer uniqueness from the aligned sequence itself
        // This reflects the actual complexity of what was aligned

        if sequence.bases.is_empty() {
            return 0.5;
        }

        // Compute Shannon entropy of the sequence
        let entropy = self.compute_sequence_entropy(&sequence.bases);

        // For DNA, max entropy is log2(4) = 2.0 bits per base
        // Normalize to [0, 1]
        let max_entropy = 2.0;
        let normalized_entropy = (entropy / max_entropy).min(1.0);

        // Higher entropy = higher uniqueness (more diverse = more unique)
        normalized_entropy
    }

    fn compute_sequence_entropy(&self, seq: &[u8]) -> f64 {
        if seq.is_empty() {
            return 0.0;
        }

        let mut counts = [0usize; 256];
        for &b in seq {
            counts[b as usize] += 1;
        }

        let len = seq.len() as f64;
        let mut entropy = 0.0;

        for &count in counts.iter() {
            if count > 0 {
                let p = count as f64 / len;
                entropy -= p * p.log2();
            }
        }

        entropy
    }

    fn detect_repeat_region(&self, position: usize) -> (bool, Option<u32>) {
        // Detect tandem repeats by looking for repeated patterns
        let window = 50;
        let start = position.saturating_sub(window / 2);
        let end = (position + window / 2).min(self.reference.len());

        let window_seq = &self.reference[start..end];

        // Check for tandem repeat patterns (2-10bp period)
        for period in 2..=10 {
            if let Some(is_repeat) = self.is_tandem_repeat_with_period(window_seq, period) {
                if is_repeat {
                    return (true, Some(period as u32));
                }
            }
        }

        (false, None)
    }

    fn is_tandem_repeat_with_period(&self, seq: &[u8], period: usize) -> Option<bool> {
        if seq.len() < period * 2 {
            return None;
        }

        let pattern = &seq[0..period];
        let mut repeat_count = 1;

        for i in 1..(seq.len() / period) {
            let chunk = &seq[i * period..(i + 1) * period];
            if chunk == pattern {
                repeat_count += 1;
            } else {
                break;
            }
        }

        Some(repeat_count >= 3)
    }

    fn detect_homopolymer(&self, position: usize) -> bool {
        // Detect homopolymer runs (4+ consecutive identical bases)
        let window = 50;
        let start = position.saturating_sub(window / 2);
        let end = (position + window / 2).min(self.reference.len());

        let window_seq = &self.reference[start..end];

        let mut current_base = window_seq[0];
        let mut run_length = 1;

        for &base in &window_seq[1..] {
            if base == current_base {
                run_length += 1;
                if run_length >= 4 {
                    return true;
                }
            } else {
                current_base = base;
                run_length = 1;
            }
        }

        false
    }

    fn estimate_snp_density(&self, _position: usize, window_size: usize, alignment_identity: f64) -> f64 {
        // Estimate SNP density based on alignment identity
        // Lower identity indicates more mismatches (higher SNP density)
        let snps_per_position = 1.0 - alignment_identity;
        (snps_per_position * window_size as f64) / window_size as f64
    }

    fn compute_base_quality_boost(&self, aligned_sequence: &Sequence) -> f64 {
        // Boost confidence if base qualities are high
        let avg_quality = aligned_sequence.average_quality();
        // Phred quality 30 = 0.999 accuracy = good threshold
        // Phred quality 20 = 0.99 accuracy = marginal
        let boost = if avg_quality >= 30.0 {
            1.05
        } else if avg_quality >= 25.0 {
            1.02
        } else if avg_quality < 15.0 {
            0.95
        } else {
            1.0
        };
        boost
    }
}

pub fn add(left: u64, right: u64) -> u64 {
    left + right
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }
}
