//! Specialized aligner for short reads with paired-end support.
//!
//! Optimizes for typical short-read sequencing data (50–150bp single reads,
//! or 2×100bp paired-end reads) using:
//! - Seed-based indexing with high k-mer specificity
//! - Paired-end mate constraints (expected insert size, orientation)
//! - Fast strategy by default (K=1 chain cap, ±150bp coverage window)
//! - Early termination once a confident match is found

use crate::executor::{AlignConfig, Strategy};

/// Configuration for paired-end read alignment.
#[derive(Debug, Clone)]
pub struct PairedEndConfig {
    /// Expected insert size (mean length between read starts on forward strands)
    pub insert_size_mean: usize,
    /// Insert size standard deviation (allows ±3σ tolerance window)
    pub insert_size_stddev: usize,
    /// If true, require both mates to map (fail if one unmapped)
    pub require_both_mates: bool,
    /// If true, require proper pair orientation (F1R2 or R1F2)
    pub require_proper_orientation: bool,
}

impl Default for PairedEndConfig {
    fn default() -> Self {
        PairedEndConfig {
            insert_size_mean: 350,
            insert_size_stddev: 35,
            require_both_mates: true,
            require_proper_orientation: true,
        }
    }
}

/// Specialized aligner for short reads.
pub struct ShortReadAligner {
    /// Strategy to use (defaults to Fast for throughput)
    pub strategy: Strategy,
    /// Configuration for paired-end alignment (None for single-end)
    pub paired_config: Option<PairedEndConfig>,
}

impl ShortReadAligner {
    /// Create a new short-read aligner with default settings (Fast strategy, single-end).
    pub fn new() -> Self {
        ShortReadAligner {
            strategy: Strategy::Fast,
            paired_config: None,
        }
    }

    /// Create a short-read aligner optimized for paired-end reads.
    pub fn paired(config: PairedEndConfig) -> Self {
        ShortReadAligner {
            strategy: Strategy::Fast,
            paired_config: Some(config),
        }
    }

    /// Get the alignment configuration for this aligner.
    pub fn alignment_config(&self) -> AlignConfig {
        AlignConfig::new(self.strategy)
    }

    /// Validate paired-end mate alignment constraints.
    ///
    /// Returns `Ok(())` if the mates satisfy the configured constraints.
    /// Returns `Err(reason)` if validation fails.
    pub fn validate_mate_pair(
        &self,
        mate1_pos: usize,
        mate1_reverse: bool,
        mate2_pos: usize,
        mate2_reverse: bool,
    ) -> Result<(), String> {
        if let Some(config) = &self.paired_config {
            if config.require_both_mates && (mate1_pos == 0 || mate2_pos == 0) {
                return Err("require_both_mates: one or both mates unmapped".to_string());
            }

            // Expected insert size window: mean ± 3σ
            let mean = config.insert_size_mean;
            let sigma = config.insert_size_stddev;
            let max_insert = mean + 3 * sigma;
            let min_insert = if mean > 3 * sigma {
                mean - 3 * sigma
            } else {
                0
            };

            let insert_size = if mate2_pos > mate1_pos {
                mate2_pos - mate1_pos
            } else {
                mate1_pos - mate2_pos
            };

            if insert_size > max_insert || insert_size < min_insert {
                return Err(format!(
                    "insert size {} outside window [{}, {}]",
                    insert_size, min_insert, max_insert
                ));
            }

            if config.require_proper_orientation {
                // Proper pair: F1R2 (first on fwd, second on rev) or R1F2 (vice versa)
                let f1r2 = !mate1_reverse && mate2_reverse;
                let r1f2 = mate1_reverse && !mate2_reverse;
                if !f1r2 && !r1f2 {
                    return Err("improper pair orientation".to_string());
                }
            }
        }

        Ok(())
    }
}

impl Default for ShortReadAligner {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_short_read_aligner_single_end() {
        let aligner = ShortReadAligner::new();
        assert_eq!(aligner.strategy, Strategy::Fast);
        assert!(aligner.paired_config.is_none());
    }

    #[test]
    fn test_short_read_aligner_paired_end() {
        let config = PairedEndConfig::default();
        let aligner = ShortReadAligner::paired(config);
        assert_eq!(aligner.strategy, Strategy::Fast);
        assert!(aligner.paired_config.is_some());
    }

    #[test]
    fn test_validate_mate_pair_proper_orientation() {
        let config = PairedEndConfig::default();
        let aligner = ShortReadAligner::paired(config);

        // Proper pair: F1R2
        let result = aligner.validate_mate_pair(100, false, 450, true);
        assert!(result.is_ok());

        // Improper pair: both forward
        let result = aligner.validate_mate_pair(100, false, 450, false);
        assert!(result.is_err());
    }

    #[test]
    fn test_validate_mate_pair_insert_size() {
        let config = PairedEndConfig {
            insert_size_mean: 350,
            insert_size_stddev: 35,
            require_both_mates: true,
            require_proper_orientation: false,
        };
        let aligner = ShortReadAligner::paired(config);

        // Valid insert size (within ±3σ = [245, 455])
        let result = aligner.validate_mate_pair(100, false, 400, true);
        assert!(result.is_ok());

        // Insert size too large
        let result = aligner.validate_mate_pair(100, false, 600, true);
        assert!(result.is_err());
    }
}
