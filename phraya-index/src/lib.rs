pub mod minimizer;

pub use minimizer::{MinimimizerSketch};
use phraya_core::types::Sequence;

/// Default k-mer length for minimizer sketching (standard for bacterial genomics)
pub const DEFAULT_K: usize = 21;

/// Default window length for minimizer sketching
pub const DEFAULT_W: usize = 11;

/// Sketch a Sequence using given k-mer and window lengths.
pub fn sketch_sequence(sequence: &Sequence, k: usize, w: usize) -> MinimimizerSketch {
    minimizer::sketch(sequence.bases(), k, w)
}

/// Sketch a Sequence using default parameters (k=21, w=11).
pub fn sketch_sequence_default(sequence: &Sequence) -> MinimimizerSketch {
    sketch_sequence(sequence, DEFAULT_K, DEFAULT_W)
}

// Also export the low-level sketch function for testing
pub use minimizer::sketch;

use std::collections::HashMap;

/// Compute k-mer uniqueness scores from multiple sketches.
///
/// For each k-mer position in the reference, counts how many sketches contain that k-mer,
/// then assigns a uniqueness score of 1.0 / occurrence_count.
///
/// # Arguments
/// * `sketches` - slice of Sketches (typically one per sequence in a population)
///
/// # Returns
/// HashMap mapping reference position → uniqueness score in range (0.0, 1.0]
pub fn compute_kmer_uniqueness(sketches: &[MinimimizerSketch]) -> HashMap<u32, f64> {
    if sketches.is_empty() {
        return HashMap::new();
    }

    // Count occurrences of each (kmer_value, position) pair across all sketches
    let mut kmer_counts: HashMap<(u64, u32), usize> = HashMap::new();

    for sketch in sketches {
        for &(kmer_val, pos) in &sketch.minimizers {
            *kmer_counts.entry((kmer_val, pos as u32)).or_insert(0) += 1;
        }
    }

    // Convert counts to uniqueness scores
    let mut uniqueness: HashMap<u32, f64> = HashMap::new();
    for ((_, pos), count) in kmer_counts {
        let score = 1.0 / count as f64;
        // Store the minimum uniqueness score if multiple k-mers at same position
        uniqueness.entry(pos)
            .and_modify(|s| *s = s.min(score))
            .or_insert(score);
    }

    uniqueness
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sketch_sequence_empty() {
        let seq = Sequence::new(vec![], None, "empty".to_string(), None);
        let sketch = sketch_sequence(&seq, 4, 2);
        assert!(sketch.is_empty());
    }

    #[test]
    fn sketch_sequence_short_sequence() {
        let seq = Sequence::new(b"ACG".to_vec(), None, "short".to_string(), None);
        let sketch = sketch_sequence(&seq, 4, 2);
        // Sequence shorter than k-mer should be handled gracefully
        assert!(sketch.len() <= 1);
    }

    #[test]
    fn sketch_sequence_with_defaults() {
        let seq = Sequence::new(b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(), None, "seq".to_string(), None);
        let sketch = sketch_sequence_default(&seq);
        assert_eq!(sketch.k, DEFAULT_K);
        assert_eq!(sketch.w, DEFAULT_W);
    }

    #[test]
    fn sketch_sequence_deterministic() {
        let seq = Sequence::new(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(), None, "seq".to_string(), None);

        let sketch1 = sketch_sequence_default(&seq);
        let sketch2 = sketch_sequence_default(&seq);
        let sketch3 = sketch_sequence_default(&seq);

        assert_eq!(sketch1.minimizers, sketch2.minimizers);
        assert_eq!(sketch2.minimizers, sketch3.minimizers);
    }

    #[test]
    fn sketch_sequence_different_sequences() {
        let seq1 = Sequence::new(b"AAAAAAAAAA".to_vec(), None, "seq1".to_string(), None);
        let seq2 = Sequence::new(b"CCCCCCCCCC".to_vec(), None, "seq2".to_string(), None);

        let sketch1 = sketch_sequence_default(&seq1);
        let sketch2 = sketch_sequence_default(&seq2);

        // Different sequences should likely produce different sketches
        // (or at least this shouldn't panic)
        assert!(sketch1.minimizers != sketch2.minimizers || (sketch1.is_empty() && sketch2.is_empty()));
    }

    #[test]
    fn sketch_sequence_with_quality_scores() {
        let seq = Sequence::new(
            b"ACGTACGTACGT".to_vec(),
            Some(vec![30, 35, 40, 38, 30, 35, 40, 38, 30, 35, 40, 38]),
            "seq_with_qual".to_string(),
            None
        );

        let sketch = sketch_sequence(&seq, 4, 2);
        assert_eq!(sketch.k, 4);
        assert_eq!(sketch.w, 2);
        // Sketch should work regardless of quality scores
    }

    #[test]
    fn sketch_sequence_with_metadata() {
        let seq = Sequence::new(
            b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
            None,
            "sequence_1".to_string(),
            Some("A test sequence".to_string())
        );

        let sketch = sketch_sequence_default(&seq);
        assert_eq!(sketch.k, DEFAULT_K);
        assert_eq!(sketch.w, DEFAULT_W);
        // Metadata should not affect sketching
    }

    #[test]
    fn compute_uniqueness_empty_sketches() {
        let sketches: Vec<MinimimizerSketch> = vec![];
        let uniqueness = compute_kmer_uniqueness(&sketches);
        assert!(uniqueness.is_empty());
    }

    #[test]
    fn compute_uniqueness_single_sketch() {
        let seq = Sequence::new(b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(), None, "seq".to_string(), None);
        let sketch = sketch_sequence(&seq, 4, 2);

        let sketches = vec![sketch.clone()];
        let uniqueness = compute_kmer_uniqueness(&sketches);

        // All k-mers in a single sketch are unique (appear once)
        for &score in uniqueness.values() {
            assert_eq!(score, 1.0);
        }
    }

    #[test]
    fn compute_uniqueness_identical_sketches() {
        let seq = Sequence::new(b"ACGTACGTACGTACGTACGT".to_vec(), None, "seq".to_string(), None);
        let sketch = sketch_sequence(&seq, 4, 2);

        let sketches = vec![sketch.clone(), sketch.clone()];
        let uniqueness = compute_kmer_uniqueness(&sketches);

        // All k-mers appear exactly twice
        for &score in uniqueness.values() {
            assert_eq!(score, 0.5);
        }
    }

    #[test]
    fn compute_uniqueness_three_identical_sketches() {
        let seq = Sequence::new(b"ACGTACGTACGTACGTACGTACGT".to_vec(), None, "seq".to_string(), None);
        let sketch = sketch_sequence(&seq, 4, 2);

        let sketches = vec![sketch.clone(), sketch.clone(), sketch.clone()];
        let uniqueness = compute_kmer_uniqueness(&sketches);

        // All k-mers appear exactly 3 times
        for &score in uniqueness.values() {
            assert!((score - 1.0/3.0).abs() < 1e-10);
        }
    }

    #[test]
    fn compute_uniqueness_mixed_sketches() {
        let seq1 = Sequence::new(b"ACGTACGTACGTACGTACGT".to_vec(), None, "seq1".to_string(), None);
        let seq2 = Sequence::new(b"TGCATGCATGCATGCATGCA".to_vec(), None, "seq2".to_string(), None);

        let sketch1 = sketch_sequence(&seq1, 4, 2);
        let sketch2 = sketch_sequence(&seq2, 4, 2);
        let sketch1_2 = sketch_sequence(&seq1, 4, 2); // seq1 appears twice

        let sketches = vec![sketch1, sketch2, sketch1_2];
        let uniqueness = compute_kmer_uniqueness(&sketches);

        // k-mers from seq1 appear twice, k-mers from seq2 appear once
        assert!(!uniqueness.is_empty());
        // All scores should be positive and <= 1.0
        for &score in uniqueness.values() {
            assert!(score > 0.0);
            assert!(score <= 1.0);
        }
    }

    #[test]
    fn compute_uniqueness_scores_in_valid_range() {
        let seq = Sequence::new(b"ACGTACGTACGTACGTACGT".to_vec(), None, "seq".to_string(), None);
        let sketch = sketch_sequence(&seq, 4, 2);

        let sketches = vec![sketch.clone(), sketch.clone(), sketch.clone(), sketch.clone(), sketch.clone()];
        let uniqueness = compute_kmer_uniqueness(&sketches);

        // All scores should be in range (0.0, 1.0]
        for &score in uniqueness.values() {
            assert!(score > 0.0, "Score must be > 0.0");
            assert!(score <= 1.0, "Score must be <= 1.0");
            // Specifically, should be 1/5 = 0.2
            assert!((score - 0.2).abs() < 1e-10);
        }
    }

    #[test]
    fn compute_uniqueness_different_sequences() {
        let seq1 = Sequence::new(b"AAAAAAAAAA".to_vec(), None, "seq1".to_string(), None);
        let seq2 = Sequence::new(b"CCCCCCCCCC".to_vec(), None, "seq2".to_string(), None);
        let seq3 = Sequence::new(b"GGGGGGGGGG".to_vec(), None, "seq3".to_string(), None);

        let sketch1 = sketch_sequence(&seq1, 3, 2);
        let sketch2 = sketch_sequence(&seq2, 3, 2);
        let sketch3 = sketch_sequence(&seq3, 3, 2);

        let sketches = vec![sketch1, sketch2, sketch3];
        let uniqueness = compute_kmer_uniqueness(&sketches);

        // Each position should have entries, all unique (appear once)
        assert!(!uniqueness.is_empty());
        for &score in uniqueness.values() {
            assert!(score > 0.0);
            assert!(score <= 1.0);
        }
    }

    #[test]
    fn compute_uniqueness_property_all_scores_sum() {
        let seq = Sequence::new(b"ACGTACGTACGTACGTACGT".to_vec(), None, "seq".to_string(), None);
        let sketch = sketch_sequence(&seq, 4, 2);

        let sketches = vec![sketch.clone(), sketch.clone()];
        let uniqueness = compute_kmer_uniqueness(&sketches);

        // Count of unique positions should match sketch length
        assert!(uniqueness.len() <= sketch.minimizers.len());
    }
}
