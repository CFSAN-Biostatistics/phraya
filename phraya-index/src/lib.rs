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
}
