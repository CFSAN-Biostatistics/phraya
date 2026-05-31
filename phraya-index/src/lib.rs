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

use std::collections::{HashMap, HashSet};

/// Compute Jaccard similarity between two sketches.
/// J(A, B) = |A ∩ B| / |A ∪ B|
fn jaccard_similarity(sketch_a: &MinimimizerSketch, sketch_b: &MinimimizerSketch) -> f64 {
    let set_a: HashSet<u64> = sketch_a.minimizers.iter().map(|&(val, _)| val).collect();
    let set_b: HashSet<u64> = sketch_b.minimizers.iter().map(|&(val, _)| val).collect();

    if set_a.is_empty() && set_b.is_empty() {
        return 1.0; // Two empty sketches are identical
    }

    let intersection = set_a.intersection(&set_b).count();
    let union = set_a.union(&set_b).count();

    if union == 0 {
        0.0
    } else {
        intersection as f64 / union as f64
    }
}

/// Select the centroid sketch (the one with median Jaccard similarity to all others).
/// Returns the index of the centroid sketch in the input slice.
///
/// # Arguments
/// * `sketches` - slice of sketches to select from
///
/// # Returns
/// Option<usize> - index of centroid, or None if empty
pub fn select_centroid(sketches: &[MinimimizerSketch]) -> Option<usize> {
    if sketches.is_empty() {
        return None;
    }

    if sketches.len() == 1 {
        return Some(0);
    }

    // Compute average Jaccard similarity for each sketch to all others
    let mut avg_similarities: Vec<f64> = Vec::new();

    for (i, sketch_i) in sketches.iter().enumerate() {
        let mut similarities = Vec::new();

        for (j, sketch_j) in sketches.iter().enumerate() {
            if i != j {
                let similarity = jaccard_similarity(sketch_i, sketch_j);
                similarities.push(similarity);
            }
        }

        // Median similarity for this sketch
        if !similarities.is_empty() {
            similarities.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let median = if similarities.len() % 2 == 0 {
                (similarities[similarities.len() / 2 - 1] + similarities[similarities.len() / 2]) / 2.0
            } else {
                similarities[similarities.len() / 2]
            };
            avg_similarities.push(median);
        } else {
            avg_similarities.push(0.0);
        }
    }

    // Find the index with the median value (middle of the distribution)
    let mut indexed_sims: Vec<(usize, f64)> = avg_similarities
        .into_iter()
        .enumerate()
        .collect();
    indexed_sims.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    let median_index = indexed_sims.len() / 2;
    Some(indexed_sims[median_index].0)
}

/// Seed from shared minimizer between query and target
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct Seed {
    pub query_pos: u32,
    pub target_pos: u32,
    pub minimizer: u64,
}

/// Find shared minimizers (seeds) between query and target sketches.
///
/// Seeds are anchor points where query and target share the same k-mer.
/// Sorted by query position.
pub fn find_seeds(query_sketch: &MinimimizerSketch, target_sketch: &MinimimizerSketch) -> Vec<Seed> {
    let shared = query_sketch.find_shared_minimizers(target_sketch);

    shared
        .into_iter()
        .map(|(minimizer, query_pos, target_pos)| Seed {
            query_pos: query_pos as u32,
            target_pos: target_pos as u32,
            minimizer,
        })
        .collect()
}

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

    #[test]
    fn jaccard_identical_sketches() {
        let seq = Sequence::new(b"ACGTACGTACGTACGTACGT".to_vec(), None, "seq".to_string(), None);
        let sketch = sketch_sequence(&seq, 4, 2);

        let j = jaccard_similarity(&sketch, &sketch);
        assert_eq!(j, 1.0);
    }

    #[test]
    fn jaccard_empty_sketches() {
        let empty_sketch = MinimimizerSketch {
            minimizers: vec![],
            k: 4,
            w: 2,
        };

        let j = jaccard_similarity(&empty_sketch, &empty_sketch);
        assert_eq!(j, 1.0);
    }

    #[test]
    fn jaccard_symmetry() {
        let seq1 = Sequence::new(b"ACGTACGTACGTACGTACGT".to_vec(), None, "seq1".to_string(), None);
        let seq2 = Sequence::new(b"TGCATGCATGCATGCATGCA".to_vec(), None, "seq2".to_string(), None);

        let sketch1 = sketch_sequence(&seq1, 4, 2);
        let sketch2 = sketch_sequence(&seq2, 4, 2);

        let j_ab = jaccard_similarity(&sketch1, &sketch2);
        let j_ba = jaccard_similarity(&sketch2, &sketch1);

        assert!((j_ab - j_ba).abs() < 1e-10);
    }

    #[test]
    fn jaccard_range() {
        let seq1 = Sequence::new(b"ACGTACGTACGTACGTACGT".to_vec(), None, "seq1".to_string(), None);
        let seq2 = Sequence::new(b"TGCATGCATGCATGCATGCA".to_vec(), None, "seq2".to_string(), None);

        let sketch1 = sketch_sequence(&seq1, 4, 2);
        let sketch2 = sketch_sequence(&seq2, 4, 2);

        let j = jaccard_similarity(&sketch1, &sketch2);
        assert!(j >= 0.0);
        assert!(j <= 1.0);
    }

    #[test]
    fn select_centroid_single_sketch() {
        let seq = Sequence::new(b"ACGTACGTACGTACGTACGT".to_vec(), None, "seq".to_string(), None);
        let sketch = sketch_sequence(&seq, 4, 2);

        let centroid = select_centroid(&[sketch]);
        assert_eq!(centroid, Some(0));
    }

    #[test]
    fn select_centroid_two_identical_sketches() {
        let seq = Sequence::new(b"ACGTACGTACGTACGTACGT".to_vec(), None, "seq".to_string(), None);
        let sketch = sketch_sequence(&seq, 4, 2);

        let centroid = select_centroid(&[sketch.clone(), sketch.clone()]);
        assert!(centroid.is_some());
        assert!(centroid.unwrap() == 0 || centroid.unwrap() == 1);
    }

    #[test]
    fn select_centroid_three_sketches() {
        let seq1 = Sequence::new(b"ACGTACGTACGTACGTACGT".to_vec(), None, "seq1".to_string(), None);
        let seq2 = Sequence::new(b"ACGTACGTACGTACGTACGT".to_vec(), None, "seq2".to_string(), None); // identical
        let seq3 = Sequence::new(b"TTTTTTTTTTTTTTTTTTTT".to_vec(), None, "seq3".to_string(), None); // different

        let sketch1 = sketch_sequence(&seq1, 4, 2);
        let sketch2 = sketch_sequence(&seq2, 4, 2);
        let sketch3 = sketch_sequence(&seq3, 4, 2);

        let centroid = select_centroid(&[sketch1, sketch2, sketch3]);
        assert!(centroid.is_some());
        // Centroid should be one of the first two (more similar to each other)
    }

    #[test]
    fn select_centroid_empty() {
        let centroid = select_centroid(&[]);
        assert_eq!(centroid, None);
    }

    #[test]
    fn select_centroid_returns_valid_index() {
        let seq1 = Sequence::new(b"ACGTACGTACGTACGTACGTACGT".to_vec(), None, "seq1".to_string(), None);
        let seq2 = Sequence::new(b"TGCATGCATGCATGCATGCATGCA".to_vec(), None, "seq2".to_string(), None);
        let seq3 = Sequence::new(b"GCTAGCTAGCTAGCTAGCTAGCTA".to_vec(), None, "seq3".to_string(), None);
        let seq4 = Sequence::new(b"ACGTACGTACGTACGTACGTACGT".to_vec(), None, "seq4".to_string(), None); // similar to seq1

        let sketch1 = sketch_sequence(&seq1, 5, 3);
        let sketch2 = sketch_sequence(&seq2, 5, 3);
        let sketch3 = sketch_sequence(&seq3, 5, 3);
        let sketch4 = sketch_sequence(&seq4, 5, 3);

        let sketches = vec![sketch1, sketch2, sketch3, sketch4];
        let centroid = select_centroid(&sketches);

        assert!(centroid.is_some());
        let idx = centroid.unwrap();
        assert!(idx < sketches.len());
    }

    #[test]
    fn select_centroid_many_identical() {
        let seq = Sequence::new(b"ACGTACGTACGTACGTACGTACGTACGT".to_vec(), None, "seq".to_string(), None);
        let sketch = sketch_sequence(&seq, 5, 3);

        let sketches = vec![sketch.clone(), sketch.clone(), sketch.clone(), sketch.clone(), sketch.clone()];
        let centroid = select_centroid(&sketches);

        assert!(centroid.is_some());
        // With all identical sketches, any index is valid
        assert!(centroid.unwrap() < sketches.len());
    }

    #[test]
    fn find_seeds_identical_sequences() {
        let seq = Sequence::new(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(), None, "seq".to_string(), None);
        let sketch1 = sketch_sequence(&seq, 5, 3);
        let sketch2 = sketch_sequence(&seq, 5, 3);

        let seeds = find_seeds(&sketch1, &sketch2);

        // Identical sequences should find shared minimizers as seeds
        // We're comparing two identical sketches, so we should find some shared minimizers
        // For identical sequences, seeds should be numerous or comprehensive
        assert!(seeds.len() > 0 || sketch1.minimizers.is_empty());
    }

    #[test]
    fn find_seeds_different_sequences() {
        let seq1 = Sequence::new(b"ACGTACGTACGTACGTACGTACGTACGT".to_vec(), None, "seq1".to_string(), None);
        let seq2 = Sequence::new(b"TGCATGCATGCATGCATGCATGCATGCA".to_vec(), None, "seq2".to_string(), None);

        let sketch1 = sketch_sequence(&seq1, 5, 3);
        let sketch2 = sketch_sequence(&seq2, 5, 3);

        let seeds = find_seeds(&sketch1, &sketch2);

        // Different sequences may have few or no seeds
        assert!(seeds.len() <= sketch1.minimizers.len());
    }

    #[test]
    fn find_seeds_empty_sketches() {
        let empty1 = MinimimizerSketch {
            minimizers: vec![],
            k: 5,
            w: 3,
        };
        let empty2 = MinimimizerSketch {
            minimizers: vec![],
            k: 5,
            w: 3,
        };

        let seeds = find_seeds(&empty1, &empty2);
        assert!(seeds.is_empty());
    }

    #[test]
    fn find_seeds_sorted_by_query_pos() {
        let seq = Sequence::new(b"ACGTACGTACGTACGTACGTACGTACGT".to_vec(), None, "seq".to_string(), None);
        let sketch = sketch_sequence(&seq, 5, 3);

        let seeds = find_seeds(&sketch, &sketch);

        // Seeds should be sorted by query position
        for window in seeds.windows(2) {
            if let [a, b] = window {
                assert!(a.query_pos <= b.query_pos);
            }
        }
    }

    #[test]
    fn find_seeds_matching_minimizers() {
        let seq = Sequence::new(b"ACGTACGTACGTACGTACGTACGTACGT".to_vec(), None, "seq".to_string(), None);
        let sketch = sketch_sequence(&seq, 5, 3);

        let seeds = find_seeds(&sketch, &sketch);

        // All seeds should have minimizers present in both sketches
        let query_minimizers: std::collections::HashSet<_> = sketch.minimizers.iter().map(|&(m, _)| m).collect();
        for seed in &seeds {
            assert!(query_minimizers.contains(&seed.minimizer));
        }
    }

    #[test]
    fn find_seeds_single_vs_empty() {
        let seq = Sequence::new(b"ACGTACGTACGTACGTACGTACGTACGT".to_vec(), None, "seq".to_string(), None);
        let sketch_nonempty = sketch_sequence(&seq, 5, 3);
        let sketch_empty = MinimimizerSketch {
            minimizers: vec![],
            k: 5,
            w: 3,
        };

        let seeds1 = find_seeds(&sketch_nonempty, &sketch_empty);
        let seeds2 = find_seeds(&sketch_empty, &sketch_nonempty);

        assert!(seeds1.is_empty());
        assert!(seeds2.is_empty());
    }

    #[test]
    fn find_seeds_seed_positions_valid() {
        let seq = Sequence::new(b"ACGTACGTACGTACGTACGTACGTACGT".to_vec(), None, "seq".to_string(), None);
        let sketch = sketch_sequence(&seq, 5, 3);

        let seeds = find_seeds(&sketch, &sketch);

        // Seed positions should be within sequence bounds
        for seed in &seeds {
            assert!(seed.query_pos < seq.len() as u32);
            assert!(seed.target_pos < seq.len() as u32);
        }
    }
}
