/// A single minimizer with its hash, position, and strand information
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub struct Minimizer {
    hash: u64,
    position: usize,
    strand: bool, // true = forward, false = reverse
}

impl Minimizer {
    fn new(hash: u64, position: usize, strand: bool) -> Self {
        Minimizer {
            hash,
            position,
            strand,
        }
    }

    pub fn hash(&self) -> u64 {
        self.hash
    }

    pub fn position(&self) -> usize {
        self.position
    }

    pub fn strand(&self) -> bool {
        self.strand
    }
}

/// Represents a shared minimizer between two sketches (a seed candidate)
#[derive(Clone, Debug)]
pub struct SharedMinimizer {
    hash: u64,
    pos_query: usize,
    pos_target: usize,
}

impl SharedMinimizer {
    pub fn new(hash: u64, pos_query: usize, pos_target: usize) -> Self {
        SharedMinimizer {
            hash,
            pos_query,
            pos_target,
        }
    }

    pub fn hash(&self) -> u64 {
        self.hash
    }

    pub fn pos_in_query(&self) -> usize {
        self.pos_query
    }

    pub fn pos_in_target(&self) -> usize {
        self.pos_target
    }
}

/// A sketch of a sequence using minimizers
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MinimimizerSketch {
    minimizers: Vec<Minimizer>,
    k: usize,
    w: usize,
}

impl MinimimizerSketch {
    pub fn new() -> Self {
        MinimimizerSketch {
            minimizers: Vec::new(),
            k: 21,
            w: 11,
        }
    }

    pub fn minimizers(&self) -> &[Minimizer] {
        &self.minimizers
    }

    pub fn k(&self) -> usize {
        self.k
    }

    pub fn w(&self) -> usize {
        self.w
    }

    /// Find shared minimizers between this sketch and another
    pub fn find_shared_minimizers(&self, other: &MinimimizerSketch) -> Vec<SharedMinimizer> {
        let mut result = Vec::new();

        // Build a set of hashes in the other sketch for fast lookup
        let other_hashes: std::collections::HashSet<u64> = other.minimizers.iter().map(|m| m.hash).collect();

        // For each minimizer in self, check if its hash exists in other
        for m in &self.minimizers {
            if other_hashes.contains(&m.hash) {
                // Find the first minimizer in other with this hash
                if let Some(other_m) = other.minimizers.iter().find(|om| om.hash == m.hash) {
                    result.push(SharedMinimizer::new(m.hash, m.position, other_m.position));
                }
            }
        }

        result
    }
}

impl Default for MinimimizerSketch {
    fn default() -> Self {
        Self::new()
    }
}

/// Convert a DNA sequence to a numeric k-mer hash using rolling 2-bit encoding
/// This is fast and deterministic
fn kmer_to_hash(kmer: &[u8]) -> u64 {
    if kmer.len() > 31 {
        return 0; // Limit to 62 bits for sign safety
    }

    let mut hash: u64 = 0;
    for &byte in kmer {
        let bits = match byte {
            b'A' => 0u8,
            b'C' => 1u8,
            b'G' => 2u8,
            b'T' => 3u8,
            _ => 0u8, // N and unknowns default to A
        };
        hash = hash << 2 | (bits as u64);
    }
    hash
}

/// Generate all k-mers and their minimizers from a sequence
/// Returns a sketch with minimizers selected by the minimizer selection algorithm
pub fn sketch(sequence: &[u8], k: usize, w: usize) -> MinimimizerSketch {
    // Validate parameters
    assert!(k > 0, "k must be greater than 0");
    assert!(w > 0, "w must be greater than 0");
    assert!(w <= k, "w must be less than or equal to k");

    let mut minimizers = Vec::new();

    // If sequence is too short, return empty sketch
    if sequence.len() < k {
        return MinimimizerSketch {
            minimizers,
            k,
            w,
        };
    }

    // Extract all k-mers and compute their hashes
    let mut kmers: Vec<(u64, usize)> = Vec::new();
    for i in 0..=(sequence.len() - k) {
        let kmer = &sequence[i..i + k];
        let hash = kmer_to_hash(kmer);
        kmers.push((hash, i));
    }

    if kmers.is_empty() {
        return MinimimizerSketch {
            minimizers,
            k,
            w,
        };
    }

    // Special case: only one k-mer
    if kmers.len() == 1 {
        let (hash, pos) = kmers[0];
        minimizers.push(Minimizer::new(hash, pos, true));
        return MinimimizerSketch {
            minimizers,
            k,
            w,
        };
    }

    // Minimizer selection using monotonic deque for O(n) performance
    let window_size = w;

    let mut deque: Vec<(u64, usize)> = Vec::with_capacity(window_size);
    let mut deque_front = 0;
    let mut last_min_kmer_idx = None;
    let mut last_forced_kmer_idx = 0;
    // Force at least one minimizer every w k-mers to ensure sequence coverage
    // This guarantees that long regions of non-minimum hashes still produce minimizers
    let force_spacing = w;

    for (kmer_idx, &(hash, seq_pos)) in kmers.iter().enumerate() {
        // Remove elements that would fall outside our window
        while deque_front < deque.len() && deque[deque_front].1 + window_size <= kmer_idx {
            deque_front += 1;
        }

        // Maintain monotonic property: remove larger or equal elements
        while deque.len() > deque_front && deque[deque.len() - 1].0 >= hash {
            deque.pop();
        }

        deque.push((hash, kmer_idx));

        // Emit minimizer when window is full
        if kmer_idx >= window_size - 1 && deque_front < deque.len() {
            let (_, min_kmer_idx) = deque[deque_front];

            // Emit if:
            // 1. The minimum k-mer position changed, OR
            // 2. We haven't emitted in force_spacing k-mers
            let should_emit = last_min_kmer_idx != Some(min_kmer_idx) ||
                             kmer_idx >= last_forced_kmer_idx + force_spacing;

            if should_emit {
                let min_pos = kmers[min_kmer_idx].1;
                let min_hash = kmers[min_kmer_idx].0;
                minimizers.push(Minimizer::new(min_hash, min_pos, true));
                last_min_kmer_idx = Some(min_kmer_idx);
                last_forced_kmer_idx = kmer_idx;
            }
        }
    }

    MinimimizerSketch {
        minimizers,
        k,
        w,
    }
}

/// Create a sketch with default parameters (k=21, w=11)
pub fn sketch_default(sequence: &[u8]) -> MinimimizerSketch {
    sketch(sequence, 21, 11)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Test: MinimimizerSketch struct exists and can be constructed
    #[test]
    fn test_minimizer_sketch_struct_creation() {
        // This should fail with a compilation error - MinimimizerSketch doesn't exist yet
        let sketch = MinimimizerSketch::new();
        assert!(sketch.minimizers().is_empty());
    }

    // Test: Basic sketch construction with explicit k and w parameters
    #[test]
    fn test_sketch_construction_basic() {
        let sequence = b"ACGTACGTACGTACGT";
        let k = 5;
        let w = 3;

        // This should fail - sketch() function doesn't exist yet
        let sketch = sketch(sequence, k, w);

        // Verify sketch was created
        assert!(!sketch.minimizers().is_empty());
    }

    // Test: Sketch construction with default parameters (k=21, w=11)
    #[test]
    fn test_sketch_construction_default_params() {
        let sequence = b"ACGTACGTACGTACGTACGTACGTACGTACGT"; // 32bp

        // This should fail - sketch_default() function doesn't exist yet
        let sketch = sketch_default(sequence);

        // Verify default parameters were used
        assert_eq!(sketch.k(), 21);
        assert_eq!(sketch.w(), 11);
    }

    // Test: Same sequence produces identical sketch
    #[test]
    fn test_same_sequence_identical_sketch() {
        let sequence = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let k = 5;
        let w = 3;

        let sketch1 = sketch(sequence, k, w);
        let sketch2 = sketch(sequence, k, w);

        // This should fail - equality check doesn't exist yet
        assert_eq!(sketch1, sketch2);
    }

    // Test: Known sequences produce expected shared minimizers
    #[test]
    fn test_shared_minimizers_known_sequences() {
        // Two sequences with overlapping k-mers
        let seq1 = b"ACGTACGTACGTACGT";
        let seq2 = b"ACGTACGTACGTACGT"; // Identical for now
        let k = 5;
        let w = 3;

        let sketch1 = sketch(seq1, k, w);
        let sketch2 = sketch(seq2, k, w);

        // This should fail - find_shared_minimizers() doesn't exist yet
        let shared = sketch1.find_shared_minimizers(&sketch2);

        // Identical sequences should share all minimizers
        assert_eq!(shared.len(), sketch1.minimizers().len());
    }

    // Test: Different sequences have fewer shared minimizers
    #[test]
    fn test_shared_minimizers_different_sequences() {
        let seq1 = b"ACGTACGTACGTACGTACGTACGT";
        let seq2 = b"TGCATGCATGCATGCATGCATGCA";
        let k = 5;
        let w = 3;

        let sketch1 = sketch(seq1, k, w);
        let sketch2 = sketch(seq2, k, w);

        let shared = sketch1.find_shared_minimizers(&sketch2);

        // Different sequences should share fewer minimizers than identical sequences
        assert!(shared.len() < sketch1.minimizers().len());
    }

    // Test: Minimizers include position information
    #[test]
    fn test_minimizers_have_positions() {
        let sequence = b"ACGTACGTACGTACGT";
        let k = 5;
        let w = 3;

        let sketch = sketch(sequence, k, w);

        // This should fail - position() method doesn't exist yet
        for minimizer in sketch.minimizers() {
            let pos = minimizer.position();
            assert!(pos < sequence.len());
        }
    }

    // Test: Minimizer selection is deterministic
    #[test]
    fn test_minimizer_selection_deterministic() {
        let sequence = b"ACGTACGTACGTACGTACGTACGT";
        let k = 7;
        let w = 5;

        // Run multiple times
        let sketch1 = sketch(sequence, k, w);
        let sketch2 = sketch(sequence, k, w);
        let sketch3 = sketch(sequence, k, w);

        // All should be identical
        assert_eq!(sketch1, sketch2);
        assert_eq!(sketch2, sketch3);
    }

    // Test: Empty sequence handling
    #[test]
    fn test_empty_sequence() {
        let sequence = b"";
        let k = 5;
        let w = 3;

        let sketch = sketch(sequence, k, w);

        // Empty sequence should produce empty sketch
        assert!(sketch.minimizers().is_empty());
    }

    // Test: Sequence shorter than k
    #[test]
    fn test_sequence_shorter_than_k() {
        let sequence = b"ACGT";
        let k = 10;
        let w = 3;

        let sketch = sketch(sequence, k, w);

        // Sequence shorter than k should produce empty sketch
        assert!(sketch.minimizers().is_empty());
    }

    // Test: Sequence exactly k length
    #[test]
    fn test_sequence_exactly_k_length() {
        let sequence = b"ACGTACGT"; // 8bp
        let k = 8;
        let w = 3;

        let sketch = sketch(sequence, k, w);

        // Should have exactly one minimizer
        assert_eq!(sketch.minimizers().len(), 1);
    }

    // Test: Invalid k parameter (k=0)
    #[test]
    #[should_panic]
    fn test_invalid_k_zero() {
        let sequence = b"ACGTACGT";
        let k = 0;
        let w = 3;

        // This should panic - invalid k parameter
        sketch(sequence, k, w);
    }

    // Test: Invalid w parameter (w=0)
    #[test]
    #[should_panic]
    fn test_invalid_w_zero() {
        let sequence = b"ACGTACGT";
        let k = 5;
        let w = 0;

        // This should panic - invalid w parameter
        sketch(sequence, k, w);
    }

    // Test: w greater than k
    #[test]
    #[should_panic]
    fn test_w_greater_than_k() {
        let sequence = b"ACGTACGTACGTACGT";
        let k = 5;
        let w = 10;

        // This should panic - w cannot be greater than k
        sketch(sequence, k, w);
    }

    // Test: Non-ACGT characters in sequence
    #[test]
    fn test_non_acgt_characters() {
        let sequence = b"ACGTNACGT"; // N is ambiguous
        let k = 5;
        let w = 3;

        // This should handle ambiguous bases gracefully
        let sketch = sketch(sequence, k, w);

        // Should either skip N-containing k-mers or handle them
        assert!(!sketch.minimizers().is_empty());
    }

    // Test: Shared minimizers are seed candidates
    #[test]
    fn test_shared_minimizers_as_seed_candidates() {
        let seq1 = b"ACGTACGTACGTACGTACGT";
        let seq2 = b"ACGTACGTACGTACGTACGT";
        let k = 7;
        let w = 4;

        let sketch1 = sketch(seq1, k, w);
        let sketch2 = sketch(seq2, k, w);

        let shared = sketch1.find_shared_minimizers(&sketch2);

        // Each shared minimizer should be usable as a seed candidate
        for seed in &shared {
            // Should have positions in both sequences
            assert!(seed.pos_in_query() < seq1.len());
            assert!(seed.pos_in_target() < seq2.len());
        }
    }

    // Test: Minimizer hash/kmer value is stored
    #[test]
    fn test_minimizer_stores_kmer_value() {
        let sequence = b"ACGTACGTACGTACGT";
        let k = 5;
        let w = 3;

        let sketch = sketch(sequence, k, w);

        // Each minimizer should have a k-mer hash or value
        for minimizer in sketch.minimizers() {
            let _hash = minimizer.hash();
            // Hash should be non-zero (or we should have explicit zero hashes)
            // Just checking it exists for now
        }
    }

    // Test: Minimizers are in order by position
    #[test]
    fn test_minimizers_ordered_by_position() {
        let sequence = b"ACGTACGTACGTACGTACGTACGT";
        let k = 7;
        let w = 5;

        let sketch = sketch(sequence, k, w);

        // Minimizers should be ordered by their position in sequence
        let positions: Vec<usize> = sketch.minimizers().iter().map(|m| m.position()).collect();

        let mut sorted_positions = positions.clone();
        sorted_positions.sort();

        assert_eq!(positions, sorted_positions);
    }

    // Test: Performance - sketch 5Mbp sequence in <1 second
    #[test]
    fn test_performance_5mbp_under_1_second() {
        use std::time::Instant;

        // Generate 5Mbp sequence
        let size = 5_000_000;
        let sequence: Vec<u8> = (0..size)
            .map(|i| match i % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            })
            .collect();

        let k = 21;
        let w = 11;

        let start = Instant::now();
        let _sketch = sketch(&sequence, k, w);
        let elapsed = start.elapsed();

        // Should complete in less than 1 second
        assert!(
            elapsed.as_secs() < 1,
            "Sketching took {:?}, expected < 1s",
            elapsed
        );
    }

    // Test: Sketch size is reasonable (not every k-mer)
    #[test]
    fn test_sketch_compression() {
        let size = 1000;
        let sequence: Vec<u8> = (0..size)
            .map(|i| match i % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            })
            .collect();

        let k = 21;
        let w = 11;

        let sketch = sketch(&sequence, k, w);

        // Number of minimizers should be much less than number of k-mers
        let num_kmers = size - k + 1;
        let num_minimizers = sketch.minimizers().len();

        // Sketch should be compressed - expect ~1/w minimizers
        assert!(num_minimizers < num_kmers / 2);
    }

    // Test: Reverse complement handling (if applicable)
    #[test]
    fn test_reverse_complement_awareness() {
        // A sequence and its reverse complement
        let seq = b"ACGTACGTACGT";
        let k = 5;
        let w = 3;

        let sketch = sketch(seq, k, w);

        // Check if minimizers include strand information
        // This will fail if strand() method doesn't exist
        for minimizer in sketch.minimizers() {
            let _strand = minimizer.strand();
            // Just checking the method exists
        }
    }

    // Test: Large window size (w close to k)
    #[test]
    fn test_large_window_size() {
        let sequence = b"ACGTACGTACGTACGTACGTACGT";
        let k = 10;
        let w = 9; // w close to k

        let sketch = sketch(sequence, k, w);

        // Should still work, but produce fewer minimizers
        assert!(!sketch.minimizers().is_empty());
    }

    // Test: Small window size (w=1)
    #[test]
    fn test_small_window_size() {
        let sequence = b"ACGTACGTACGTACGTACGTACGT";
        let k = 5;
        let w = 1;

        let sketch = sketch(sequence, k, w);

        // w=1 means every position gets a minimizer (most dense)
        assert!(!sketch.minimizers().is_empty());
    }

    // Test: Partial overlap between sequences
    #[test]
    fn test_partial_overlap_shared_minimizers() {
        // Two sequences with partial overlap
        let seq1 = b"ACGTACGTACGTACGTACGTACGT"; // prefix
        let seq2 = b"ACGTACGTACGTTTTTTTTTTTTT"; // shared prefix, different suffix
        let k = 7;
        let w = 5;

        let sketch1 = sketch(seq1, k, w);
        let sketch2 = sketch(seq2, k, w);

        let shared = sketch1.find_shared_minimizers(&sketch2);

        // Should find some shared minimizers from the overlapping region
        assert!(shared.len() > 0);
        // But not all minimizers should be shared
        assert!(shared.len() < sketch1.minimizers().len());
    }

    // Test: No shared minimizers between completely different sequences
    #[test]
    fn test_no_shared_minimizers_different_sequences() {
        let seq1 = b"AAAAAAAAAAAAAAAAAAAAAA";
        let seq2 = b"CCCCCCCCCCCCCCCCCCCCCC";
        let k = 5;
        let w = 3;

        let sketch1 = sketch(seq1, k, w);
        let sketch2 = sketch(seq2, k, w);

        let shared = sketch1.find_shared_minimizers(&sketch2);

        // Completely different sequences should share no minimizers
        assert_eq!(shared.len(), 0);
    }
}
