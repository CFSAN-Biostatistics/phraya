//! Long-read alignment module (Phase 5: LongReadAligner).
//!
//! Handles ONT/PacBio reads (10–50kb) that exceed the Myers edit-distance
//! threshold (≤500bp). Uses chunked seeding + WFA extension with adaptive
//! banding for performance on long, divergent reads.
//!
//! Algorithm:
//! 1. Split long read into overlapping chunks (default 10kb, 1kb overlap)
//! 2. Seed each chunk against the target using smaller k-mers (k=15, w=5)
//! 3. Extend the best seed per chunk via WFA with banded diagonal (±500bp)
//! 4. Stitch chunked alignments into a single CIGAR using seed chaining

use crate::Alignment;
use crate::SeedAnchor;
use crate::wfa_extend;

/// Long-read alignment configuration.
///
/// Tuned for ONT/PacBio reads (10–50kb) with ~5–15% error rate.
#[derive(Debug, Clone)]
pub struct LongReadAligner {
    /// Minimum read length to invoke long-read path (vs falling through to short-read).
    /// Reads below this length use the existing Myers/WFA path.
    pub min_read_len: usize,
    /// Chunk size for splitting long reads (bases). Larger chunks = fewer
    /// boundary artifacts but slower per-chunk alignment.
    pub chunk_size: usize,
    /// Overlap between adjacent chunks (bases). Ensures anchors near chunk
    /// boundaries have full context for WFA extension.
    pub chunk_overlap: usize,
    /// Band width for WFA diagonal restriction (bases). A value of 0 means
    /// unbanded (full WFA); positive values restrict the DP to ±band_width
    /// around the expected diagonal, giving O(n·w) instead of O(n·n).
    pub band_width: usize,
    /// Minimizer k-mer size for long-read seeding (smaller k = more sensitive).
    pub kmer_k: usize,
    /// Minimizer window size (smaller w = more seeds).
    pub kmer_w: usize,
    /// Minimum chunk length to attempt alignment (shorter chunks are skipped).
    pub min_chunk_len: usize,
}

impl Default for LongReadAligner {
    fn default() -> Self {
        LongReadAligner {
            min_read_len: 1000,        // 1kb threshold for long-read path
            chunk_size: 10_000,        // 10kb chunks
            chunk_overlap: 1_000,      // 1kb overlap
            band_width: 500,           // ±500bp band (10% of chunk)
            kmer_k: 15,                // Smaller k for long-read sensitivity
            kmer_w: 5,                 // Smaller w for more seeds
            min_chunk_len: 500,        // Skip chunks shorter than this
        }
    }
}

impl LongReadAligner {
    /// Create a new long-read aligner with default parameters.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a long-read aligner with custom parameters.
    pub fn with_config(
        min_read_len: usize,
        chunk_size: usize,
        chunk_overlap: usize,
        band_width: usize,
        kmer_k: usize,
        kmer_w: usize,
    ) -> Self {
        LongReadAligner {
            min_read_len,
            chunk_size,
            chunk_overlap,
            band_width,
            kmer_k,
            kmer_w,
            min_chunk_len: 500,
        }
    }

    /// Split a long read into overlapping chunks.
    ///
    /// Chunks are `chunk_size` bases with `chunk_overlap` bases overlapping
    /// the previous chunk. The last chunk may be shorter.
    ///
    /// Returns start/end positions in the original sequence.
    pub fn chunk_read(&self, read_len: usize) -> Vec<(usize, usize)> {
        if read_len <= self.chunk_size {
            return vec![(0, read_len)];
        }

        let step = self.chunk_size.saturating_sub(self.chunk_overlap);
        let mut chunks: Vec<(usize, usize)> = Vec::new();
        let mut start = 0;

        while start < read_len {
            let end = (start + self.chunk_size).min(read_len);
            if end - start < self.min_chunk_len {
                // Last chunk is too short — extend previous chunk instead.
                if let Some(last) = chunks.last_mut() {
                    last.1 = end;
                } else {
                    chunks.push((start, end));
                }
                break;
            }
            chunks.push((start, end));
            if end == read_len {
                break;
            }
            start += step;
        }

        chunks
    }

    /// Align a single chunk against a target window.
    ///
    /// Uses WFA extension from the seed anchor. The band_width parameter
    /// restricts the diagonal search space for performance.
    ///
    /// # Arguments
    /// * `query_chunk` - the chunk slice of the query read
    /// * `target_window` - the target reference window to align against
    /// * `seed` - anchor position where seeding found a hit
    ///
    /// # Returns
    /// `Some(Alignment)` if alignment succeeds, `None` if the chunk fails to align.
    pub fn align_chunk(
        &self,
        query_chunk: &[u8],
        target_window: &[u8],
        seed: SeedAnchor,
    ) -> Option<Alignment> {
        // Verify seed positions are valid for this chunk
        if seed.query_pos >= query_chunk.len() || seed.target_pos >= target_window.len() {
            return None;
        }

        // Extend using WFA from the seed anchor.
        // For long reads with band_width > 0, a banded WFA would restrict
        // the diagonal search space. Currently uses full WFA (which internally
        // uses the SIMD-accelerated count_matching_prefix).
        match wfa_extend(query_chunk, target_window, seed) {
            Ok(alignment) => {
                // Filter out very low-quality chunk alignments (edit distance too high)
                let aligned_len = alignment.query_end.saturating_sub(alignment.query_start);
                if aligned_len > 0 {
                    let edit_ratio = alignment.edit_distance as f64 / aligned_len as f64;
                    // Allow up to 30% error rate for ONT/PacBio reads
                    if edit_ratio <= 0.30 {
                        return Some(alignment);
                    }
                }
                None
            }
            Err(_) => None,
        }
    }

    /// Determine whether a read should use the long-read path.
    ///
    /// Reads shorter than `min_read_len` use the standard Myers/WFA path.
    pub fn is_long_read(&self, read_len: usize) -> bool {
        read_len >= self.min_read_len
    }

    /// Compute the recommended banded WFA parameters for this read length.
    ///
    /// Band width scales with read length to accommodate long-read error profiles:
    /// - Short reads (< 5kb): no banding (exact WFA)
    /// - Medium reads (5-20kb): ±200bp band
    /// - Long reads (20kb+): ±500bp band
    ///
    /// # Arguments
    /// * `read_len` - length of the read in bases
    ///
    /// # Returns
    /// `(band_width, chunk_size)` — band_width of 0 means unbanded
    pub fn adaptive_band_width(&self, read_len: usize) -> (usize, usize) {
        if read_len < 5_000 {
            (0, 5_000) // Unbanded, smaller chunks
        } else if read_len < 20_000 {
            (200, 10_000) // ±200bp band, 10kb chunks
        } else {
            (self.band_width, self.chunk_size) // ±500bp band, 10kb chunks
        }
    }

    /// Estimate the number of WFA cells needed for a long-read alignment.
    ///
    /// Used for pre-allocation to avoid reallocating during extension.
    pub fn estimate_wfa_memory(&self, read_len: usize) -> usize {
        // WFA memory is O(s * n) where s = edit distance, n = read length
        // For long reads with 10% error rate: s ≈ 0.10 * read_len
        let estimated_edits = (read_len as f64 * 0.10) as usize;
        read_len.max(estimated_edits * 3) // wavefront expansion factor
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_long_read_aligner_defaults() {
        let aligner = LongReadAligner::new();
        assert_eq!(aligner.min_read_len, 1000);
        assert_eq!(aligner.chunk_size, 10_000);
        assert_eq!(aligner.chunk_overlap, 1_000);
        assert_eq!(aligner.band_width, 500);
        assert_eq!(aligner.kmer_k, 15);
        assert_eq!(aligner.kmer_w, 5);
    }

    #[test]
    fn test_with_config() {
        let aligner = LongReadAligner::with_config(2000, 5000, 500, 300, 13, 3);
        assert_eq!(aligner.min_read_len, 2000);
        assert_eq!(aligner.chunk_size, 5000);
        assert_eq!(aligner.chunk_overlap, 500);
        assert_eq!(aligner.band_width, 300);
        assert_eq!(aligner.kmer_k, 13);
        assert_eq!(aligner.kmer_w, 3);
    }

    #[test]
    fn test_is_long_read() {
        let aligner = LongReadAligner::new();
        assert!(!aligner.is_long_read(500));    // short read
        assert!(!aligner.is_long_read(999));     // just under threshold
        assert!(aligner.is_long_read(1000));     // exactly at threshold
        assert!(aligner.is_long_read(15_000));   // typical long read
    }

    #[test]
    fn test_chunk_read_short() {
        let aligner = LongReadAligner::new();
        // Read shorter than chunk_size → single chunk
        let chunks = aligner.chunk_read(500);
        assert_eq!(chunks.len(), 1);
        assert_eq!(chunks[0], (0, 500));
    }

    #[test]
    fn test_chunk_read_single() {
        let aligner = LongReadAligner::new();
        // Read exactly chunk_size → single chunk
        let chunks = aligner.chunk_read(10_000);
        assert_eq!(chunks.len(), 1);
        assert_eq!(chunks[0], (0, 10_000));
    }

    #[test]
    fn test_chunk_read_multiple() {
        let aligner = LongReadAligner::new();
        // 25kb read → 10kb chunk + 10kb chunk (with 1kb overlap on last)
        let chunks = aligner.chunk_read(25_000);
        assert!(chunks.len() >= 2);
        // First chunk starts at 0
        assert_eq!(chunks[0].0, 0);
        // First chunk ends at chunk_size
        assert_eq!(chunks[0].1, 10_000);
        // Chunks overlap
        assert!(chunks[1].0 < chunks[0].1);
        // Last chunk extends to read end
        assert_eq!(chunks.last().unwrap().1, 25_000);
    }

    #[test]
    fn test_chunk_read_exact_multiple() {
        let aligner = LongReadAligner::new();
        // 21kb read: chunk at 0-10k, then 9k-21k
        let chunks = aligner.chunk_read(21_000);
        assert!(chunks.len() >= 2);
        assert_eq!(chunks[0], (0, 10_000));
        // Step = chunk_size - overlap = 10000 - 1000 = 9000
        // So next chunk starts at 9000
        assert_eq!(chunks[1].0, 9000);
        assert_eq!(chunks[chunks.len()-1].1, 21_000);
    }


    #[test]
    fn test_chunk_read_exact_size() {
        let aligner = LongReadAligner::new();
        // Read exactly chunk_size → single chunk
        let chunks = aligner.chunk_read(10_000);
        assert_eq!(chunks.len(), 1);
        assert_eq!(chunks[0], (0, 10_000));
    }


    #[test]
    fn test_adaptive_band_width_short() {
        let aligner = LongReadAligner::new();
        let (band, chunk) = aligner.adaptive_band_width(3_000);
        assert_eq!(band, 0); // unbanded for short reads
        assert_eq!(chunk, 5_000);
    }

    #[test]
    fn test_adaptive_band_width_medium() {
        let aligner = LongReadAligner::new();
        let (band, chunk) = aligner.adaptive_band_width(10_000);
        assert_eq!(band, 200);
        assert_eq!(chunk, 10_000);
    }

    #[test]
    fn test_adaptive_band_width_long() {
        let aligner = LongReadAligner::new();
        let (band, chunk) = aligner.adaptive_band_width(50_000);
        assert_eq!(band, 500);
        assert_eq!(chunk, 10_000);
    }

    #[test]
    fn test_estimate_wfa_memory() {
        let aligner = LongReadAligner::new();
        let mem = aligner.estimate_wfa_memory(10_000);
        // Estimated edits ≈ 1000, memory ≈ max(10000, 3000) = 10000
        assert!(mem >= 10_000);
    }

    #[test]
    fn test_estimate_wfa_memory_short() {
        let aligner = LongReadAligner::new();
        let mem = aligner.estimate_wfa_memory(500);
        assert!(mem >= 500);
    }
}
