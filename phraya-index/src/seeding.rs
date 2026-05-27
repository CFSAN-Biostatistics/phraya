// Sketch-based seeding for pairwise alignment
// Finds shared minimizers between two sketches and returns seed positions for WFA extension

#[cfg(test)]
mod tests {
    use super::*;

    /// Test that find_seeds_from_sketches exists and returns Vec<Seed>
    #[test]
    fn test_find_seeds_from_sketches_function_exists() {
        // This should fail with: cannot find function `find_seeds_from_sketches` in this scope
        let query = MinimimizerSketch::default();
        let target = MinimimizerSketch::default();
        let seeds = find_seeds_from_sketches(&query, &target);
        assert_eq!(seeds.len(), 0);
    }

    /// Test Seed struct has required fields
    #[test]
    fn test_seed_struct_fields() {
        // This should fail with: cannot find struct `Seed` in this scope
        let seed = Seed {
            query_pos: 100,
            target_pos: 200,
            minimizer: 0x123456789ABCDEF0,
        };
        assert_eq!(seed.query_pos, 100);
        assert_eq!(seed.target_pos, 200);
        assert_eq!(seed.minimizer, 0x123456789ABCDEF0);
    }

    /// Test that identical sequences produce seeds at matching positions
    #[test]
    fn test_identical_sequences_produce_colinear_seeds() {
        // Create a simple sequence and sketch it twice
        let sequence = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

        // Sketch with k=21, w=11 (default from issue #16)
        let sketch1 = MinimimizerSketch::sketch(sequence, 21, 11);
        let sketch2 = MinimimizerSketch::sketch(sequence, 21, 11);

        let seeds = find_seeds_from_sketches(&sketch1, &sketch2);

        // Should have seeds, and they should be diagonal (query_pos ≈ target_pos)
        assert!(
            !seeds.is_empty(),
            "Identical sequences should produce seeds"
        );

        for seed in &seeds {
            assert_eq!(
                seed.query_pos, seed.target_pos,
                "Identical sequences should have diagonal seeds"
            );
        }
    }

    /// Test that seeds are found for sequences with known shared k-mers
    #[test]
    fn test_seeds_for_known_shared_kmers() {
        // Two sequences with a known shared region
        let query = b"AAAAAAAAAAAAAAAAAAAAAACGTACGTACGTACGTACGTACGTBBBBBBBBBBBBBBBBBBBBBB";
        let target = b"CCCCCCCCCCCCCCCCCCCCCCGTACGTACGTACGTACGTACGTDDDDDDDDDDDDDDDDDDDDDD";

        let query_sketch = MinimimizerSketch::sketch(query, 21, 11);
        let target_sketch = MinimimizerSketch::sketch(target, 21, 11);

        let seeds = find_seeds_from_sketches(&query_sketch, &target_sketch);

        // Should find seeds in the shared "CGTACGTACGTACGTACGT" region
        assert!(!seeds.is_empty(), "Should find seeds in shared region");

        // Verify seeds are in the expected region (positions 21-44 in query, 21-44 in target)
        for seed in &seeds {
            assert!(
                seed.query_pos >= 21 && seed.query_pos <= 44,
                "Query seed position {} should be in shared region [21, 44]",
                seed.query_pos
            );
            assert!(
                seed.target_pos >= 21 && seed.target_pos <= 44,
                "Target seed position {} should be in shared region [21, 44]",
                seed.target_pos
            );
        }
    }

    /// Test that no seeds are found for completely different sequences
    #[test]
    fn test_no_seeds_for_different_sequences() {
        let query = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let target = b"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";

        let query_sketch = MinimimizerSketch::sketch(query, 21, 11);
        let target_sketch = MinimimizerSketch::sketch(target, 21, 11);

        let seeds = find_seeds_from_sketches(&query_sketch, &target_sketch);

        // Should find no seeds for completely different sequences
        assert_eq!(
            seeds.len(),
            0,
            "Should find no seeds for different sequences"
        );
    }

    /// Test chain filtering: isolated seeds should be discarded
    #[test]
    fn test_chain_filtering_discards_isolated_seeds() {
        // Create sequences where there's a small shared region far from a larger shared region
        // The isolated seed should be filtered out
        let query = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
                      TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\
                      GGCAGGCAGGCAGGCAGGCAGGCA\
                      AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\
                      CGATCGATCGATCGATCGATCGAT";

        let target = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
                       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\
                       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\
                       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\
                       CGATCGATCGATCGATCGATCGAT";

        let query_sketch = MinimimizerSketch::sketch(query, 21, 11);
        let target_sketch = MinimimizerSketch::sketch(target, 21, 11);

        let seeds = find_seeds_from_sketches(&query_sketch, &target_sketch);

        // Should have seeds from the two main colinear regions
        assert!(!seeds.is_empty(), "Should find seeds in colinear chains");

        // Seeds should form chains - check that seeds are reasonably close to each other
        // Sort by query position and verify gaps aren't too large
        let mut sorted_seeds = seeds.clone();
        sorted_seeds.sort_by_key(|s| s.query_pos);

        if sorted_seeds.len() > 1 {
            for window in sorted_seeds.windows(2) {
                let gap = window[1].query_pos.saturating_sub(window[0].query_pos);
                // Isolated seeds should be >50bp away from any chain
                // Chain seeds should be <50bp apart
                assert!(
                    gap < 50,
                    "Seeds should form chains without large gaps (gap: {})",
                    gap
                );
            }
        }
    }

    /// Test chain filtering: colinear chains should be kept
    #[test]
    fn test_chain_filtering_keeps_colinear_chains() {
        // Create sequences with a long shared region that should form a clear chain
        let shared = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let query = [
            b"TTTTTTTTTTTTTTTTTTTTTTTT".as_slice(),
            shared,
            b"AAAAAAAAAAAAAAAAAAAAAA",
        ]
        .concat();
        let target = [
            b"GGGGGGGGGGGGGGGGGGGGGGGG".as_slice(),
            shared,
            b"CCCCCCCCCCCCCCCCCCCCCC",
        ]
        .concat();

        let query_sketch = MinimimizerSketch::sketch(&query, 21, 11);
        let target_sketch = MinimimizerSketch::sketch(&target, 21, 11);

        let seeds = find_seeds_from_sketches(&query_sketch, &target_sketch);

        assert!(!seeds.is_empty(), "Should find seeds in shared region");

        // Verify seeds form a colinear chain (monotonically increasing in both query and target)
        let mut sorted_seeds = seeds.clone();
        sorted_seeds.sort_by_key(|s| s.query_pos);

        for window in sorted_seeds.windows(2) {
            assert!(
                window[1].query_pos > window[0].query_pos,
                "Seeds should be monotonically increasing in query"
            );
            assert!(
                window[1].target_pos > window[0].target_pos,
                "Seeds should be colinear (monotonically increasing in target)"
            );
        }
    }

    /// Test that seeds enable successful WFA extension
    #[test]
    fn test_seeds_enable_wfa_extension() {
        // Create two similar but not identical sequences
        let query = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let target = b"ACGTACGTACGTACGGACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        //                            ^ single SNP at position 16

        let query_sketch = MinimimizerSketch::sketch(query, 21, 11);
        let target_sketch = MinimimizerSketch::sketch(target, 21, 11);

        let seeds = find_seeds_from_sketches(&query_sketch, &target_sketch);

        assert!(!seeds.is_empty(), "Should find seeds despite single SNP");

        // Seeds should bracket the SNP position
        // We expect seeds before and/or after position 16
        let has_seed_before_snp = seeds.iter().any(|s| s.query_pos < 16);
        let has_seed_after_snp = seeds.iter().any(|s| s.query_pos > 16 + 21); // after the k-mer containing SNP

        assert!(
            has_seed_before_snp || has_seed_after_snp,
            "Seeds should provide anchor points for WFA extension around variant"
        );

        // Verify that the seeds have matching minimizers (this is the contract)
        for seed in &seeds {
            // This test will pass once the implementation ensures seeds have matching minimizers
            // For now, we just document the expectation
            assert!(
                seed.minimizer != 0,
                "Seed minimizer should be non-zero (represents the shared k-mer)"
            );
        }
    }

    /// Test seeds with a longer sequence representing a realistic bacterial genome segment
    #[test]
    fn test_seeds_for_bacterial_genome_segment() {
        // Simulate a 10Kbp segment with 99.5% identity (5 SNPs per kb)
        // For simplicity, we'll create sequences programmatically
        let mut query = Vec::with_capacity(10_000);
        let mut target = Vec::with_capacity(10_000);

        // Create a repetitive but not uniform pattern
        let bases = [b'A', b'C', b'G', b'T'];
        for i in 0..10_000 {
            let base = bases[i % 4];
            query.push(base);
            // Introduce a SNP every 200bp
            if i % 200 == 100 {
                target.push(bases[(i % 4 + 1) % 4]); // different base
            } else {
                target.push(base);
            }
        }

        let query_sketch = MinimimizerSketch::sketch(&query, 21, 11);
        let target_sketch = MinimimizerSketch::sketch(&target, 21, 11);

        let seeds = find_seeds_from_sketches(&query_sketch, &target_sketch);

        // Should find many seeds across the 10Kbp region
        assert!(
            seeds.len() >= 10,
            "Should find at least 10 seeds in 10Kbp with 99.5% identity, found {}",
            seeds.len()
        );

        // Seeds should span the full length of the sequence
        let min_query_pos = seeds.iter().map(|s| s.query_pos).min().unwrap_or(0);
        let max_query_pos = seeds.iter().map(|s| s.query_pos).max().unwrap_or(0);

        assert!(
            max_query_pos - min_query_pos > 8_000,
            "Seeds should span most of the sequence (span: {})",
            max_query_pos - min_query_pos
        );
    }

    /// Test MinimimizerSketch integration (assumes issue #16 implementation)
    #[test]
    fn test_minimizer_sketch_integration() {
        let sequence = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT";

        // Should be able to create a sketch with default parameters
        let sketch = MinimimizerSketch::sketch(sequence, 21, 11);

        // Sketch should have some minimizers
        assert!(
            sketch.minimizer_count() > 0,
            "Sketch should contain minimizers"
        );

        // Should be able to use sketch in find_seeds_from_sketches
        let sketch2 = MinimimizerSketch::sketch(sequence, 21, 11);
        let seeds = find_seeds_from_sketches(&sketch, &sketch2);

        assert!(
            !seeds.is_empty(),
            "Should find seeds from identical sketches"
        );
    }

    /// Performance test: seed finding for two 5Mbp sketches should complete in <100ms
    #[test]
    fn test_performance_5mbp_sketches() {
        // Create two 5Mbp sequences with 95% identity (realistic for bacterial genomes)
        let mut query = Vec::with_capacity(5_000_000);
        let mut target = Vec::with_capacity(5_000_000);

        let bases = [b'A', b'C', b'G', b'T'];
        for i in 0..5_000_000 {
            let base = bases[(i / 1000) % 4]; // creates 1Kbp blocks of each base
            query.push(base);
            // 5% divergence: change every 20th base
            if i % 20 == 0 {
                target.push(bases[((i / 1000) % 4 + 1) % 4]);
            } else {
                target.push(base);
            }
        }

        let query_sketch = MinimimizerSketch::sketch(&query, 21, 11);
        let target_sketch = MinimimizerSketch::sketch(&target, 21, 11);

        // Time the seed finding
        let start = std::time::Instant::now();
        let seeds = find_seeds_from_sketches(&query_sketch, &target_sketch);
        let elapsed = start.elapsed();

        // Should complete in <100ms
        assert!(
            elapsed.as_millis() < 100,
            "Seed finding for 5Mbp sketches took {}ms, should be <100ms",
            elapsed.as_millis()
        );

        // Should find a reasonable number of seeds
        assert!(
            seeds.len() > 100,
            "Should find many seeds in 5Mbp comparison, found {}",
            seeds.len()
        );
    }

    /// Test that seeds are sorted in a useful order (by query position)
    #[test]
    fn test_seeds_ordering() {
        let query = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let target = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

        let query_sketch = MinimimizerSketch::sketch(query, 21, 11);
        let target_sketch = MinimimizerSketch::sketch(target, 21, 11);

        let seeds = find_seeds_from_sketches(&query_sketch, &target_sketch);

        // Seeds should be sorted by query position for easier WFA extension
        for window in seeds.windows(2) {
            assert!(
                window[0].query_pos <= window[1].query_pos,
                "Seeds should be sorted by query position"
            );
        }
    }

    /// Test edge case: empty sketches
    #[test]
    fn test_empty_sketches() {
        let empty_query = MinimimizerSketch::default();
        let empty_target = MinimimizerSketch::default();

        let seeds = find_seeds_from_sketches(&empty_query, &empty_target);

        assert_eq!(seeds.len(), 0, "Empty sketches should produce no seeds");
    }

    /// Test edge case: one empty sketch
    #[test]
    fn test_one_empty_sketch() {
        let sequence = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let sketch = MinimimizerSketch::sketch(sequence, 21, 11);
        let empty = MinimimizerSketch::default();

        let seeds1 = find_seeds_from_sketches(&sketch, &empty);
        let seeds2 = find_seeds_from_sketches(&empty, &sketch);

        assert_eq!(seeds1.len(), 0, "Empty target should produce no seeds");
        assert_eq!(seeds2.len(), 0, "Empty query should produce no seeds");
    }

    /// Test that Seed struct implements common traits
    #[test]
    fn test_seed_traits() {
        let seed1 = Seed {
            query_pos: 100,
            target_pos: 200,
            minimizer: 0x123,
        };
        let seed2 = Seed {
            query_pos: 100,
            target_pos: 200,
            minimizer: 0x123,
        };
        let seed3 = Seed {
            query_pos: 150,
            target_pos: 250,
            minimizer: 0x456,
        };

        // Should implement Clone
        let _cloned = seed1.clone();

        // Should implement Debug
        let _debug = format!("{:?}", seed1);

        // Should implement PartialEq
        assert_eq!(seed1, seed2);
        assert_ne!(seed1, seed3);
    }

    /// Test chain filtering with inverted/non-colinear seeds
    #[test]
    fn test_chain_filtering_removes_inversions() {
        // This is a more advanced test for chain filtering
        // Manually create a scenario where we'd expect seeds but some are non-colinear
        // In real implementation, inversions would be filtered

        // For now, this test documents the expectation that the chain filtering
        // should remove seeds that don't maintain colinear order
        let query = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let target = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

        let query_sketch = MinimimizerSketch::sketch(query, 21, 11);
        let target_sketch = MinimimizerSketch::sketch(target, 21, 11);

        let seeds = find_seeds_from_sketches(&query_sketch, &target_sketch);

        // All seeds should maintain colinear order
        let mut sorted_seeds = seeds.clone();
        sorted_seeds.sort_by_key(|s| s.query_pos);

        let mut prev_target_pos = 0;
        for seed in sorted_seeds {
            assert!(
                seed.target_pos >= prev_target_pos,
                "Seeds should be colinear (no inversions)"
            );
            prev_target_pos = seed.target_pos;
        }
    }
}
