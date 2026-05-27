pub mod seeding;

pub use seeding::{MinimimizerSketch, Seed, find_seeds_from_sketches};

#[cfg(test)]
mod tests {
    use super::*;

    /// Test that find_seeds_from_sketches exists and returns Vec<Seed>
    #[test]
    fn test_find_seeds_from_sketches_function_exists() {
        let query = MinimimizerSketch::default();
        let target = MinimimizerSketch::default();
        let seeds = find_seeds_from_sketches(&query, &target);
        assert_eq!(seeds.len(), 0);
    }

    /// Test Seed struct has required fields
    #[test]
    fn test_seed_struct_fields() {
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
        let sequence = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let sketch1 = MinimimizerSketch::sketch(sequence, 21, 11);
        let sketch2 = MinimimizerSketch::sketch(sequence, 21, 11);
        let seeds = find_seeds_from_sketches(&sketch1, &sketch2);
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
        let query = b"AAAAAAAAAAAAAAAAAAAAAACGTACGTACGTACGTACGTACGTBBBBBBBBBBBBBBBBBBBBBB";
        let target = b"CCCCCCCCCCCCCCCCCCCCCCGTACGTACGTACGTACGTACGTDDDDDDDDDDDDDDDDDDDDDD";
        let query_sketch = MinimimizerSketch::sketch(query, 21, 11);
        let target_sketch = MinimimizerSketch::sketch(target, 21, 11);
        let seeds = find_seeds_from_sketches(&query_sketch, &target_sketch);
        assert!(!seeds.is_empty(), "Should find seeds in shared region");
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
        assert_eq!(
            seeds.len(),
            0,
            "Should find no seeds for different sequences"
        );
    }

    /// Test chain filtering: isolated seeds should be discarded
    #[test]
    fn test_chain_filtering_discards_isolated_seeds() {
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
        assert!(!seeds.is_empty(), "Should find seeds in colinear chains");

        let mut sorted_seeds = seeds.clone();
        sorted_seeds.sort_by_key(|s| s.query_pos);

        if sorted_seeds.len() > 1 {
            for window in sorted_seeds.windows(2) {
                let gap = window[1].query_pos.saturating_sub(window[0].query_pos);
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
        let query = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let target = b"ACGTACGTACGTACGGACGTACGTACGTACGTACGTACGTACGTACGTACGT";

        let query_sketch = MinimimizerSketch::sketch(query, 21, 11);
        let target_sketch = MinimimizerSketch::sketch(target, 21, 11);
        let seeds = find_seeds_from_sketches(&query_sketch, &target_sketch);
        assert!(!seeds.is_empty(), "Should find seeds despite single SNP");

        let has_seed_before_snp = seeds.iter().any(|s| s.query_pos < 16);
        let has_seed_after_snp = seeds.iter().any(|s| s.query_pos > 16 + 21);

        assert!(
            has_seed_before_snp || has_seed_after_snp,
            "Seeds should provide anchor points for WFA extension around variant"
        );

        for seed in &seeds {
            assert!(
                seed.minimizer != 0,
                "Seed minimizer should be non-zero (represents the shared k-mer)"
            );
        }
    }

    /// Test seeds with a longer sequence representing a realistic bacterial genome segment
    #[test]
    fn test_seeds_for_bacterial_genome_segment() {
        let mut query = Vec::with_capacity(10_000);
        let mut target = Vec::with_capacity(10_000);

        let bases = [b'A', b'C', b'G', b'T'];
        for i in 0..10_000 {
            let base = bases[i % 4];
            query.push(base);
            if i % 200 == 100 {
                target.push(bases[(i % 4 + 1) % 4]);
            } else {
                target.push(base);
            }
        }

        let query_sketch = MinimimizerSketch::sketch(&query, 21, 11);
        let target_sketch = MinimimizerSketch::sketch(&target, 21, 11);
        let seeds = find_seeds_from_sketches(&query_sketch, &target_sketch);
        assert!(
            seeds.len() >= 10,
            "Should find at least 10 seeds in 10Kbp with 99.5% identity, found {}",
            seeds.len()
        );

        let min_query_pos = seeds.iter().map(|s| s.query_pos).min().unwrap_or(0);
        let max_query_pos = seeds.iter().map(|s| s.query_pos).max().unwrap_or(0);

        assert!(
            max_query_pos - min_query_pos > 8_000,
            "Seeds should span most of the sequence (span: {})",
            max_query_pos - min_query_pos
        );
    }

    /// Test MinimimizerSketch integration
    #[test]
    fn test_minimizer_sketch_integration() {
        let sequence = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let sketch = MinimimizerSketch::sketch(sequence, 21, 11);
        assert!(
            sketch.minimizer_count() > 0,
            "Sketch should contain minimizers"
        );

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
        let mut query = Vec::with_capacity(5_000_000);
        let mut target = Vec::with_capacity(5_000_000);

        let bases = [b'A', b'C', b'G', b'T'];
        for i in 0..5_000_000 {
            let base = bases[(i / 1000) % 4];
            query.push(base);
            if i % 20 == 0 {
                target.push(bases[((i / 1000) % 4 + 1) % 4]);
            } else {
                target.push(base);
            }
        }

        let query_sketch = MinimimizerSketch::sketch(&query, 21, 11);
        let target_sketch = MinimimizerSketch::sketch(&target, 21, 11);

        let start = std::time::Instant::now();
        let seeds = find_seeds_from_sketches(&query_sketch, &target_sketch);
        let elapsed = start.elapsed();

        assert!(
            elapsed.as_millis() < 100,
            "Seed finding for 5Mbp sketches took {}ms, should be <100ms",
            elapsed.as_millis()
        );

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
        let sequence = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT";
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

        let _cloned = seed1.clone();
        let _debug = format!("{:?}", seed1);
        assert_eq!(seed1, seed2);
        assert_ne!(seed1, seed3);
    }

    /// Test chain filtering with inverted/non-colinear seeds
    #[test]
    fn test_chain_filtering_removes_inversions() {
        let query = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let target = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

        let query_sketch = MinimimizerSketch::sketch(query, 21, 11);
        let target_sketch = MinimimizerSketch::sketch(target, 21, 11);
        let seeds = find_seeds_from_sketches(&query_sketch, &target_sketch);

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
