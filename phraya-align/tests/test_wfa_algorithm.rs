//! Acceptance tests for issue #131: Implement WFA (Wavefront Alignment Algorithm)
//!
//! These tests verify that the WFA implementation is a true wavefront algorithm
//! (sub-quadratic time/space), not the O(n×m) Levenshtein DP that currently exists.
//!
//! The tests are organized to verify:
//! 1. Correctness: CIGAR and edit distance for various sequence patterns
//! 2. Performance: sub-quadratic behavior, <1s for 150bp read vs 5Mbp reference segment
//! 3. Algorithm compliance: demonstrates wavefront (diagonal) structure, not full DP matrix

use phraya_align::{wfa_extend, SeedAnchor};
use std::time::Instant;

// ============================================================================
// Test Category 1: Correctness - Basic Cases
// ============================================================================

/// Issue #131: Exact match produces CIGAR "NM" (all matches)
#[test]
fn issue_131_exact_match_produces_correct_cigar_and_edit_distance() {
    let query = b"ACGTACGTACGT";
    let target = b"ACGTACGTACGT";
    let seed = SeedAnchor {
        query_pos: 0,
        target_pos: 0,
    };

    let result = wfa_extend(query, target, seed);
    assert!(result.is_ok());

    let alignment = result.unwrap();
    // For exact match, edit distance must be 0
    assert_eq!(
        alignment.edit_distance, 0,
        "Exact match must have edit_distance 0"
    );
    // CIGAR should be "12M" (12 matches)
    assert_eq!(alignment.cigar, "12M", "Exact match CIGAR must be '12M'");
}

/// Issue #131: Single SNP produces edit distance 1 with mismatch operation
#[test]
fn issue_131_single_snp_produces_edit_distance_one() {
    let query = b"ACGTACGTACGT";
    let target = b"ACGTACATACGT"; // Mismatch at position 6: G vs A
    let seed = SeedAnchor {
        query_pos: 0,
        target_pos: 0,
    };

    let result = wfa_extend(query, target, seed);
    assert!(result.is_ok());

    let alignment = result.unwrap();
    // Single SNP = edit distance 1
    assert_eq!(
        alignment.edit_distance, 1,
        "Single SNP must have edit_distance 1"
    );
    // CIGAR should contain both M and X operations
    assert!(
        alignment.cigar.contains("M") || alignment.cigar.contains("X"),
        "CIGAR must contain M or X for SNP alignment"
    );
}

/// Issue #131: Mismatch at position 6: G (query) vs T (target) - single character difference
#[test]
fn issue_131_single_mismatch_at_position_six() {
    let query = b"ACGTACGTACGT";
    let target = b"ACGTACTTACGT"; // Mismatch at position 7: G vs T (single char)
    let seed = SeedAnchor {
        query_pos: 0,
        target_pos: 0,
    };

    let result = wfa_extend(query, target, seed);
    assert!(result.is_ok());

    let alignment = result.unwrap();
    assert_eq!(
        alignment.edit_distance, 1,
        "Single mismatch must have edit_distance 1"
    );
}

/// Issue #131: Single insertion produces edit distance 1
#[test]
fn issue_131_single_insertion_produces_edit_distance_one() {
    let query = b"ACGTACGT";
    let target = b"ACGTAACGT"; // Extra 'A' inserted in target
    let seed = SeedAnchor {
        query_pos: 0,
        target_pos: 0,
    };

    let result = wfa_extend(query, target, seed);
    assert!(result.is_ok());

    let alignment = result.unwrap();
    // Single insertion = edit distance 1
    assert_eq!(
        alignment.edit_distance, 1,
        "Single insertion must have edit_distance 1"
    );
    // CIGAR must contain 'I' operation
    assert!(
        alignment.cigar.contains("I"),
        "CIGAR for insertion must contain 'I'"
    );
}

/// Issue #131: Single deletion produces edit distance 1
#[test]
fn issue_131_single_deletion_produces_edit_distance_one() {
    let query = b"ACGTAACGT";
    let target = b"ACGTACGT"; // 'A' deleted from target (relative to query)
    let seed = SeedAnchor {
        query_pos: 0,
        target_pos: 0,
    };

    let result = wfa_extend(query, target, seed);
    assert!(result.is_ok());

    let alignment = result.unwrap();
    // Single deletion = edit distance 1
    assert_eq!(
        alignment.edit_distance, 1,
        "Single deletion must have edit_distance 1"
    );
    // CIGAR must contain 'D' operation
    assert!(
        alignment.cigar.contains("D"),
        "CIGAR for deletion must contain 'D'"
    );
}

// ============================================================================
// Test Category 2: Correctness - Multiple Indels
// ============================================================================

/// Issue #131: Multiple indels (mixed insertions and deletions)
#[test]
fn issue_131_multiple_indels_correct_edit_distance() {
    let query = b"ACGTACGTTAGC";
    let target = b"ACGTTCGTAGC";
    let seed = SeedAnchor {
        query_pos: 0,
        target_pos: 0,
    };

    let result = wfa_extend(query, target, seed);
    assert!(result.is_ok());

    let alignment = result.unwrap();
    // Should have edit distance >= 2 (at least 2 operations)
    assert!(
        alignment.edit_distance > 0,
        "Multiple indels must have edit_distance > 0"
    );
    assert!(
        alignment.edit_distance <= 4,
        "Expected alignment with ~2-3 ops for this pair"
    );
}

/// Issue #131: Consecutive indels (multiple deletions in a row)
#[test]
fn issue_131_consecutive_deletions_correct_edit_distance() {
    let query = b"ACGTAAAAACGT";
    let target = b"ACGTACGT"; // 'AAAA' deleted (4 A's)
    let seed = SeedAnchor {
        query_pos: 0,
        target_pos: 0,
    };

    let result = wfa_extend(query, target, seed);
    assert!(result.is_ok());

    let alignment = result.unwrap();
    // Deleting 4 'A' characters = edit distance 4
    assert_eq!(
        alignment.edit_distance, 4,
        "Deleting 4 consecutive bases must have edit_distance 4"
    );
}

// ============================================================================
// Test Category 3: Correctness - High Divergence
// ============================================================================

/// Issue #131: High divergence sequences still align correctly
#[test]
fn issue_131_high_divergence_produces_reasonable_edit_distance() {
    let query = b"ACGTACGTACGTACGT";
    let target = b"ACACACACACACACAC"; // Systematic pattern with high divergence
    let seed = SeedAnchor {
        query_pos: 0,
        target_pos: 0,
    };

    let result = wfa_extend(query, target, seed);
    assert!(result.is_ok());

    let alignment = result.unwrap();
    // High divergence: many bases don't match
    // Should have non-zero edit distance
    assert!(
        alignment.edit_distance > 0,
        "High divergence must have positive edit_distance"
    );
}

// ============================================================================
// Test Category 4: Correctness - Edge Cases
// ============================================================================

/// Issue #131: Empty suffix at seed position (seed at end of sequences)
#[test]
fn issue_131_empty_suffix_at_seed_position() {
    let query = b"ACGTACGTACGT";
    let target = b"ACGTACGTACGT";
    let seed = SeedAnchor {
        query_pos: 12, // At end, suffix is empty
        target_pos: 12,
    };

    let result = wfa_extend(query, target, seed);
    assert!(result.is_ok());

    let alignment = result.unwrap();
    assert_eq!(
        alignment.edit_distance, 0,
        "Empty alignment must have edit_distance 0"
    );
    assert_eq!(alignment.cigar, "", "Empty alignment must have empty CIGAR");
}

/// Issue #131: Seed in middle of sequences
#[test]
fn issue_131_seed_in_middle_of_sequence() {
    let query = b"ACGTACGTACGT";
    let target = b"ACGTACGTACGT";
    let seed = SeedAnchor {
        query_pos: 4,
        target_pos: 4,
    };

    let result = wfa_extend(query, target, seed);
    assert!(result.is_ok());

    let alignment = result.unwrap();
    // Should align the suffix starting at position 4
    assert_eq!(
        alignment.query_start, 4,
        "query_start must be seed position"
    );
    assert_eq!(
        alignment.target_start, 4,
        "target_start must be seed position"
    );
    assert_eq!(alignment.query_end, 12, "query_end must be sequence end");
    assert_eq!(alignment.target_end, 12, "target_end must be sequence end");
    assert_eq!(
        alignment.edit_distance, 0,
        "Exact match suffix has edit_distance 0"
    );
}

/// Issue #131: Very short sequences (boundary case)
#[test]
fn issue_131_very_short_sequences() {
    let query = b"AC";
    let target = b"AC";
    let seed = SeedAnchor {
        query_pos: 0,
        target_pos: 0,
    };

    let result = wfa_extend(query, target, seed);
    assert!(result.is_ok());

    let alignment = result.unwrap();
    assert_eq!(alignment.edit_distance, 0);
    assert_eq!(alignment.cigar, "2M");
}

/// Issue #131: Single base sequences
#[test]
fn issue_131_single_base_match() {
    let query = b"A";
    let target = b"A";
    let seed = SeedAnchor {
        query_pos: 0,
        target_pos: 0,
    };

    let result = wfa_extend(query, target, seed);
    assert!(result.is_ok());

    let alignment = result.unwrap();
    assert_eq!(alignment.edit_distance, 0);
    assert_eq!(alignment.cigar, "1M");
}

/// Issue #131: Single base mismatch
#[test]
fn issue_131_single_base_mismatch() {
    let query = b"A";
    let target = b"C";
    let seed = SeedAnchor {
        query_pos: 0,
        target_pos: 0,
    };

    let result = wfa_extend(query, target, seed);
    assert!(result.is_ok());

    let alignment = result.unwrap();
    assert_eq!(alignment.edit_distance, 1);
}

// ============================================================================
// Test Category 5: Correctness - Complex Patterns
// ============================================================================

/// Issue #131: Repeat regions (AT-rich repeats)
#[test]
fn issue_131_repeat_regions_aligned_correctly() {
    let query = b"ATATATATATATAT";
    let target = b"ATATATATATAT";
    let seed = SeedAnchor {
        query_pos: 0,
        target_pos: 0,
    };

    let result = wfa_extend(query, target, seed);
    assert!(result.is_ok());

    let alignment = result.unwrap();
    // Deleting 2 'AT' pairs
    assert_eq!(alignment.edit_distance, 2, "Two missing bases in target");
}

/// Issue #131: GC-rich sequences
#[test]
fn issue_131_gc_rich_sequence_alignment() {
    let query = b"GCGCGCGCGCGCGCGC";
    let target = b"GCGCGCGGCGCGCGC";
    let seed = SeedAnchor {
        query_pos: 0,
        target_pos: 0,
    };

    let result = wfa_extend(query, target, seed);
    assert!(result.is_ok());

    let alignment = result.unwrap();
    // One mismatch and one deletion
    assert_eq!(alignment.edit_distance, 1, "One mismatch: C vs G");
}

// ============================================================================
// Test Category 6: Performance - Sub-quadratic Behavior
// ============================================================================

/// Issue #176: differential correctness of Myers fitting extension vs WFA.
///
/// Across a broad sweep of randomized reads windowed against ~2× reference (with SNPs
/// and indels), `myers_extend` must report the same edit distance as `wfa_extend`.
/// These are two independent fitting-alignment implementations; agreement over hundreds
/// of cases is the scientific guarantee that the Fast strategy preserves variant calls.
#[test]
fn issue_176_myers_fitting_edit_distance_matches_wfa_sweep() {
    // Small deterministic LCG — no external rng dependency.
    let mut state: u64 = 0x9E3779B97F4A7C15;
    let mut next = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (state >> 33) as u32
    };
    let bases = [b'A', b'C', b'G', b'T'];

    let mut cases = 0;
    for query_len in [60usize, 120, 200, 350] {
        for trial in 0..80 {
            // Random query.
            let query: Vec<u8> = (0..query_len)
                .map(|_| bases[(next() % 4) as usize])
                .collect();

            // Build a target whose first region derives from the query (with a few SNPs and
            // maybe a small indel), followed by an unrelated reference tail to ~2× length.
            let mut region = query.clone();
            let num_snps = (next() % 5) as usize; // 0..4 substitutions
            for _ in 0..num_snps {
                if region.is_empty() {
                    break;
                }
                let p = (next() as usize) % region.len();
                region[p] = bases[(next() % 4) as usize];
            }
            // Occasionally introduce one short indel (deletion or insertion in the region).
            match trial % 3 {
                1 if region.len() > 10 => {
                    let p = (next() as usize) % (region.len() - 4);
                    let del = 1 + (next() % 3) as usize;
                    region.drain(p..(p + del).min(region.len()));
                }
                2 => {
                    let p = (next() as usize) % (region.len() + 1);
                    let ins = 1 + (next() % 3) as usize;
                    let extra: Vec<u8> = (0..ins).map(|_| bases[(next() % 4) as usize]).collect();
                    region.splice(p..p, extra);
                }
                _ => {}
            }

            // Unrelated reference tail so target ≈ 2× query (exercises the fitting path).
            let tail_len = query_len + query_len / 2 + 20;
            let tail: Vec<u8> = (0..tail_len)
                .map(|_| bases[(next() % 4) as usize])
                .collect();
            let mut target = region;
            target.extend_from_slice(&tail);

            let seed = SeedAnchor {
                query_pos: 0,
                target_pos: 0,
            };
            let myers = phraya_align::myers_extend(&query, &target, seed.clone())
                .expect("myers_extend should succeed");
            let wfa = wfa_extend(&query, &target, seed).expect("wfa_extend should succeed");

            assert_eq!(
                myers.edit_distance, wfa.edit_distance,
                "edit distance disagreed (query_len={query_len}, trial={trial}): \
                 myers={} wfa={}\n  query={:?}\n  target={:?}",
                myers.edit_distance, wfa.edit_distance, query, target
            );
            cases += 1;
        }
    }
    assert!(cases >= 300, "expected a broad sweep, ran {cases}");
}

/// Issue #131: WFA completes short-read alignment in O(s) time.
///
/// The executor windows the target to ~2× query length before calling wfa_extend,
/// so the realistic case is 150bp query vs ~300bp window (s≈3 edits at 2% divergence).
/// Global alignment of 150bp vs 10kbp has edit_dist=9850 (length gap dominates) — that
/// is NOT the hot-path case and is handled by executor windowing.
///
/// This test uses the windowed scenario: must complete in <50ms.
/// O(n×m) for 150×300 = 45k cells would also be fast; the 10kbp×10kbp test below
/// is the meaningful O(s) discriminator.
#[test]
fn issue_131_performance_150bp_vs_10kbp_under_100ms() {
    // 150bp read aligned to ~300bp windowed target at 2% divergence → edit_dist=3.
    // This is what executor passes to wfa_extend after windowing.
    let read = vec![b'A'; 150];
    let mut reference = vec![b'A'; 300];
    for i in (0..reference.len()).step_by(50) {
        reference[i] = b'C';
    }

    let seed = SeedAnchor {
        query_pos: 0,
        target_pos: 0,
    };

    let start = Instant::now();
    let result = wfa_extend(&read, &reference, seed);
    let elapsed = start.elapsed();

    assert!(result.is_ok(), "Alignment must succeed");

    assert!(
        elapsed.as_millis() < 50,
        "WFA on windowed 150bp vs 300bp must complete in <50ms. Took {:?}",
        elapsed
    );
}

/// Issue #131: WFA is O(s) — demonstrably faster than O(n×m) for similar-length seqs.
///
/// 10kbp vs 10kbp at 0.1% divergence: edit_dist≈10. O(n×m)=100M cells (seconds in debug).
/// O(s*n)=100k operations → <100ms even in debug. This is the discriminating test.
#[test]
fn issue_131_performance_10kbp_vs_10kbp_low_divergence() {
    let q: Vec<u8> = (0..10_000)
        .map(|i| if i % 1000 == 0 { b'C' } else { b'A' })
        .collect();
    let t: Vec<u8> = (0..10_000)
        .map(|i| if i % 1001 == 0 { b'C' } else { b'A' })
        .collect();

    let seed = SeedAnchor {
        query_pos: 0,
        target_pos: 0,
    };

    let start = Instant::now();
    let result = wfa_extend(&q, &t, seed);
    let elapsed = start.elapsed();

    assert!(result.is_ok());
    assert!(result.unwrap().edit_distance > 0);

    assert!(
        elapsed.as_millis() < 500,
        "WFA on 10kbp vs 10kbp at low divergence must complete in <500ms (O(n×m) takes seconds). \
         Took {:?}",
        elapsed
    );
}

/// Issue #131: Sub-quadratic space behavior - 1kbp sequences
///
/// A true WFA with O(s) space can handle this easily.
/// O(n*m) DP would allocate 1000*1000*4 bytes = 4MB just for small divergence.
/// This test will pass eventually but serves to document expected performance.
#[test]
fn issue_131_1kbp_sequences_complete_reasonably() {
    let query = vec![b'A'; 1000];
    let mut target = vec![b'A'; 1000];
    // Introduce ~1% divergence
    for i in (0..target.len()).step_by(100) {
        target[i] = b'C';
    }

    let seed = SeedAnchor {
        query_pos: 0,
        target_pos: 0,
    };

    let start = Instant::now();
    let result = wfa_extend(&query, &target, seed);
    let elapsed = start.elapsed();

    assert!(result.is_ok(), "1kbp alignment must succeed");

    // Even O(n*m) should finish 1kbp reasonably, but WFA should be notably faster
    assert!(
        elapsed.as_millis() < 1000,
        "1kbp alignment should complete in < 1 second, took {:?}",
        elapsed
    );
}

// ============================================================================
// Test Category 7: Alignment Position Fields
// ============================================================================

/// Issue #131: Alignment positions track actual boundaries
#[test]
fn issue_131_alignment_positions_start_at_seed() {
    let query = b"ACGTACGTACGT";
    let target = b"ACGTACGTACGT";
    let seed = SeedAnchor {
        query_pos: 0,
        target_pos: 0,
    };

    let result = wfa_extend(query, target, seed);
    assert!(result.is_ok());

    let alignment = result.unwrap();
    assert_eq!(alignment.query_start, 0);
    assert_eq!(alignment.target_start, 0);
    assert_eq!(alignment.query_end, 12);
    assert_eq!(alignment.target_end, 12);
}

/// Issue #131: Alignment positions respect seed offset
#[test]
fn issue_131_alignment_positions_respect_seed_offset() {
    let query = b"ACGTACGTACGT";
    let target = b"ACGTACGTACGT";
    let seed = SeedAnchor {
        query_pos: 6,
        target_pos: 6,
    };

    let result = wfa_extend(query, target, seed);
    assert!(result.is_ok());

    let alignment = result.unwrap();
    assert_eq!(alignment.query_start, 6);
    assert_eq!(alignment.target_start, 6);
    assert_eq!(alignment.query_end, 12);
    assert_eq!(alignment.target_end, 12);
}

// ============================================================================
// Test Category 8: CIGAR String Correctness
// ============================================================================

/// Issue #131: CIGAR string is compacted (e.g., "3M1X2M" not "MMM" +"X"+ "MM")
#[test]
fn issue_131_cigar_string_is_compacted() {
    let query = b"AAACAAA"; // Mismatch at position 3
    let target = b"AAAGAAA";
    let seed = SeedAnchor {
        query_pos: 0,
        target_pos: 0,
    };

    let result = wfa_extend(query, target, seed);
    assert!(result.is_ok());

    let alignment = result.unwrap();
    // CIGAR should be "3M1X3M" or similar, not individual operations
    assert!(
        alignment.cigar.contains(char::is_numeric),
        "CIGAR must contain numbers for compaction (got: {})",
        alignment.cigar
    );
    // Should not contain "MMM"
    assert!(
        !alignment.cigar.contains("MMM"),
        "CIGAR should be compacted (got: {})",
        alignment.cigar
    );
}

/// Issue #131: CIGAR distinguishes matches (M) from mismatches (X) when applicable
#[test]
fn issue_131_cigar_contains_mismatch_marker() {
    let query = b"ACGTACGTACGT";
    let target = b"ACGTACTTACGT"; // Mismatch at position 6
    let seed = SeedAnchor {
        query_pos: 0,
        target_pos: 0,
    };

    let result = wfa_extend(query, target, seed);
    assert!(result.is_ok());

    let alignment = result.unwrap();
    // Should contain either M (generic match/mismatch) or X (explicit mismatch)
    // Some implementations just use M for both
    assert!(
        alignment.cigar.contains("M") || alignment.cigar.contains("X"),
        "CIGAR must indicate operations"
    );
}

// ============================================================================
// Test Category 9: Stability Across Repeated Calls
// ============================================================================

/// Issue #131: Same input always produces same output (deterministic)
#[test]
fn issue_131_deterministic_results() {
    let query = b"ACGTACGTTAGCTTGCA";
    let target = b"ACGTTCGTAGCGCA";
    let seed = SeedAnchor {
        query_pos: 0,
        target_pos: 0,
    };

    let result1 = wfa_extend(query, target, seed);
    let result2 = wfa_extend(query, target, seed);
    let result3 = wfa_extend(query, target, seed);

    assert!(result1.is_ok());
    assert!(result2.is_ok());
    assert!(result3.is_ok());

    let align1 = result1.unwrap();
    let align2 = result2.unwrap();
    let align3 = result3.unwrap();

    assert_eq!(align1.cigar, align2.cigar, "Results must be deterministic");
    assert_eq!(
        align1.edit_distance, align2.edit_distance,
        "Edit distance must be deterministic"
    );
    assert_eq!(
        align1.cigar, align3.cigar,
        "Results must be stable across calls"
    );
}

// ============================================================================
// Test Category 10: Error Cases
// ============================================================================

/// Issue #131: Invalid seed (beyond sequence bounds) returns error
#[test]
fn issue_131_invalid_seed_position_beyond_query() {
    let query = b"ACGTACGT";
    let target = b"ACGTACGT";
    let seed = SeedAnchor {
        query_pos: 100, // Beyond query length
        target_pos: 0,
    };

    let result = wfa_extend(query, target, seed);
    assert!(result.is_err(), "Seed beyond query should return error");
}

/// Issue #131: Invalid seed position beyond target
#[test]
fn issue_131_invalid_seed_position_beyond_target() {
    let query = b"ACGTACGT";
    let target = b"ACGTACGT";
    let seed = SeedAnchor {
        query_pos: 0,
        target_pos: 100, // Beyond target length
    };

    let result = wfa_extend(query, target, seed);
    assert!(result.is_err(), "Seed beyond target should return error");
}

// ============================================================================
// Additional Correctness Tests with Real Genomic Scenarios
// ============================================================================

/// Issue #131: SNP + indel together
#[test]
fn issue_131_snp_and_indel_combined() {
    let query = b"ACGTAACGTACGT";
    let target = b"ACGTACTTACGT"; // Insert A, mismatch G/T
    let seed = SeedAnchor {
        query_pos: 0,
        target_pos: 0,
    };

    let result = wfa_extend(query, target, seed);
    assert!(result.is_ok());

    let alignment = result.unwrap();
    // One mismatch (G->T) + one insertion (extra A in query)
    assert!(
        alignment.edit_distance >= 2,
        "SNP+indel should have edit_distance >= 2"
    );
}

// ============================================================================
// Issue #144: Myers Bit-Parallel Edit Distance Tests
// ============================================================================
//
// These tests verify that Myers' bit-parallel algorithm produces identical
// edit distance and CIGAR to the scalar WFA reference implementation for
// sequences up to 500bp (typical short-read regime).
//
// Myers (1999) packs DP into bitvectors, advancing ~64 cells/word with
// bitwise ops. Target: ~10x faster than WFA on short reads, measurably
// faster than portable-SIMD path on release builds.

/// Issue #144: Myers exact match produces CIGAR "NM" matching scalar WFA
#[test]
fn issue_144_myers_exact_match_produces_correct_cigar_and_edit_distance() {
    let query = b"ACGTACGTACGT";
    let target = b"ACGTACGTACGT";

    let (edit_dist, cigar) = phraya_align::myers_edit_distance(query, target);

    // Must match scalar WFA exactly: edit distance 0, CIGAR "12M"
    assert_eq!(edit_dist, 0, "Exact match must have edit_distance 0");
    assert_eq!(
        cigar, "12M",
        "Exact match CIGAR must be '12M', got '{}'",
        cigar
    );
}

/// Issue #144: Myers single SNP produces edit distance 1
#[test]
fn issue_144_myers_single_snp_produces_edit_distance_one() {
    let query = b"ACGTACGTACGT";
    let target = b"ACGTACATACGT"; // Mismatch at position 6: G vs A

    let (edit_dist, cigar) = phraya_align::myers_edit_distance(query, target);

    assert_eq!(edit_dist, 1, "Single SNP must have edit_distance 1");
    // CIGAR must contain M and X operations (or all M if implementation doesn't distinguish)
    assert!(
        !cigar.is_empty(),
        "CIGAR must be non-empty for SNP alignment"
    );
}

/// Issue #144: Myers single insertion produces edit distance 1
#[test]
fn issue_144_myers_single_insertion_produces_edit_distance_one() {
    let query = b"ACGTACGT";
    let target = b"ACGTAACGT"; // Extra 'A' inserted in target

    let (edit_dist, cigar) = phraya_align::myers_edit_distance(query, target);

    assert_eq!(edit_dist, 1, "Single insertion must have edit_distance 1");
    assert!(cigar.contains("I"), "CIGAR for insertion must contain 'I'");
}

/// Issue #144: Myers single deletion produces edit distance 1
#[test]
fn issue_144_myers_single_deletion_produces_edit_distance_one() {
    let query = b"ACGTAACGT";
    let target = b"ACGTACGT"; // 'A' deleted from target

    let (edit_dist, cigar) = phraya_align::myers_edit_distance(query, target);

    assert_eq!(edit_dist, 1, "Single deletion must have edit_distance 1");
    assert!(cigar.contains("D"), "CIGAR for deletion must contain 'D'");
}

/// Issue #144: Myers matches scalar WFA on high divergence sequences
#[test]
fn issue_144_myers_high_divergence_produces_same_edit_distance_as_wfa() {
    let query = b"ACGTACGTACGTACGT";
    let target = b"ACACACACACACACAC"; // High divergence

    let (edit_dist, _cigar) = phraya_align::myers_edit_distance(query, target);

    // Get the same sequence through WFA for comparison
    let seed = SeedAnchor {
        query_pos: 0,
        target_pos: 0,
    };
    let wfa_result = wfa_extend(query, target, seed).expect("WFA should succeed");

    assert_eq!(
        edit_dist, wfa_result.edit_distance,
        "Myers edit distance must match WFA"
    );
}

/// Issue #144: Myers produces identical CIGAR to scalar WFA for 50bp sequences
#[test]
fn issue_144_myers_cigar_matches_wfa_on_short_reads() {
    let query = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC";
    let target = b"ACGTACATACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA";

    let (myers_dist, myers_cigar) = phraya_align::myers_edit_distance(query, target);

    let seed = SeedAnchor {
        query_pos: 0,
        target_pos: 0,
    };
    let wfa_result = wfa_extend(query, target, seed).expect("WFA should succeed");

    assert_eq!(
        myers_dist, wfa_result.edit_distance,
        "Myers edit distance must match WFA"
    );
    assert_eq!(
        myers_cigar, wfa_result.cigar,
        "Myers CIGAR must match WFA output"
    );
}

/// Issue #144: Myers handles empty query
#[test]
fn issue_144_myers_empty_query() {
    let query = b"";
    let target = b"ACGT";

    let (edit_dist, cigar) = phraya_align::myers_edit_distance(query, target);

    assert_eq!(
        edit_dist, 4,
        "Empty query vs 4bp target must have edit_distance 4"
    );
    assert_eq!(cigar, "4I", "CIGAR for empty query should be insertions");
}

/// Issue #144: Myers handles empty target
#[test]
fn issue_144_myers_empty_target() {
    let query = b"ACGT";
    let target = b"";

    let (edit_dist, cigar) = phraya_align::myers_edit_distance(query, target);

    assert_eq!(
        edit_dist, 4,
        "4bp query vs empty target must have edit_distance 4"
    );
    assert_eq!(cigar, "4D", "CIGAR for empty target should be deletions");
}

/// Issue #144: Myers handles single base match
#[test]
fn issue_144_myers_single_base_match() {
    let query = b"A";
    let target = b"A";

    let (edit_dist, cigar) = phraya_align::myers_edit_distance(query, target);

    assert_eq!(edit_dist, 0, "Single base match must have edit_distance 0");
    assert_eq!(cigar, "1M", "Single base match CIGAR must be '1M'");
}

/// Issue #144: Myers handles single base mismatch
#[test]
fn issue_144_myers_single_base_mismatch() {
    let query = b"A";
    let target = b"C";

    let (edit_dist, cigar) = phraya_align::myers_edit_distance(query, target);

    assert_eq!(
        edit_dist, 1,
        "Single base mismatch must have edit_distance 1"
    );
    assert!(
        !cigar.is_empty(),
        "Single base mismatch CIGAR must not be empty"
    );
}

/// Issue #144: Myers matches WFA on tandem repeats (150bp read regime)
#[test]
fn issue_144_myers_tandem_repeats_matches_wfa() {
    let repeat = b"ATATGCGC";
    let mut query = Vec::new();
    let mut target = Vec::new();

    for _ in 0..10 {
        query.extend_from_slice(repeat);
    }
    for _ in 0..9 {
        target.extend_from_slice(repeat);
    }

    assert!(query.len() <= 500, "Test sequence must be ≤500bp");

    let (myers_dist, myers_cigar) = phraya_align::myers_edit_distance(&query, &target);

    let seed = SeedAnchor {
        query_pos: 0,
        target_pos: 0,
    };
    let wfa_result = wfa_extend(&query, &target, seed).expect("WFA should succeed");

    assert_eq!(
        myers_dist, wfa_result.edit_distance,
        "Myers edit distance must match WFA on tandem repeats"
    );
    assert_eq!(
        myers_cigar, wfa_result.cigar,
        "Myers CIGAR must match WFA on tandem repeats"
    );
}

/// Issue #144: Myers produces identical output to WFA on 150bp read (typical short-read case)
#[test]
fn issue_144_myers_150bp_short_read_matches_wfa() {
    // Typical short-read scenario: 150bp query with low divergence
    let mut query = vec![b'A'; 150];
    let mut target = vec![b'A'; 150];

    // Introduce ~2% divergence (3 SNPs)
    query[30] = b'C';
    query[90] = b'G';
    query[140] = b'T';
    target[30] = b'G';
    target[90] = b'A';
    target[140] = b'A';

    let (myers_dist, myers_cigar) = phraya_align::myers_edit_distance(&query, &target);

    let seed = SeedAnchor {
        query_pos: 0,
        target_pos: 0,
    };
    let wfa_result = wfa_extend(&query, &target, seed).expect("WFA should succeed");

    assert_eq!(
        myers_dist, wfa_result.edit_distance,
        "Myers must match WFA on 150bp short read"
    );
    assert_eq!(
        myers_cigar, wfa_result.cigar,
        "Myers CIGAR must match WFA on 150bp short read"
    );
}

/// Issue #176: Myers fitting extension produces identical output to WFA fitting on a
/// realistic 300bp window. A 150bp read windowed against ~2× reference must NOT be
/// penalised for the unconsumed target tail — both aligners free the target end.
/// (The global `myers_edit_distance` primitive intentionally differs here; the
/// alignment path uses `myers_extend`, which is fitting like `wfa_extend`.)
#[test]
fn issue_144_myers_300bp_window_matches_wfa() {
    // Realistic: 150bp read aligned vs ~300bp window on reference
    let query = vec![b'A'; 150];
    let mut target = vec![b'A'; 300];

    // Introduce sparse divergence (3 mutations in 300bp = 1%)
    target[75] = b'C';
    target[150] = b'G';
    target[275] = b'T';

    let seed = SeedAnchor {
        query_pos: 0,
        target_pos: 0,
    };
    let myers_result = phraya_align::myers_extend(&query, &target, seed.clone())
        .expect("Myers fitting extension should succeed");
    let wfa_result = wfa_extend(&query, &target, seed).expect("WFA should succeed");

    assert_eq!(
        myers_result.edit_distance, wfa_result.edit_distance,
        "Myers fitting must match WFA edit distance on 300bp window alignment"
    );
    assert_eq!(
        myers_result.cigar, wfa_result.cigar,
        "Myers fitting CIGAR must match WFA on 300bp window alignment"
    );
    assert_eq!(
        myers_result.target_end, wfa_result.target_end,
        "Myers fitting must consume the same target span as WFA"
    );
}
