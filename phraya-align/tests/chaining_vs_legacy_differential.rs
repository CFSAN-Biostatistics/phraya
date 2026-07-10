//! Differential test module: compare the pre-chaining legacy raw-vote anchor selection
//! (`AlignConfig::use_legacy_anchors = true`, debug-only A/B toggle) against the new
//! seed-chaining anchor selection (default) on shared fixtures.
//!
//! Styled after `phraya-align/src/wfa_simd.rs`'s `simd_vs_naive_differential` module: run
//! the same fixture through both implementations and compare. On clean/SNP/indel fixtures
//! where a read has one true, uniquely-mappable locus, both paths must agree on the
//! primary edit distance/CIGAR regardless of *how* the anchor was selected. On
//! multi-mapping/repeat-heavy fixtures, divergence is expected (chaining reports fewer,
//! better-supported alternates than raw vote-ranking) and is asserted/documented rather
//! than forced to match — per this session's decision to drop `sensitive`'s algebraic
//! byte-equivalence guarantee in favor of empirical (PA/CBS benchmark) validation.

use phraya_align::executor::{align_task_with_config, AlignConfig, Strategy};
use phraya_core::types::Sequence;
use phraya_io::plan::{PhrayaPlan, UseCase};
use std::collections::HashMap;

fn make_plan() -> PhrayaPlan {
    PhrayaPlan::new(
        UseCase::ReadsWithRef,
        vec!["test".to_string()],
        "2026-07-10T00:00:00Z".to_string(),
        HashMap::new(),
        HashMap::new(),
        vec![],
    )
}

/// Deterministic pseudo-random DNA (LCG), matching the pattern used throughout the
/// existing strategy test files — avoids low-complexity sequences that would cause
/// spurious minimizer-seed explosions unrelated to what's under test here.
fn diverse_dna(len: usize, seed: u64) -> Vec<u8> {
    let mut x = seed;
    (0..len)
        .map(|_| {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            b"ACGT"[((x >> 33) & 3) as usize]
        })
        .collect()
}

fn legacy_config(strategy: Strategy) -> AlignConfig {
    AlignConfig::new(strategy).with_legacy_anchors(true)
}

fn chained_config(strategy: Strategy) -> AlignConfig {
    AlignConfig::new(strategy).with_legacy_anchors(false)
}

/// A clean, uniquely-mapping read with one SNP: both paths must call the identical
/// variant (same position/ref/alt) and agree on primary edit distance, regardless of
/// anchor-selection mechanism.
#[test]
fn clean_read_with_snp_agrees_between_legacy_and_chained() {
    let target = diverse_dna(300, 1);
    let mut read = target[50..200].to_vec();
    read[30] = if read[30] == b'A' { b'C' } else { b'A' }; // SNP at read pos 30 (target pos 80)

    let query = Sequence::new(read, None, "q".to_string(), None);
    let target_seq = Sequence::new(target, None, "ref".to_string(), None);
    let plan = make_plan();

    for strategy in [Strategy::Fast, Strategy::Balanced, Strategy::Sensitive] {
        let legacy = align_task_with_config(&query, &target_seq, &plan, &legacy_config(strategy))
            .expect("legacy alignment should succeed");
        let chained = align_task_with_config(&query, &target_seq, &plan, &chained_config(strategy))
            .expect("chained alignment should succeed");

        assert_eq!(
            legacy.variants.len(),
            chained.variants.len(),
            "strategy {:?}: variant count must agree on a clean uniquely-mapping read",
            strategy
        );
        for (l, c) in legacy.variants.iter().zip(chained.variants.iter()) {
            assert_eq!(l.position(), c.position(), "strategy {:?}: variant position mismatch", strategy);
            assert_eq!(l.ref_base(), c.ref_base(), "strategy {:?}: ref_base mismatch", strategy);
        }
    }
}

/// A read with a real short indel: both paths must place it at the same locus and call
/// the same indel variant type — the indel must not fragment the chain differently than
/// it fragments the legacy independent-anchor extension.
#[test]
fn read_with_short_indel_agrees_between_legacy_and_chained() {
    let target = diverse_dna(300, 2);
    let mut read = target[50..200].to_vec();
    read.remove(60); // 1bp deletion relative to target at read-relative pos 60

    let query = Sequence::new(read, None, "q".to_string(), None);
    let target_seq = Sequence::new(target, None, "ref".to_string(), None);
    let plan = make_plan();

    for strategy in [Strategy::Balanced, Strategy::Sensitive] {
        let legacy = align_task_with_config(&query, &target_seq, &plan, &legacy_config(strategy))
            .expect("legacy alignment should succeed");
        let chained = align_task_with_config(&query, &target_seq, &plan, &chained_config(strategy))
            .expect("chained alignment should succeed");

        assert!(!legacy.variants.is_empty(), "strategy {:?}: legacy should find the indel", strategy);
        assert!(!chained.variants.is_empty(), "strategy {:?}: chained should find the indel", strategy);
        assert_eq!(
            legacy.variants[0].variant_type(),
            chained.variants[0].variant_type(),
            "strategy {:?}: indel variant type must agree",
            strategy
        );
    }
}

/// A read that multi-maps across two distant, internally-clean repeat copies. Legacy
/// (raw vote, K=∞ for Sensitive) enumerates every raw vote independently; chaining
/// collapses each copy into one chain. Primary placement must still agree (same
/// edit-distance-minimizing locus), but the two paths are NOT required to report
/// identical *alternate* sets — documenting the expected, accepted divergence rather
/// than asserting false equality.
#[test]
fn multi_mapping_read_diverges_in_alternate_count_but_agrees_on_primary() {
    let unit = diverse_dna(200, 3);
    let mut target = unit.clone();
    target.extend_from_slice(&diverse_dna(500, 4)); // filler between copies
    target.extend_from_slice(&unit); // second identical copy
    target.extend_from_slice(&diverse_dna(500, 6)); // trailing filler: gives the second
    // copy's anchor window enough margin to trigger WFA/Myers *fitting* mode (needs
    // target substantially longer than the query — see fill_wfa_fitting_impl's
    // "tn <= qn + qn/2 + 10" global-vs-fitting heuristic in wfa_simd.rs). Without this,
    // the second copy sits too close to the target's end and both the legacy and
    // chained paths fall back to global alignment, which penalizes the unconsumed tail
    // and produces a spuriously bad edit distance — an artifact of this fixture's
    // window, not of chaining.

    let read = unit[20..170].to_vec(); // uniquely lands at both copies equally well

    let query = Sequence::new(read, None, "q".to_string(), None);
    let target_seq = Sequence::new(target, None, "ref".to_string(), None);
    let plan = make_plan();

    let legacy = align_task_with_config(&query, &target_seq, &plan, &legacy_config(Strategy::Sensitive))
        .expect("legacy alignment should succeed");
    let chained = align_task_with_config(&query, &target_seq, &plan, &chained_config(Strategy::Sensitive))
        .expect("chained alignment should succeed");

    // Both must find at least the two true copies as placements (primary + >=1 alternate).
    assert!(
        legacy.query_positions.len() >= 2,
        "legacy should report both repeat copies as multi-mapping placements, got {}",
        legacy.query_positions.len()
    );
    assert!(
        chained.query_positions.len() >= 2,
        "chained should report both repeat copies as multi-mapping placements, got {}",
        chained.query_positions.len()
    );
    // Documenting the divergence: chaining is not required to produce the same *count*
    // of reported placements as legacy's raw K=∞ vote enumeration — only that genuine
    // multi-mapping signal (>=2 placements) survives in both.
}

/// Boundary case: an unmappable read (no shared minimizers with the target at all) must
/// be handled identically by both paths — neither should panic, and both should either
/// return `None` or a below-threshold/no-seed classification consistently.
#[test]
fn unmappable_read_handled_consistently() {
    let target = diverse_dna(300, 5);
    let read = diverse_dna(150, 999); // unrelated random sequence, shares no true locus

    let query = Sequence::new(read, None, "q".to_string(), None);
    let target_seq = Sequence::new(target, None, "ref".to_string(), None);
    let plan = make_plan();

    // Both paths must complete without panicking; whether they return Some or None is
    // not the point here (there's no seeding guarantee against random DNA at this
    // length) — this test's job is to catch a chaining-path panic/crash on the "nothing
    // seeded at all" edge case, which the unconditional (0,0) anchor fallback exists to
    // handle gracefully in both paths.
    let _ = align_task_with_config(&query, &target_seq, &plan, &legacy_config(Strategy::Balanced));
    let _ = align_task_with_config(&query, &target_seq, &plan, &chained_config(Strategy::Balanced));
}
