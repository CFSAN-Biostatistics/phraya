//! End-to-end correctness tests for the seed-chaining anchor selection (ADR-0012), at
//! the full `align_task_with_config` pipeline level — complementing `chaining.rs`'s unit
//! tests, which only exercise `chain_seeds` in isolation.
//!
//! These fixtures originated as a differential comparison against the pre-chaining
//! legacy raw-vote anchor selection during the chaining redesign's validation phase.
//! Now that the legacy path has been fully removed (validated on the full HPC benchmark
//! ladder — see `~/idiot-binfie/benchmarking-5.md`), they assert known-correct outcomes
//! directly instead of comparing two implementations.

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

/// A clean, uniquely-mapping read with one SNP: every strategy must call the correct
/// variant (same target position, same reference base) regardless of anchor cap K.
#[test]
fn clean_read_with_snp_calls_correct_variant_across_strategies() {
    let target = diverse_dna(300, 1);
    let mut read = target[50..200].to_vec();
    read[30] = if read[30] == b'A' { b'C' } else { b'A' }; // SNP at read pos 30 (target pos 80)
    let expected_ref_base = target[80];

    let query = Sequence::new(read, None, "q".to_string(), None);
    let target_seq = Sequence::new(target, None, "ref".to_string(), None);
    let plan = make_plan();

    for strategy in [Strategy::Fast, Strategy::Balanced, Strategy::Sensitive] {
        let result = align_task_with_config(&query, &target_seq, &plan, &AlignConfig::new(strategy))
            .expect("alignment should succeed");

        assert_eq!(
            result.variants.len(),
            1,
            "strategy {:?}: clean uniquely-mapping read must call exactly 1 variant",
            strategy
        );
        assert_eq!(
            result.variants[0].position(),
            80,
            "strategy {:?}: variant must be at target position 80",
            strategy
        );
        assert_eq!(
            result.variants[0].ref_base(),
            expected_ref_base,
            "strategy {:?}: ref_base must match the true target base",
            strategy
        );
    }
}

/// A read with a real short indel: chaining must place the seeds before and after the
/// gap into one chain (not fragment across two), so the indel is called correctly.
#[test]
fn read_with_short_indel_calls_correct_variant_type() {
    use phraya_core::types::VariantType;

    let target = diverse_dna(300, 2);
    let mut read = target[50..200].to_vec();
    read.remove(60); // 1bp deletion relative to target at read-relative pos 60

    let query = Sequence::new(read, None, "q".to_string(), None);
    let target_seq = Sequence::new(target, None, "ref".to_string(), None);
    let plan = make_plan();

    for strategy in [Strategy::Balanced, Strategy::Sensitive] {
        let result = align_task_with_config(&query, &target_seq, &plan, &AlignConfig::new(strategy))
            .expect("alignment should succeed");

        assert!(!result.variants.is_empty(), "strategy {:?}: should find the indel", strategy);
        assert_eq!(
            result.variants[0].variant_type(),
            VariantType::Deletion,
            "strategy {:?}: a 1bp deletion in the read must be called as Deletion",
            strategy
        );
    }
}

/// A read that multi-maps across two distant, internally-clean repeat copies: chaining
/// must collapse each copy into its own chain and report both as placements (genuine
/// multi-mapping signal preserved), not merge them into one or lose the second entirely.
#[test]
fn multi_mapping_read_reports_both_repeat_copies() {
    let unit = diverse_dna(200, 3);
    let mut target = unit.clone();
    target.extend_from_slice(&diverse_dna(500, 4)); // filler between copies
    target.extend_from_slice(&unit); // second identical copy
    target.extend_from_slice(&diverse_dna(500, 6)); // trailing filler: gives the second
    // copy's anchor window enough margin to trigger WFA/Myers *fitting* mode (needs
    // target substantially longer than the query — see fill_wfa_fitting_impl's
    // "tn <= qn + qn/2 + 10" global-vs-fitting heuristic in wfa_simd.rs). Without this,
    // the second copy sits too close to the target's end and extension falls back to
    // global alignment, penalizing the unconsumed tail with a spuriously bad edit
    // distance — an artifact of a too-short fixture, not of chaining.

    let read = unit[20..170].to_vec(); // uniquely lands at both copies equally well

    let query = Sequence::new(read, None, "q".to_string(), None);
    let target_seq = Sequence::new(target, None, "ref".to_string(), None);
    let plan = make_plan();

    let result = align_task_with_config(&query, &target_seq, &plan, &AlignConfig::new(Strategy::Sensitive))
        .expect("alignment should succeed");

    assert!(
        result.query_positions.len() >= 2,
        "sensitive must report both repeat copies as multi-mapping placements, got {}",
        result.query_positions.len()
    );
}

/// Boundary case: an unmappable read (no shared minimizers with the target at all) must
/// not panic — the unconditional (0,0) fallback anchor exists to handle this gracefully.
#[test]
fn unmappable_read_does_not_panic() {
    let target = diverse_dna(300, 5);
    let read = diverse_dna(150, 999); // unrelated random sequence, shares no true locus

    let query = Sequence::new(read, None, "q".to_string(), None);
    let target_seq = Sequence::new(target, None, "ref".to_string(), None);
    let plan = make_plan();

    // No seeding guarantee against random DNA at this length — the only assertion this
    // test makes is that alignment completes without panicking on the "nothing seeded
    // at all" edge case.
    let _ = align_task_with_config(&query, &target_seq, &plan, &AlignConfig::new(Strategy::Balanced));
}
