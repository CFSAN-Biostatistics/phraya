//! The `(0,0)` fallback anchor is scoped to comparable-length query/target pairs.
//!
//! `anchors_from_chains` falls back to a `(0,0)` anchor when no chain survives, so that two
//! comparable-length sequences that share no minimizer can still align (issue #146, and the
//! Case 4 contigs-only comparisons that depend on it). Applied unconditionally, that same
//! fallback force-aligns every read against every reference it does *not* come from: the
//! alignment can only cover `target[0..2*query_len)`, scores ~50-60% identity, and — because
//! only the cross-space sidecar filters on identity, not the `.phraya` — deposits
//! ~`query_len/2` fabricated `VariantObservation`s in a reference the read never came from.
//!
//! These tests pin both halves of the resulting scope rule: no seedless alignment when the
//! target dwarfs the query, but the fallback preserved when the two are comparable.
//!
//! `Alphabet::Protein` is scoped further still: the fallback never applies at all, comparable
//! length or not. A proteome's candidate proteins are *naturally* comparable in length to a
//! query (that's the normal shape of a database, not a homology signal the way it is for a
//! handful of DNA contigs), so the DNA ratio rule would force a real extension attempt
//! against every candidate with no seed — the combinatorial cost issue #146's rule was never
//! meant to pay in a reference-palette protein search.

use phraya_align::executor::{
    align_read, fallback_anchor_applies, AlignConfig, Strategy, TargetContext,
};
use phraya_core::types::{
    sketch_sequence_alphabet, sketch_sequence_default, Alphabet, Sequence, DEFAULT_K_PROTEIN,
    DEFAULT_W_PROTEIN,
};
use phraya_io::plan::{PhrayaPlan, UseCase};
use std::collections::HashSet;

fn make_plan() -> PhrayaPlan {
    PhrayaPlan::new(
        UseCase::ReadsWithRef,
        vec!["test".to_string()],
        "2026-01-01T00:00:00Z".to_string(),
        std::collections::HashMap::new(),
        std::collections::HashMap::new(),
        vec![],
    )
}

fn diverse_dna(len: usize, seed: u64) -> Vec<u8> {
    let mut x = seed;
    (0..len)
        .map(|_| {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"ACGT"[((x >> 33) & 3) as usize]
        })
        .collect()
}

/// Assert the fixture really is seedless, so a test that expects "no placement" is not
/// passing for the wrong reason.
fn assert_no_shared_minimizer(read: &Sequence, target: &Sequence) {
    let read_vals: HashSet<u64> = sketch_sequence_default(read)
        .minimizers
        .iter()
        .map(|&(v, _)| v)
        .collect();
    let shared = sketch_sequence_default(target)
        .minimizers
        .iter()
        .filter(|(v, _)| read_vals.contains(v))
        .count();
    assert_eq!(
        shared, 0,
        "fixture precondition: read and target must share no minimizer"
    );
}

fn diverse_protein(len: usize, seed: u64) -> Vec<u8> {
    const AA: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
    let mut x = seed;
    (0..len)
        .map(|_| {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            AA[((x >> 33) as usize) % AA.len()]
        })
        .collect()
}

/// Protein counterpart of [`assert_no_shared_minimizer`], sketched at protein k/w.
fn assert_no_shared_minimizer_protein(query: &Sequence, target: &Sequence) {
    let query_vals: HashSet<u64> =
        sketch_sequence_alphabet(query, DEFAULT_K_PROTEIN, DEFAULT_W_PROTEIN, Alphabet::Protein)
            .minimizers
            .iter()
            .map(|&(v, _)| v)
            .collect();
    let shared =
        sketch_sequence_alphabet(target, DEFAULT_K_PROTEIN, DEFAULT_W_PROTEIN, Alphabet::Protein)
            .minimizers
            .iter()
            .filter(|(v, _)| query_vals.contains(v))
            .count();
    assert_eq!(
        shared, 0,
        "fixture precondition: protein query and target must share no minimizer"
    );
}

/// A read that shares no seed with a target far longer than itself must not place at all.
///
/// Before the scope rule this returned `Some`, carrying a ~59%-identity alignment pinned at
/// target offset 0 and ~62 fabricated variants.
#[test]
fn seedless_read_does_not_place_against_a_much_longer_target() {
    let target = Sequence::new(diverse_dna(20_000, 7), None, "ref".to_string(), None);
    let read = Sequence::new(diverse_dna(150, 99_991), None, "read".to_string(), None);
    assert_no_shared_minimizer(&read, &target);

    let plan = make_plan();
    let config = AlignConfig::new(Strategy::Balanced);
    let ctx = TargetContext::build(&target, &plan, config.strategy);

    assert!(
        align_read(&ctx, &read, &plan, &config, None).is_none(),
        "a read sharing no seed with a 20kb target must not be force-aligned at offset 0"
    );
}

/// The contig-vs-contig fallback (issue #146 / Case 4) survives: comparable lengths still
/// attempt a seedless alignment, because there a real homolog genuinely can share no
/// minimizer.
#[test]
fn seedless_query_still_attempts_alignment_against_comparable_length_target() {
    let target = Sequence::new(diverse_dna(1_200, 13), None, "contig_a".to_string(), None);
    let query = Sequence::new(diverse_dna(150, 99_991), None, "contig_b".to_string(), None);
    assert_no_shared_minimizer(&query, &target);

    let plan = make_plan();
    let config = AlignConfig::new(Strategy::Balanced);
    let ctx = TargetContext::build(&target, &plan, config.strategy);

    assert!(
        align_read(&ctx, &query, &plan, &config, None).is_some(),
        "comparable-length sequences must still reach the (0,0) fallback anchor"
    );
}

/// The predicate reference-palette alignment shares with `extend_chains` to decide a space
/// is unreachable. Pinned directly because the two callers must agree: if the palette
/// prefilter skipped a space the aligner would still have placed reads in, output would
/// silently go missing.
#[test]
fn fallback_scope_boundary_is_ten_times_query_length() {
    assert!(fallback_anchor_applies(150, 1_500, Alphabet::Dna), "exactly 10x still applies");
    assert!(!fallback_anchor_applies(150, 1_501, Alphabet::Dna), "past 10x does not apply");
    assert!(fallback_anchor_applies(150, 150, Alphabet::Dna), "equal lengths apply");
    assert!(
        fallback_anchor_applies(0, 0, Alphabet::Dna),
        "degenerate lengths must not panic or divide"
    );
}

/// Protein disables the fallback entirely, at every ratio the DNA test above applies it at —
/// including equal lengths, the case a proteome's candidates hit by construction.
#[test]
fn fallback_scope_is_disabled_entirely_for_protein() {
    assert!(!fallback_anchor_applies(150, 150, Alphabet::Protein));
    assert!(!fallback_anchor_applies(150, 1_500, Alphabet::Protein));
    assert!(!fallback_anchor_applies(0, 0, Alphabet::Protein));
}

/// The `align_read`-level counterpart of the predicate test above: a protein query and a
/// same-length candidate protein sharing zero minimizers must not force-align, unlike the
/// DNA `seedless_query_still_attempts_alignment_against_comparable_length_target` case.
/// Without this, reference-palette protein search would attempt a real Myers/WFA extension
/// against every candidate in a proteome regardless of seed evidence.
#[test]
fn seedless_protein_query_does_not_force_align_even_at_comparable_length() {
    let target = Sequence::new(diverse_protein(300, 13), None, "protein_a".to_string(), None);
    let query = Sequence::new(diverse_protein(300, 99_991), None, "protein_b".to_string(), None);
    assert_no_shared_minimizer_protein(&query, &target);

    let mut plan = make_plan();
    plan.alphabet = Alphabet::Protein;
    let config = AlignConfig::new(Strategy::Balanced);
    let ctx = TargetContext::build(&target, &plan, config.strategy);

    assert!(
        align_read(&ctx, &query, &plan, &config, None).is_none(),
        "comparable-length protein sequences with zero shared minimizers must not force-align"
    );
}
