/// Issue #183: perf(align): score-bounded branch-and-bound alternate extension (ADR-0007)
///
/// RED acceptance tests against the real seams the implementation must fill in
/// (`phraya_align::score_bound_max_s`, `wfa_extend_capped`, `extend_alternates_bounded` —
/// all currently `todo!()` stubs in `phraya-align/src/executor.rs`). Every test here calls
/// production code directly; none re-derive the formula or behavior locally, so they fail
/// for the right reason (unimplemented!) rather than passing vacuously against a
/// self-referential check.

use phraya_align::{score_bound_max_s, wfa_extend_capped, extend_alternates_bounded};
use phraya_align::{wfa_extend_naive, SeedAnchor};

/// A 150bp perfect primary (d_best=0) must yield max_s = floor(0.05*150) = 7.
#[test]
fn issue_183_max_s_formula_150bp_perfect_primary() {
    assert_eq!(score_bound_max_s(150, 0), 7);
}

/// Safety invariant: max_s >= d_best for every query length / d_best pair (never prunes a
/// potential new primary). Checked against the real function, not a duplicate formula.
#[test]
fn issue_183_max_s_never_below_d_best() {
    for query_len in [50, 100, 150, 200, 500] {
        for d_best in [0, 1, 5, 10, 20, 50] {
            if d_best > query_len {
                continue;
            }
            let max_s = score_bound_max_s(query_len, d_best);
            assert!(
                max_s >= d_best,
                "max_s ({max_s}) must be >= d_best ({d_best}) at query_len={query_len}"
            );
        }
    }
}

/// Boundary: worst-case primary (d_best == query_len) must allow max_s == query_len (no
/// premature pruning when the primary itself is terrible).
#[test]
fn issue_183_max_s_worst_case_primary_equals_query_len() {
    assert_eq!(score_bound_max_s(100, 100), 100);
}

/// Boundary: max_s must scale with d_best per the 0.95 coefficient — a better primary
/// tightens the cap for alternates.
#[test]
fn issue_183_max_s_tightens_as_d_best_improves() {
    let query_len = 150;
    let max_s_at_d0 = score_bound_max_s(query_len, 0);
    let max_s_at_d10 = score_bound_max_s(query_len, 10);
    let max_s_at_d20 = score_bound_max_s(query_len, 20);
    assert!(max_s_at_d0 < max_s_at_d10, "cap must grow with d_best");
    assert!(max_s_at_d10 < max_s_at_d20, "cap must grow with d_best");
}

/// Build a query/target pair with a known, ground-truth edit distance (measured via the
/// existing, already-shipped `wfa_extend_naive` — independent of anything #183 adds), so
/// the capped-extension tests below have a real oracle rather than a hand-computed guess.
fn ground_truth_edit_distance(query: &[u8], target_window: &[u8]) -> usize {
    wfa_extend_naive(query, target_window, SeedAnchor { query_pos: 0, target_pos: 0 })
        .expect("naive WFA must succeed on a well-formed test fixture")
        .edit_distance
}

/// An alternate whose true edit distance is within the cap must be retained (Ok), with the
/// edit distance matching the independently-measured ground truth.
#[test]
fn issue_183_wfa_extend_capped_retains_alignment_within_cap() {
    let query = vec![b'A'; 150];
    let mut target = query.clone();
    // 3 substitutions -> ground-truth edit distance is small and known to be <= a generous cap.
    for i in [10, 50, 90] {
        target[i] = b'T';
    }
    let d_true = ground_truth_edit_distance(&query, &target);
    assert!(d_true <= 3, "fixture should have a small ground-truth edit distance, got {d_true}");

    let seed = SeedAnchor { query_pos: 0, target_pos: 0 };
    let result = wfa_extend_capped(&query, &target, seed, 10);
    let aln = result.expect("alignment within cap must be retained, not abandoned");
    assert_eq!(
        aln.edit_distance, d_true,
        "capped extension must report the same edit distance as the uncapped ground truth"
    );
}

/// An alternate whose true edit distance exceeds the cap must be abandoned (Err) rather
/// than silently truncated or reported with a wrong distance.
#[test]
fn issue_183_wfa_extend_capped_abandons_alignment_beyond_cap() {
    let query = vec![b'A'; 150];
    let mut target = query.clone();
    // 30 substitutions spread through the sequence -> ground-truth edit distance is large.
    for i in (0..150).step_by(5) {
        target[i] = b'T';
    }
    let d_true = ground_truth_edit_distance(&query, &target);
    assert!(d_true > 10, "fixture should have a large ground-truth edit distance, got {d_true}");

    let seed = SeedAnchor { query_pos: 0, target_pos: 0 };
    // Cap well below the true distance -> must abandon.
    let result = wfa_extend_capped(&query, &target, seed, 5);
    assert!(
        result.is_err(),
        "alignment whose true edit distance ({d_true}) exceeds the cap (5) must be abandoned"
    );
}

/// Core acceptance criterion: given a primary edit distance and a mix of alternates (some
/// within the score bound, some beyond it), `extend_alternates_bounded` must return exactly
/// the alternates within the bound, each with the correct (ground-truth-matching) edit
/// distance — the ones beyond the bound must be silently dropped, not errored out to the
/// caller or included with a wrong distance.
#[test]
fn issue_183_extend_alternates_bounded_drops_only_alternates_beyond_bound() {
    let query = vec![b'A'; 150];
    // primary is a perfect match (d_best = 0) -> max_s = score_bound_max_s(150, 0) = 7
    let primary_edit_distance = 0;
    let expected_max_s = score_bound_max_s(150, 0);
    assert_eq!(expected_max_s, 7);

    // Build one long target containing multiple windows (each 150bp) at distinct offsets,
    // one retained (needs 7 substitutions - editable at exactly the cap) and one dropped
    // (needs far more than 7).
    let window_len = 150;
    let mut target = Vec::new();

    // Window A (offset 0): retained -- 7 substitutions, right at the cap.
    let mut window_a = query.clone();
    for i in [5, 20, 40, 60, 80, 100, 130] {
        window_a[i] = b'T';
    }
    let d_a = ground_truth_edit_distance(&query, &window_a);
    assert!(d_a <= expected_max_s, "window A ground truth ({d_a}) must be within the cap ({expected_max_s})");
    target.extend_from_slice(&window_a);

    // Window B (offset window_len): dropped -- heavily diverged, far past the cap.
    let mut window_b = query.clone();
    for i in (0..150).step_by(3) {
        window_b[i] = b'T';
    }
    let d_b = ground_truth_edit_distance(&query, &window_b);
    assert!(d_b > expected_max_s, "window B ground truth ({d_b}) must exceed the cap ({expected_max_s})");
    target.extend_from_slice(&window_b);

    let alternates = vec![
        SeedAnchor { query_pos: 0, target_pos: 0 },
        SeedAnchor { query_pos: 0, target_pos: window_len },
    ];

    let retained = extend_alternates_bounded(&query, &target, primary_edit_distance, &alternates);

    assert_eq!(
        retained.len(),
        1,
        "exactly one of the two alternates should survive the score bound"
    );
    assert_eq!(
        retained[0].target_start, 0,
        "the retained alternate must be window A (offset 0), not the abandoned window B"
    );
    assert_eq!(
        retained[0].edit_distance, d_a,
        "the retained alternate's edit distance must match the ground truth"
    );
}

/// Bound-tightening acceptance criterion: when an early alternate turns out to beat the
/// primary, the bound used for *later* alternates must tighten accordingly — an alternate
/// that would have been retained under the original (looser) bound, but not under the
/// tightened one, must be dropped.
#[test]
fn issue_183_extend_alternates_bounded_tightens_after_finding_a_better_alternate() {
    let query = vec![b'A'; 150];
    let primary_edit_distance = 20; // deliberately poor primary
    let loose_max_s = score_bound_max_s(150, 20);
    assert!(loose_max_s > 20, "sanity: a poor primary should allow a generous initial cap");

    let window_len = 150;
    let mut target = Vec::new();

    // Window A (processed first): a much better alternate (d ~ 0), which should become the
    // new incumbent and tighten the bound for subsequent alternates.
    let window_a = query.clone(); // 0 substitutions -> d_a == 0
    let d_a = ground_truth_edit_distance(&query, &window_a);
    assert_eq!(d_a, 0);
    target.extend_from_slice(&window_a);

    // Window B (processed second): needs more edits than score_bound_max_s(150, d_a) allows,
    // but fewer than the original loose_max_s -- only abandoned if the bound actually
    // tightened after window A was processed.
    let tightened_max_s = score_bound_max_s(150, 0);
    let mut window_b = query.clone();
    // Substitution count strictly between tightened_max_s and loose_max_s.
    let mid = (tightened_max_s + loose_max_s) / 2 + 1;
    assert!(mid > tightened_max_s && mid <= loose_max_s, "test fixture setup invariant");
    for i in 0..mid.min(window_len) {
        window_b[i] = b'T';
    }
    let d_b = ground_truth_edit_distance(&query, &window_b);
    assert!(
        d_b > tightened_max_s,
        "window B ground truth ({d_b}) must exceed the tightened cap ({tightened_max_s}) for this test to be meaningful"
    );
    target.extend_from_slice(&window_b);

    let alternates = vec![
        SeedAnchor { query_pos: 0, target_pos: 0 },
        SeedAnchor { query_pos: 0, target_pos: window_len },
    ];

    let retained = extend_alternates_bounded(&query, &target, primary_edit_distance, &alternates);

    // Only window A (the new incumbent) should survive; window B must be dropped because
    // the bound tightened to reflect window A's distance before window B was processed.
    assert_eq!(
        retained.len(),
        1,
        "window B must be dropped once the bound tightens after window A is found"
    );
    assert_eq!(retained[0].target_start, 0, "surviving alternate must be window A");
}

/// End-to-end guardrail (not a substitute for the unit tests above): the *existing*
/// alignment path must still detect a real SNP once the branch-and-bound wiring lands in
/// `align_read` -- output must be unchanged by a pure performance optimization. This test
/// passes today (before #183) and must continue to pass after; its purpose is to be part of
/// the regression net during implementation, not to prove the new feature exists.
#[test]
fn issue_183_end_to_end_variant_output_unaffected_by_bounding() {
    use phraya_align::executor::{align_task_with_config, AlignConfig};
    use phraya_core::types::Sequence;
    use phraya_io::plan::{PhrayaPlan, UseCase};
    use std::collections::HashMap;

    let mut query_bases = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let mut target_bases = query_bases.clone();
    query_bases[25] = b'T';
    target_bases[25] = b'C';
    query_bases.resize(150, b'A');
    target_bases.resize(150, b'A');

    let query = Sequence::new(query_bases, None, "query".to_string(), None);
    let target = Sequence::new(target_bases, None, "target".to_string(), None);
    let plan = PhrayaPlan::new(
        UseCase::ReadsWithRef,
        vec!["test".to_string()],
        "2026-06-01T00:00:00Z".to_string(),
        HashMap::new(),
        HashMap::new(),
        vec![],
    );

    let result = align_task_with_config(&query, &target, &plan, &AlignConfig::balanced());
    let aln = result.expect("alignment should succeed");
    assert!(
        aln.variants.iter().any(|v| v.position() == 25),
        "SNP at position 25 must still be reported"
    );
}
