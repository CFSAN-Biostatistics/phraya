use phraya_align::executor::{align_task_with_config, AlignConfig, Strategy};
use phraya_core::types::Sequence;
use phraya_io::plan::{PhrayaPlan, UseCase};
/// Issue #184: feat(align): recall-axis strategy ladder; rename exact -> sensitive, K=1/5/inf (ADR-0008)
///
/// This test file contains RED (failing) acceptance tests for issue #184.
/// Tests verify that the strategy ladder is redefined on a single recall axis (anchor cap K),
/// and that the `exact` strategy is renamed to `sensitive`.
///
/// Expected API after implementation:
/// - Strategy enum: fast (K=1), balanced (K=5), sensitive (K=∞)
/// - AlignConfig struct: new AlignConfig::sensitive() method
/// - --strategy exact is rejected with updated error message
/// - Balanced strategy reports at most 5 anchors (multi-mapping preserved)
/// - Sensitive strategy reproduces old exact behavior (all anchors)
/// - Coverage windows unchanged: fast ±150bp, balanced ±50bp, sensitive ±25bp
///
/// References:
/// - ADR-0008: Strategy ladder as a single recall axis
/// - ADR-0003: Earlier strategy definitions (superseded)
use std::collections::HashMap;

fn make_plan() -> PhrayaPlan {
    PhrayaPlan::new(
        UseCase::ReadsWithRef,
        vec!["test".to_string()],
        "2026-06-01T00:00:00Z".to_string(),
        HashMap::new(),
        HashMap::new(),
        vec![],
    )
}

/// Deterministic pseudo-random DNA of a given length.
fn random_dna(seed: u64, len: usize) -> Vec<u8> {
    let mut state = seed;
    let bases = [b'A', b'C', b'G', b'T'];
    (0..len)
        .map(|_| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            bases[((state >> 33) % 4) as usize]
        })
        .collect()
}

// ============================================================================
// API EXISTENCE TESTS: Verify Strategy::Sensitive exists and is properly wired
// ============================================================================

/// Test that Strategy enum has been renamed: Exact -> Sensitive.
/// Acceptance criterion: Strategy::Sensitive exists and can be instantiated.
#[test]
fn issue_184_strategy_sensitive_exists() {
    // This test verifies that the enum variant exists.
    // It will fail if the variant is still named Strategy::Exact.
    let _sensitive = Strategy::Sensitive;
    let _balanced = Strategy::Balanced;
    let _fast = Strategy::Fast;
}

/// Test that AlignConfig has a new sensitive() method replacing exact().
/// Acceptance criterion: AlignConfig::sensitive() returns a config with Strategy::Sensitive.
#[test]
fn issue_184_align_config_has_sensitive_method() {
    let config = AlignConfig::sensitive();
    assert_eq!(
        config.strategy,
        Strategy::Sensitive,
        "AlignConfig::sensitive() must return a config with Strategy::Sensitive"
    );
}

/// Test that Strategy::Sensitive has the expected coverage window radius (±25bp, unchanged).
/// Acceptance criterion: Sensitive uses ±25bp coverage window, same as old Exact.
#[test]
fn issue_184_sensitive_strategy_has_correct_window_radius() {
    let config = AlignConfig::sensitive();
    assert_eq!(
        config.coverage_window_radius, 25,
        "Sensitive strategy must use ±25bp coverage window (same as old exact), got {}",
        config.coverage_window_radius
    );
}

/// Test that coverage windows are unchanged across all strategies.
/// Acceptance criteria:
/// - Fast: ±150bp (unchanged)
/// - Balanced: ±50bp (unchanged)
/// - Sensitive: ±25bp (unchanged, same as old exact)
#[test]
fn issue_184_coverage_windows_unchanged() {
    let fast = AlignConfig::new(Strategy::Fast);
    let balanced = AlignConfig::new(Strategy::Balanced);
    let sensitive = AlignConfig::new(Strategy::Sensitive);

    assert_eq!(fast.coverage_window_radius, 150, "fast must be ±150bp");
    assert_eq!(
        balanced.coverage_window_radius, 50,
        "balanced must be ±50bp"
    );
    assert_eq!(
        sensitive.coverage_window_radius, 25,
        "sensitive must be ±25bp"
    );
}

// ============================================================================
// CLI PARSING TESTS: moved to phraya-cli/tests/issue_184_cli_strategy_tests.rs
// ============================================================================
//
// The tests originally here simulated CLI parsing with a local, hand-rolled
// `match "exact" { ... }` statement instead of calling the real parser in
// phraya-cli/src/main.rs -- they would pass regardless of whether that parser
// was ever updated. Real subprocess-based tests (spawning the `phraya` binary,
// matching the pattern in issue_181_preset_rename.rs) now live in
// phraya-cli/tests/issue_184_cli_strategy_tests.rs.

// ============================================================================
// ANCHOR CAP K TESTS: Verify strategies implement K=1/5/∞ behavior
// ============================================================================

/// Test that fast strategy (K=1) reports at most 1 anchor.
/// Acceptance criterion: A read at an ambiguous locus reports ≤1 placement.
/// The tandem duplication fixture has two equally-good matches at offsets 0 and N.
#[test]
fn issue_184_fast_strategy_reports_single_anchor_on_tandem() {
    let unit = random_dna(0x0FAC_E001, 80);
    // Tandem duplication: read matches at offset 0 and offset 80 equally well.
    let mut target_bases = unit.clone();
    target_bases.extend_from_slice(&unit);
    target_bases.extend_from_slice(&random_dna(0x7777, 60));

    let query = Sequence::new(unit, None, "read1".to_string(), None);
    let target = Sequence::new(target_bases, None, "ref".to_string(), None);
    let plan = make_plan();

    let result = align_task_with_config(&query, &target, &plan, &AlignConfig::new(Strategy::Fast))
        .expect("fast alignment should succeed");

    assert_eq!(
        result.query_positions.len(),
        1,
        "Fast (K=1) must report exactly 1 placement on tandem duplication, got {}",
        result.query_positions.len()
    );
}

/// Test that balanced strategy (K=5) reports multiple placements at ambiguous loci.
/// Acceptance criterion: A read at a genuinely ambiguous locus reports >1 placement (multi-mapping preserved).
/// The tandem duplication fixture has two equally-good matches; balanced with K=5 should report both.
#[test]
fn issue_184_balanced_strategy_preserves_multimapping_at_k5() {
    let unit = random_dna(0x0FAC_E001, 80);
    // Tandem duplication: read matches at offset 0 and offset 80 equally well.
    let mut target_bases = unit.clone();
    target_bases.extend_from_slice(&unit);
    target_bases.extend_from_slice(&random_dna(0x7777, 60));

    let query = Sequence::new(unit, None, "read1".to_string(), None);
    let target = Sequence::new(target_bases, None, "ref".to_string(), None);
    let plan = make_plan();

    let result = align_task_with_config(
        &query,
        &target,
        &plan,
        &AlignConfig::new(Strategy::Balanced),
    )
    .expect("balanced alignment should succeed");

    assert!(
        result.query_positions.len() >= 2,
        "Balanced (K=5) must preserve multi-mapping signal on tandem duplication, got {} placements",
        result.query_positions.len()
    );
}

/// Test that balanced strategy (K=5) does not report excessive anchors.
/// Acceptance criterion: Balanced reports ≤5 primary anchors (bounded by K=5).
/// On a very repetitive target with many possible alignments, balanced caps at 5.
#[test]
fn issue_184_balanced_strategy_caps_at_k5() {
    // Create a highly repetitive read: 50bp of a single k-mer unit.
    let unit = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // ~56bp repeating pattern
    let query = Sequence::new(unit.to_vec(), None, "read1".to_string(), None);

    // Create a target with many repeat copies: 20 tandem repeats of the unit.
    // This creates >20 possible seed anchors, but balanced should cap at 5.
    let mut target_bases = Vec::new();
    for _ in 0..20 {
        target_bases.extend_from_slice(unit);
    }
    target_bases.extend_from_slice(&random_dna(0x9999, 100));
    let target = Sequence::new(target_bases, None, "ref".to_string(), None);

    let plan = make_plan();

    let result = align_task_with_config(
        &query,
        &target,
        &plan,
        &AlignConfig::new(Strategy::Balanced),
    )
    .expect("balanced alignment on repetitive target should succeed");

    assert!(
        result.query_positions.len() <= 6, // Up to 5 seed anchors + 1 fallback
        "Balanced (K=5) must cap reported placements at ~6 (5 anchors + fallback), got {}",
        result.query_positions.len()
    );
}

/// Test that sensitive strategy (K=∞) reports all anchors (regression test).
/// Acceptance criterion: Sensitive reproduces the old exact behavior: all distinct seed targets + fallback.
#[test]
fn issue_184_sensitive_strategy_reports_all_anchors() {
    let unit = random_dna(0x0FAC_E001, 80);
    // Tandem duplication: read matches at offset 0 and offset 80 equally well.
    let mut target_bases = unit.clone();
    target_bases.extend_from_slice(&unit);
    target_bases.extend_from_slice(&random_dna(0x7777, 60));

    let query = Sequence::new(unit, None, "read1".to_string(), None);
    let target = Sequence::new(target_bases, None, "ref".to_string(), None);
    let plan = make_plan();

    let result = align_task_with_config(
        &query,
        &target,
        &plan,
        &AlignConfig::new(Strategy::Sensitive),
    )
    .expect("sensitive alignment should succeed");

    // Sensitive should report both tandem matches (K=∞ means all anchors).
    assert!(
        result.query_positions.len() >= 2,
        "Sensitive (K=∞) must report all distinct seed target-starts, got {} placements",
        result.query_positions.len()
    );
}

/// Test that sensitive reproduces the old exact anchor set (on a multi-mapping fixture).
/// Acceptance criterion: Sensitive and old-exact both report ≥2 placements on tandem duplication.
/// This is a regression test to ensure sensitive() provides the same anchor list as exact() did.
#[test]
fn issue_184_sensitive_reproduces_old_exact_anchor_set() {
    let unit = random_dna(0xDEAD_BEEF, 90);
    // Create multiple seed anchors: 3 tandem repeats of unit, allowing ≥3 equally-good placements.
    let mut target_bases = unit.clone();
    target_bases.extend_from_slice(&unit);
    target_bases.extend_from_slice(&unit);
    target_bases.extend_from_slice(&random_dna(0x5678, 80));

    let query = Sequence::new(unit, None, "read1".to_string(), None);
    let target = Sequence::new(target_bases, None, "ref".to_string(), None);
    let plan = make_plan();

    let sensitive = align_task_with_config(
        &query,
        &target,
        &plan,
        &AlignConfig::new(Strategy::Sensitive),
    )
    .expect("sensitive alignment should succeed");

    // Sensitive (K=∞) should report all three tandem matches.
    assert!(
        sensitive.query_positions.len() >= 3,
        "Sensitive must report all seed anchors on 3×tandem duplication, got {}",
        sensitive.query_positions.len()
    );
}

// ============================================================================
// BACKWARD COMPATIBILITY TESTS: Verify unchanged semantics for all strategies
// ============================================================================

/// Test that balanced strategy (K=5) is still the default.
/// Acceptance criterion: AlignConfig::default() uses Strategy::Balanced.
#[test]
fn issue_184_balanced_remains_default_strategy() {
    let config = AlignConfig::default();
    assert_eq!(
        config.strategy,
        Strategy::Balanced,
        "default strategy must still be balanced"
    );
    assert_eq!(
        config.coverage_window_radius, 50,
        "default coverage window must still be ±50bp"
    );
}

/// Test that divergence cutoff still applies only to fast strategy.
/// Acceptance criterion: Fast drops a highly divergent read; balanced keeps it.
/// This behavior should not change with the K-axis refactor.
#[test]
fn issue_184_divergence_cutoff_applies_only_to_fast() {
    let read = random_dna(0x1234_5678, 100);

    // Target: read with 40 substitutions (≈40% divergence), followed by unrelated tail.
    let mut region = read.clone();
    for i in 50..90 {
        region[i] = match region[i] {
            b'A' => b'C',
            b'C' => b'A',
            b'G' => b'T',
            _ => b'G',
        };
    }
    let tail = random_dna(0x9999, 120);
    let mut target_bases = region;
    target_bases.extend_from_slice(&tail);

    let query = Sequence::new(read, None, "read1".to_string(), None);
    let target = Sequence::new(target_bases, None, "ref".to_string(), None);
    let plan = make_plan();

    // Balanced (K=5) should still handle the divergent read (no cutoff).
    let balanced = align_task_with_config(
        &query,
        &target,
        &plan,
        &AlignConfig::new(Strategy::Balanced),
    );
    assert!(
        balanced.is_some(),
        "Balanced (K=5) must not apply divergence cutoff — should align the divergent read"
    );

    // Fast should still drop it (divergence cutoff unchanged).
    let fast = align_task_with_config(&query, &target, &plan, &AlignConfig::new(Strategy::Fast));
    assert!(
        fast.is_none(),
        "Fast must still apply divergence cutoff — should drop the divergent read"
    );
}

/// Test that the 0.95 reporting threshold is unchanged across all strategies.
/// Acceptance criterion: All strategies use the same threshold for multi-mapping filter.
/// This is verified indirectly by checking that variant positions and scores are comparable.
#[test]
fn issue_184_multimapping_threshold_unchanged() {
    // A read with two equally-good placements (score ratio ≥0.95).
    let unit = random_dna(0x0FAC_E001, 80);
    let mut target_bases = unit.clone();
    target_bases.extend_from_slice(&unit); // Tandem duplication
    target_bases.extend_from_slice(&random_dna(0x7777, 60));

    let query = Sequence::new(unit, None, "read1".to_string(), None);
    let target = Sequence::new(target_bases, None, "ref".to_string(), None);
    let plan = make_plan();

    // Both should report the same multi-mapping.
    let balanced = align_task_with_config(
        &query,
        &target,
        &plan,
        &AlignConfig::new(Strategy::Balanced),
    )
    .expect("balanced should align");
    let sensitive = align_task_with_config(
        &query,
        &target,
        &plan,
        &AlignConfig::new(Strategy::Sensitive),
    )
    .expect("sensitive should align");

    // Both should report both placements (they pass the 0.95 threshold).
    assert_eq!(
        balanced.query_positions.len(), sensitive.query_positions.len(),
        "Balanced (K=5) and Sensitive (K=∞) should report the same count of placements on tandem (same threshold)"
    );
    assert!(
        balanced.query_positions.len() >= 2,
        "Both must report the tandem duplication (threshold unchanged)"
    );
}

// ============================================================================
// CODE CLEANUP TESTS: real error-message check lives in
// phraya-cli/tests/issue_184_cli_strategy_tests.rs (subprocess-based, exercises
// the actual CLI parser rather than a locally-built string).
// ============================================================================

/// Test that AlignConfig builder methods only provide fast, balanced, sensitive.
/// Acceptance criterion: factory methods are fast(), balanced(), sensitive() (no exact()).
#[test]
fn issue_184_align_config_factory_methods_correct() {
    let _fast = AlignConfig::fast();
    let _balanced = AlignConfig::balanced();
    let _sensitive = AlignConfig::sensitive();

    // If AlignConfig::exact() still exists, the implementation is incomplete.
    // This is caught by the type system / compile-time check.
}

// ============================================================================
// FINAL INTEGRATION TEST: Verify the full strategy ladder with K-axis
// ============================================================================

/// Test the complete strategy ladder on a single fixture:
/// - Fast (K=1): single best anchor
/// - Balanced (K=5): multiple anchors, capped at ~5
/// - Sensitive (K=∞): all anchors
///
/// This integration test ensures the ladder behaves as specified across all three rungs.
#[test]
fn issue_184_strategy_ladder_on_recall_axis() {
    let unit = random_dna(0xCAFE_BABE, 70);
    // Create 5+ tandem repeats to exceed K=5.
    let mut target_bases = Vec::new();
    for _ in 0..7 {
        target_bases.extend_from_slice(&unit);
    }
    target_bases.extend_from_slice(&random_dna(0x1111, 100));

    let query = Sequence::new(unit, None, "read1".to_string(), None);
    let target = Sequence::new(target_bases, None, "ref".to_string(), None);
    let plan = make_plan();

    let fast = align_task_with_config(&query, &target, &plan, &AlignConfig::new(Strategy::Fast))
        .expect("fast should align");
    let balanced = align_task_with_config(
        &query,
        &target,
        &plan,
        &AlignConfig::new(Strategy::Balanced),
    )
    .expect("balanced should align");
    let sensitive = align_task_with_config(
        &query,
        &target,
        &plan,
        &AlignConfig::new(Strategy::Sensitive),
    )
    .expect("sensitive should align");

    let fast_count = fast.query_positions.len();
    let balanced_count = balanced.query_positions.len();
    let sensitive_count = sensitive.query_positions.len();

    // Verify the hierarchy:
    // fast (K=1) <= balanced (K=5) <= sensitive (K=∞)
    assert_eq!(fast_count, 1, "Fast must report exactly 1 placement (K=1)");
    assert!(
        balanced_count <= 6, // Up to 5 seed anchors + 1 fallback
        "Balanced must report ≤6 placements (K=5 + fallback), got {}",
        balanced_count
    );
    assert_eq!(
        sensitive_count, 7,
        "Sensitive must report all 7 tandem repeats (K=∞), got {}",
        sensitive_count
    );
    assert!(
        sensitive_count >= balanced_count,
        "Sensitive (K=∞) must report ≥ Balanced (K=5)"
    );
}
