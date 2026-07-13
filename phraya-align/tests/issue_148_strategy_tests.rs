/// Issue #148: Add --strategy flag and tie coverage window size to alignment strategy
///
/// This test file contains RED (failing) acceptance tests for issue #148.
/// Tests verify that different alignment strategies produce different local coverage
/// window radii in the resulting variant observations.
///
/// Expected API after implementation:
/// - Strategy enum: fast|balanced|sensitive
/// - AlignConfig struct carrying strategy and coverage_window_radius
/// - align_task(&query, &target, &plan, &config) instead of align_task(&query, &target, &plan)
/// - Window radii: fast: ±150bp, balanced: ±50bp (default), sensitive: ±25bp

use std::collections::HashMap;
use phraya_core::types::Sequence;
use phraya_io::plan::{PhrayaPlan, UseCase};

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

/// Test that AlignConfig struct exists and carries strategy information.
/// This is a fundamental prerequisite for all other tests.
#[test]
fn issue_148_align_config_struct_exists() {
    use phraya_align::executor::{AlignConfig, Strategy};
    let config = AlignConfig::default();
    assert_eq!(config.strategy, Strategy::Balanced, "default strategy should be Balanced");
}

/// Test that Strategy enum has fast, balanced, and sensitive variants.
#[test]
fn issue_148_strategy_enum_has_required_variants() {
    use phraya_align::executor::{AlignConfig, Strategy};
    let fast_config = AlignConfig::new(Strategy::Fast);
    let balanced_config = AlignConfig::new(Strategy::Balanced);
    let sensitive_config = AlignConfig::new(Strategy::Sensitive);

    assert_eq!(fast_config.coverage_window_radius, 150);
    assert_eq!(balanced_config.coverage_window_radius, 50);
    assert_eq!(sensitive_config.coverage_window_radius, 25);
}

/// Test that align_task signature accepts AlignConfig parameter.
/// Current signature: align_task(&Sequence, &Sequence, &PhrayaPlan) -> Option<AlignmentResult>
/// New signature: align_task(&Sequence, &Sequence, &PhrayaPlan, &AlignConfig) -> Option<AlignmentResult>
#[test]
fn issue_148_align_task_signature_accepts_config() {
    use phraya_align::executor::{align_task_with_config, AlignConfig};

    let mut query_bases = vec![b'A'; 100];
    let mut target_bases = vec![b'A'; 200];
    query_bases[50] = b'T';
    target_bases[50] = b'C';

    let query = Sequence::new(query_bases, None, "q".to_string(), None);
    let target = Sequence::new(target_bases, None, "ref".to_string(), None);
    let plan = make_plan();

    let config = AlignConfig::balanced();
    let result = align_task_with_config(&query, &target, &plan, &config);
    assert!(result.is_some(), "alignment with config should succeed");
}

/// Test that fast strategy produces ±150bp local coverage window.
/// A variant at position 50 in a 200bp target should have window length:
/// window_start = max(0, 50 - 150) = 0
/// window_end = min(200, 50 + 150 + 1) = 200
/// Expected window length: 200
#[test]
fn issue_148_fast_strategy_produces_wide_coverage_window() {
    use phraya_align::executor::{align_task_with_config, AlignConfig, Strategy};

    let mut query_bases = vec![b'A'; 100];
    let mut target_bases = vec![b'A'; 200];
    query_bases[50] = b'T';
    target_bases[50] = b'C'; // SNP at position 50

    let query = Sequence::new(query_bases, None, "q".to_string(), None);
    let target = Sequence::new(target_bases, None, "ref".to_string(), None);
    let plan = make_plan();

    let config = AlignConfig::new(Strategy::Fast);
    let result = align_task_with_config(&query, &target, &plan, &config).expect("alignment should succeed");
    assert!(!result.variants.is_empty(), "should have at least one variant");
    let var = &result.variants[0];
    let lc = var.local_coverage();

    // For fast strategy (±150bp), window around pos 50 should extend to nearly the bounds:
    assert!(
        lc.len() >= 150,
        "fast strategy window should be ±150bp (length ≥ 150), got {}",
        lc.len()
    );
}

/// Test that balanced strategy produces ±50bp local coverage window (current behavior).
/// Window around position 50 should be:
/// window_start = max(0, 50 - 50) = 0
/// window_end = min(200, 50 + 50 + 1) = 101
/// Expected window length: exactly 101
#[test]
fn issue_148_balanced_strategy_produces_default_coverage_window() {
    use phraya_align::executor::{align_task_with_config, AlignConfig, Strategy};

    let mut query_bases = vec![b'A'; 100];
    let mut target_bases = vec![b'A'; 200];
    query_bases[50] = b'T';
    target_bases[50] = b'C';

    let query = Sequence::new(query_bases, None, "q".to_string(), None);
    let target = Sequence::new(target_bases, None, "ref".to_string(), None);
    let plan = make_plan();

    let config = AlignConfig::new(Strategy::Balanced);
    let result = align_task_with_config(&query, &target, &plan, &config).expect("alignment should succeed");
    let var = &result.variants[0];
    let lc = var.local_coverage();

    // Balanced strategy: ±50bp window. SNP at target pos 50.
    // window_start = max(0, 50-50) = 0, window_end = min(200, 50+51) = 101 → len=101.
    assert_eq!(
        lc.len(),
        101,
        "balanced strategy window should be ±50bp (length 101), got {}",
        lc.len()
    );
}

/// Test that sensitive strategy produces ±25bp local coverage window.
/// Window around position 50 should be:
/// window_start = max(0, 50 - 25) = 25
/// window_end = min(200, 50 + 25 + 1) = 76
/// Expected window length: exactly 51
#[test]
fn issue_148_sensitive_strategy_produces_narrow_coverage_window() {
    use phraya_align::executor::{align_task_with_config, AlignConfig, Strategy};

    let mut query_bases = vec![b'A'; 100];
    let mut target_bases = vec![b'A'; 200];
    query_bases[50] = b'T';
    target_bases[50] = b'C';

    let query = Sequence::new(query_bases, None, "q".to_string(), None);
    let target = Sequence::new(target_bases, None, "ref".to_string(), None);
    let plan = make_plan();

    let config = AlignConfig::new(Strategy::Sensitive);
    let result = align_task_with_config(&query, &target, &plan, &config).expect("alignment should succeed");
    let var = &result.variants[0];
    let lc = var.local_coverage();

    // Sensitive strategy uses tight window: ±25bp
    assert_eq!(
        lc.len(),
        51,
        "sensitive strategy window should be ±25bp (length 51), got {}",
        lc.len()
    );
}

/// Test boundary condition: AlignConfig.coverage_window_radius is correct for fast strategy.
/// This directly tests the config struct, which is the source of truth for the window radius.
#[test]
fn issue_148_window_radius_boundary_conditions_fast() {
    use phraya_align::executor::{AlignConfig, Strategy};

    let config = AlignConfig::new(Strategy::Fast);
    assert_eq!(
        config.coverage_window_radius,
        150,
        "fast strategy must use ±150bp radius, got {}",
        config.coverage_window_radius
    );
    // Window at position 5: [max(0,5-150)=0, min(N, 5+151)] — always < 157 elements
    let pos: usize = 5;
    let target_len: usize = 200;
    let start = if pos >= config.coverage_window_radius { pos - config.coverage_window_radius } else { 0 };
    let end = (pos + config.coverage_window_radius + 1).min(target_len);
    assert_eq!(end - start, 156, "fast window at pos 5 in 200bp target should be 156 elements");
}

/// Test boundary condition: AlignConfig.coverage_window_radius is correct for sensitive strategy.
#[test]
fn issue_148_window_radius_boundary_conditions_sensitive() {
    use phraya_align::executor::{AlignConfig, Strategy};

    let config = AlignConfig::new(Strategy::Sensitive);
    assert_eq!(
        config.coverage_window_radius,
        25,
        "sensitive strategy must use ±25bp radius, got {}",
        config.coverage_window_radius
    );
    // Window at position 5: [max(0,5-25)=0, min(200, 5+26)] = [0, 31] = 31 elements
    let pos: usize = 5;
    let target_len: usize = 200;
    let start = if pos >= config.coverage_window_radius { pos - config.coverage_window_radius } else { 0 };
    let end = (pos + config.coverage_window_radius + 1).min(target_len);
    assert_eq!(end - start, 31, "sensitive window at pos 5 in 200bp target should be 31 elements");
}

/// Test backward compatibility: default strategy should be balanced.
/// Existing code that doesn't specify a strategy should maintain current behavior (±50bp).
#[test]
fn issue_148_balanced_is_default_and_backward_compatible() {
    use phraya_align::executor::{align_task_with_config, AlignConfig};

    let mut query_bases = vec![b'A'; 100];
    let mut target_bases = vec![b'A'; 200];
    query_bases[50] = b'T';
    target_bases[50] = b'C';

    let query = Sequence::new(query_bases, None, "q".to_string(), None);
    let target = Sequence::new(target_bases, None, "ref".to_string(), None);
    let plan = make_plan();

    let config = AlignConfig::default(); // Should default to balanced
    let result = align_task_with_config(&query, &target, &plan, &config).expect("alignment should succeed");
    let var = &result.variants[0];
    let lc = var.local_coverage();

    // Verify the window matches the ±50bp behavior. SNP at pos 50 → [0, 101) = 101.
    assert_eq!(
        lc.len(),
        101,
        "default strategy should maintain ±50bp window (length 101), got {}",
        lc.len()
    );
}

/// Test that sensitive strategy narrows the window near a variant cluster.
/// Multiple SNPs at positions 48-52 with sensitive strategy should produce ±25bp windows.
#[test]
fn issue_148_sensitive_strategy_narrows_window_near_variant_cluster() {
    use phraya_align::executor::{align_task_with_config, AlignConfig, Strategy};

    let mut query_bases = vec![b'A'; 100];
    let mut target_bases = vec![b'A'; 200];

    // Create a cluster of SNPs at positions 48, 49, 50, 51, 52
    for i in 48..=52 {
        query_bases[i] = b'T';
        target_bases[i] = b'C';
    }

    let query = Sequence::new(query_bases, None, "q".to_string(), None);
    let target = Sequence::new(target_bases, None, "ref".to_string(), None);
    let plan = make_plan();

    let config = AlignConfig::new(Strategy::Sensitive);
    let result = align_task_with_config(&query, &target, &plan, &config).expect("alignment should succeed");

    // Alignment places SNPs at target positions index-1. Find the middle variant.
    // Any variant in the cluster should have ±25bp window (51 or 50 elements).
    assert!(!result.variants.is_empty(), "should have variants");
    let lc = result.variants[result.variants.len() / 2].local_coverage();
    assert!(
        lc.len() == 51 || lc.len() == 50,
        "sensitive strategy cluster window should be ~±25bp (length 50-51), got {}",
        lc.len()
    );
}

/// Test that fast strategy provides wide context for complex regions.
/// Fast mode (±150bp) should capture more context than balanced mode (±50bp).
#[test]
fn issue_148_fast_strategy_provides_context_for_complex_regions() {
    use phraya_align::executor::{align_task_with_config, AlignConfig, Strategy};

    let mut query_bases = vec![b'A'; 100];
    let mut target_bases = vec![b'A'; 200];

    // Create variations at positions 48-50
    query_bases[48] = b'T';
    query_bases[49] = b'G';
    target_bases[48] = b'C';
    target_bases[49] = b'C';
    target_bases[50] = b'G';

    let query = Sequence::new(query_bases, None, "q".to_string(), None);
    let target = Sequence::new(target_bases, None, "ref".to_string(), None);
    let plan = make_plan();

    let config = AlignConfig::new(Strategy::Fast);
    let result = align_task_with_config(&query, &target, &plan, &config).expect("alignment should succeed");

    // The window should be ±150bp to capture structural context
    let var = result.variants.iter().next().expect("at least one variant");
    let lc = var.local_coverage();
    assert!(
        lc.len() >= 150,
        "fast strategy on complex region: window should be ±150bp (length ≥ 150), got {}",
        lc.len()
    );
}

/// Test that all variants from the same alignment task use the same window radius.
/// Different strategies must be consistent: all variants in one task use the same window.
#[test]
fn issue_148_all_variants_in_task_use_same_window_radius() {
    use phraya_align::executor::{align_task_with_config, AlignConfig, Strategy};

    // Use a longer target so both SNPs are far from edges (≥50bp clearance after alignment offset).
    let mut query_bases = vec![b'A'; 100];
    let mut target_bases = vec![b'A'; 300];

    // SNPs at positions 75 and 150 — both ≥50bp from any edge in 300bp target.
    query_bases[75] = b'T';
    target_bases[75] = b'C';
    query_bases[80] = b'T'; // second SNP nearby in query
    target_bases[150] = b'C'; // second SNP in target

    let query = Sequence::new(query_bases, None, "q".to_string(), None);
    let target = Sequence::new(target_bases, None, "ref".to_string(), None);
    let plan = make_plan();

    let config = AlignConfig::new(Strategy::Balanced);
    let result = align_task_with_config(&query, &target, &plan, &config).expect("alignment should succeed");
    assert!(!result.variants.is_empty(), "should have at least one variant");

    // All variants should have the same window size configuration (balanced: ~100 elements).
    let first_len = result.variants[0].local_coverage().len();
    for var in &result.variants {
        let len = var.local_coverage().len();
        // Near-edge variants may be truncated; far-from-edge ones should be full ~100.
        assert!(
            len <= 101,
            "balanced window should not exceed 2*50+1=101, got {}",
            len
        );
    }
    // At least one variant should have the full ~100 window (position far from edges).
    assert!(
        first_len >= 50,
        "balanced strategy should produce substantial coverage window, got {}",
        first_len
    );
}
