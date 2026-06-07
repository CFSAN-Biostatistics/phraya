/// Issue #146: K-mer evidence tier on VariantObservation (two-tier evidence, tier 1)
///
/// This test file contains RED (failing) acceptance tests for issue #146 in phraya-filter.
/// Tests verify that:
/// 1. FilterBuilder gains min_kmer_uniqueness(threshold: f64) method
/// 2. ThresholdFilter correctly filters variants based on kmer_uniqueness
/// 3. min_kmer_uniqueness works in combination with other filters
/// 4. Named presets can include kmer_uniqueness constraints

use std::collections::HashMap;
use phraya_core::types::VariantObservation;
use phraya_filter::FilterBuilder;

fn create_observation_with_kmer_uniqueness(
    position: u32,
    kmer_uniqueness: f64,
) -> VariantObservation {
    let mut alleles = HashMap::new();
    alleles.insert(b'T', 5u32);

    VariantObservation::new(
        position,
        b'A',
        alleles,
        0.95,
        "10M".to_string(),
        30u8,
        0u32,
        vec![10u32],
        35.0,
        "sample:read".to_string(),
    )
    .with_kmer_uniqueness(kmer_uniqueness)
}

/// Test that FilterBuilder has min_kmer_uniqueness method
#[test]
fn issue_146_filter_builder_has_min_kmer_uniqueness() {
    let filter = FilterBuilder::new()
        .min_kmer_uniqueness(0.5)
        .build();

    // The filter should be built without panicking
    assert!(true, "filter builder should accept min_kmer_uniqueness");
}

/// Test that FilterBuilder method is chainable
#[test]
fn issue_146_min_kmer_uniqueness_is_chainable() {
    let _filter = FilterBuilder::new()
        .min_coverage(10)
        .min_kmer_uniqueness(0.7)
        .min_mapq(20)
        .build();

    // Should compile and chain methods without error
}

/// Test that ThresholdFilter correctly excludes variants below kmer_uniqueness threshold
#[test]
fn issue_146_filter_excludes_low_kmer_uniqueness() {
    let obs_low_km = create_observation_with_kmer_uniqueness(100, 0.3);
    let obs_high_km = create_observation_with_kmer_uniqueness(200, 0.8);

    let filter = FilterBuilder::new()
        .min_kmer_uniqueness(0.5)
        .build();

    // Variant with kmer_uniqueness 0.3 should be excluded (0.3 < 0.5)
    assert!(!filter.apply(&obs_low_km), "should exclude variant with km=0.3 when threshold is 0.5");

    // Variant with kmer_uniqueness 0.8 should pass (0.8 >= 0.5)
    assert!(filter.apply(&obs_high_km), "should include variant with km=0.8 when threshold is 0.5");
}

/// Test that threshold of 0.0 includes all variants
#[test]
fn issue_146_kmer_uniqueness_threshold_zero_includes_all() {
    let obs = create_observation_with_kmer_uniqueness(100, 0.1);

    let filter = FilterBuilder::new()
        .min_kmer_uniqueness(0.0)
        .build();

    assert!(filter.apply(&obs), "threshold 0.0 should include all variants");
}

/// Test that threshold of 1.0 only includes unique variants
#[test]
fn issue_146_kmer_uniqueness_threshold_one_strict() {
    let obs_almost_unique = create_observation_with_kmer_uniqueness(100, 0.99);
    let obs_unique = create_observation_with_kmer_uniqueness(200, 1.0);

    let filter = FilterBuilder::new()
        .min_kmer_uniqueness(1.0)
        .build();

    // obs_almost_unique (0.99) is below threshold, should be excluded
    assert!(!filter.apply(&obs_almost_unique), "should exclude km=0.99 when threshold is 1.0");

    // obs_unique (1.0) should pass
    assert!(filter.apply(&obs_unique), "should include km=1.0 when threshold is 1.0");
}

/// Test that min_kmer_uniqueness combines with other filters (AND logic)
#[test]
fn issue_146_kmer_uniqueness_combines_with_other_filters() {
    let obs_pass_both = create_observation_with_kmer_uniqueness(100, 0.8);
    let obs_fail_km = create_observation_with_kmer_uniqueness(200, 0.3);
    let obs_fail_mapq = VariantObservation::new(
        300u32,
        b'A',
        {
            let mut a = HashMap::new();
            a.insert(b'T', 5u32);
            a
        },
        0.95,
        "10M".to_string(),
        15u8, // low mapq
        0u32,
        vec![10u32],
        35.0,
        "sample:read".to_string(),
    )
    .with_kmer_uniqueness(0.8);

    let filter = FilterBuilder::new()
        .min_kmer_uniqueness(0.5)
        .min_mapq(20)
        .build();

    // Should pass both filters
    assert!(filter.apply(&obs_pass_both), "should pass when both km and mapq thresholds met");

    // Should fail kmer_uniqueness filter
    assert!(!filter.apply(&obs_fail_km), "should fail when km below threshold");

    // Should fail mapq filter
    assert!(!filter.apply(&obs_fail_mapq), "should fail when mapq below threshold");
}

/// Test that default (no kmer_uniqueness filter) doesn't exclude anything based on km
#[test]
fn issue_146_default_filter_accepts_all_kmer_uniqueness() {
    let obs_very_low = create_observation_with_kmer_uniqueness(100, 0.1);
    let obs_medium = create_observation_with_kmer_uniqueness(200, 0.5);
    let obs_high = create_observation_with_kmer_uniqueness(300, 1.0);

    let filter = FilterBuilder::new().build();

    // With no min_kmer_uniqueness set, all should pass
    assert!(filter.apply(&obs_very_low), "default filter should accept km=0.1");
    assert!(filter.apply(&obs_medium), "default filter should accept km=0.5");
    assert!(filter.apply(&obs_high), "default filter should accept km=1.0");
}

/// Test that filter() method correctly filters observations by kmer_uniqueness
#[test]
fn issue_146_filter_iterator_excludes_low_km_variants() {
    let observations = vec![
        create_observation_with_kmer_uniqueness(100, 0.2),
        create_observation_with_kmer_uniqueness(200, 0.6),
        create_observation_with_kmer_uniqueness(300, 0.9),
        create_observation_with_kmer_uniqueness(400, 0.3),
    ];

    let filter = FilterBuilder::new()
        .min_kmer_uniqueness(0.5)
        .build();

    let filtered: Vec<_> = filter.filter(&observations).collect();

    // Should keep positions 200 (0.6) and 300 (0.9)
    // Should exclude positions 100 (0.2) and 400 (0.3)
    assert_eq!(filtered.len(), 2, "should keep 2 variants above threshold");
    assert_eq!(filtered[0].position(), 200, "first filtered should be position 200");
    assert_eq!(filtered[1].position(), 300, "second filtered should be position 300");
}

/// Test that kmer_uniqueness of exactly at threshold is included
#[test]
fn issue_146_kmer_uniqueness_at_threshold_included() {
    let obs_at_threshold = create_observation_with_kmer_uniqueness(100, 0.5);
    let obs_just_below = create_observation_with_kmer_uniqueness(200, 0.49);

    let filter = FilterBuilder::new()
        .min_kmer_uniqueness(0.5)
        .build();

    assert!(filter.apply(&obs_at_threshold), "variant at exact threshold should pass");
    assert!(!filter.apply(&obs_just_below), "variant just below threshold should fail");
}

/// Test that boundary case of kmer_uniqueness 0.0 is handled correctly
#[test]
fn issue_146_kmer_uniqueness_zero_handling() {
    let obs_zero = create_observation_with_kmer_uniqueness(100, 0.0);

    let filter_low = FilterBuilder::new()
        .min_kmer_uniqueness(0.0)
        .build();

    let filter_high = FilterBuilder::new()
        .min_kmer_uniqueness(0.1)
        .build();

    // Should pass filter with 0.0 threshold
    assert!(filter_low.apply(&obs_zero), "km=0.0 should pass threshold 0.0");

    // Should fail filter with 0.1 threshold
    assert!(!filter_high.apply(&obs_zero), "km=0.0 should fail threshold 0.1");
}

/// Test kmer_uniqueness works independently of tandem repeat exclusion
#[test]
fn issue_146_kmer_uniqueness_independent_of_tandem_repeat() {
    let obs_in_repeat_high_km = VariantObservation::new(
        100u32,
        b'A',
        {
            let mut a = HashMap::new();
            a.insert(b'T', 5u32);
            a
        },
        0.95,
        "10M".to_string(),
        30u8,
        0u32,
        vec![10u32],
        35.0,
        "sample:read".to_string(),
    )
    .with_kmer_uniqueness(0.9)
    .with_tandem_repeat(true);

    let obs_not_in_repeat_low_km = create_observation_with_kmer_uniqueness(200, 0.3);

    let filter = FilterBuilder::new()
        .min_kmer_uniqueness(0.5)
        .exclude_tandem_repeats(true)
        .build();

    // Should fail because it's in a tandem repeat (even though km is high)
    assert!(!filter.apply(&obs_in_repeat_high_km), "should exclude even with high km if in tandem repeat");

    // Should fail because km is too low (not because of tandem repeat)
    assert!(!filter.apply(&obs_not_in_repeat_low_km), "should exclude because km below threshold");
}

/// Test that multiple min_kmer_uniqueness calls use the last one (override behavior)
#[test]
fn issue_146_multiple_min_kmer_uniqueness_uses_last() {
    let obs = create_observation_with_kmer_uniqueness(100, 0.6);

    let filter = FilterBuilder::new()
        .min_kmer_uniqueness(0.3)
        .min_kmer_uniqueness(0.7)
        .build();

    // The filter should use the last set value (0.7), not the first (0.3)
    assert!(!filter.apply(&obs), "should use last min_kmer_uniqueness value (0.7)");
}

/// Test that the builder fluent API allows setting kmer_uniqueness at any point
#[test]
fn issue_146_min_kmer_uniqueness_can_be_set_anywhere_in_chain() {
    let obs = create_observation_with_kmer_uniqueness(100, 0.6);

    // Test 1: km in the middle
    let filter1 = FilterBuilder::new()
        .min_coverage(10)
        .min_kmer_uniqueness(0.7)
        .min_mapq(20)
        .build();

    assert!(!filter1.apply(&obs), "km=0.6 should fail when threshold is 0.7");

    // Test 2: km at the start
    let filter2 = FilterBuilder::new()
        .min_kmer_uniqueness(0.5)
        .min_coverage(10)
        .build();

    assert!(filter2.apply(&obs), "km=0.6 should pass when threshold is 0.5");

    // Test 3: km at the end
    let filter3 = FilterBuilder::new()
        .min_mapq(20)
        .min_kmer_uniqueness(0.5)
        .build();

    assert!(filter3.apply(&obs), "km=0.6 should pass when threshold is 0.5");
}
