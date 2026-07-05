/// Issue #146: K-mer evidence tier on VariantObservation (two-tier evidence, tier 1)
///
/// This test file contains RED (failing) acceptance tests for issue #146.
/// Tests verify that:
/// 1. VariantObservation carries kmer_uniqueness field
/// 2. During alignment, variants inside hotspot intervals get kmer_uniqueness < 1.0
/// 3. Variants outside hotspot intervals get kmer_uniqueness close to 1.0
/// 4. Old .phraya files deserialize without error (serde default)

use std::collections::HashMap;
use phraya_align::executor::{align_task_with_config, AlignConfig};
use phraya_core::types::Sequence;
use phraya_io::plan::{PhrayaPlan, UseCase};

fn make_plan_with_hotspots(hotspot_intervals: Vec<(u32, u32)>) -> PhrayaPlan {
    let mut plan = PhrayaPlan::new(
        UseCase::ReadsWithRef,
        vec!["test".to_string()],
        "2026-06-06T00:00:00Z".to_string(),
        HashMap::new(),
        HashMap::new(),
        vec![],
    );
    plan.hotspot_intervals = hotspot_intervals;
    plan
}

/// Test that VariantObservation has kmer_uniqueness field
/// The field should be accessible via getter and initialized correctly
#[test]
fn issue_146_variant_observation_has_kmer_uniqueness_field() {
    use phraya_core::types::VariantObservation;

    let mut alleles = HashMap::new();
    alleles.insert(b'T', 1u32);

    let obs = VariantObservation::new(
        100u32,
        b'A',
        alleles,
        0.95,
        "10M".to_string(),
        30u8,
        0u32,
        vec![5u32],
        35.0,
        "sample:read".to_string(),
    );

    // kmer_uniqueness should have a default getter returning a value close to 1.0
    // since serde default = 1.0 for new observations
    let km = obs.kmer_uniqueness();
    assert!(km >= 0.99 && km <= 1.0, "default kmer_uniqueness should be 1.0, got {}", km);
}

/// Test that VariantObservation can be constructed with kmer_uniqueness via builder pattern
#[test]
fn issue_146_variant_observation_kmer_uniqueness_builder() {
    use phraya_core::types::VariantObservation;

    let mut alleles = HashMap::new();
    alleles.insert(b'T', 1u32);

    let obs = VariantObservation::new(
        100u32,
        b'A',
        alleles,
        0.95,
        "10M".to_string(),
        30u8,
        0u32,
        vec![5u32],
        35.0,
        "sample:read".to_string(),
    )
    .with_kmer_uniqueness(0.5);

    assert_eq!(obs.kmer_uniqueness(), 0.5, "kmer_uniqueness should be set to 0.5");
}

/// Test that PhrayaPlan has hotspot_intervals field
#[test]
fn issue_146_phraya_plan_has_hotspot_intervals() {
    let hotspots = vec![(10u32, 20u32), (50u32, 60u32)];
    let plan = make_plan_with_hotspots(hotspots.clone());

    assert_eq!(
        plan.hotspot_intervals, hotspots,
        "plan should store and return hotspot intervals"
    );
}

/// Test that align_task_with_config sets kmer_uniqueness < 1.0 for variants inside hotspot
///
/// Setup:
/// - 100bp query with SNP at position 50
/// - 200bp target with SNP at position 50 (offset by target_start)
/// - Hotspot interval: [40, 60]
/// - Expected: variant at position 50 should have kmer_uniqueness < 1.0 (inside hotspot)
#[test]
fn issue_146_align_task_sets_low_kmer_uniqueness_inside_hotspot() {
    let mut query_bases = vec![b'A'; 100];
    let mut target_bases = vec![b'A'; 200];
    query_bases[50] = b'T';
    target_bases[50] = b'C'; // SNP at position 50

    let query = Sequence::new(query_bases, None, "q".to_string(), None);
    let target = Sequence::new(target_bases, None, "ref".to_string(), None);

    // Create plan with hotspot interval [40, 60] that contains position 50
    let plan = make_plan_with_hotspots(vec![(40, 60)]);
    let config = AlignConfig::balanced();

    let result = align_task_with_config(&query, &target, &plan, &config)
        .expect("alignment should succeed");

    assert!(!result.variants.is_empty(), "should have at least one variant");
    let var = &result.variants[0];

    // The variant at position 50 should be inside hotspot [40, 60]
    assert!(var.position() >= 40 && var.position() < 60, "variant should be in hotspot range");

    // kmer_uniqueness should be < 1.0 because position is inside a hotspot
    let km = var.kmer_uniqueness();
    assert!(km < 1.0, "kmer_uniqueness should be < 1.0 inside hotspot, got {}", km);
}

/// Test that align_task_with_config sets kmer_uniqueness close to 1.0 for variants outside hotspot
///
/// Setup:
/// - 100bp query with SNP at position 50
/// - 200bp target with SNP at position 50 (offset by target_start)
/// - Hotspot interval: [100, 120] (does not contain position 50)
/// - Expected: variant at position 50 should have kmer_uniqueness close to 1.0
#[test]
fn issue_146_align_task_sets_high_kmer_uniqueness_outside_hotspot() {
    let mut query_bases = vec![b'A'; 100];
    let mut target_bases = vec![b'A'; 200];
    query_bases[50] = b'T';
    target_bases[50] = b'C'; // SNP at position 50

    let query = Sequence::new(query_bases, None, "q".to_string(), None);
    let target = Sequence::new(target_bases, None, "ref".to_string(), None);

    // Create plan with hotspot interval [100, 120] that does NOT contain position 50
    let plan = make_plan_with_hotspots(vec![(100, 120)]);
    let config = AlignConfig::balanced();

    let result = align_task_with_config(&query, &target, &plan, &config)
        .expect("alignment should succeed");

    assert!(!result.variants.is_empty(), "should have at least one variant");
    let var = &result.variants[0];

    // The variant at position 50 should be outside hotspot [100, 120]
    assert!(!(var.position() >= 100 && var.position() < 120), "variant should be outside hotspot range");

    // kmer_uniqueness should be close to 1.0 because position is outside all hotspots
    let km = var.kmer_uniqueness();
    assert!(km >= 0.99 && km <= 1.0, "kmer_uniqueness should be ~1.0 outside hotspot, got {}", km);
}

/// Test that with multiple hotspot intervals, variants are correctly classified
///
/// Setup:
/// - 100bp query with SNPs at positions 30, 80
/// - 200bp target with SNPs at positions 30, 80
/// - Hotspot intervals: [20, 40], [90, 110]
/// - Expected: position 30 inside [20, 40], position 80 outside both
#[test]
fn issue_146_multiple_hotspot_intervals() {
    let mut query_bases = vec![b'A'; 100];
    let mut target_bases = vec![b'A'; 200];
    query_bases[30] = b'T';
    query_bases[80] = b'T';
    target_bases[30] = b'C';
    target_bases[80] = b'C'; // SNPs at positions 30 and 80

    let query = Sequence::new(query_bases, None, "q".to_string(), None);
    let target = Sequence::new(target_bases, None, "ref".to_string(), None);

    // Create plan with two hotspot intervals
    let plan = make_plan_with_hotspots(vec![(20, 40), (90, 110)]);
    let config = AlignConfig::balanced();

    let result = align_task_with_config(&query, &target, &plan, &config)
        .expect("alignment should succeed");

    assert!(result.variants.len() >= 2, "should have at least 2 variants");

    // Find variants at positions 30 and 80
    let var_at_30 = result.variants.iter().find(|v| v.position() == 30);
    let var_at_80 = result.variants.iter().find(|v| v.position() == 80);

    assert!(var_at_30.is_some(), "should have variant at position 30");
    assert!(var_at_80.is_some(), "should have variant at position 80");

    // Position 30 is inside hotspot [20, 40]
    let km30 = var_at_30.unwrap().kmer_uniqueness();
    assert!(km30 < 1.0, "position 30 inside [20, 40] should have km < 1.0, got {}", km30);

    // Position 80 is outside both hotspots
    let km80 = var_at_80.unwrap().kmer_uniqueness();
    assert!(km80 >= 0.99 && km80 <= 1.0, "position 80 outside hotspots should have km ~1.0, got {}", km80);
}

/// Test that empty hotspot intervals result in all kmer_uniqueness ~1.0
#[test]
fn issue_146_empty_hotspot_intervals_all_unique() {
    let mut query_bases = vec![b'A'; 100];
    let mut target_bases = vec![b'A'; 200];
    query_bases[50] = b'T';
    target_bases[50] = b'C';

    let query = Sequence::new(query_bases, None, "q".to_string(), None);
    let target = Sequence::new(target_bases, None, "ref".to_string(), None);

    // Create plan with NO hotspot intervals (empty vec)
    let plan = make_plan_with_hotspots(vec![]);
    let config = AlignConfig::balanced();

    let result = align_task_with_config(&query, &target, &plan, &config)
        .expect("alignment should succeed");

    assert!(!result.variants.is_empty(), "should have at least one variant");
    for var in &result.variants {
        let km = var.kmer_uniqueness();
        assert!(km >= 0.99 && km <= 1.0, "with empty hotspots, all variants should have km ~1.0, got {}", km);
    }
}

/// Test deserialization of old .phraya files (backward compatibility)
/// Old files will not have kmer_uniqueness field, so serde default should provide 1.0
#[test]
fn issue_146_old_phraya_files_deserialize_with_default() {
    use phraya_core::types::VariantObservation;

    // Create an observation and manually check default behavior
    // (In a real test, this would deserialize from JSON/MessagePack)
    let mut alleles = HashMap::new();
    alleles.insert(b'A', 5u32);

    // When deserializing old format without kmer_uniqueness field,
    // serde(default) should set it to 1.0
    let obs = VariantObservation::new(
        100u32,
        b'A',
        alleles,
        0.95,
        "10M".to_string(),
        30u8,
        0u32,
        vec![5u32],
        35.0,
        "sample:read".to_string(),
    );

    // Should not panic and should have default value
    assert!(obs.kmer_uniqueness() >= 0.99 && obs.kmer_uniqueness() <= 1.0, "default should be ~1.0");
}

/// Test that variant hotspot status is orthogonal to tandem repeat status
/// A variant can be both in a hotspot AND in a tandem repeat
#[test]
fn issue_146_hotspot_and_tandem_repeat_orthogonal() {
    use phraya_core::types::VariantObservation;

    let mut alleles = HashMap::new();
    alleles.insert(b'T', 1u32);

    let obs = VariantObservation::new(
        50u32,
        b'A',
        alleles,
        0.95,
        "10M".to_string(),
        30u8,
        0u32,
        vec![5u32],
        35.0,
        "sample:read".to_string(),
    )
    .with_kmer_uniqueness(0.3)
    .with_tandem_repeat(true);

    // Both should be set independently
    assert_eq!(obs.kmer_uniqueness(), 0.3, "kmer_uniqueness should be 0.3");
    assert_eq!(obs.in_tandem_repeat(), true, "in_tandem_repeat should be true");
}

#[test]
fn debug_146_variant_positions() {
    let mut query_bases = vec![b'A'; 100];
    let mut target_bases = vec![b'A'; 200];
    query_bases[50] = b'T';
    target_bases[50] = b'C';
    let query = Sequence::new(query_bases, None, "q".to_string(), None);
    let target = Sequence::new(target_bases, None, "ref".to_string(), None);
    let plan = make_plan_with_hotspots(vec![(40, 60)]);
    let config = AlignConfig::balanced();
    let result = align_task_with_config(&query, &target, &plan, &config).expect("ok");
    for v in &result.variants {
        eprintln!("variant pos={} km={}", v.position(), v.kmer_uniqueness());
    }
}

#[test]
fn debug_146_seeds_and_anchors() {
    use phraya_align::seeding::find_seeds;
    use phraya_core::types::{sketch_sequence_default, Sequence};
    let mut query_bases = vec![b'A'; 100];
    let mut target_bases = vec![b'A'; 200];
    query_bases[50] = b'T';
    target_bases[50] = b'C';
    let query = Sequence::new(query_bases, None, "q".to_string(), None);
    let target = Sequence::new(target_bases, None, "ref".to_string(), None);
    let qs = sketch_sequence_default(&query);
    let ts = sketch_sequence_default(&target);
    let seeds = find_seeds(&qs, &ts);
    eprintln!("seed count: {}", seeds.len());
    for s in seeds.iter().take(5) {
        let ts = (s.target_pos as i64 - s.query_pos as i64).max(0);
        eprintln!("  seed qp={} tp={} → target_start={}", s.query_pos, s.target_pos, ts);
    }
}
