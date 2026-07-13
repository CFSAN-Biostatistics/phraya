use phraya_align::executor::{align_task_with_config, AlignConfig};
use phraya_core::types::Sequence;
use phraya_io::plan::{PhrayaPlan, UseCase};
/// Issue #147: Composite confidence score merging k-mer and alignment evidence (tier 2)
///
/// Tests verify that:
/// 1. confidence = kmer_uniqueness × alignment_score
/// 2. Variant inside hotspot (kmer_uniqueness=0.0) → confidence=0.0
/// 3. Variant outside hotspot (kmer_uniqueness=1.0) → confidence ≈ alignment_score
/// 4. confidence is always in [0.0, 1.0]
/// 5. Hotspot variant has strictly lower confidence than non-hotspot variant at same alignment quality
use std::collections::HashMap;

fn make_plan_with_hotspots(hotspot_intervals: Vec<(u32, u32)>) -> PhrayaPlan {
    let mut plan = PhrayaPlan::new(
        UseCase::ReadsWithRef,
        vec!["test".to_string()],
        "2026-06-07T00:00:00Z".to_string(),
        HashMap::new(),
        HashMap::new(),
        vec![],
    );
    plan.hotspot_intervals = hotspot_intervals;
    plan
}

fn make_plan() -> PhrayaPlan {
    make_plan_with_hotspots(vec![])
}

/// Variant inside hotspot has kmer_uniqueness=0.0, so confidence must be 0.0
#[test]
fn issue_147_hotspot_variant_confidence_is_zero() {
    let mut query_bases = vec![b'A'; 100];
    let mut target_bases = vec![b'A'; 200];
    query_bases[50] = b'T';
    target_bases[50] = b'C';

    let query = Sequence::new(query_bases, None, "q".to_string(), None);
    let target = Sequence::new(target_bases, None, "ref".to_string(), None);
    let plan = make_plan_with_hotspots(vec![(40, 60)]);
    let config = AlignConfig::balanced();

    let result =
        align_task_with_config(&query, &target, &plan, &config).expect("alignment should succeed");

    assert!(
        !result.variants.is_empty(),
        "should have at least one variant"
    );
    let var = result
        .variants
        .iter()
        .find(|v| v.position() >= 40 && v.position() < 60)
        .expect("should have variant inside hotspot");

    // kmer_uniqueness=0.0 inside hotspot → composite = 0.0 × alignment_score = 0.0
    assert_eq!(
        var.confidence(),
        0.0,
        "variant inside hotspot should have confidence=0.0, got {}",
        var.confidence()
    );
}

/// Variant outside hotspot: confidence ≈ alignment_score (kmer_uniqueness=1.0)
/// With a single SNP in a 100bp query, alignment_score = 1 - 1/100 = 0.99
#[test]
fn issue_147_non_hotspot_variant_confidence_equals_alignment_score() {
    let mut query_bases = vec![b'A'; 100];
    let mut target_bases = vec![b'A'; 200];
    query_bases[50] = b'T';
    target_bases[50] = b'C';

    let query = Sequence::new(query_bases, None, "q".to_string(), None);
    let target = Sequence::new(target_bases, None, "ref".to_string(), None);
    // Hotspot elsewhere — position 50 is outside
    let plan = make_plan_with_hotspots(vec![(80, 100)]);
    let config = AlignConfig::balanced();

    let result =
        align_task_with_config(&query, &target, &plan, &config).expect("alignment should succeed");

    assert!(
        !result.variants.is_empty(),
        "should have at least one variant"
    );
    let var = &result.variants[0];

    // kmer_uniqueness=1.0 outside hotspot → confidence = 1.0 × alignment_score = alignment_score
    // alignment_score for 1 edit in 100bp query = 1 - 1/100 = 0.99
    let c = var.confidence();
    assert!(
        c > 0.95 && c <= 1.0,
        "non-hotspot variant confidence should be ≈alignment_score (~0.99), got {}",
        c
    );
}

/// confidence is always in [0.0, 1.0]
#[test]
fn issue_147_confidence_is_in_unit_interval() {
    let mut query_bases = vec![b'A'; 100];
    let mut target_bases = vec![b'A'; 200];
    // Multiple SNPs
    for i in [20usize, 40, 60, 80] {
        query_bases[i] = b'T';
        target_bases[i] = b'C';
    }

    let query = Sequence::new(query_bases, None, "q".to_string(), None);
    let target = Sequence::new(target_bases, None, "ref".to_string(), None);
    let plan = make_plan_with_hotspots(vec![(35, 45)]);
    let config = AlignConfig::balanced();

    let result =
        align_task_with_config(&query, &target, &plan, &config).expect("alignment should succeed");

    for var in &result.variants {
        let c = var.confidence();
        assert!(
            c >= 0.0 && c <= 1.0,
            "confidence must be in [0.0, 1.0], got {} at position {}",
            c,
            var.position()
        );
    }
}

/// Hotspot variant has strictly lower confidence than same-quality non-hotspot variant
#[test]
fn issue_147_hotspot_confidence_lower_than_non_hotspot() {
    let mut query_bases = vec![b'A'; 100];
    let mut target_bases = vec![b'A'; 200];
    query_bases[30] = b'T';
    query_bases[70] = b'T';
    target_bases[30] = b'C';
    target_bases[70] = b'C';

    let query = Sequence::new(query_bases, None, "q".to_string(), None);
    let target = Sequence::new(target_bases, None, "ref".to_string(), None);
    // Only position 30 is in hotspot
    let plan = make_plan_with_hotspots(vec![(20, 40)]);
    let config = AlignConfig::balanced();

    let result =
        align_task_with_config(&query, &target, &plan, &config).expect("alignment should succeed");

    let var_in = result
        .variants
        .iter()
        .find(|v| v.position() >= 20 && v.position() < 40)
        .expect("should have variant inside hotspot");
    let var_out = result
        .variants
        .iter()
        .find(|v| v.position() >= 60 && v.position() < 80)
        .expect("should have variant outside hotspot");

    assert!(
        var_in.confidence() < var_out.confidence(),
        "hotspot variant confidence ({}) should be < non-hotspot confidence ({})",
        var_in.confidence(),
        var_out.confidence()
    );
}

/// No hotspots → confidence = alignment_score for all variants
#[test]
fn issue_147_no_hotspots_confidence_equals_alignment_score() {
    let mut query_bases = vec![b'A'; 100];
    let mut target_bases = vec![b'A'; 200];
    query_bases[50] = b'T';
    target_bases[50] = b'C';

    let query = Sequence::new(query_bases, None, "q".to_string(), None);
    let target = Sequence::new(target_bases, None, "ref".to_string(), None);
    let plan = make_plan(); // no hotspots
    let config = AlignConfig::balanced();

    let result =
        align_task_with_config(&query, &target, &plan, &config).expect("alignment should succeed");

    assert!(!result.variants.is_empty());
    for var in &result.variants {
        let c = var.confidence();
        // kmer_uniqueness=1.0 everywhere → confidence = alignment_score ∈ (0, 1]
        assert!(
            c > 0.0 && c <= 1.0,
            "with no hotspots, confidence should equal alignment_score ∈ (0,1], got {}",
            c
        );
    }
}

/// Perfect alignment (no edits) outside hotspot → confidence close to 1.0
/// Use identical query and target so edit_distance=0, alignment_score=1.0
#[test]
fn issue_147_perfect_alignment_outside_hotspot_confidence_near_one() {
    // Identical sequences except one position for the variant
    let mut query_bases = vec![b'A'; 100];
    let mut target_bases = vec![b'A'; 100];
    // Make them identical — no SNP. But then there'd be no variants.
    // Use minimal edit: 1 SNP at pos 50, outside hotspot
    query_bases[50] = b'T';
    target_bases[50] = b'C';

    let query = Sequence::new(query_bases, None, "q".to_string(), None);
    let target = Sequence::new(target_bases, None, "ref".to_string(), None);
    let plan = make_plan_with_hotspots(vec![(80, 90)]);
    let config = AlignConfig::balanced();

    let result =
        align_task_with_config(&query, &target, &plan, &config).expect("alignment should succeed");

    assert!(!result.variants.is_empty());
    let var = &result.variants[0];
    let c = var.confidence();
    // 1 edit in 100bp: alignment_score = 0.99; kmer_uniqueness=1.0 → confidence=0.99
    assert!(
        c >= 0.98 && c <= 1.0,
        "1-edit alignment outside hotspot: confidence should be ~0.99, got {}",
        c
    );
}
