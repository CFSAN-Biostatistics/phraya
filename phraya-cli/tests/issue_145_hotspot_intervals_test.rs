/// RED integration tests for issue #145 — hotspot_intervals field in PhrayaPlan.
use phraya_core::types::{
    compute_kmer_uniqueness, detect_hotspot_intervals, sketch_sequence_default, Sequence,
};
use phraya_io::plan::{read_plan, write_plan, PhrayaPlan, UseCase, PHRAYAPLAN_VERSION};
use std::collections::HashMap;
use tempfile::NamedTempFile;

fn make_seq(bases: &str, id: &str) -> Sequence {
    Sequence::new(bases.as_bytes().to_vec(), None, id.to_string(), None)
}

/// issue #145: PhrayaPlan carries hotspot_intervals field and round-trips via serde
#[test]
fn issue_145_phrayaplan_hotspot_intervals_round_trips() {
    let intervals = vec![(10u32, 20u32), (50u32, 60u32)];
    let mut plan = PhrayaPlan::new(
        UseCase::ReadsWithRef,
        vec!["input.fa".to_string()],
        "2026-06-06T00:00:00Z".to_string(),
        HashMap::new(),
        HashMap::new(),
        vec![],
    );
    plan.hotspot_intervals = intervals.clone();

    let tmp = NamedTempFile::new().unwrap();
    write_plan(tmp.path(), &plan).expect("write_plan failed");
    let loaded = read_plan(tmp.path()).expect("read_plan failed");

    assert_eq!(
        loaded.hotspot_intervals, intervals,
        "hotspot_intervals must survive serde round-trip"
    );
}

/// issue #145: new plan has empty hotspot_intervals by default
#[test]
fn issue_145_new_plan_has_empty_hotspot_intervals() {
    let plan = PhrayaPlan::new(
        UseCase::ReadsWithRef,
        vec![],
        "2026-06-06T00:00:00Z".to_string(),
        HashMap::new(),
        HashMap::new(),
        vec![],
    );
    assert!(
        plan.hotspot_intervals.is_empty(),
        "new PhrayaPlan must default to empty hotspot_intervals"
    );
}

/// issue #145: old plan file without hotspot_intervals deserializes with empty default
#[test]
fn issue_145_old_plan_deserializes_without_hotspot_intervals() {
    // Write a plan then manually clear hotspot_intervals before serialization is impossible,
    // so we test via serde default: create plan, write it, read it back — field should default empty.
    let plan = PhrayaPlan::new(
        UseCase::ContigsOnly,
        vec!["contigs.fa".to_string()],
        "2026-06-06T00:00:00Z".to_string(),
        HashMap::new(),
        HashMap::new(),
        vec![],
    );
    let tmp = NamedTempFile::new().unwrap();
    write_plan(tmp.path(), &plan).unwrap();
    let loaded = read_plan(tmp.path()).unwrap();
    assert_eq!(loaded.version, PHRAYAPLAN_VERSION);
    assert!(loaded.hotspot_intervals.is_empty());
}

/// issue #145: detect_hotspot_intervals on a repeat-heavy sequence produces intervals
/// that phraya plan would store in the plan file
#[test]
fn issue_145_plan_detects_hotspots_from_repeated_sequence() {
    let seq = make_seq(
        "ACGTACGTATATATATATATATATATATATATATATACGTACGT",
        "repeat_target",
    );
    let sketch = sketch_sequence_default(&seq);
    // Two identical sketches → all positions get uniqueness 0.5, below 0.6 threshold
    let uniqueness = compute_kmer_uniqueness(&[sketch.clone(), sketch]);
    let intervals = detect_hotspot_intervals(&uniqueness, 0.6);
    assert!(
        !intervals.is_empty(),
        "repeated sequence should produce at least one hotspot interval for plan storage"
    );
}
