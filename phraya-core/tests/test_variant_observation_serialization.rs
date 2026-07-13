use phraya_core::types::{MateInfo, VariantObservation};
use std::collections::HashMap;

#[test]
fn variant_observation_with_mate_info_serializes() {
    let mate_info = MateInfo::new("mate/2".to_string(), true, 450, true, false, true);

    let mut alleles = HashMap::new();
    alleles.insert(b'A', 10);

    let obs = VariantObservation::new(
        100,
        b'A',
        alleles,
        0.95,
        "10M".to_string(),
        60,
        0,
        vec![10],
        35.0,
        "sample:read".to_string(),
    )
    .with_mate_info(mate_info);

    // Serialize to JSON
    let json = serde_json::to_string(&obs).unwrap();

    // Should contain mate_info fields
    assert!(json.contains("mate_info"), "JSON should contain mate_info");
    assert!(json.contains("mate/2"), "JSON should contain mate_id");
    assert!(
        json.contains("insert_size"),
        "JSON should contain insert_size"
    );
    assert!(
        json.contains("450"),
        "JSON should contain insert_size value"
    );
}

#[test]
fn variant_observation_with_mate_info_deserializes() {
    let mate_info = MateInfo::new("mate/2".to_string(), true, 450, true, false, true);

    let mut alleles = HashMap::new();
    alleles.insert(b'A', 10);

    let obs = VariantObservation::new(
        100,
        b'A',
        alleles.clone(),
        0.95,
        "10M".to_string(),
        60,
        0,
        vec![10],
        35.0,
        "sample:read".to_string(),
    )
    .with_mate_info(mate_info.clone());

    // Round-trip through JSON
    let json = serde_json::to_string(&obs).unwrap();
    let deserialized: VariantObservation = serde_json::from_str(&json).unwrap();

    // Verify mate_info preserved
    assert!(deserialized.mate_info().is_some());
    let deser_mate = deserialized.mate_info().unwrap();
    assert_eq!(deser_mate.mate_id, "mate/2");
    assert_eq!(deser_mate.insert_size, 450);
    assert!(deser_mate.proper_pair);
    assert!(deser_mate.is_first_in_pair);
    assert!(!deser_mate.is_second_in_pair);
    assert!(deser_mate.mate_mapped);
}

#[test]
fn variant_observation_without_mate_info_serializes() {
    let mut alleles = HashMap::new();
    alleles.insert(b'A', 10);

    let obs = VariantObservation::new(
        100,
        b'A',
        alleles,
        0.95,
        "10M".to_string(),
        60,
        0,
        vec![10],
        35.0,
        "sample:read".to_string(),
    );

    // Serialize to JSON
    let json = serde_json::to_string(&obs).unwrap();

    // mate_info should be omitted (skip_serializing_if)
    assert!(
        !json.contains("mate_info"),
        "JSON should not contain mate_info when None"
    );
}

#[test]
fn variant_observation_without_mate_info_deserializes() {
    let mut alleles = HashMap::new();
    alleles.insert(b'A', 10);

    let obs = VariantObservation::new(
        100,
        b'A',
        alleles,
        0.95,
        "10M".to_string(),
        60,
        0,
        vec![10],
        35.0,
        "sample:read".to_string(),
    );

    // Round-trip through JSON
    let json = serde_json::to_string(&obs).unwrap();
    let deserialized: VariantObservation = serde_json::from_str(&json).unwrap();

    // mate_info should be None
    assert!(deserialized.mate_info().is_none());
}

#[test]
fn backward_compatibility_old_json_without_mate_info() {
    // Simulate old .phraya file JSON (before mate_info was added)
    let old_json = r#"{
        "position": 100,
        "ref_base": 65,
        "all_alleles": {"65": 10},
        "confidence": 0.95,
        "cigar": "10M",
        "mapq": 60,
        "edit_distance": 0,
        "local_coverage": [10],
        "avg_base_quality": 35.0,
        "provenance": "sample:read"
    }"#;

    // Should deserialize successfully with mate_info = None
    let deserialized: VariantObservation = serde_json::from_str(old_json).unwrap();

    assert_eq!(deserialized.position(), 100);
    assert!(
        deserialized.mate_info().is_none(),
        "Old JSON should deserialize with mate_info = None"
    );
}

#[test]
fn mate_info_with_zero_insert_size() {
    let mate_info = MateInfo::new(
        "mate/2".to_string(),
        false, // Not proper pair (unmapped)
        0,     // Zero insert size
        true,
        false,
        false,
    );

    let mut alleles = HashMap::new();
    alleles.insert(b'A', 10);

    let obs = VariantObservation::new(
        100,
        b'A',
        alleles,
        0.95,
        "10M".to_string(),
        60,
        0,
        vec![10],
        35.0,
        "sample:read".to_string(),
    )
    .with_mate_info(mate_info);

    // Round-trip
    let json = serde_json::to_string(&obs).unwrap();
    let deserialized: VariantObservation = serde_json::from_str(&json).unwrap();

    let deser_mate = deserialized.mate_info().unwrap();
    assert_eq!(deser_mate.insert_size, 0);
    assert!(!deser_mate.proper_pair);
    assert!(!deser_mate.mate_mapped);
}

#[test]
fn mate_info_with_negative_insert_size() {
    let mate_info = MateInfo::new(
        "mate/2".to_string(),
        true,
        -450, // Negative = mate upstream
        true,
        false,
        true,
    );

    let mut alleles = HashMap::new();
    alleles.insert(b'A', 10);

    let obs = VariantObservation::new(
        100,
        b'A',
        alleles,
        0.95,
        "10M".to_string(),
        60,
        0,
        vec![10],
        35.0,
        "sample:read".to_string(),
    )
    .with_mate_info(mate_info);

    // Round-trip
    let json = serde_json::to_string(&obs).unwrap();
    let deserialized: VariantObservation = serde_json::from_str(&json).unwrap();

    let deser_mate = deserialized.mate_info().unwrap();
    assert_eq!(deser_mate.insert_size, -450);
}
