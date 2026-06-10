use phraya_core::types::MateInfo;

#[test]
fn mate_info_creation() {
    let mate_info = MateInfo::new(
        "read123/2".to_string(),
        true,  // proper_pair
        450,   // insert_size
        true,  // is_first_in_pair
        false, // is_second_in_pair
        true,  // mate_mapped
    );

    assert_eq!(mate_info.mate_id, "read123/2");
    assert!(mate_info.proper_pair);
    assert_eq!(mate_info.insert_size, 450);
    assert!(mate_info.is_first_in_pair);
    assert!(!mate_info.is_second_in_pair);
    assert!(mate_info.mate_mapped);
}

#[test]
fn is_discordant_detects_outliers_beyond_threshold() {
    let mate_info = MateInfo::new(
        "read1/2".to_string(),
        true,
        800, // Large insert size
        true,
        false,
        true,
    );

    // Mean=400, StdDev=50, 3σ=150 → threshold at 400+150=550
    // 800 > 550 → discordant
    assert!(mate_info.is_discordant(400, 50, 3.0));
}

#[test]
fn is_discordant_within_threshold() {
    let mate_info = MateInfo::new(
        "read2/2".to_string(),
        true,
        450, // Close to mean
        true,
        false,
        true,
    );

    // Mean=400, StdDev=50, 3σ=150
    // |450 - 400| = 50 < 150 → concordant
    assert!(!mate_info.is_discordant(400, 50, 3.0));
}

#[test]
fn is_discordant_symmetric_around_mean() {
    let mate_info_above = MateInfo::new("r1/2".to_string(), true, 600, true, false, true);
    let mate_info_below = MateInfo::new("r2/2".to_string(), true, 200, true, false, true);

    let mean = 400;
    let std_dev = 50;

    // Both are 200 away from mean (beyond 3σ=150)
    assert!(mate_info_above.is_discordant(mean, std_dev, 3.0));
    assert!(mate_info_below.is_discordant(mean, std_dev, 3.0));
}

#[test]
fn is_discordant_handles_zero_insert_size() {
    let mate_info = MateInfo::new(
        "unmapped/2".to_string(),
        false,
        0, // Unmapped or unpaired
        true,
        false,
        false,
    );

    // TLEN=0 should not be flagged as discordant
    assert!(!mate_info.is_discordant(400, 50, 3.0));
}

#[test]
fn is_discordant_custom_sigma_threshold() {
    let mate_info = MateInfo::new("read/2".to_string(), true, 650, true, false, true);

    // Mean=400, StdDev=50
    // At 3σ: 400 ± 150 = [250, 550] → 650 is discordant
    assert!(mate_info.is_discordant(400, 50, 3.0));

    // At 5σ: 400 ± 250 = [150, 650] → 650 is NOT discordant (right at boundary)
    assert!(!mate_info.is_discordant(400, 50, 5.0));
}

#[test]
fn is_discordant_negative_insert_size() {
    // Negative TLEN indicates mate upstream (reverse orientation)
    let mate_info = MateInfo::new("read/2".to_string(), true, -450, true, false, true);

    // Should use absolute value: |-450| = 450
    // Mean=400, StdDev=50, 3σ=150
    // |450 - 400| = 50 < 150 → concordant
    assert!(!mate_info.is_discordant(400, 50, 3.0));
}
