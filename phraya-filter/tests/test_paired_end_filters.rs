use phraya_core::types::{MateInfo, VariantObservation};
use phraya_filter::FilterBuilder;
use phraya_io::plan::InsertSizeDistribution;
use std::collections::HashMap;

/// Helper to create test VariantObservation
fn create_test_obs(insert_size: i32, proper_pair: bool, mate_mapped: bool) -> VariantObservation {
    let mate_info = MateInfo::new(
        "mate/2".to_string(),
        proper_pair,
        insert_size,
        true,  // is_first_in_pair
        false, // is_second_in_pair
        mate_mapped,
    );

    let mut alleles = HashMap::new();
    alleles.insert(b'A', 10);

    VariantObservation::new(
        100,             // position
        b'A',            // ref_base
        alleles,
        0.95,            // confidence
        "10M".to_string(), // cigar
        60,              // mapq
        0,               // edit_distance
        vec![10],        // local_coverage
        35.0,            // avg_base_quality
        "sample:read".to_string(), // provenance
    )
    .with_mate_info(mate_info)
    .with_pair_counts(1, if proper_pair { 1 } else { 0 })
}

#[test]
fn filter_requires_proper_pairs() {
    let filter = FilterBuilder::new()
        .require_proper_pairs(1.0)
        .build();

    let obs_proper = create_test_obs(400, true, true);
    let obs_improper = create_test_obs(400, false, true);

    assert!(filter.apply(&obs_proper), "Should pass proper pair");
    assert!(!filter.apply(&obs_improper), "Should fail improper pair");
}

#[test]
fn filter_excludes_discordant_pairs() {
    let dist = InsertSizeDistribution {
        mean: 400,
        std_dev: 50,
        orientation: "FR".to_string(),
        sample_size: 1000,
    };

    let filter = FilterBuilder::new()
        .exclude_discordant_pairs(true)
        .with_insert_distribution(dist)
        .build();

    // Concordant: 450 within 400 ± 150 (3σ)
    let obs_concordant = create_test_obs(450, true, true);
    assert!(filter.apply(&obs_concordant), "Should pass concordant pair");

    // Discordant: 800 > 400+150=550
    let obs_discordant = create_test_obs(800, true, true);
    assert!(!filter.apply(&obs_discordant), "Should fail discordant pair");
}

#[test]
fn filter_excludes_discordant_without_distribution_passes() {
    // If no distribution provided, discordant filter has no effect
    let filter = FilterBuilder::new()
        .exclude_discordant_pairs(true)
        .build();

    let obs = create_test_obs(800, true, true);
    assert!(filter.apply(&obs), "Should pass without distribution");
}

#[test]
fn filter_insert_size_range() {
    let filter = FilterBuilder::new()
        .min_insert_size(300)
        .max_insert_size(600)
        .build();

    let obs_within = create_test_obs(450, true, true);
    assert!(filter.apply(&obs_within), "450 within [300, 600]");

    let obs_below = create_test_obs(250, true, true);
    assert!(!filter.apply(&obs_below), "250 < 300");

    let obs_above = create_test_obs(700, true, true);
    assert!(!filter.apply(&obs_above), "700 > 600");
}

#[test]
fn filter_insert_size_uses_absolute_value() {
    let filter = FilterBuilder::new()
        .min_insert_size(300)
        .max_insert_size(600)
        .build();

    // Negative insert size (mate upstream)
    let obs_negative = create_test_obs(-450, true, true);
    assert!(filter.apply(&obs_negative), "|-450| = 450 within [300, 600]");

    let obs_negative_above = create_test_obs(-700, true, true);
    assert!(!filter.apply(&obs_negative_above), "|-700| = 700 > 600");
}

#[test]
fn filter_requires_both_mates_mapped() {
    let filter = FilterBuilder::new()
        .require_both_mates_mapped(true)
        .build();

    let obs_both_mapped = create_test_obs(400, true, true);
    assert!(filter.apply(&obs_both_mapped), "Both mates mapped");

    let obs_mate_unmapped = create_test_obs(400, true, false);
    assert!(!filter.apply(&obs_mate_unmapped), "Mate unmapped");
}

#[test]
fn filter_skips_paired_checks_for_unpaired_reads() {
    let filter = FilterBuilder::new()
        .min_coverage(5) // Only non-paired filter
        .build();

    // Observation without mate_info (FASTQ input)
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

    // Should pass (no paired filters active)
    assert!(filter.apply(&obs), "Unpaired read should pass");
}

#[test]
fn filter_rejects_unpaired_when_proper_pairs_required() {
    let filter = FilterBuilder::new()
        .require_proper_pairs(1.0)
        .build();

    // Observation without mate_info
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

    // Should fail (proper pairs required but no mate_info)
    assert!(!filter.apply(&obs), "Should reject unpaired when proper_pairs required");
}

#[test]
fn filter_combined_proper_and_discordant() {
    let dist = InsertSizeDistribution {
        mean: 400,
        std_dev: 50,
        orientation: "FR".to_string(),
        sample_size: 1000,
    };

    let filter = FilterBuilder::new()
        .require_proper_pairs(1.0)
        .exclude_discordant_pairs(true)
        .with_insert_distribution(dist)
        .build();

    // Pass: proper + concordant
    let obs_pass = create_test_obs(450, true, true);
    assert!(filter.apply(&obs_pass));

    // Fail: improper (even if concordant insert size)
    let obs_improper = create_test_obs(450, false, true);
    assert!(!filter.apply(&obs_improper));

    // Fail: proper but discordant
    let obs_discordant = create_test_obs(800, true, true);
    assert!(!filter.apply(&obs_discordant));

    // Fail: both improper and discordant
    let obs_both_bad = create_test_obs(800, false, true);
    assert!(!filter.apply(&obs_both_bad));
}

#[test]
fn filter_custom_sigma_threshold() {
    let dist = InsertSizeDistribution {
        mean: 400,
        std_dev: 50,
        orientation: "FR".to_string(),
        sample_size: 1000,
    };

    // Default 3σ: 400 ± 150 = [250, 550]
    let filter_3sigma = FilterBuilder::new()
        .exclude_discordant_pairs(true)
        .discordant_sigma_threshold(3.0)
        .with_insert_distribution(dist.clone())
        .build();

    let obs_600 = create_test_obs(600, true, true);
    assert!(!filter_3sigma.apply(&obs_600), "600 > 550 at 3σ");

    // 5σ: 400 ± 250 = [150, 650]
    let filter_5sigma = FilterBuilder::new()
        .exclude_discordant_pairs(true)
        .discordant_sigma_threshold(5.0)
        .with_insert_distribution(dist)
        .build();

    assert!(filter_5sigma.apply(&obs_600), "600 < 650 at 5σ");
}

#[test]
fn filter_all_paired_constraints_together() {
    let dist = InsertSizeDistribution {
        mean: 400,
        std_dev: 50,
        orientation: "FR".to_string(),
        sample_size: 1000,
    };

    let filter = FilterBuilder::new()
        .require_proper_pairs(1.0)
        .exclude_discordant_pairs(true)
        .min_insert_size(300)
        .max_insert_size(600)
        .require_both_mates_mapped(true)
        .with_insert_distribution(dist)
        .build();

    // Perfect read: proper, concordant, within range, both mapped
    let obs_perfect = create_test_obs(450, true, true);
    assert!(filter.apply(&obs_perfect), "Perfect read should pass");

    // Fail proper pair
    let obs_fail_proper = create_test_obs(450, false, true);
    assert!(!filter.apply(&obs_fail_proper));

    // Fail discordant (beyond 3σ)
    let obs_fail_discordant = create_test_obs(700, true, true);
    assert!(!filter.apply(&obs_fail_discordant));

    // Fail min insert size
    let obs_fail_min = create_test_obs(250, true, true);
    assert!(!filter.apply(&obs_fail_min));

    // Fail mate unmapped
    let obs_fail_mate = create_test_obs(450, true, false);
    assert!(!filter.apply(&obs_fail_mate));
}

#[test]
fn filter_with_coverage_and_paired_constraints() {
    let dist = InsertSizeDistribution {
        mean: 400,
        std_dev: 50,
        orientation: "FR".to_string(),
        sample_size: 1000,
    };

    let filter = FilterBuilder::new()
        .min_coverage(5)
        .min_mapq(20)
        .require_proper_pairs(1.0)
        .exclude_discordant_pairs(true)
        .with_insert_distribution(dist)
        .build();

    // Pass all filters
    let obs_pass = create_test_obs(450, true, true);
    assert!(filter.apply(&obs_pass));

    // Fail coverage (local_coverage[0] < 5)
    let mate_info = MateInfo::new("mate/2".to_string(), true, 450, true, false, true);
    let mut alleles = HashMap::new();
    alleles.insert(b'A', 10);
    let obs_low_cov = VariantObservation::new(
        100, b'A', alleles, 0.95, "10M".to_string(), 60, 0,
        vec![3], // Low coverage
        35.0, "sample:read".to_string(),
    ).with_mate_info(mate_info);

    assert!(!filter.apply(&obs_low_cov), "Should fail low coverage");
}
