use phraya_io::plan::InsertSizeDistribution;

#[test]
fn from_bam_proper_pairs_calculates_mean_and_stddev() {
    // Need at least 100 samples
    let mut tlens = Vec::new();
    for i in 0..100 {
        tlens.push(400 + (i % 10) * 5); // Range 400-445
    }
    let dist = InsertSizeDistribution::from_bam_proper_pairs(&tlens).unwrap();

    // Mean should be ~420
    assert!(dist.mean >= 415 && dist.mean <= 425, "Mean was {}", dist.mean);

    // StdDev should be > 0
    assert!(dist.std_dev > 0, "StdDev was {}", dist.std_dev);

    assert_eq!(dist.sample_size, 100);
    assert_eq!(dist.orientation, "FR");
}

#[test]
fn from_bam_proper_pairs_requires_minimum_samples() {
    // Less than 100 samples
    let tlens = vec![400, 420, 380];
    assert!(InsertSizeDistribution::from_bam_proper_pairs(&tlens).is_none());
}

#[test]
fn from_bam_proper_pairs_accepts_exactly_100_samples() {
    let tlens: Vec<i32> = (0..100).map(|i| 400 + (i % 20) - 10).collect();
    let dist = InsertSizeDistribution::from_bam_proper_pairs(&tlens);

    assert!(dist.is_some());
    let dist = dist.unwrap();
    assert_eq!(dist.sample_size, 100);
}

#[test]
fn from_bam_proper_pairs_handles_large_dataset() {
    let tlens: Vec<i32> = (0..10000).map(|i| 450 + (i % 100) - 50).collect();
    let dist = InsertSizeDistribution::from_bam_proper_pairs(&tlens).unwrap();

    assert_eq!(dist.sample_size, 10000);
    assert!(dist.mean >= 440 && dist.mean <= 460, "Mean was {}", dist.mean);
}

#[test]
fn from_bam_proper_pairs_handles_uniform_insert_sizes() {
    // All identical insert sizes
    let tlens = vec![500; 200];
    let dist = InsertSizeDistribution::from_bam_proper_pairs(&tlens).unwrap();

    assert_eq!(dist.mean, 500);
    assert_eq!(dist.std_dev, 0); // No variance
    assert_eq!(dist.sample_size, 200);
}

#[test]
fn from_bam_proper_pairs_handles_negative_tlens() {
    // Negative TLEN = mate is upstream (reverse orientation)
    let mut tlens = Vec::new();
    for i in 0..100 {
        tlens.push(-(400 + (i % 10) * 5)); // Range -400 to -445
    }
    let dist = InsertSizeDistribution::from_bam_proper_pairs(&tlens);

    // Should still compute (using signed values)
    assert!(dist.is_some());
    let dist = dist.unwrap();
    assert!(dist.mean < 0, "Mean should be negative: {}", dist.mean);
}

#[test]
fn from_bam_proper_pairs_mixed_positive_negative() {
    // Mixed orientations (unusual but possible)
    let mut tlens = Vec::new();
    for i in 0..50 {
        tlens.push(400 + i * 2); // Positive
    }
    for i in 0..50 {
        tlens.push(-(400 + i * 2)); // Negative
    }

    let dist = InsertSizeDistribution::from_bam_proper_pairs(&tlens);
    assert!(dist.is_some());
    let dist = dist.unwrap();

    // Mean should be near 0 (balanced positive/negative)
    assert!(dist.mean.abs() < 50, "Mean should be near 0, was {}", dist.mean);
}
