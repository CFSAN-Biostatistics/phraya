use phraya_core::types::VariantObservation;
use phraya_io::queries::QueryIndex;

pub mod extractors;
pub mod tsv;
pub mod vcf;

pub use extractors::{extract_allele_frequency, extract_cigar_ops, extract_multi_map_fraction};

/// Named filter presets.
///
/// Presets return a pre-configured `FilterBuilder` as a starting point.
/// Any individual threshold set afterward overrides the preset value.
///
/// - **strict**: high-confidence calls only. Good for clinical/typing use.
/// - **tolerant**: catches low-frequency variants at the cost of more noise.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FilterPreset {
    Strict,
    Tolerant,
}

impl FilterPreset {
    /// Return a `FilterBuilder` seeded with this preset's defaults.
    pub fn builder(self) -> FilterBuilder {
        match self {
            FilterPreset::Strict => FilterBuilder::new()
                .min_coverage(10)
                .min_mapq(30)
                .min_allele_frequency(0.10)
                .exclude_tandem_repeats(true),
            FilterPreset::Tolerant => FilterBuilder::new()
                .min_coverage(3)
                .min_mapq(20)
                .min_allele_frequency(0.02),
        }
    }
}

/// Threshold-based filter for VariantObservations.
///
/// Constructed via the builder pattern: `ThresholdFilter::new().min_coverage(10).build()`.
/// `FilterBuilder` is a type alias for this type — callers using `FilterBuilder` continue to work.
#[derive(Debug, Clone)]
pub struct ThresholdFilter {
    min_coverage: Option<u32>,
    max_coverage: Option<u32>,
    min_mapq: Option<u8>,
    max_mapq: Option<u8>,
    min_base_quality: Option<f64>,
    min_allele_frequency: Option<f64>,
    min_kmer_uniqueness: Option<f64>,
    exclude_tandem_repeats: bool,
    exclude_discordant_pairs: bool,
    discordant_sigma_threshold: f64,
    require_proper_pairs: Option<f64>,
    min_insert_size: Option<i32>,
    max_insert_size: Option<i32>,
    require_both_mates_mapped: bool,
    insert_distribution: Option<phraya_io::plan::InsertSizeDistribution>,
}

/// Backward-compatible alias. Prefer `ThresholdFilter` directly.
pub type FilterBuilder = ThresholdFilter;

impl ThresholdFilter {
    pub fn new() -> Self {
        ThresholdFilter {
            min_coverage: None,
            max_coverage: None,
            min_mapq: None,
            max_mapq: None,
            min_base_quality: None,
            min_allele_frequency: None,
            min_kmer_uniqueness: None,
            exclude_tandem_repeats: false,
            exclude_discordant_pairs: false,
            discordant_sigma_threshold: 3.0,
            require_proper_pairs: None,
            min_insert_size: None,
            max_insert_size: None,
            require_both_mates_mapped: false,
            insert_distribution: None,
        }
    }

    /// No-op: returns `self`. Exists for backward compatibility with `.build()` call sites.
    pub fn build(self) -> ThresholdFilter {
        self
    }

    pub fn min_coverage(mut self, threshold: u32) -> Self {
        self.min_coverage = Some(threshold);
        self
    }

    pub fn max_coverage(mut self, threshold: u32) -> Self {
        self.max_coverage = Some(threshold);
        self
    }

    pub fn min_mapq(mut self, threshold: u8) -> Self {
        self.min_mapq = Some(threshold);
        self
    }

    pub fn max_mapq(mut self, threshold: u8) -> Self {
        self.max_mapq = Some(threshold);
        self
    }

    pub fn min_base_quality(mut self, threshold: f64) -> Self {
        self.min_base_quality = Some(threshold);
        self
    }

    pub fn min_allele_frequency(mut self, threshold: f64) -> Self {
        self.min_allele_frequency = Some(threshold);
        self
    }

    pub fn min_kmer_uniqueness(mut self, threshold: f64) -> Self {
        self.min_kmer_uniqueness = Some(threshold);
        self
    }

    pub fn exclude_tandem_repeats(mut self, value: bool) -> Self {
        self.exclude_tandem_repeats = value;
        self
    }

    pub fn exclude_discordant_pairs(mut self, value: bool) -> Self {
        self.exclude_discordant_pairs = value;
        self
    }

    pub fn discordant_sigma_threshold(mut self, sigma: f64) -> Self {
        self.discordant_sigma_threshold = sigma;
        self
    }

    pub fn require_proper_pairs(mut self, min_fraction: f64) -> Self {
        self.require_proper_pairs = Some(min_fraction);
        self
    }

    pub fn min_insert_size(mut self, size: i32) -> Self {
        self.min_insert_size = Some(size);
        self
    }

    pub fn max_insert_size(mut self, size: i32) -> Self {
        self.max_insert_size = Some(size);
        self
    }

    pub fn require_both_mates_mapped(mut self, value: bool) -> Self {
        self.require_both_mates_mapped = value;
        self
    }

    pub fn with_insert_distribution(mut self, dist: phraya_io::plan::InsertSizeDistribution) -> Self {
        self.insert_distribution = Some(dist);
        self
    }
}

impl Default for ThresholdFilter {
    fn default() -> Self {
        Self::new()
    }
}

impl ThresholdFilter {
    /// Apply filter to an observation
    pub fn apply(&self, obs: &VariantObservation) -> bool {
        if let Some(min) = self.min_coverage {
            match obs.coverage_at_variant() {
                None => return false,
                Some(cov) if cov < min => return false,
                _ => {}
            }
        }

        if let Some(max) = self.max_coverage {
            if let Some(cov) = obs.coverage_at_variant() {
                if cov > max {
                    return false;
                }
            }
        }

        // Check MAPQ
        if let Some(min) = self.min_mapq {
            if obs.mapq() < min {
                return false;
            }
        }

        if let Some(max) = self.max_mapq {
            if obs.mapq() > max {
                return false;
            }
        }

        // Check base quality
        if let Some(min) = self.min_base_quality {
            if obs.avg_base_quality() < min {
                return false;
            }
        }

        // Check allele frequency
        if let Some(min_freq) = self.min_allele_frequency {
            let total: u32 = obs.all_alleles().values().sum();
            if total == 0 {
                return false;
            }

            let _ref_count = obs.all_alleles().get(&obs.ref_base()).copied().unwrap_or(0);

            // Only fail if ref base doesn't meet frequency (check all alleles)
            let max_alt_freq = obs
                .all_alleles()
                .iter()
                .filter(|(&base, _)| base != obs.ref_base())
                .map(|(_, &count)| count as f64 / total as f64)
                .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
                .unwrap_or(0.0);

            if max_alt_freq < min_freq {
                return false;
            }
        }

        // Check k-mer uniqueness
        if let Some(min) = self.min_kmer_uniqueness {
            if obs.kmer_uniqueness() < min {
                return false;
            }
        }

        // Exclude tandem repeat variants
        if self.exclude_tandem_repeats && obs.in_tandem_repeat() {
            return false;
        }

        // Filter 1: Require minimum proper-pair fraction
        if let Some(min_fraction) = self.require_proper_pairs {
            match obs.proper_pair_fraction() {
                None => return false, // no paired reads → reject
                Some(frac) if frac < min_fraction => return false,
                _ => {}
            }
        }

        // Filter 2: Exclude discordant pairs — uses mean_insert_size (merge-stable) compared
        // against the plan's insert-size distribution at the configured sigma threshold.
        if self.exclude_discordant_pairs {
            if let (Some(mean_ins), Some(ref dist)) = (obs.mean_insert_size(), &self.insert_distribution) {
                let threshold = dist.mean as f64 + self.discordant_sigma_threshold * dist.std_dev as f64;
                if mean_ins > threshold {
                    return false;
                }
            }
            // No mean insert size or no distribution → cannot determine discordance → pass
        }

        // Filter 3: Insert size range — uses mean_insert_size (merge-stable).
        if self.min_insert_size.is_some() || self.max_insert_size.is_some() {
            if let Some(mean_ins) = obs.mean_insert_size() {
                let mean_ins_i32 = mean_ins.round() as i32;
                if let Some(min) = self.min_insert_size {
                    if mean_ins_i32 < min {
                        return false;
                    }
                }
                if let Some(max) = self.max_insert_size {
                    if mean_ins_i32 > max {
                        return false;
                    }
                }
            }
            // If no insert size data available, pass (insert size filter is not applicable)
        }

        // Filter 4: Both mates mapped — falls back to mate_info (pre-merge only).
        if let Some(mate_info) = obs.mate_info() {
            if self.require_both_mates_mapped && !mate_info.mate_mapped {
                return false;
            }
        } else if self.require_both_mates_mapped {
            return false;
        }

        true
    }

    /// Filter observations, returning an iterator
    pub fn filter<'a>(
        &'a self,
        observations: &'a [VariantObservation],
    ) -> impl Iterator<Item = &'a VariantObservation> {
        observations.iter().filter(move |obs| self.apply(obs))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;
    use phraya_io;

    fn create_observation(
        position: u32,
        mapq: u8,
        coverage: u32,
        base_quality: f64,
    ) -> VariantObservation {
        let mut alleles = HashMap::new();
        alleles.insert(b'A', coverage);

        VariantObservation::new(
            position,
            b'A',
            alleles,
            0.95,
            "10M".to_string(),
            mapq,
            0,
            vec![coverage],
            base_quality,
            "test:read".to_string(),
        )
    }

    #[test]
    fn min_coverage_filter() {
        let obs_low = create_observation(100, 60, 5, 35.0);
        let obs_high = create_observation(100, 60, 15, 35.0);

        let filter = FilterBuilder::new().min_coverage(10).build();

        assert!(!filter.apply(&obs_low));
        assert!(filter.apply(&obs_high));
    }

    #[test]
    fn min_coverage_rejects_observation_with_no_coverage_data() {
        let obs = VariantObservation::new(
            100,
            b'A',
            HashMap::from([(b'A', 10)]),
            0.95,
            "10M".to_string(),
            60,
            0,
            vec![], // no local_coverage window data
            35.0,
            "test:read".to_string(),
        );

        let filter = FilterBuilder::new().min_coverage(10).build();
        assert!(!filter.apply(&obs));
    }

    #[test]
    fn max_coverage_filter() {
        let obs_low = create_observation(100, 60, 5, 35.0);
        let obs_high = create_observation(100, 60, 25, 35.0);

        let filter = FilterBuilder::new().max_coverage(20).build();

        assert!(filter.apply(&obs_low));
        assert!(!filter.apply(&obs_high));
    }

    #[test]
    fn max_coverage_passes_observation_with_no_coverage_data() {
        // max_coverage only rejects when coverage data IS available and exceeds
        // the threshold; missing coverage data is not itself a rejection reason.
        let obs = VariantObservation::new(
            100,
            b'A',
            HashMap::from([(b'A', 10)]),
            0.95,
            "10M".to_string(),
            60,
            0,
            vec![], // no local_coverage window data
            35.0,
            "test:read".to_string(),
        );

        let filter = FilterBuilder::new().max_coverage(20).build();
        assert!(filter.apply(&obs));
    }

    #[test]
    fn min_mapq_filter() {
        let obs_low = create_observation(100, 30, 10, 35.0);
        let obs_high = create_observation(100, 50, 10, 35.0);

        let filter = FilterBuilder::new().min_mapq(40).build();

        assert!(!filter.apply(&obs_low));
        assert!(filter.apply(&obs_high));
    }

    #[test]
    fn max_mapq_filter() {
        let obs_low = create_observation(100, 30, 10, 35.0);
        let obs_high = create_observation(100, 50, 10, 35.0);

        let filter = FilterBuilder::new().max_mapq(40).build();

        assert!(filter.apply(&obs_low));
        assert!(!filter.apply(&obs_high));
    }

    #[test]
    fn min_base_quality_filter() {
        let obs_low = create_observation(100, 60, 10, 25.0);
        let obs_high = create_observation(100, 60, 10, 35.0);

        let filter = FilterBuilder::new().min_base_quality(30.0).build();

        assert!(!filter.apply(&obs_low));
        assert!(filter.apply(&obs_high));
    }

    #[test]
    fn boundary_values() {
        let obs_at_boundary = create_observation(100, 60, 10, 30.0);

        let filter = FilterBuilder::new().min_coverage(10).build();
        assert!(filter.apply(&obs_at_boundary)); // At boundary passes

        let filter = FilterBuilder::new().min_coverage(11).build();
        assert!(!filter.apply(&obs_at_boundary)); // Just above boundary fails
    }

    #[test]
    fn chaining_filters() {
        let obs = create_observation(100, 60, 25, 35.0);

        let filter1 = FilterBuilder::new().min_coverage(10).build();
        assert!(filter1.apply(&obs));

        let filter2 = FilterBuilder::new().min_coverage(20).build();
        assert!(filter2.apply(&obs));

        // Stricter filter
        let filter3 = FilterBuilder::new().min_coverage(30).build();
        assert!(!filter3.apply(&obs));

        // Composition: filter(min=10) then filter(min=20) should equal filter(min=20)
        let observations = vec![obs.clone()];
        let filtered1: Vec<_> = observations.iter().filter(|o| filter1.apply(o)).collect();
        let filtered2: Vec<_> = filtered1.iter().filter(|o| filter2.apply(o)).collect();

        let filtered_direct: Vec<_> = observations.iter().filter(|o| filter2.apply(o)).collect();

        assert_eq!(filtered2.len(), filtered_direct.len());
    }

    #[test]
    fn filter_composition_monotonic() {
        let observations = vec![
            create_observation(100, 60, 5, 35.0),
            create_observation(100, 60, 15, 35.0),
            create_observation(100, 60, 25, 35.0),
        ];

        let filter1 = FilterBuilder::new().min_coverage(10).build();
        let filter2 = FilterBuilder::new().min_coverage(20).build();

        let count1: usize = observations.iter().filter(|o| filter1.apply(o)).count();
        let count2: usize = observations.iter().filter(|o| filter2.apply(o)).count();

        // Stricter filter should have fewer results
        assert!(count2 <= count1);
    }

    #[test]
    fn no_filters_passes_all() {
        let obs = create_observation(100, 60, 10, 35.0);
        let filter = FilterBuilder::new().build();

        assert!(filter.apply(&obs));
    }

    #[test]
    fn threshold_filter_default_matches_new() {
        let obs = create_observation(100, 60, 10, 35.0);
        let filter = ThresholdFilter::default();

        assert!(filter.apply(&obs));
    }

    #[test]
    fn multiple_thresholds() {
        let obs = create_observation(100, 50, 15, 33.0);

        let filter = FilterBuilder::new()
            .min_coverage(10)
            .min_mapq(40)
            .min_base_quality(30.0)
            .build();

        assert!(filter.apply(&obs));

        // Fail one threshold
        let strict_filter = FilterBuilder::new()
            .min_coverage(10)
            .min_mapq(60) // This fails
            .min_base_quality(30.0)
            .build();

        assert!(!strict_filter.apply(&obs));
    }

    #[test]
    fn iterator_filter() {
        let observations = vec![
            create_observation(100, 60, 5, 35.0),
            create_observation(100, 60, 15, 35.0),
            create_observation(100, 60, 25, 35.0),
        ];

        let filter = FilterBuilder::new().min_coverage(10).build();
        let filtered: Vec<_> = filter.filter(&observations).collect();

        assert_eq!(filtered.len(), 2); // Only 15 and 25 pass
    }

    #[test]
    fn empty_observations() {
        let observations: Vec<VariantObservation> = vec![];
        let filter = FilterBuilder::new().min_coverage(10).build();
        let filtered: Vec<_> = filter.filter(&observations).collect();

        assert_eq!(filtered.len(), 0);
    }

    #[test]
    fn min_allele_frequency() {
        let mut alleles = HashMap::new();
        alleles.insert(b'A', 90); // 90% frequency
        alleles.insert(b'T', 10); // 10% frequency

        let obs = VariantObservation::new(
            100,
            b'A',
            alleles,
            0.95,
            "10M".to_string(),
            60,
            0,
            vec![100],
            35.0,
            "test:read".to_string(),
        );

        let filter = FilterBuilder::new().min_allele_frequency(0.05).build();
        assert!(filter.apply(&obs)); // 10% passes 5% threshold

        let strict_filter = FilterBuilder::new().min_allele_frequency(0.15).build();
        assert!(!strict_filter.apply(&obs)); // 10% fails 15% threshold
    }

    #[test]
    fn min_allele_frequency_rejects_zero_total_alleles() {
        let obs = VariantObservation::new(
            100,
            b'A',
            HashMap::new(), // no alleles at all → total count is 0
            0.95,
            "10M".to_string(),
            60,
            0,
            vec![100],
            35.0,
            "test:read".to_string(),
        );

        let filter = FilterBuilder::new().min_allele_frequency(0.05).build();
        assert!(!filter.apply(&obs));
    }

    #[test]
    fn require_both_mates_mapped_rejects_when_mate_unmapped() {
        use phraya_core::types::MateInfo;

        let obs = create_observation(100, 60, 10, 35.0).with_mate_info(MateInfo::new(
            "read/2".to_string(),
            true,
            300,
            true,
            false,
            false, // mate not mapped
        ));

        let filter = FilterBuilder::new().require_both_mates_mapped(true).build();
        assert!(!filter.apply(&obs));
    }

    #[test]
    fn require_both_mates_mapped_accepts_when_mate_mapped() {
        use phraya_core::types::MateInfo;

        let obs = create_observation(100, 60, 10, 35.0).with_mate_info(MateInfo::new(
            "read/2".to_string(),
            true,
            300,
            true,
            false,
            true, // mate mapped
        ));

        let filter = FilterBuilder::new().require_both_mates_mapped(true).build();
        assert!(filter.apply(&obs));
    }

    #[test]
    fn require_both_mates_mapped_rejects_when_no_mate_info() {
        let obs = create_observation(100, 60, 10, 35.0);

        let filter = FilterBuilder::new().require_both_mates_mapped(true).build();
        assert!(!filter.apply(&obs));
    }

    #[test]
    fn require_both_mates_mapped_false_ignores_missing_mate_info() {
        let obs = create_observation(100, 60, 10, 35.0);

        let filter = FilterBuilder::new().require_both_mates_mapped(false).build();
        assert!(filter.apply(&obs));
    }

    #[test]
    fn exclude_tandem_repeats_filter() {
        let mut alleles = HashMap::new();
        alleles.insert(b'T', 10u32);

        let obs_in_repeat = VariantObservation::new(
            100,
            b'A',
            alleles.clone(),
            0.95,
            "10M".to_string(),
            60,
            0,
            vec![10],
            35.0,
            "test:read".to_string(),
        )
        .with_tandem_repeat(true);

        let obs_not_in_repeat = VariantObservation::new(
            200,
            b'A',
            alleles,
            0.95,
            "10M".to_string(),
            60,
            0,
            vec![10],
            35.0,
            "test:read".to_string(),
        )
        .with_tandem_repeat(false);

        let filter = FilterBuilder::new().exclude_tandem_repeats(true).build();

        assert!(
            !filter.apply(&obs_in_repeat),
            "variant in tandem repeat must be excluded"
        );
        assert!(
            filter.apply(&obs_not_in_repeat),
            "variant outside tandem repeat must pass"
        );

        let permissive = FilterBuilder::new().exclude_tandem_repeats(false).build();
        assert!(
            permissive.apply(&obs_in_repeat),
            "when exclude_tandem_repeats=false, repeat variants must pass"
        );
    }

    #[test]
    fn strict_preset_rejects_low_quality() {
        // obs with coverage=5, mapq=25, allele_freq=0.05 — below all strict thresholds
        let mut alleles = HashMap::new();
        alleles.insert(b'T', 5u32); // 5/100 = 5% alt freq — below 10%
        alleles.insert(b'A', 95u32);
        let low_quality = VariantObservation::new(
            100, b'A', alleles, 0.5, "10M".to_string(), 25, 1, vec![100], 30.0, "s:r".to_string(),
        );

        let filter = FilterPreset::Strict.builder().build();
        assert!(!filter.apply(&low_quality), "strict must reject low-quality obs");
    }

    #[test]
    fn strict_preset_passes_high_quality() {
        let mut alleles = HashMap::new();
        alleles.insert(b'T', 20u32); // 20/100 = 20% alt freq — above 10%
        alleles.insert(b'A', 80u32);
        let high_quality = VariantObservation::new(
            100, b'A', alleles, 0.95, "10M".to_string(), 40, 0, vec![100], 35.0, "s:r".to_string(),
        );

        let filter = FilterPreset::Strict.builder().build();
        assert!(filter.apply(&high_quality), "strict must pass high-quality obs");
    }

    #[test]
    fn tolerant_preset_passes_low_coverage_variant() {
        // coverage=3, mapq=20, allele_freq=0.03 — just at tolerant thresholds
        let mut alleles = HashMap::new();
        alleles.insert(b'T', 3u32);
        alleles.insert(b'A', 97u32); // 3% alt freq
        let obs = VariantObservation::new(
            100, b'A', alleles, 0.9, "10M".to_string(), 20, 1, vec![3], 30.0, "s:r".to_string(),
        );

        let filter = FilterPreset::Tolerant.builder().build();
        assert!(filter.apply(&obs), "tolerant must pass low-coverage variants");
    }

    #[test]
    fn tolerant_preset_rejects_below_its_thresholds() {
        // mapq=19 — below tolerant's min_mapq=20
        let obs = create_observation(100, 19, 5, 30.0);
        let filter = FilterPreset::Tolerant.builder().build();
        assert!(!filter.apply(&obs), "tolerant must still reject mapq<20");
    }

    #[test]
    fn preset_can_be_overridden_by_explicit_threshold() {
        // Strict preset has min_coverage=10. Override to min_coverage=5.
        // Use an obs that passes all other strict thresholds (mapq=40, allele_freq=20%).
        let mut alleles = HashMap::new();
        alleles.insert(b'T', 2u32); // 2/10 = 20% alt freq → passes strict's 10% threshold
        alleles.insert(b'A', 8u32);
        // local_coverage[0]=7 → below strict's min_coverage=10, above override of 5
        let obs = VariantObservation::new(
            100, b'A', alleles, 0.95, "10M".to_string(), 40, 0, vec![7], 35.0, "s:r".to_string(),
        );

        let default_filter = FilterPreset::Strict.builder().build();
        assert!(!default_filter.apply(&obs), "local_coverage[0]=7 fails strict default (min=10)");

        let overridden = FilterPreset::Strict.builder().min_coverage(5).build();
        assert!(overridden.apply(&obs), "local_coverage[0]=7 passes after override to min_coverage=5");
    }

    #[test]
    fn strict_excludes_tandem_repeats_by_default() {
        let mut alleles = HashMap::new();
        alleles.insert(b'T', 20u32);
        alleles.insert(b'A', 80u32);
        let obs = VariantObservation::new(
            100, b'A', alleles, 0.95, "10M".to_string(), 40, 0, vec![100], 35.0, "s:r".to_string(),
        ).with_tandem_repeat(true);

        let filter = FilterPreset::Strict.builder().build();
        assert!(!filter.apply(&obs), "strict must exclude tandem repeat variants");
    }

    #[test]
    fn tolerant_does_not_exclude_tandem_repeats_by_default() {
        let mut alleles = HashMap::new();
        alleles.insert(b'T', 5u32);
        alleles.insert(b'A', 95u32);
        let obs = VariantObservation::new(
            100, b'A', alleles, 0.9, "10M".to_string(), 25, 0, vec![10], 30.0, "s:r".to_string(),
        ).with_tandem_repeat(true);

        let filter = FilterPreset::Tolerant.builder().build();
        assert!(filter.apply(&obs), "tolerant must not exclude tandem repeat variants");
    }

    /// Issue #181: strict preset exists and rejects low-quality observations
    #[test]
    fn issue_181_strict_preset_rejects_low_quality() {
        // obs with coverage=5, mapq=25, allele_freq=0.05 — below all strict thresholds
        let mut alleles = HashMap::new();
        alleles.insert(b'T', 5u32); // 5/100 = 5% alt freq — below 10%
        alleles.insert(b'A', 95u32);
        let low_quality = VariantObservation::new(
            100, b'A', alleles, 0.5, "10M".to_string(), 25, 1, vec![100], 30.0, "s:r".to_string(),
        );

        let filter = FilterPreset::Strict.builder().build();
        assert!(!filter.apply(&low_quality), "strict must reject low-quality obs");
    }

    /// Issue #181: strict preset exists and passes high-quality observations
    #[test]
    fn issue_181_strict_preset_passes_high_quality() {
        let mut alleles = HashMap::new();
        alleles.insert(b'T', 20u32); // 20/100 = 20% alt freq — above 10%
        alleles.insert(b'A', 80u32);
        let high_quality = VariantObservation::new(
            100, b'A', alleles, 0.95, "10M".to_string(), 40, 0, vec![100], 35.0, "s:r".to_string(),
        );

        let filter = FilterPreset::Strict.builder().build();
        assert!(filter.apply(&high_quality), "strict must pass high-quality obs");
    }

    /// Issue #181: tolerant preset exists and passes low-coverage variants
    #[test]
    fn issue_181_tolerant_preset_passes_low_coverage_variant() {
        // coverage=3, mapq=20, allele_freq=0.03 — just at tolerant thresholds
        let mut alleles = HashMap::new();
        alleles.insert(b'T', 3u32);
        alleles.insert(b'A', 97u32); // 3% alt freq
        let obs = VariantObservation::new(
            100, b'A', alleles, 0.9, "10M".to_string(), 20, 1, vec![3], 30.0, "s:r".to_string(),
        );

        let filter = FilterPreset::Tolerant.builder().build();
        assert!(filter.apply(&obs), "tolerant must pass low-coverage variants");
    }

    /// Issue #181: tolerant preset rejects below its thresholds
    #[test]
    fn issue_181_tolerant_preset_rejects_below_its_thresholds() {
        // mapq=19 — below tolerant's min_mapq=20
        let obs = create_observation(100, 19, 5, 30.0);
        let filter = FilterPreset::Tolerant.builder().build();
        assert!(!filter.apply(&obs), "tolerant must still reject mapq<20");
    }

    /// Issue #181: strict preset excludes tandem repeats by default
    #[test]
    fn issue_181_strict_excludes_tandem_repeats_by_default() {
        let mut alleles = HashMap::new();
        alleles.insert(b'T', 20u32);
        alleles.insert(b'A', 80u32);
        let obs = VariantObservation::new(
            100, b'A', alleles, 0.95, "10M".to_string(), 40, 0, vec![100], 35.0, "s:r".to_string(),
        ).with_tandem_repeat(true);

        let filter = FilterPreset::Strict.builder().build();
        assert!(!filter.apply(&obs), "strict must exclude tandem repeat variants");
    }

    /// Issue #181: tolerant does not exclude tandem repeats by default
    #[test]
    fn issue_181_tolerant_does_not_exclude_tandem_repeats_by_default() {
        let mut alleles = HashMap::new();
        alleles.insert(b'T', 5u32);
        alleles.insert(b'A', 95u32);
        let obs = VariantObservation::new(
            100, b'A', alleles, 0.9, "10M".to_string(), 25, 0, vec![10], 30.0, "s:r".to_string(),
        ).with_tandem_repeat(true);

        let filter = FilterPreset::Tolerant.builder().build();
        assert!(filter.apply(&obs), "tolerant must not exclude tandem repeat variants");
    }

    /// Issue #181: strict preset threshold values match old conservative
    #[test]
    fn issue_181_strict_has_correct_thresholds() {
        // Verify strict preset has: min_coverage=10, min_mapq=30, min_allele_frequency=0.10, exclude_tandem_repeats=true
        let mut alleles = HashMap::new();
        alleles.insert(b'T', 10u32); // 10/100 = 10% alt freq — exactly at threshold
        alleles.insert(b'A', 90u32);
        let obs_at_boundary = VariantObservation::new(
            100, b'A', alleles, 0.95, "10M".to_string(), 30, 0, vec![10], 35.0, "s:r".to_string(),
        );

        let filter = FilterPreset::Strict.builder().build();
        assert!(filter.apply(&obs_at_boundary), "strict must pass obs at coverage=10, mapq=30, allele_freq=10%");

        // Just below coverage
        let mut alleles2 = HashMap::new();
        alleles2.insert(b'T', 10u32);
        alleles2.insert(b'A', 90u32);
        let obs_below_cov = VariantObservation::new(
            100, b'A', alleles2, 0.95, "10M".to_string(), 30, 0, vec![9], 35.0, "s:r".to_string(),
        );
        assert!(!filter.apply(&obs_below_cov), "strict must reject coverage<10");
    }

    /// Issue #181: tolerant preset threshold values match old sensitive
    #[test]
    fn issue_181_tolerant_has_correct_thresholds() {
        // Verify tolerant preset has: min_coverage=3, min_mapq=20, min_allele_frequency=0.02, no exclude_tandem_repeats
        let mut alleles = HashMap::new();
        alleles.insert(b'T', 2u32); // 2/100 = 2% alt freq — exactly at threshold
        alleles.insert(b'A', 98u32);
        let obs_at_boundary = VariantObservation::new(
            100, b'A', alleles, 0.9, "10M".to_string(), 20, 0, vec![3], 30.0, "s:r".to_string(),
        );

        let filter = FilterPreset::Tolerant.builder().build();
        assert!(filter.apply(&obs_at_boundary), "tolerant must pass obs at coverage=3, mapq=20, allele_freq=2%");

        // Just below allele frequency
        let mut alleles2 = HashMap::new();
        alleles2.insert(b'T', 1u32); // 1/100 = 1% alt freq — below 2%
        alleles2.insert(b'A', 99u32);
        let obs_below_freq = VariantObservation::new(
            100, b'A', alleles2, 0.9, "10M".to_string(), 20, 0, vec![3], 30.0, "s:r".to_string(),
        );
        assert!(!filter.apply(&obs_below_freq), "tolerant must reject allele_freq<2%");
    }

    /// Issue #181: FilterPreset::Strict and FilterPreset::Tolerant variants exist
    #[test]
    fn issue_181_enum_variants_exist() {
        // This test verifies that the enum variants can be constructed
        let _strict = FilterPreset::Strict;
        let _tolerant = FilterPreset::Tolerant;

        // Test that they can be used in pattern matching
        let presets = vec![FilterPreset::Strict, FilterPreset::Tolerant];
        for preset in presets {
            let _builder = preset.builder();
        }
    }

    fn obs_with_insert_stats(insert_size_sum: i64, insert_size_count: u32, total_paired: u32) -> VariantObservation {
        create_observation(100, 60, 10, 35.0)
            .with_pair_counts(total_paired, total_paired) // all properly paired for simplicity
            .with_insert_stats(insert_size_sum, insert_size_count)
    }

    fn test_dist() -> phraya_io::plan::InsertSizeDistribution {
        phraya_io::plan::InsertSizeDistribution {
            mean: 400,
            std_dev: 50,
            orientation: "FR".to_string(),
            sample_size: 1000,
        }
    }

    #[test]
    fn min_insert_size_filter_passes_when_mean_above_threshold() {
        // mean = 400 / 1 = 400 → passes min=200
        let obs = obs_with_insert_stats(400, 1, 1);
        let filter = FilterBuilder::new().min_insert_size(200).build();
        assert!(filter.apply(&obs));
    }

    #[test]
    fn min_insert_size_filter_rejects_when_mean_below_threshold() {
        // mean = 100 / 1 = 100 → fails min=200
        let obs = obs_with_insert_stats(100, 1, 1);
        let filter = FilterBuilder::new().min_insert_size(200).build();
        assert!(!filter.apply(&obs));
    }

    #[test]
    fn max_insert_size_filter_passes_when_mean_below_threshold() {
        // mean = 300 / 1 = 300 → passes max=500
        let obs = obs_with_insert_stats(300, 1, 1);
        let filter = FilterBuilder::new().max_insert_size(500).build();
        assert!(filter.apply(&obs));
    }

    #[test]
    fn max_insert_size_filter_rejects_when_mean_above_threshold() {
        // mean = 800 / 1 = 800 → fails max=500
        let obs = obs_with_insert_stats(800, 1, 1);
        let filter = FilterBuilder::new().max_insert_size(500).build();
        assert!(!filter.apply(&obs));
    }

    #[test]
    fn insert_size_filter_passes_when_no_paired_reads() {
        // No insert data (count=0) → filter is not applicable, pass through
        let obs = obs_with_insert_stats(0, 0, 0);
        let filter = FilterBuilder::new().min_insert_size(200).max_insert_size(500).build();
        assert!(filter.apply(&obs), "no paired reads → insert size filter should not apply");
    }

    #[test]
    fn insert_size_filter_works_after_merge_with_aggregate_stats() {
        // Two reads merged: inserts 300+500=800, count=2 → mean=400 → passes [300,500]
        let obs = obs_with_insert_stats(800, 2, 2);
        let filter = FilterBuilder::new().min_insert_size(300).max_insert_size(500).build();
        assert!(filter.apply(&obs), "mean=400 passes [300,500]");

        // Two reads: 50+150=200, count=2 → mean=100 → fails min=300
        let obs_fail = obs_with_insert_stats(200, 2, 2);
        assert!(!filter.apply(&obs_fail), "mean=100 fails min=300");
    }

    #[test]
    fn exclude_discordant_pairs_passes_when_mean_within_distribution() {
        // mean=450, dist=400±50×3=550 → 450 < 550 → concordant → pass
        let obs = obs_with_insert_stats(450, 1, 1);
        let filter = FilterBuilder::new()
            .exclude_discordant_pairs(true)
            .with_insert_distribution(test_dist())
            .build();
        assert!(filter.apply(&obs), "mean=450 is within 3σ of 400");
    }

    #[test]
    fn exclude_discordant_pairs_rejects_when_mean_beyond_distribution() {
        // mean=800, dist=400±50×3=550 → 800 > 550 → discordant → reject
        let obs = obs_with_insert_stats(800, 1, 1);
        let filter = FilterBuilder::new()
            .exclude_discordant_pairs(true)
            .with_insert_distribution(test_dist())
            .build();
        assert!(!filter.apply(&obs), "mean=800 exceeds 3σ threshold of 550");
    }

    #[test]
    fn exclude_discordant_pairs_passes_when_no_paired_reads() {
        // No insert data → cannot determine discordance → pass
        let obs = obs_with_insert_stats(0, 0, 0);
        let filter = FilterBuilder::new()
            .exclude_discordant_pairs(true)
            .with_insert_distribution(test_dist())
            .build();
        assert!(filter.apply(&obs), "no paired reads → discordant filter should not apply");
    }

    #[test]
    fn exclude_discordant_pairs_passes_without_distribution() {
        // No distribution provided → cannot determine discordance → pass
        let obs = obs_with_insert_stats(800, 1, 1);
        let filter = FilterBuilder::new().exclude_discordant_pairs(true).build();
        assert!(filter.apply(&obs), "no distribution → discordant filter has no effect");
    }

    #[test]
    fn discordant_filter_respects_sigma_threshold() {
        // dist=400±50, obs mean=600; at 3σ threshold=550 → 600>550 → discordant
        let obs = obs_with_insert_stats(600, 1, 1);
        let filter_3s = FilterBuilder::new()
            .exclude_discordant_pairs(true)
            .discordant_sigma_threshold(3.0)
            .with_insert_distribution(test_dist())
            .build();
        assert!(!filter_3s.apply(&obs), "mean=600 > 550 at 3σ");

        // at 5σ threshold=650 → 600<650 → concordant
        let filter_5s = FilterBuilder::new()
            .exclude_discordant_pairs(true)
            .discordant_sigma_threshold(5.0)
            .with_insert_distribution(test_dist())
            .build();
        assert!(filter_5s.apply(&obs), "mean=600 < 650 at 5σ");
    }

    #[test]
    fn discordant_filter_works_after_merge_with_aggregate_mean() {
        // 10 reads merged: 8 concordant (insert~400) + 2 extreme (insert~1000)
        // sum = 8*400 + 2*1000 = 3200+2000 = 5200, count=10 → mean=520 > 550? No, 520<550 → pass
        let obs_pass = obs_with_insert_stats(5200, 10, 10);
        let filter = FilterBuilder::new()
            .exclude_discordant_pairs(true)
            .with_insert_distribution(test_dist())
            .build();
        assert!(filter.apply(&obs_pass), "mean=520 is within 3σ=550 post-merge");

        // 2 extreme reads: sum=2000, count=2 → mean=1000 > 550 → reject
        let obs_fail = obs_with_insert_stats(2000, 2, 2);
        assert!(!filter.apply(&obs_fail), "mean=1000 exceeds 3σ=550 post-merge");
    }
}

/// Internal AST representation for expression filters.
#[derive(Debug, Clone)]
enum Expr {
    Or(Box<Expr>, Box<Expr>),
    And(Box<Expr>, Box<Expr>),
    Not(Box<Expr>),
    Comparison {
        field: String,
        op: CompOp,
        value: f64,
    },
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum CompOp {
    GreaterEq,
    Greater,
    LessEq,
    Less,
    Equal,
    NotEqual,
}

/// Expression-based filter for VariantObservations.
///
/// Parses and evaluates boolean expressions over VariantObservation fields.
/// Example: `"coverage >= 10 && mapq > 30"`
///
/// Supported fields: coverage, mapq, allele_frequency, base_quality, confidence,
/// kmer_uniqueness, edit_distance, in_tandem_repeat
///
/// Supported operators: >=, >, <=, <, ==, !=, &&, ||, !, parentheses
#[derive(Debug, Clone)]
pub struct ExprFilter {
    ast: Expr,
}

impl ExprFilter {
    /// Create a new expression filter from a string expression.
    /// Returns an error if the expression is malformed or references unknown fields.
    pub fn new(expr: &str) -> Result<Self, ExprParseError> {
        let parser = Parser::new(expr);
        let ast = parser.parse()?;
        Ok(ExprFilter { ast })
    }

    /// Apply the filter to an observation.
    pub fn apply(&self, obs: &VariantObservation) -> bool {
        evaluate_expr(&self.ast, obs)
    }

    /// Filter observations, returning an iterator
    pub fn filter<'a>(
        &'a self,
        observations: &'a [VariantObservation],
    ) -> impl Iterator<Item = &'a VariantObservation> {
        observations.iter().filter(move |obs| self.apply(obs))
    }
}

/// Parse error for expression parsing
#[derive(Debug, Clone, PartialEq)]
pub enum ExprParseError {
    UnknownField(String),
    MalformedExpression(String),
    UnmatchedParen,
    MissingOperand,
}

impl std::fmt::Display for ExprParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ExprParseError::UnknownField(name) => {
                write!(f, "unknown field: '{}'", name)
            }
            ExprParseError::MalformedExpression(msg) => {
                write!(f, "malformed expression: {}", msg)
            }
            ExprParseError::UnmatchedParen => {
                write!(f, "unmatched parenthesis")
            }
            ExprParseError::MissingOperand => {
                write!(f, "missing operand")
            }
        }
    }
}

impl std::error::Error for ExprParseError {}

/// Extract a field value from a VariantObservation
fn extract_field(field: &str, obs: &VariantObservation) -> Result<f64, ExprParseError> {
    match field {
        "coverage" => Ok(obs.coverage_at_variant().unwrap_or(0) as f64),
        "mapq" => Ok(obs.mapq() as f64),
        "allele_frequency" => {
            let total: u32 = obs.all_alleles().values().sum();
            if total == 0 {
                Ok(0.0)
            } else {
                let max_alt_freq = obs
                    .all_alleles()
                    .iter()
                    .filter(|(&base, _)| base != obs.ref_base())
                    .map(|(_, &count)| count as f64 / total as f64)
                    .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
                    .unwrap_or(0.0);
                Ok(max_alt_freq)
            }
        }
        "base_quality" => Ok(obs.avg_base_quality()),
        "confidence" => Ok(obs.confidence()),
        "kmer_uniqueness" => Ok(obs.kmer_uniqueness()),
        "edit_distance" => Ok(obs.edit_distance() as f64),
        "in_tandem_repeat" => Ok(if obs.in_tandem_repeat() { 1.0 } else { 0.0 }),
        _ => Err(ExprParseError::UnknownField(field.to_string())),
    }
}

/// Evaluate an expression AST against a VariantObservation
fn evaluate_expr(expr: &Expr, obs: &VariantObservation) -> bool {
    match expr {
        Expr::Or(left, right) => evaluate_expr(left, obs) || evaluate_expr(right, obs),
        Expr::And(left, right) => evaluate_expr(left, obs) && evaluate_expr(right, obs),
        Expr::Not(inner) => !evaluate_expr(inner, obs),
        Expr::Comparison { field, op, value } => {
            match extract_field(field, obs) {
                Ok(field_value) => compare_values(field_value, *op, *value),
                Err(_) => false, // Shouldn't happen if parser is correct
            }
        }
    }
}

/// Compare two values using an operator
fn compare_values(left: f64, op: CompOp, right: f64) -> bool {
    match op {
        CompOp::GreaterEq => (left - right).abs() < f64::EPSILON || left > right,
        CompOp::Greater => left > right,
        CompOp::LessEq => (left - right).abs() < f64::EPSILON || left < right,
        CompOp::Less => left < right,
        CompOp::Equal => (left - right).abs() < f64::EPSILON,
        CompOp::NotEqual => (left - right).abs() >= f64::EPSILON,
    }
}

/// Recursive-descent parser for filter expressions
struct Parser {
    input: Vec<char>,
    pos: usize,
}

impl Parser {
    fn new(input: &str) -> Self {
        Parser {
            input: input.chars().collect(),
            pos: 0,
        }
    }

    fn parse(mut self) -> Result<Expr, ExprParseError> {
        let expr = self.parse_or()?;
        self.skip_whitespace();
        if self.pos < self.input.len() {
            return Err(ExprParseError::MalformedExpression(
                "unexpected characters after expression".to_string(),
            ));
        }
        Ok(expr)
    }

    fn parse_or(&mut self) -> Result<Expr, ExprParseError> {
        let mut left = self.parse_and()?;
        self.skip_whitespace();
        while self.match_token("||") {
            self.skip_whitespace();
            let right = self.parse_and()?;
            left = Expr::Or(Box::new(left), Box::new(right));
            self.skip_whitespace();
        }
        Ok(left)
    }

    fn parse_and(&mut self) -> Result<Expr, ExprParseError> {
        let mut left = self.parse_not()?;
        self.skip_whitespace();
        while self.match_token("&&") {
            self.skip_whitespace();
            let right = self.parse_not()?;
            left = Expr::And(Box::new(left), Box::new(right));
            self.skip_whitespace();
        }
        Ok(left)
    }

    fn parse_not(&mut self) -> Result<Expr, ExprParseError> {
        self.skip_whitespace();
        if self.match_token("!") {
            self.skip_whitespace();
            let expr = self.parse_not()?;
            Ok(Expr::Not(Box::new(expr)))
        } else {
            self.parse_atom()
        }
    }

    fn parse_atom(&mut self) -> Result<Expr, ExprParseError> {
        self.skip_whitespace();
        if self.match_token("(") {
            self.skip_whitespace();
            let expr = self.parse_or()?;
            self.skip_whitespace();
            if !self.match_token(")") {
                return Err(ExprParseError::UnmatchedParen);
            }
            Ok(expr)
        } else {
            self.parse_comparison()
        }
    }

    fn parse_comparison(&mut self) -> Result<Expr, ExprParseError> {
        self.skip_whitespace();
        let field = self.parse_field()?;
        self.skip_whitespace();
        let op = self.parse_op()?;
        self.skip_whitespace();
        let value = self.parse_number()?;
        Ok(Expr::Comparison { field, op, value })
    }

    fn parse_field(&mut self) -> Result<String, ExprParseError> {
        self.skip_whitespace();
        let start = self.pos;
        while self.pos < self.input.len() && (self.input[self.pos].is_alphanumeric() || self.input[self.pos] == '_') {
            self.pos += 1;
        }
        if start == self.pos {
            return Err(ExprParseError::MissingOperand);
        }
        let field: String = self.input[start..self.pos].iter().collect();

        // Validate field name
        match field.as_str() {
            "coverage" | "mapq" | "allele_frequency" | "base_quality" | "confidence" |
            "kmer_uniqueness" | "edit_distance" | "in_tandem_repeat" => Ok(field),
            _ => Err(ExprParseError::UnknownField(field)),
        }
    }

    fn parse_op(&mut self) -> Result<CompOp, ExprParseError> {
        self.skip_whitespace();
        if self.match_token(">=") {
            Ok(CompOp::GreaterEq)
        } else if self.match_token("<=") {
            Ok(CompOp::LessEq)
        } else if self.match_token("==") {
            Ok(CompOp::Equal)
        } else if self.match_token("!=") {
            Ok(CompOp::NotEqual)
        } else if self.match_token(">") {
            Ok(CompOp::Greater)
        } else if self.match_token("<") {
            Ok(CompOp::Less)
        } else {
            Err(ExprParseError::MalformedExpression(
                "expected comparison operator".to_string(),
            ))
        }
    }

    fn parse_number(&mut self) -> Result<f64, ExprParseError> {
        self.skip_whitespace();
        let start = self.pos;

        // Optional sign
        if self.pos < self.input.len() && (self.input[self.pos] == '+' || self.input[self.pos] == '-') {
            self.pos += 1;
        }

        // Integer part
        if self.pos >= self.input.len() || !self.input[self.pos].is_ascii_digit() {
            return Err(ExprParseError::MissingOperand);
        }

        while self.pos < self.input.len() && self.input[self.pos].is_ascii_digit() {
            self.pos += 1;
        }

        // Decimal part
        if self.pos < self.input.len() && self.input[self.pos] == '.' {
            self.pos += 1;
            while self.pos < self.input.len() && self.input[self.pos].is_ascii_digit() {
                self.pos += 1;
            }
        }

        let num_str: String = self.input[start..self.pos].iter().collect();
        num_str.parse::<f64>()
            .map_err(|_| ExprParseError::MalformedExpression("invalid number".to_string()))
    }

    fn match_token(&mut self, token: &str) -> bool {
        let remaining = &self.input[self.pos..];
        if remaining.len() < token.len() {
            return false;
        }
        let check: String = remaining.iter().take(token.len()).collect();
        if check == token {
            self.pos += token.len();
            true
        } else {
            false
        }
    }

    fn skip_whitespace(&mut self) {
        while self.pos < self.input.len() && self.input[self.pos].is_whitespace() {
            self.pos += 1;
        }
    }
}

#[cfg(test)]
mod expr_filter_tests {
    use super::*;
    use std::collections::HashMap;

    fn create_obs(
        position: u32,
        mapq: u8,
        coverage: u32,
        base_quality: f64,
        confidence: f64,
        edit_distance: u32,
        in_tandem_repeat: bool,
    ) -> VariantObservation {
        // coverage = alt allele count; total depth = 100, allele_frequency = coverage/100
        let mut alleles = HashMap::new();
        alleles.insert(b'A', 100u32.saturating_sub(coverage)); // ref
        alleles.insert(b'T', coverage); // alt

        let mut obs = VariantObservation::new(
            position,
            b'A',
            alleles,
            confidence,
            "10M".to_string(),
            mapq,
            edit_distance,
            vec![coverage],
            base_quality,
            "test:read".to_string(),
        );

        if in_tandem_repeat {
            obs = obs.with_tandem_repeat(true);
        }

        obs
    }

    /// Issue #150: Single coverage comparison
    #[test]
    fn issue_150_expr_single_coverage_gte() {
        let filter = ExprFilter::new("coverage >= 10").expect("valid expr");
        let obs_pass = create_obs(100, 60, 15, 35.0, 0.95, 0, false);
        let obs_fail = create_obs(100, 60, 5, 35.0, 0.95, 0, false);

        assert!(
            filter.apply(&obs_pass),
            "coverage >= 10 should pass obs with coverage=15"
        );
        assert!(
            !filter.apply(&obs_fail),
            "coverage >= 10 should fail obs with coverage=5"
        );
    }

    /// Issue #150: AND operator combines two conditions
    #[test]
    fn issue_150_expr_and_operator() {
        let filter = ExprFilter::new("mapq > 30 && allele_frequency >= 0.1")
            .expect("valid expr");

        let obs_both_pass = create_obs(100, 40, 20, 35.0, 0.95, 0, false);
        let obs_mapq_fails = create_obs(100, 20, 20, 35.0, 0.95, 0, false);
        let obs_freq_fails = create_obs(100, 40, 5, 35.0, 0.95, 0, false);

        assert!(
            filter.apply(&obs_both_pass),
            "mapq > 30 && allele_frequency >= 0.1 should pass when both conditions true"
        );
        assert!(
            !filter.apply(&obs_mapq_fails),
            "mapq > 30 && allele_frequency >= 0.1 should fail when mapq condition fails"
        );
        assert!(
            !filter.apply(&obs_freq_fails),
            "mapq > 30 && allele_frequency >= 0.1 should fail when allele_frequency condition fails"
        );
    }

    /// Issue #150: OR operator
    #[test]
    fn issue_150_expr_or_operator() {
        let filter = ExprFilter::new("coverage >= 10 || mapq > 40").expect("valid expr");

        let obs_coverage_pass = create_obs(100, 20, 15, 35.0, 0.95, 0, false);
        let obs_mapq_pass = create_obs(100, 50, 5, 35.0, 0.95, 0, false);
        let obs_both_fail = create_obs(100, 20, 5, 35.0, 0.95, 0, false);

        assert!(
            filter.apply(&obs_coverage_pass),
            "coverage >= 10 || mapq > 40 should pass when coverage condition true"
        );
        assert!(
            filter.apply(&obs_mapq_pass),
            "coverage >= 10 || mapq > 40 should pass when mapq condition true"
        );
        assert!(
            !filter.apply(&obs_both_fail),
            "coverage >= 10 || mapq > 40 should fail when both conditions false"
        );
    }

    /// Issue #150: Unknown field name produces parse error with field name in message
    #[test]
    fn issue_150_expr_unknown_field() {
        let result = ExprFilter::new("unknown_field >= 10");

        assert!(
            result.is_err(),
            "expression with unknown field should produce error"
        );
        let err = result.unwrap_err();
        assert!(
            err.to_string().contains("unknown_field"),
            "error message should mention the unknown field name"
        );
    }

    /// Issue #150: Malformed expression — unclosed paren
    #[test]
    fn issue_150_expr_unclosed_paren() {
        let result = ExprFilter::new("(coverage >= 10");

        assert!(
            result.is_err(),
            "expression with unclosed paren should produce error"
        );
        assert!(result.unwrap_err().to_string().contains("parenthesis"));
    }

    #[test]
    fn expr_parse_error_display_variants() {
        assert!(ExprParseError::UnmatchedParen.to_string().contains("parenthesis"));
        assert!(ExprParseError::MissingOperand.to_string().contains("missing operand"));
        assert!(ExprParseError::MalformedExpression("bad token".to_string())
            .to_string()
            .contains("bad token"));
        assert!(ExprParseError::UnknownField("foo".to_string())
            .to_string()
            .contains("foo"));
    }

    /// Issue #150: Malformed expression — missing operand
    #[test]
    fn issue_150_expr_missing_operand() {
        let result = ExprFilter::new("coverage >= ");

        assert!(
            result.is_err(),
            "expression with missing operand should produce error"
        );
    }

    #[test]
    fn expr_trailing_garbage_after_expression_is_malformed() {
        let result = ExprFilter::new("coverage >= 10 extra");
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("unexpected characters"));
    }

    #[test]
    fn expr_missing_field_name_is_missing_operand() {
        let result = ExprFilter::new(">= 10");
        assert!(matches!(result, Err(ExprParseError::MissingOperand)));
    }

    #[test]
    fn expr_unrecognized_operator_is_malformed() {
        let result = ExprFilter::new("coverage ~ 10");
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("comparison operator"));
    }

    #[test]
    fn expr_negative_number_literal_parses() {
        let filter = ExprFilter::new("edit_distance >= -5").expect("valid expr");
        let obs = create_obs(100, 40, 15, 35.0, 0.95, 2, false);
        assert!(filter.apply(&obs), "edit_distance=2 should satisfy >= -5");
    }

    #[test]
    fn expr_allele_frequency_zero_when_no_alleles() {
        let obs = VariantObservation::new(
            100,
            b'A',
            HashMap::new(),
            0.95,
            "10M".to_string(),
            40,
            0,
            vec![10],
            35.0,
            "test:read".to_string(),
        );
        let filter = ExprFilter::new("allele_frequency == 0").expect("valid expr");
        assert!(filter.apply(&obs));
    }

    #[test]
    fn extract_field_unknown_field_errs_directly() {
        // extract_field's own field validation is a defensive fallback behind the
        // parser's identical validation in parse_field; exercise it directly since
        // ExprFilter::new can never reach it through normal parsing.
        let obs = create_obs(100, 40, 15, 35.0, 0.95, 0, false);
        let result = super::extract_field("not_a_real_field", &obs);
        assert!(matches!(result, Err(ExprParseError::UnknownField(_))));
    }

    #[test]
    fn evaluate_expr_false_when_field_extraction_fails() {
        // Mirrors extract_field's own defensive fallback: evaluate_expr treats an
        // extraction error as a non-match rather than panicking.
        let obs = create_obs(100, 40, 15, 35.0, 0.95, 0, false);
        let expr = super::Expr::Comparison {
            field: "not_a_real_field".to_string(),
            op: super::CompOp::GreaterEq,
            value: 0.0,
        };
        assert!(!super::evaluate_expr(&expr, &obs));
    }

    /// Issue #150: Both expr and threshold filters can be composed
    #[test]
    fn issue_150_expr_and_threshold_together() {
        let expr_filter = ExprFilter::new("coverage >= 10").expect("valid expr");
        let threshold_filter = FilterBuilder::new().min_mapq(30).build();

        let obs = create_obs(100, 40, 15, 35.0, 0.95, 0, false);

        // Both filters should pass
        assert!(
            expr_filter.apply(&obs),
            "expr filter should pass obs with coverage=15"
        );
        assert!(
            threshold_filter.apply(&obs),
            "threshold filter should pass obs with mapq=40"
        );
    }

    /// Issue #150: Parentheses work correctly
    #[test]
    fn issue_150_expr_parentheses() {
        let filter = ExprFilter::new("(coverage >= 10 && mapq > 30) || mapq > 50")
            .expect("valid expr");

        let obs1 = create_obs(100, 20, 15, 35.0, 0.95, 0, false); // coverage ok, mapq not
        let obs2 = create_obs(100, 60, 5, 35.0, 0.95, 0, false); // high mapq saves it

        assert!(
            !filter.apply(&obs1),
            "coverage ok but mapq not ok, and overall mapq not > 50"
        );
        assert!(
            filter.apply(&obs2),
            "mapq > 50 should pass despite low coverage"
        );
    }

    /// Issue #150: base_quality field comparison
    #[test]
    fn issue_150_expr_base_quality() {
        let filter = ExprFilter::new("base_quality >= 30.0").expect("valid expr");

        let obs_pass = create_obs(100, 60, 10, 35.0, 0.95, 0, false);
        let obs_fail = create_obs(100, 60, 10, 25.0, 0.95, 0, false);

        assert!(filter.apply(&obs_pass), "base_quality >= 30.0 should pass");
        assert!(
            !filter.apply(&obs_fail),
            "base_quality >= 30.0 should fail"
        );
    }

    /// Issue #150: confidence field comparison
    #[test]
    fn issue_150_expr_confidence() {
        let filter = ExprFilter::new("confidence >= 0.9").expect("valid expr");

        let obs_pass = create_obs(100, 60, 10, 35.0, 0.95, 0, false);
        let obs_fail = create_obs(100, 60, 10, 35.0, 0.85, 0, false);

        assert!(filter.apply(&obs_pass), "confidence >= 0.9 should pass");
        assert!(
            !filter.apply(&obs_fail),
            "confidence >= 0.9 should fail"
        );
    }

    /// Issue #150: edit_distance field comparison
    #[test]
    fn issue_150_expr_edit_distance() {
        let filter = ExprFilter::new("edit_distance <= 5").expect("valid expr");

        let obs_pass = create_obs(100, 60, 10, 35.0, 0.95, 3, false);
        let obs_fail = create_obs(100, 60, 10, 35.0, 0.95, 10, false);

        assert!(
            filter.apply(&obs_pass),
            "edit_distance <= 5 should pass"
        );
        assert!(
            !filter.apply(&obs_fail),
            "edit_distance <= 5 should fail"
        );
    }

    /// Issue #150: in_tandem_repeat boolean field
    #[test]
    fn issue_150_expr_in_tandem_repeat() {
        let filter = ExprFilter::new("in_tandem_repeat == 0").expect("valid expr");

        let obs_pass = create_obs(100, 60, 10, 35.0, 0.95, 0, false);
        let obs_fail = create_obs(100, 60, 10, 35.0, 0.95, 0, true);

        assert!(
            filter.apply(&obs_pass),
            "in_tandem_repeat == 0 should pass for non-repeat variant"
        );
        assert!(
            !filter.apply(&obs_fail),
            "in_tandem_repeat == 0 should fail for repeat variant"
        );
    }

    /// Issue #150: NOT operator
    #[test]
    fn issue_150_expr_not_operator() {
        let filter = ExprFilter::new("!(mapq < 30)").expect("valid expr");

        let obs_pass = create_obs(100, 40, 10, 35.0, 0.95, 0, false);
        let obs_fail = create_obs(100, 20, 10, 35.0, 0.95, 0, false);

        assert!(
            filter.apply(&obs_pass),
            "!(mapq < 30) should pass for mapq >= 30"
        );
        assert!(
            !filter.apply(&obs_fail),
            "!(mapq < 30) should fail for mapq < 30"
        );
    }

    /// Issue #150: Complex expression with multiple operators
    #[test]
    fn issue_150_expr_complex() {
        let filter = ExprFilter::new("(coverage >= 10 && mapq > 20) || base_quality >= 35.0")
            .expect("valid expr");

        let obs1 = create_obs(100, 25, 15, 30.0, 0.95, 0, false); // coverage & mapq ok
        let obs2 = create_obs(100, 15, 5, 36.0, 0.95, 0, false); // base_quality ok
        let obs3 = create_obs(100, 15, 5, 25.0, 0.95, 0, false); // all fail

        assert!(
            filter.apply(&obs1),
            "should pass with coverage & mapq condition"
        );
        assert!(filter.apply(&obs2), "should pass with base_quality");
        assert!(!filter.apply(&obs3), "should fail when all conditions false");
    }

    /// Issue #150: Filter iterator
    #[test]
    fn issue_150_expr_filter_iterator() {
        let filter = ExprFilter::new("coverage >= 10").expect("valid expr");
        let observations = vec![
            create_obs(100, 60, 5, 35.0, 0.95, 0, false),
            create_obs(100, 60, 15, 35.0, 0.95, 0, false),
            create_obs(100, 60, 25, 35.0, 0.95, 0, false),
        ];

        let filtered: Vec<_> = filter.filter(&observations).collect();

        assert_eq!(
            filtered.len(),
            2,
            "filter iterator should return only obs with coverage >= 10"
        );
        assert_eq!(filtered[0].position(), 100);
        assert_eq!(filtered[1].position(), 100);
    }

    /// Issue #150: Equality operator
    #[test]
    fn issue_150_expr_equality() {
        let filter = ExprFilter::new("mapq == 60").expect("valid expr");

        let obs_equal = create_obs(100, 60, 10, 35.0, 0.95, 0, false);
        let obs_not_equal = create_obs(100, 50, 10, 35.0, 0.95, 0, false);

        assert!(
            filter.apply(&obs_equal),
            "mapq == 60 should pass for mapq=60"
        );
        assert!(
            !filter.apply(&obs_not_equal),
            "mapq == 60 should fail for mapq=50"
        );
    }

    /// ExpressionFilter and ThresholdFilter must agree on kmer_uniqueness.
    /// Previously extract_field("kmer_uniqueness") always returned 1.0, so
    /// `kmer_uniqueness < 0.5` always evaluated false while the threshold
    /// filter correctly rejected the same observation.
    #[test]
    fn expr_kmer_uniqueness_agrees_with_threshold_filter() {
        // Observation with kmer_uniqueness = 0.2 (low — non-unique region)
        let low_uniq_obs = {
            let mut alleles = HashMap::new();
            alleles.insert(b'A', 80u32);
            alleles.insert(b'T', 20u32);
            VariantObservation::new(
                100, b'A', alleles, 0.9, "10M".to_string(), 60, 0,
                vec![20], 35.0, "test:read".to_string(),
            ).with_kmer_uniqueness(0.2)
        };

        let expr_filter  = ExprFilter::new("kmer_uniqueness >= 0.5").expect("valid expr");
        let thresh_filter = FilterBuilder::new().min_kmer_uniqueness(0.5).build();

        // Both filters must agree: 0.2 < 0.5, so reject.
        assert!(!expr_filter.apply(&low_uniq_obs),  "expr filter must reject low kmer_uniqueness");
        assert!(!thresh_filter.apply(&low_uniq_obs), "threshold filter must reject low kmer_uniqueness");

        // High uniqueness: both must pass.
        let high_uniq_obs = {
            let mut alleles = HashMap::new();
            alleles.insert(b'A', 80u32);
            alleles.insert(b'T', 20u32);
            VariantObservation::new(
                200, b'A', alleles, 0.9, "10M".to_string(), 60, 0,
                vec![20], 35.0, "test:read".to_string(),
            ).with_kmer_uniqueness(0.8)
        };
        assert!(expr_filter.apply(&high_uniq_obs),  "expr filter must pass high kmer_uniqueness");
        assert!(thresh_filter.apply(&high_uniq_obs), "threshold filter must pass high kmer_uniqueness");
    }

    /// Issue #150: Inequality operator
    #[test]
    fn issue_150_expr_inequality() {
        let filter = ExprFilter::new("mapq != 30").expect("valid expr");

        let obs_pass1 = create_obs(100, 60, 10, 35.0, 0.95, 0, false);
        let obs_pass2 = create_obs(100, 20, 10, 35.0, 0.95, 0, false);
        let obs_fail = create_obs(100, 30, 10, 35.0, 0.95, 0, false);

        assert!(filter.apply(&obs_pass1), "mapq != 30 should pass for mapq=60");
        assert!(
            filter.apply(&obs_pass2),
            "mapq != 30 should pass for mapq=20"
        );
        assert!(
            !filter.apply(&obs_fail),
            "mapq != 30 should fail for mapq=30"
        );
    }
}
