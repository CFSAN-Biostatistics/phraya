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
/// - **conservative**: high-confidence calls only. Good for clinical/typing use.
/// - **sensitive**: catches low-frequency variants at the cost of more noise.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FilterPreset {
    Conservative,
    Sensitive,
}

impl FilterPreset {
    /// Return a `FilterBuilder` seeded with this preset's defaults.
    pub fn builder(self) -> FilterBuilder {
        match self {
            FilterPreset::Conservative => FilterBuilder::new()
                .min_coverage(10)
                .min_mapq(30)
                .min_allele_frequency(0.10)
                .exclude_tandem_repeats(true),
            FilterPreset::Sensitive => FilterBuilder::new()
                .min_coverage(3)
                .min_mapq(20)
                .min_allele_frequency(0.02),
        }
    }
}

/// Threshold-based filter configuration
#[derive(Debug, Clone)]
pub struct FilterBuilder {
    min_coverage: Option<u32>,
    max_coverage: Option<u32>,
    min_mapq: Option<u8>,
    max_mapq: Option<u8>,
    min_base_quality: Option<f64>,
    min_allele_frequency: Option<f64>,
    exclude_tandem_repeats: bool,
}

impl FilterBuilder {
    /// Create a new filter builder with no filters
    pub fn new() -> Self {
        FilterBuilder {
            min_coverage: None,
            max_coverage: None,
            min_mapq: None,
            max_mapq: None,
            min_base_quality: None,
            min_allele_frequency: None,
            exclude_tandem_repeats: false,
        }
    }

    /// Set minimum coverage threshold
    pub fn min_coverage(mut self, threshold: u32) -> Self {
        self.min_coverage = Some(threshold);
        self
    }

    /// Set maximum coverage threshold
    pub fn max_coverage(mut self, threshold: u32) -> Self {
        self.max_coverage = Some(threshold);
        self
    }

    /// Set minimum MAPQ threshold
    pub fn min_mapq(mut self, threshold: u8) -> Self {
        self.min_mapq = Some(threshold);
        self
    }

    /// Set maximum MAPQ threshold
    pub fn max_mapq(mut self, threshold: u8) -> Self {
        self.max_mapq = Some(threshold);
        self
    }

    /// Set minimum base quality threshold
    pub fn min_base_quality(mut self, threshold: f64) -> Self {
        self.min_base_quality = Some(threshold);
        self
    }

    /// Set minimum allele frequency threshold
    pub fn min_allele_frequency(mut self, threshold: f64) -> Self {
        self.min_allele_frequency = Some(threshold);
        self
    }

    /// Exclude variants in tandem repeat regions.
    pub fn exclude_tandem_repeats(mut self, value: bool) -> Self {
        self.exclude_tandem_repeats = value;
        self
    }

    /// Build the filter
    pub fn build(self) -> ThresholdFilter {
        ThresholdFilter {
            min_coverage: self.min_coverage,
            max_coverage: self.max_coverage,
            min_mapq: self.min_mapq,
            max_mapq: self.max_mapq,
            min_base_quality: self.min_base_quality,
            min_allele_frequency: self.min_allele_frequency,
            exclude_tandem_repeats: self.exclude_tandem_repeats,
        }
    }
}

impl Default for FilterBuilder {
    fn default() -> Self {
        Self::new()
    }
}

/// Threshold-based filter for VariantObservations
#[derive(Debug, Clone)]
pub struct ThresholdFilter {
    min_coverage: Option<u32>,
    max_coverage: Option<u32>,
    min_mapq: Option<u8>,
    max_mapq: Option<u8>,
    min_base_quality: Option<f64>,
    min_allele_frequency: Option<f64>,
    exclude_tandem_repeats: bool,
}

impl ThresholdFilter {
    /// Apply filter to an observation
    pub fn apply(&self, obs: &VariantObservation) -> bool {
        // Check coverage (local_coverage[0] as proxy for coverage)
        if let Some(min) = self.min_coverage {
            if obs.local_coverage().is_empty() || obs.local_coverage()[0] < min {
                return false;
            }
        }

        if let Some(max) = self.max_coverage {
            if !obs.local_coverage().is_empty() && obs.local_coverage()[0] > max {
                return false;
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

        // Exclude tandem repeat variants
        if self.exclude_tandem_repeats && obs.in_tandem_repeat() {
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
    fn max_coverage_filter() {
        let obs_low = create_observation(100, 60, 5, 35.0);
        let obs_high = create_observation(100, 60, 25, 35.0);

        let filter = FilterBuilder::new().max_coverage(20).build();

        assert!(filter.apply(&obs_low));
        assert!(!filter.apply(&obs_high));
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
    fn conservative_preset_rejects_low_quality() {
        // obs with coverage=5, mapq=25, allele_freq=0.05 — below all conservative thresholds
        let mut alleles = HashMap::new();
        alleles.insert(b'T', 5u32); // 5/100 = 5% alt freq — below 10%
        alleles.insert(b'A', 95u32);
        let low_quality = VariantObservation::new(
            100, b'A', alleles, 0.5, "10M".to_string(), 25, 1, vec![100], 30.0, "s:r".to_string(),
        );

        let filter = FilterPreset::Conservative.builder().build();
        assert!(!filter.apply(&low_quality), "conservative must reject low-quality obs");
    }

    #[test]
    fn conservative_preset_passes_high_quality() {
        let mut alleles = HashMap::new();
        alleles.insert(b'T', 20u32); // 20/100 = 20% alt freq — above 10%
        alleles.insert(b'A', 80u32);
        let high_quality = VariantObservation::new(
            100, b'A', alleles, 0.95, "10M".to_string(), 40, 0, vec![100], 35.0, "s:r".to_string(),
        );

        let filter = FilterPreset::Conservative.builder().build();
        assert!(filter.apply(&high_quality), "conservative must pass high-quality obs");
    }

    #[test]
    fn sensitive_preset_passes_low_coverage_variant() {
        // coverage=3, mapq=20, allele_freq=0.03 — just at sensitive thresholds
        let mut alleles = HashMap::new();
        alleles.insert(b'T', 3u32);
        alleles.insert(b'A', 97u32); // 3% alt freq
        let obs = VariantObservation::new(
            100, b'A', alleles, 0.9, "10M".to_string(), 20, 1, vec![3], 30.0, "s:r".to_string(),
        );

        let filter = FilterPreset::Sensitive.builder().build();
        assert!(filter.apply(&obs), "sensitive must pass low-coverage variants");
    }

    #[test]
    fn sensitive_preset_rejects_below_its_thresholds() {
        // mapq=19 — below sensitive's min_mapq=20
        let obs = create_observation(100, 19, 5, 30.0);
        let filter = FilterPreset::Sensitive.builder().build();
        assert!(!filter.apply(&obs), "sensitive must still reject mapq<20");
    }

    #[test]
    fn preset_can_be_overridden_by_explicit_threshold() {
        // Conservative preset has min_coverage=10. Override to min_coverage=5.
        // Use an obs that passes all other conservative thresholds (mapq=40, allele_freq=20%).
        let mut alleles = HashMap::new();
        alleles.insert(b'T', 2u32); // 2/10 = 20% alt freq → passes conservative's 10% threshold
        alleles.insert(b'A', 8u32);
        // local_coverage[0]=7 → below conservative's min_coverage=10, above override of 5
        let obs = VariantObservation::new(
            100, b'A', alleles, 0.95, "10M".to_string(), 40, 0, vec![7], 35.0, "s:r".to_string(),
        );

        let default_filter = FilterPreset::Conservative.builder().build();
        assert!(!default_filter.apply(&obs), "local_coverage[0]=7 fails conservative default (min=10)");

        let overridden = FilterPreset::Conservative.builder().min_coverage(5).build();
        assert!(overridden.apply(&obs), "local_coverage[0]=7 passes after override to min_coverage=5");
    }

    #[test]
    fn conservative_excludes_tandem_repeats_by_default() {
        let mut alleles = HashMap::new();
        alleles.insert(b'T', 20u32);
        alleles.insert(b'A', 80u32);
        let obs = VariantObservation::new(
            100, b'A', alleles, 0.95, "10M".to_string(), 40, 0, vec![100], 35.0, "s:r".to_string(),
        ).with_tandem_repeat(true);

        let filter = FilterPreset::Conservative.builder().build();
        assert!(!filter.apply(&obs), "conservative must exclude tandem repeat variants");
    }

    #[test]
    fn sensitive_does_not_exclude_tandem_repeats_by_default() {
        let mut alleles = HashMap::new();
        alleles.insert(b'T', 5u32);
        alleles.insert(b'A', 95u32);
        let obs = VariantObservation::new(
            100, b'A', alleles, 0.9, "10M".to_string(), 25, 0, vec![10], 30.0, "s:r".to_string(),
        ).with_tandem_repeat(true);

        let filter = FilterPreset::Sensitive.builder().build();
        assert!(filter.apply(&obs), "sensitive must not exclude tandem repeat variants");
    }
}
