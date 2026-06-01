use phraya_core::types::VariantObservation;
use std::collections::HashMap;

pub mod tsv;
pub mod vcf;

/// Threshold-based filter configuration
#[derive(Debug, Clone)]
pub struct FilterBuilder {
    min_coverage: Option<u32>,
    max_coverage: Option<u32>,
    min_mapq: Option<u8>,
    max_mapq: Option<u8>,
    min_base_quality: Option<f64>,
    min_allele_frequency: Option<f64>,
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

    /// Build the filter
    pub fn build(self) -> ThresholdFilter {
        ThresholdFilter {
            min_coverage: self.min_coverage,
            max_coverage: self.max_coverage,
            min_mapq: self.min_mapq,
            max_mapq: self.max_mapq,
            min_base_quality: self.min_base_quality,
            min_allele_frequency: self.min_allele_frequency,
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

            let ref_count = obs.all_alleles().get(&obs.ref_base()).copied().unwrap_or(0);
            let ref_freq = ref_count as f64 / total as f64;

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
}
