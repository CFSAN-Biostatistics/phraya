/// RED acceptance tests for issue #145 — variation hotspot estimation at plan time.
/// All tests call `detect_hotspot_intervals` which is `unimplemented!()` in production code.

#[cfg(test)]
mod tests {
    use crate::types::{
        compute_kmer_uniqueness, detect_hotspot_intervals, sketch_sequence_default, Sequence,
    };
    use std::collections::HashMap;

    fn make_seq(bases: &str, id: &str) -> Sequence {
        Sequence::new(bases.as_bytes().to_vec(), None, id.to_string(), None)
    }

    /// issue #145: unique sequence with no repeats → no intervals emitted
    #[test]
    fn issue_145_unique_sequence_no_hotspots() {
        let seq = make_seq("ACGTACGTACGTACGTACGTACGTACGTACGT", "unique");
        let uniqueness = compute_kmer_uniqueness(&[sketch_sequence_default(&seq)]);
        let intervals = detect_hotspot_intervals(&uniqueness, 0.5);
        assert!(
            intervals.is_empty(),
            "unique sequence should produce no hotspot intervals"
        );
    }

    /// issue #145: synthetic tandem repeat → at least one interval detected
    #[test]
    fn issue_145_tandem_repeat_produces_hotspot_interval() {
        // Two identical sketches → every shared minimizer has uniqueness = 0.5 < threshold
        let seq = make_seq(
            "ACGTACGTACGTATATATATATATATATATATATATACGTACGTACGT",
            "repeat",
        );
        let sketch = sketch_sequence_default(&seq);
        let uniqueness = compute_kmer_uniqueness(&[sketch.clone(), sketch]);
        let intervals = detect_hotspot_intervals(&uniqueness, 0.6);
        assert!(
            !intervals.is_empty(),
            "repeated sequence should produce at least one hotspot interval"
        );
    }

    /// issue #145: empty uniqueness map → empty intervals
    #[test]
    fn issue_145_empty_uniqueness_map_returns_empty() {
        let intervals = detect_hotspot_intervals(&HashMap::new(), 0.5);
        assert!(intervals.is_empty());
    }

    /// issue #145: all positions below threshold → one spanning interval
    #[test]
    fn issue_145_all_positions_below_threshold_yields_one_interval() {
        let mut map = HashMap::new();
        map.insert(0u32, 0.1);
        map.insert(1u32, 0.2);
        map.insert(2u32, 0.3);
        let intervals = detect_hotspot_intervals(&map, 0.5);
        assert_eq!(intervals.len(), 1, "contiguous low-uniqueness run → one interval");
        assert_eq!(intervals[0], (0, 2));
    }

    /// issue #145: gap between low-uniqueness positions → two separate intervals
    #[test]
    fn issue_145_gap_separates_into_two_intervals() {
        let mut map = HashMap::new();
        map.insert(0u32, 0.1);
        map.insert(1u32, 0.1);
        map.insert(5u32, 0.1);
        map.insert(6u32, 0.1);
        let intervals = detect_hotspot_intervals(&map, 0.5);
        assert_eq!(intervals.len(), 2, "gap should separate into two intervals");
    }

    /// issue #145: adjacent positions (consecutive) merge into one interval
    #[test]
    fn issue_145_adjacent_positions_merge_into_one_interval() {
        let mut map = HashMap::new();
        map.insert(10u32, 0.1);
        map.insert(11u32, 0.2);
        map.insert(12u32, 0.3);
        let intervals = detect_hotspot_intervals(&map, 0.5);
        assert_eq!(intervals.len(), 1);
        assert_eq!(intervals[0], (10, 12));
    }

    /// issue #145: positions above threshold are excluded
    #[test]
    fn issue_145_positions_above_threshold_excluded() {
        let mut map = HashMap::new();
        map.insert(0u32, 0.1);
        map.insert(1u32, 0.9); // above threshold
        map.insert(2u32, 0.1);
        let intervals = detect_hotspot_intervals(&map, 0.5);
        // positions 0 and 2 are not adjacent, so two intervals
        assert_eq!(intervals.len(), 2);
    }

    /// issue #145: intervals are returned sorted by start position
    #[test]
    fn issue_145_intervals_sorted_by_start() {
        let mut map = HashMap::new();
        map.insert(100u32, 0.1);
        map.insert(0u32, 0.1);
        map.insert(50u32, 0.1);
        let intervals = detect_hotspot_intervals(&map, 0.5);
        for w in intervals.windows(2) {
            assert!(w[0].0 < w[1].0, "intervals must be sorted");
        }
    }
}
