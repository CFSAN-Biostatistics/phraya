/// RED Acceptance Tests for Issue #130: Implement real local coverage ±50bp window
///
/// These tests verify that:
/// 1. VariantObservation.local_coverage contains real read depth values, not placeholder [1]
/// 2. Local coverage reflects the actual ±50bp window around variant sites
/// 3. Coverage values match the expected window size for each variant position

#[cfg(test)]
mod issue_130_local_coverage_tests {
    use crate::executor::align_task;
    use phraya_core::types::{Sequence, VariantObservation};
    use phraya_io::plan::PhrayaPlan;
    use std::collections::HashMap;

    /// Helper: create a plan for testing
    fn create_test_plan() -> PhrayaPlan {
        PhrayaPlan::new(
            phraya_io::plan::UseCase::ReadsWithRef,
            vec!["test".to_string()],
            "2026-06-03T00:00:00Z".to_string(),
            vec![],
            HashMap::new(),
            vec![],
        )
    }

    /// Test that local_coverage contains the ±50bp window, not just [1]
    /// When a variant is found at position P, local_coverage should contain
    /// per-position read depth for positions [max(0, P-50), min(len, P+50)].
    #[test]
    fn issue_130_local_coverage_window_at_middle_position() {
        // Sequence long enough that window fits entirely (position 100)
        let query = Sequence::new(
            b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
            None,
            "query1".to_string(),
            None,
        );
        let mut target_bases = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
        target_bases[100] = b'T'; // Create variant at position 100
        let target = Sequence::new(target_bases, None, "ref".to_string(), None);
        let plan = create_test_plan();

        let result = align_task(&query, &target, &plan);
        assert!(result.is_some());

        let result = result.unwrap();
        // Find variant at position 100
        let variant = result.variants.iter().find(|v| v.position() == 100);
        assert!(variant.is_some(), "Should find variant at position 100");

        let var = variant.unwrap();
        let local_coverage = var.local_coverage();

        // At position 100 with ±50bp window: window is [50, 150]
        // That's 101 positions (50 through 150 inclusive)
        let expected_window_size = 101;
        assert_eq!(
            local_coverage.len(),
            expected_window_size,
            "At position 100, ±50bp window should contain {} values, got {}",
            expected_window_size,
            local_coverage.len()
        );

        // Should NOT be the placeholder [1]
        assert!(
            !(local_coverage.len() == 1 && local_coverage[0] == 1),
            "local_coverage should be the ±50bp window, not placeholder vec![1]"
        );
    }

    /// Test local_coverage window at start of sequence (position 0)
    /// Window at position 0 is [max(0, -50), min(len, 50)] = [0, 50] = 51 values
    #[test]
    fn issue_130_local_coverage_window_at_start() {
        let query = Sequence::new(
            b"TACTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
            None,
            "query1".to_string(),
            None,
        );
        let mut target_bases = b"TACTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
        target_bases[0] = b'A'; // Variant at position 0
        let target = Sequence::new(target_bases, None, "ref".to_string(), None);
        let plan = create_test_plan();

        let result = align_task(&query, &target, &plan);
        assert!(result.is_some());

        let result = result.unwrap();
        let variant = result.variants.iter().find(|v| v.position() == 0);
        assert!(variant.is_some());

        let var = variant.unwrap();
        let local_coverage = var.local_coverage();

        // At position 0: window is [0, 50] = 51 values
        let expected_window_size = 51;
        assert_eq!(
            local_coverage.len(),
            expected_window_size,
            "At position 0, ±50bp window should be 51 values"
        );
    }

    /// Test local_coverage window at end of sequence
    /// For a 150bp sequence with variant at position 149, window is [99, 149] = 51 values
    #[test]
    fn issue_130_local_coverage_window_at_end() {
        let seq_len = 150usize;
        let query = Sequence::new(
            b"A".to_vec().repeat(seq_len),
            None,
            "query1".to_string(),
            None,
        );
        let mut target_bases = b"A".to_vec().repeat(seq_len);
        target_bases[seq_len - 1] = b'T'; // Variant at last position
        let target = Sequence::new(target_bases, None, "ref".to_string(), None);
        let plan = create_test_plan();

        let result = align_task(&query, &target, &plan);
        assert!(result.is_some());

        let result = result.unwrap();
        let variant = result
            .variants
            .iter()
            .find(|v| v.position() == (seq_len - 1) as u32);
        assert!(variant.is_some());

        let var = variant.unwrap();
        let local_coverage = var.local_coverage();

        // At position 149 in 150bp sequence: window is [99, 149] = 51 values
        let expected_window_size = 51;
        assert_eq!(
            local_coverage.len(),
            expected_window_size,
            "At end position, ±50bp window should be 51 values"
        );
    }

    /// Test that local_coverage is NOT just the hardcoded placeholder [1]
    /// The old implementation hardcoded vec![1] in executor.rs line 49
    #[test]
    fn issue_130_local_coverage_not_placeholder() {
        let query = Sequence::new(
            b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
            None,
            "query1".to_string(),
            None,
        );
        let mut target_bases = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
        target_bases[100] = b'T'; // Variant at position 100
        let target = Sequence::new(target_bases, None, "ref".to_string(), None);
        let plan = create_test_plan();

        let result = align_task(&query, &target, &plan);
        assert!(result.is_some());

        let result = result.unwrap();
        assert!(
            !result.variants.is_empty(),
            "Should have found at least one variant"
        );

        for var in &result.variants {
            let local_coverage = var.local_coverage();

            // MUST NOT be the old placeholder vec![1]
            // Real coverage should be a window, typically 51-101 values
            let is_placeholder = local_coverage.len() == 1 && local_coverage[0] == 1;
            assert!(
                !is_placeholder,
                "local_coverage must not be placeholder vec![1] \
                 at position {}, got {:?}",
                var.position(),
                local_coverage
            );
        }
    }

    /// Test that local_coverage reflects per-position depth, not count of observations
    /// Issue: "Coverage is count-of-reads-with-variant-at-position, not actual local read depth"
    /// After fix: each value in local_coverage should represent actual depth at that position
    #[test]
    fn issue_130_local_coverage_per_position_depth() {
        let query = Sequence::new(
            b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
            None,
            "query1".to_string(),
            None,
        );
        let mut target_bases = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
        target_bases[75] = b'T'; // Create variant
        let target = Sequence::new(target_bases, None, "ref".to_string(), None);
        let plan = create_test_plan();

        let result = align_task(&query, &target, &plan);
        assert!(result.is_some());

        let result = result.unwrap();
        let variant = result.variants.iter().find(|v| v.position() == 75);
        assert!(variant.is_some());

        let var = variant.unwrap();
        let local_coverage = var.local_coverage();

        // The window must reflect per-position depth
        // For a single aligned read, depth should be 1 at all positions in window
        // But it must be a vector of positions, not a single [1] placeholder

        // Window at position 75: [25, 125] = 101 values
        assert_eq!(
            local_coverage.len(),
            101,
            "Window should have 101 positions, one for each base in ±50bp window"
        );

        // Each value should be a depth (typically 1 for single read)
        // The key point: we have per-position values, not a single observation count
        for (i, &depth) in local_coverage.iter().enumerate() {
            assert!(
                depth > 0,
                "Coverage at window position {} should be > 0, got {}",
                i,
                depth
            );
        }
    }

    /// Test multiple variants have correct window sizes
    /// Each variant at different position should get its own properly-sized window
    #[test]
    fn issue_130_multiple_variants_correct_windows() {
        let query = Sequence::new(
            b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
            None,
            "query1".to_string(),
            None,
        );
        let mut target_bases = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
        target_bases[30] = b'T'; // Variant at position 30
        target_bases[100] = b'T'; // Variant at position 100
        let target = Sequence::new(target_bases, None, "ref".to_string(), None);
        let plan = create_test_plan();

        let result = align_task(&query, &target, &plan);
        assert!(result.is_some());

        let result = result.unwrap();
        assert_eq!(result.variants.len(), 2, "Should have 2 variants");

        for var in &result.variants {
            let pos = var.position();
            let local_coverage = var.local_coverage();

            // Each should have properly-sized window
            let window_start = (pos as i32 - 50).max(0) as u32;
            let window_end = std::cmp::min(pos + 51, query.len() as u32);
            let expected_size = (window_end - window_start) as usize;

            assert_eq!(
                local_coverage.len(),
                expected_size,
                "Variant at position {} should have window size {}, got {}",
                pos,
                expected_size,
                local_coverage.len()
            );

            // Should NOT be placeholder
            assert!(
                !(local_coverage.len() == 1 && local_coverage[0] == 1),
                "Variant at position {} should not have placeholder vec![1]",
                pos
            );
        }
    }

    /// Test that local_coverage[0] is the center position depth
    /// After fix, local_coverage[0] should be the depth at or near the variant position,
    /// not "count of observations at this position" (which would be 1 before merge)
    #[test]
    fn issue_130_coverage_center_meaningful() {
        let query = Sequence::new(
            b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
            None,
            "query1".to_string(),
            None,
        );
        let mut target_bases = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
        target_bases[90] = b'T'; // Variant at position 90
        let target = Sequence::new(target_bases, None, "ref".to_string(), None);
        let plan = create_test_plan();

        let result = align_task(&query, &target, &plan);
        assert!(result.is_some());

        let result = result.unwrap();
        let variant = result.variants.iter().find(|v| v.position() == 90);
        assert!(variant.is_some());

        let var = variant.unwrap();
        let local_coverage = var.local_coverage();

        // Window at position 90: [40, 140]
        // local_coverage[0] should be depth at position 40 (window start)
        // local_coverage[50] should be depth at position 90 (variant center)

        assert!(local_coverage.len() > 50, "Window should be large enough");

        // The first element and center elements should all be reasonable depth values
        // For single read: should be 1
        // But crucially: they should exist as separate values per position
        // NOT be a single [1] placeholder

        for &depth in local_coverage {
            assert!(depth > 0, "All coverage values must be > 0");
        }
    }
}
