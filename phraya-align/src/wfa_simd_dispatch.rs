/// Runtime dispatch testing for multiversion SSE4.2 integration.
///
/// This module tests that the multiversion crate correctly dispatches between
/// SSE4.2 and naive implementations based on runtime CPU feature detection.

#[cfg(test)]
mod tests {
    use crate::{wfa_extend, wfa_extend_naive, wfa_extend_simd, SeedAnchor};

    // Test will fail: get_active_dispatch_target function does not exist yet
    #[test]
    fn test_multiversion_dispatch_selects_correct_target() {
        // This test verifies that multiversion correctly identifies the best
        // implementation for the current CPU.

        let dispatch_target = crate::wfa_simd::get_active_dispatch_target();

        // Should be either "sse42" or "naive" depending on CPU
        assert!(
            dispatch_target == "sse42" || dispatch_target == "naive",
            "Dispatch target must be 'sse42' or 'naive', got: {}",
            dispatch_target
        );
    }

    // Test will fail: is_sse42_available function does not exist yet
    #[test]
    fn test_sse42_feature_detection() {
        // This test verifies that SSE4.2 feature detection works correctly.

        let sse42_available = crate::wfa_simd::is_sse42_available();

        // Detection should be deterministic and match CPU capabilities
        assert!(
            sse42_available == crate::wfa_simd::is_sse42_available(),
            "SSE4.2 detection must be consistent across calls"
        );
    }

    // Test will fail: wfa_extend function does not exist yet
    #[test]
    fn test_dispatched_function_matches_manual_selection() {
        // This test verifies that the multiversion-dispatched function
        // produces the same result as manually selecting the implementation.

        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let dispatched = wfa_extend(query, target, seed.clone());

        // Manually select based on feature detection
        let manual = if crate::wfa_simd::is_sse42_available() {
            wfa_extend_simd(query, target, seed)
        } else {
            wfa_extend_naive(query, target, seed)
        };

        assert!(dispatched.is_ok());
        assert!(manual.is_ok());

        let dispatched_result = dispatched.unwrap();
        let manual_result = manual.unwrap();

        assert_eq!(dispatched_result.cigar, manual_result.cigar);
        assert_eq!(dispatched_result.edit_distance, manual_result.edit_distance);
    }

    // Test will fail: multiversion attribute does not exist yet
    #[test]
    #[cfg(target_arch = "x86_64")]
    fn test_x86_64_has_multiple_implementations() {
        // On x86_64, we should have both naive and SSE4.2 implementations compiled in.

        let implementations = crate::wfa_simd::get_compiled_implementations();

        assert!(
            implementations.contains(&"naive"),
            "x86_64 builds must include naive implementation"
        );

        assert!(
            implementations.contains(&"sse42"),
            "x86_64 builds must include SSE4.2 implementation"
        );

        assert!(
            implementations.len() >= 2,
            "x86_64 builds should have at least 2 implementations (naive + sse42)"
        );
    }

    // Test will fail: multiversion fallback does not exist yet
    #[test]
    #[cfg(not(target_arch = "x86_64"))]
    fn test_non_x86_64_uses_naive_only() {
        // On non-x86_64 architectures, only the naive implementation should be available.

        let implementations = crate::wfa_simd::get_compiled_implementations();

        assert!(
            implementations.contains(&"naive"),
            "Non-x86_64 builds must include naive implementation"
        );

        assert!(
            !implementations.contains(&"sse42"),
            "Non-x86_64 builds should not include SSE4.2 implementation"
        );

        assert_eq!(
            implementations.len(),
            1,
            "Non-x86_64 builds should have exactly 1 implementation (naive)"
        );
    }

    // Test will fail: dispatch overhead tracking does not exist yet
    #[test]
    fn test_dispatch_has_minimal_overhead() {
        // This test verifies that the runtime dispatch overhead is negligible.
        // Multiversion should use static dispatch after first call.

        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        // First call - may include dispatch overhead
        let _ = wfa_extend(query, target, seed.clone());

        // Subsequent calls should have no dispatch overhead
        let iterations = 1000;
        let start = std::time::Instant::now();

        for _ in 0..iterations {
            let _ = wfa_extend(query, target, seed.clone());
        }

        let elapsed = start.elapsed();
        let avg_time = elapsed / iterations;

        // Dispatch overhead should be negligible (< 1 microsecond per call)
        assert!(
            avg_time.as_nanos() < 1_000_000, // 1ms is very generous
            "Dispatch overhead too high: {:?} average per call",
            avg_time
        );
    }

    // Test will fail: force_implementation function does not exist yet
    #[test]
    #[cfg(target_arch = "x86_64")]
    fn test_can_force_naive_implementation() {
        // This test verifies that we can explicitly force the naive implementation
        // even on SSE4.2-capable CPUs (useful for testing and validation).

        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        // Force naive implementation
        let naive_result =
            crate::wfa_simd::force_implementation("naive", query, target, seed.clone());

        assert!(naive_result.is_ok());
    }

    // Test will fail: force_implementation function does not exist yet
    #[test]
    #[cfg(target_arch = "x86_64")]
    fn test_can_force_sse42_implementation() {
        // This test verifies that we can explicitly force the SSE4.2 implementation
        // if available (useful for testing and benchmarking).

        if !crate::wfa_simd::is_sse42_available() {
            // Skip if SSE4.2 not available
            return;
        }

        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        // Force SSE4.2 implementation
        let simd_result =
            crate::wfa_simd::force_implementation("sse42", query, target, seed.clone());

        assert!(simd_result.is_ok());
    }

    // Test will fail: implementation verification does not exist yet
    #[test]
    fn test_dispatched_implementation_logs_selection() {
        // This test verifies that the selected implementation is logged or
        // can be queried for debugging purposes.

        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let _ = wfa_extend(query, target, seed);

        let selected_impl = crate::wfa_simd::get_last_selected_implementation();

        assert!(
            !selected_impl.is_empty(),
            "Should be able to query which implementation was selected"
        );

        assert!(
            selected_impl == "naive" || selected_impl == "sse42",
            "Selected implementation must be valid: {}",
            selected_impl
        );
    }

    // Test will fail: multiversion macro usage does not exist yet
    #[test]
    fn test_multiversion_macro_properly_applied() {
        // This test verifies that the multiversion macro is correctly applied
        // to the wfa_extend function with appropriate targets.

        let has_multiversion = crate::wfa_simd::has_multiversion_attribute();

        assert!(
            has_multiversion,
            "wfa_extend must use multiversion attribute for runtime dispatch"
        );
    }

    // Test will fail: CPUID access does not exist yet
    #[test]
    #[cfg(target_arch = "x86_64")]
    fn test_cpuid_detection_works() {
        // This test verifies that CPUID detection is working correctly
        // and can identify SSE4.2 support.

        let cpuid_result = crate::wfa_simd::query_cpuid_features();

        // Should have basic feature information
        assert!(cpuid_result.contains_key("sse42"));

        // Result should match is_sse42_available()
        assert_eq!(cpuid_result["sse42"], crate::wfa_simd::is_sse42_available());
    }

    // Test will fail: fallback mechanism does not exist yet
    #[test]
    fn test_fallback_on_unsupported_feature() {
        // This test verifies that attempting to use an unsupported feature
        // gracefully falls back to naive implementation.

        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        // Should not panic, should fall back to naive
        let result = wfa_extend(query, target, seed);

        assert!(
            result.is_ok(),
            "Should gracefully fall back to naive implementation"
        );
    }
}
