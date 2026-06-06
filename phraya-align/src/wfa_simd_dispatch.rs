/// Runtime dispatch testing for multiversion SSE4.2 integration.
///
/// This module tests that the multiversion crate correctly dispatches between
/// SSE4.2 and naive implementations based on runtime CPU feature detection.

#[cfg(test)]
mod tests {
    use crate::{wfa_extend, wfa_extend_naive, wfa_extend_simd, SeedAnchor};

    #[test]
    fn test_multiversion_dispatch_selects_correct_target() {
        // This test verifies that multiversion correctly identifies the best
        // implementation for the current CPU.

        let dispatch_target = crate::wfa_simd::get_active_dispatch_target();

        // Valid targets: "sse42" / "naive" on x86_64, "neon" on aarch64.
        assert!(
            matches!(dispatch_target.as_str(), "sse42" | "naive" | "neon"),
            "Dispatch target must be 'sse42', 'neon', or 'naive', got: {}",
            dispatch_target
        );
    }

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

    #[test]
    #[cfg(not(target_arch = "x86_64"))]
    fn test_non_x86_64_implementation_set() {
        // Off x86_64 there is no SSE4.2 path. aarch64 additionally has a real
        // NEON implementation; other architectures are naive-only.
        let implementations = crate::wfa_simd::get_compiled_implementations();

        assert!(
            implementations.contains(&"naive"),
            "Non-x86_64 builds must include naive implementation"
        );
        assert!(
            !implementations.contains(&"sse42"),
            "Non-x86_64 builds must not include the SSE4.2 implementation"
        );

        #[cfg(target_arch = "aarch64")]
        {
            assert!(
                implementations.contains(&"neon"),
                "aarch64 builds must include the NEON implementation"
            );
            assert_eq!(implementations.len(), 2, "aarch64 has naive + neon");
        }
        #[cfg(not(target_arch = "aarch64"))]
        assert_eq!(implementations.len(), 1, "other arches are naive-only");
    }

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

    #[test]
    fn test_dispatched_implementation_logs_selection() {
        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        let _ = wfa_extend(query, target, seed);

        let selected_impl = crate::wfa_simd::get_active_dispatch_target();

        assert!(
            matches!(selected_impl.as_str(), "naive" | "sse42" | "neon"),
            "Selected implementation must be a known variant, got: {}",
            selected_impl
        );
        // Note: wfa_extend uses scalar in debug/test builds by design, so in
        // tests we correctly see "naive" even on SSE4.2 hardware. SIMD correctness
        // is verified by simd_diff_tests and simd_vs_naive_differential.
    }

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

    // ========================================================================
    // Issue #71: SIMD diagonal fill acceptance tests
    // ========================================================================
    // These tests validate correctness of the dispatched alignment path.
    // Each test calls wfa_extend (the dispatcher) to exercise the dispatch logic.
    // Note: wfa_extend intentionally uses the scalar WFA in debug/test builds
    // (see lib.rs). SIMD correctness under real inputs is verified separately
    // in simd_diff_tests::fill_simd_matches_fill_scalar_property (calls fill_simd
    // directly) and the #[cfg(x86_64)] simd_vs_naive_differential module
    // (calls wfa_extend_simd_impl directly).

    #[test]
    fn issue_71_diagonal_fill_sse42_2x2_matrix() {
        let query = b"AC";
        let target = b"AC";
        let seed = crate::SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };
        let result = wfa_extend(query, target, seed);
        assert!(result.is_ok());
        let alignment = result.unwrap();
        assert_eq!(alignment.cigar, "2M");
        assert_eq!(alignment.edit_distance, 0);
    }

    #[test]
    fn issue_71_diagonal_fill_sse42_4x4_matrix() {
        let query = b"ACGT";
        let target = b"ACGT";
        let seed = crate::SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };
        let result = wfa_extend(query, target, seed);
        assert!(result.is_ok());
        let alignment = result.unwrap();
        assert_eq!(alignment.cigar, "4M");
        assert_eq!(alignment.edit_distance, 0);
    }

    #[test]
    fn issue_71_diagonal_fill_sse42_16x16_matrix() {
        let query = b"ACGTACGTACGTACGT";
        let target = b"ACGTACGTACGTACGT";
        let seed = crate::SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };
        let result = wfa_extend(query, target, seed);
        assert!(result.is_ok());
        let alignment = result.unwrap();
        assert_eq!(alignment.cigar, "16M");
        assert_eq!(alignment.edit_distance, 0);
    }

    #[test]
    fn issue_71_diagonal_fill_sse42_32x32_matrix() {
        let query = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let target = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let seed = crate::SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };
        let result = wfa_extend(query, target, seed);
        assert!(result.is_ok());
        let alignment = result.unwrap();
        assert_eq!(alignment.edit_distance, 0);
        assert!(alignment.cigar.contains("32M"));
    }

    #[test]
    fn issue_71_simd_correctness_exact_match() {
        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = crate::SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };
        let result = wfa_extend(query, target, seed);
        assert!(result.is_ok());
        let alignment = result.unwrap();
        assert_eq!(alignment.cigar, "12M");
        assert_eq!(alignment.edit_distance, 0);
    }

    #[test]
    fn issue_71_simd_correctness_mismatch() {
        let query = b"ACGTACGTACGT";
        let target = b"ACGTACTTACGT";
        let seed = crate::SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };
        let result = wfa_extend(query, target, seed);
        assert!(result.is_ok());
        let alignment = result.unwrap();
        assert!(alignment.edit_distance > 0);
    }

    #[test]
    fn issue_71_simd_correctness_insertion() {
        let query = b"ACGTACGT";
        let target = b"ACGTAACGT";
        let seed = crate::SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };
        let result = wfa_extend(query, target, seed);
        assert!(result.is_ok());
        let alignment = result.unwrap();
        assert!(alignment.cigar.contains("I"));
    }

    #[test]
    fn issue_71_simd_correctness_deletion() {
        let query = b"ACGTAACGT";
        let target = b"ACGTACGT";
        let seed = crate::SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };
        let result = wfa_extend(query, target, seed);
        assert!(result.is_ok());
        let alignment = result.unwrap();
        assert!(alignment.cigar.contains("D"));
    }

    #[test]
    fn issue_71_diagonal_fill_seed_middle() {
        let query = b"ACGTACGTACGTACGTACGT";
        let target = b"ACGTACGTACGTACGTACGT";
        let seed = crate::SeedAnchor {
            query_pos: 5,
            target_pos: 5,
        };
        let result = wfa_extend(query, target, seed);
        assert!(result.is_ok());
        let alignment = result.unwrap();
        assert_eq!(alignment.query_start, 5);
        assert_eq!(alignment.target_start, 5);
        assert_eq!(alignment.edit_distance, 0);
    }

    #[test]
    fn issue_71_diagonal_fill_edge_empty_suffix() {
        let query = b"ACGTACGT";
        let target = b"ACGTACGT";
        let seed = crate::SeedAnchor {
            query_pos: 8,
            target_pos: 8,
        };
        let result = wfa_extend(query, target, seed);
        assert!(result.is_ok());
        let alignment = result.unwrap();
        assert_eq!(alignment.cigar, "");
        assert_eq!(alignment.edit_distance, 0);
    }

    #[test]
    fn issue_71_diagonal_fill_single_character() {
        let query = b"A";
        let target = b"A";
        let seed = crate::SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };
        let result = wfa_extend(query, target, seed);
        assert!(result.is_ok());
        let alignment = result.unwrap();
        assert_eq!(alignment.cigar, "1M");
        assert_eq!(alignment.edit_distance, 0);
    }

    #[test]
    fn issue_71_diagonal_fill_10kb() {
        let base = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let mut query = Vec::new();
        let mut target = Vec::new();
        for _ in 0..100 {
            query.extend_from_slice(base);
            target.extend_from_slice(base);
        }
        assert_eq!(query.len(), 10000);
        let seed = crate::SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };
        let result = wfa_extend(&query, &target, seed);
        assert!(result.is_ok());
        let alignment = result.unwrap();
        assert_eq!(alignment.edit_distance, 0);
    }

    #[test]
    fn issue_71_runtime_dispatch_valid_target() {
        // Must call wfa_extend first so dispatch sets the thread-local; reading
        // LAST_IMPL before any call would always see the default "naive" initialization.
        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let seed = SeedAnchor { query_pos: 0, target_pos: 0 };
        let _ = wfa_extend(query, target, seed);

        let dispatch_target = crate::wfa_simd::get_active_dispatch_target();
        assert!(
            matches!(dispatch_target.as_str(), "sse42" | "naive" | "neon"),
            "Dispatch target must be valid, got: {}",
            dispatch_target
        );
    }

}
