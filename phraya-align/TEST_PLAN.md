# Test Plan for Issue #13: SSE4.2 SIMD for WFA diagonal fill

## Overview

This test plan covers all acceptance criteria for implementing SSE4.2 SIMD optimization
for WFA diagonal fill with runtime dispatch via the `multiversion` crate.

## Acceptance Criteria Coverage

### 1. SSE4.2 implementation of WFA diagonal fill using `core::arch::x86_64` intrinsics

**Test Files**: `src/wfa_simd.rs`

**Tests**:
- `test_simd_exact_match` - Verifies SIMD handles exact sequence matches
- `test_simd_single_mismatch` - Verifies SIMD handles single mismatches
- `test_simd_insertion` - Verifies SIMD handles insertions
- `test_simd_deletion` - Verifies SIMD handles deletions
- `test_simd_complex_alignment` - Verifies SIMD handles complex mixed events

These tests verify that the SSE4.2 implementation produces correct alignment results
for various sequence patterns.

### 2. Uses `multiversion` crate for runtime CPUID dispatch (SSE4.2 vs fallback to naive)

**Test Files**: `src/wfa_simd_dispatch.rs`

**Tests**:
- `test_multiversion_dispatch_selects_correct_target` - Verifies dispatch selects correct implementation
- `test_sse42_feature_detection` - Verifies SSE4.2 detection works
- `test_dispatched_function_matches_manual_selection` - Verifies dispatch matches manual selection
- `test_x86_64_has_multiple_implementations` - Verifies both implementations compiled on x86_64
- `test_non_x86_64_uses_naive_only` - Verifies only naive on non-x86 platforms
- `test_dispatch_has_minimal_overhead` - Verifies dispatch overhead is negligible
- `test_can_force_naive_implementation` - Verifies can explicitly use naive version
- `test_can_force_sse42_implementation` - Verifies can explicitly use SSE4.2 version
- `test_dispatched_implementation_logs_selection` - Verifies dispatch logging
- `test_multiversion_macro_properly_applied` - Verifies multiversion attribute present
- `test_cpuid_detection_works` - Verifies CPUID works correctly
- `test_fallback_on_unsupported_feature` - Verifies graceful fallback

These tests verify that the `multiversion` crate correctly performs runtime CPU feature
detection and dispatches to the appropriate implementation.

### 3. Correctness tests: SSE4.2 results match naive WFA on 20+ test cases

**Test Files**: `src/wfa_simd.rs`

**Tests** (20+ test cases comparing SIMD to naive):
1. `test_simd_matches_naive_exact` - Exact match
2. `test_simd_matches_naive_mismatch` - Single mismatch
3. `test_simd_matches_naive_insertion` - Single insertion
4. `test_simd_matches_naive_deletion` - Single deletion
5. `test_simd_matches_naive_long_sequence` - Long sequences (40bp)
6. `test_simd_matches_naive_high_divergence` - High divergence sequences
7. `test_simd_matches_naive_short_sequences` - Short sequences (4bp)
8. `test_simd_matches_naive_mid_seed` - Seed in middle of sequence
9. `test_simd_matches_naive_multiple_indels` - Multiple indels
10. `test_simd_matches_naive_consecutive_indels` - Consecutive indels
11. `test_simd_matches_naive_complex_pattern_1` - Complex pattern 1
12. `test_simd_matches_naive_complex_pattern_2` - Complex pattern 2
13. `test_simd_matches_naive_repeat_regions` - Repeat regions (AT repeats)
14. `test_simd_matches_naive_gc_rich` - GC-rich sequences
15. `test_simd_matches_naive_at_rich` - AT-rich sequences
16. `test_simd_matches_naive_edge_case_empty_prefix` - Edge case: start at beginning
17. `test_simd_matches_naive_edge_case_near_end` - Edge case: seed near end
18. `test_simd_matches_naive_random_sequence_1` - Random sequence 1
19. `test_simd_matches_naive_random_sequence_2` - Random sequence 2
20. `test_simd_matches_naive_random_sequence_3` - Random sequence 3

Each test compares the CIGAR string and alignment score from both implementations
to ensure they produce identical results.

### 4. Benchmark: ≥1.5× speedup on SSE4.2 CPU vs naive for 10kb alignment

**Test Files**: `benches/wfa_simd_bench.rs`

**Benchmarks**:
- `bench_10kb_alignment_naive` - Baseline naive implementation on 10kb
- `bench_10kb_alignment_simd` - SSE4.2 implementation on 10kb
- `bench_comparison` - Direct comparison of both on 10kb
- `bench_varying_sizes` - Comparison across sizes (100bp, 500bp, 1kb, 5kb, 10kb)
- `bench_varying_divergence` - Comparison across divergence levels (1%, 5%, 10%, 20%)

These benchmarks use criterion to measure performance. The `bench_comparison` benchmark
will directly show the speedup ratio. The target is ≥1.5× speedup on 10kb alignments.

### 5. Safety: document invariants for unsafe SIMD code

**Test Files**: `src/wfa_simd_safety.rs`

**Tests**:
- `test_safety_invariants_documented` - Verifies safety invariants are documented
- `test_simd_memory_safety_proof` - Verifies memory safety proofs exist
- `test_simd_alignment_requirements` - Verifies alignment requirements documented
- `test_simd_bounds_checking` - Verifies bounds checking is documented
- `test_simd_runtime_feature_detection` - Verifies feature detection documented
- `test_all_unsafe_blocks_documented` - Verifies all unsafe blocks have SAFETY comments
- `test_sse42_intrinsics_have_safety_comments` - Verifies intrinsics have safety docs
- `test_no_undefined_behavior_in_simd` - Verifies no UB scenarios
- `test_safety_docs_include_invariant_examples` - Verifies examples in docs
- `test_all_used_intrinsics_documented` - Verifies all intrinsics documented

These tests enforce that:
1. Module-level safety documentation exists
2. Every unsafe block has a SAFETY comment
3. All SSE4.2 intrinsics are documented with purpose and safety requirements
4. Memory safety invariants are explicit
5. No undefined behavior scenarios exist

### 6. Compiles and runs on non-x86 (falls back to naive automatically)

**Test Files**: `src/wfa_simd.rs`, `src/wfa_simd_dispatch.rs`

**Tests**:
- `test_compiles_and_runs_on_non_x86` - Verifies compilation on non-x86_64
- `test_arm64_fallback` - Explicitly tests ARM64 fallback
- `test_non_x86_64_uses_naive_only` - Verifies only naive on non-x86

These tests use `#[cfg(not(target_arch = "x86_64"))]` and `#[cfg(target_arch = "aarch64")]`
to verify the code compiles and runs on non-x86 platforms by falling back to the naive
implementation.

## Test Execution Status

**Expected Result**: ALL TESTS FAIL (RED phase)

This is Test-Driven Development. The tests are written first and should all fail because:
- The SSE4.2 implementation does not exist yet
- The multiversion dispatch does not exist yet
- The naive WFA baseline (from issue #5) does not exist yet
- Safety documentation is not complete yet

The implementation phase will make these tests pass (GREEN phase).

## Running Tests

### Unit Tests
```bash
cargo test --package phraya-align
```

### Benchmarks
```bash
cargo bench --package phraya-align --bench wfa_simd_bench
```

### Platform-Specific Tests
```bash
# On x86_64
cargo test --package phraya-align --lib -- test_x86_64

# On ARM64
cargo test --package phraya-align --lib -- test_arm64
```

## Dependencies Added

- `multiversion = "0.7"` - Runtime CPU feature dispatch
- `criterion = "0.5"` (dev) - Benchmarking framework

## Files Created

1. `src/wfa_simd.rs` - Core SSE4.2 implementation tests (27 tests)
2. `src/wfa_simd_safety.rs` - Safety documentation tests (10 tests)
3. `src/wfa_simd_dispatch.rs` - Runtime dispatch tests (11 tests)
4. `benches/wfa_simd_bench.rs` - Performance benchmarks (5 benchmarks)

Total: **48 unit tests + 5 benchmarks**

## Next Steps for Implementation

1. Implement naive WFA baseline (issue #5)
2. Implement SSE4.2 intrinsics for diagonal fill
3. Add multiversion attribute for runtime dispatch
4. Document all safety invariants
5. Optimize to achieve ≥1.5× speedup
6. Verify all tests pass (GREEN phase)
