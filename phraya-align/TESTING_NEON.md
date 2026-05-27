# Testing NEON SIMD WFA Implementation

## Overview

Issue #14 requires NEON SIMD optimizations for WFA diagonal fill on ARM64 platforms. This document describes how to verify the tests are in RED phase (failing) and how they should be tested on ARM64 hardware.

## Test Status: RED Phase

All NEON tests are written to **FAIL** until the implementation is complete. The tests currently call `unimplemented!()` which will panic when executed on ARM64 hardware.

## Platform-Specific Testing

### On x86_64 (Development Machine)

The NEON tests are conditionally compiled **only** on `target_arch = "aarch64"`. On x86_64:

```bash
cargo test --package phraya-align
```

**Expected behavior:** Only the `test_neon_noop_on_non_arm64` test runs and passes. All ARM64-specific tests are not compiled.

### On ARM64 (Graviton, Apple Silicon, etc.)

To verify RED phase on ARM64 hardware:

```bash
cargo test --package phraya-align --lib
```

**Expected behavior:** All 22 NEON tests should **FAIL** with:

```
thread 'wfa_neon::tests::test_neon_exact_match' panicked at phraya-align/src/wfa_neon.rs:37:5:
not yet implemented: NEON WFA diagonal fill not yet implemented
```

## Test Coverage

The test suite includes 22 comprehensive tests covering:

### Correctness Tests (Tests 1-15)
1. Exact match sequences
2. Single mismatch
3. All mismatches (stress test)
4. Long sequence (10kb) - benchmark target
5. Short sequence (boundary case)
6. Aligned to NEON vector width (16 bytes)
7. Unaligned to NEON vector width (17 bytes)
8. Empty sequences (edge case)
9. Non-ACGT characters
10. Negative scores in previous wavefront
11. High-scoring previous wavefront
12. GC-rich sequences (bacterial genomics)
13. AT-rich sequences
14. Homopolymer runs
15. Mixed events (multiple mismatches)

### Cross-Architecture Test (Test 16)
16. Verify no-op on non-ARM64 platforms (x86_64)

### NEON vs Naive Correctness (Tests 17-21)
17-21. Five test cases comparing NEON output to naive WFA implementation
    - Requires issue #5 (naive WFA) to be implemented first
    - Currently have placeholder assertions that will be enhanced

### Performance Hint (Test 22)
22. Documents ≥1.5× speedup expectation (real benchmarks in criterion)

## Benchmarks

Performance benchmarks are in `benches/wfa_neon_benchmark.rs` using criterion:

```bash
# On ARM64 hardware only
cargo bench --package phraya-align
```

**Acceptance criteria:** NEON implementation must achieve ≥1.5× speedup vs naive on 10kb sequences.

## Implementation Requirements

When implementing the NEON optimizations (issue #14), the implementation agent must:

1. Replace the `unimplemented!()` call in `wfa_diagonal_fill_neon()`
2. Use `core::arch::aarch64` intrinsics
3. Document safety invariants for all unsafe blocks
4. Ensure all 22 tests PASS
5. Achieve ≥1.5× speedup in benchmarks

## Verifying GREEN Phase

After implementation:

```bash
# On ARM64 hardware
cargo test --package phraya-align --lib
# All tests should PASS

cargo bench --package phraya-align
# Should show ≥1.5× speedup for NEON vs naive on 10kb sequences
```

## Cross-Compilation Testing

To verify compilation (not execution) for ARM64 from x86_64:

```bash
rustup target add aarch64-unknown-linux-gnu
cargo check --package phraya-align --target aarch64-unknown-linux-gnu
```

This ensures the NEON code compiles correctly even when developing on x86_64.

## CI/CD Considerations

The CI pipeline should include:

1. **x86_64 runner:** Verify non-ARM64 tests pass, NEON tests not compiled
2. **ARM64 runner (Graviton/Apple Silicon):** Run full test suite + benchmarks
3. **Cross-compilation check:** Verify ARM64 target compiles on x86_64

## Dependencies

- Issue #5: Naive WFA implementation must be complete for correctness comparison tests
- ARM64 hardware access for full test verification
