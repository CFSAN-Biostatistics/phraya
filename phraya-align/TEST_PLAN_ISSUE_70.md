# Test Plan for Issue #70: WFA extension (scalar)

## Overview

This test plan validates the acceptance criteria for implementing a scalar (non-SIMD) WFA
(Wavefront Alignment) extension algorithm in phraya-align. This establishes the correctness
baseline before SIMD optimization.

## Current Status

**RED PHASE (Expected)**: All tests fail with compilation errors because the API specified
in the acceptance criteria does not match the current implementation.

## Blocking Dependency

Issue #59 (error types) - **RESOLVED**. Error types exist in `phraya-align/src/lib.rs`.

## Acceptance Criteria Coverage

### 1. Function Signature: `wfa_extend(query: &[u8], target: &[u8], seed_pos: usize) → Alignment`

**Test File**: `tests/test_wfa_scalar_acceptance.rs`

**Tests**:
- `test_wfa_extend_signature_accepts_seed_pos` - Verifies function accepts usize seed_pos
- `test_seed_pos_at_start` - Verifies seed_pos=0 behavior
- `test_seed_pos_in_middle` - Verifies seed_pos in middle of sequence
- `test_seed_pos_near_end` - Verifies seed_pos near end of sequence

**Current Failure**: Function expects `SeedAnchor` struct, not `usize`.

### 2. Alignment Struct Fields

**Test File**: `tests/test_wfa_scalar_acceptance.rs`

**Required Fields**:
- `cigar: String` - CIGAR string representation
- `edit_distance: usize` - Count of edits (mismatches + insertions + deletions)
- `query_start: usize` - Starting position in query
- `query_end: usize` - Ending position in query
- `target_start: usize` - Starting position in target
- `target_end: usize` - Ending position in target

**Tests**:
- `test_alignment_has_edit_distance_field` - Verifies edit_distance field exists
- `test_alignment_has_query_positions` - Verifies query_start/query_end exist
- `test_alignment_has_target_positions` - Verifies target_start/target_end exist

**Current Failure**: Alignment struct has `score: i32` instead of `edit_distance: usize`,
and lacks position fields.

### 3. CIGAR Format: M (match/mismatch), I (insertion), D (deletion)

**Test File**: `tests/test_wfa_scalar_acceptance.rs`

**Tests**:
- `test_cigar_format_valid` - Verifies only M, I, D (and optionally =, X) in CIGAR
- `test_cigar_length_consistency` - Verifies CIGAR accounts for all positions

### 4. Edit Distance: count of mismatches + insertions + deletions

**Test File**: `tests/test_wfa_scalar_acceptance.rs`

**Tests**:
- `test_edit_distance_is_non_negative` - Verifies edit_distance >= 0
- `test_edit_distance_matches_cigar_operations` - Verifies consistency with CIGAR
- `test_edit_distance_bounded_by_sequence_length` - Verifies <= max(query_len, target_len)

### 5. Exact match (100bp) → CIGAR "100M", edit_dist 0

**Test File**: `tests/test_wfa_scalar_acceptance.rs`

**Tests**:
- `test_exact_match_100bp` - 100bp perfect match
- `test_exact_match_produces_zero_edit_distance` - Verifies 0 edit distance for exact matches

### 6. Single mismatch at pos 50 → CIGAR includes mismatch operation

**Test File**: `tests/test_wfa_scalar_acceptance.rs`

**Tests**:
- `test_single_mismatch_at_position_50` - Mismatch at position 50 in 100bp sequence
- `test_single_mismatch_early_position` - Mismatch near start
- `test_single_mismatch_late_position` - Mismatch near end

### 7. Single insertion → CIGAR includes 'I'

**Test File**: `tests/test_wfa_scalar_acceptance.rs`

**Tests**:
- `test_single_insertion` - Single insertion in middle
- `test_insertion_at_start` - Insertion at start
- `test_insertion_at_end` - Insertion at end
- `test_multiple_insertions` - Multiple insertions

### 8. Single deletion → CIGAR includes 'D'

**Test File**: `tests/test_wfa_scalar_acceptance.rs`

**Tests**:
- `test_single_deletion` - Single deletion in middle
- `test_deletion_at_start` - Deletion at start
- `test_deletion_at_end` - Deletion at end
- `test_multiple_deletions` - Multiple deletions

### 9. Complex alignment (mix of M/I/D) → verify CIGAR + edit_dist correct

**Test File**: `tests/test_wfa_scalar_acceptance.rs`

**Tests**:
- `test_complex_alignment_mixed_operations` - Mix of operations
- `test_complex_alignment_consecutive_indels` - Multiple consecutive indels
- `test_complex_alignment_alternating_events` - Alternating events
- `test_complex_alignment_high_divergence` - Highly divergent sequences

### 10. Empty sequences → handle gracefully

**Test File**: `tests/test_wfa_scalar_acceptance.rs`

**Tests**:
- `test_empty_query_sequence` - Empty query, non-empty target
- `test_empty_target_sequence` - Non-empty query, empty target
- `test_both_sequences_empty` - Both empty

### 11. Benchmark: align 10kb sequences (measure time)

**Test File**: `benches/bench_wfa_scalar.rs`

**Benchmarks**:
- `bench_10kb_exact_match` - 10kb exact match
- `bench_10kb_with_1pct_divergence` - 10kb with 1% divergence
- `bench_10kb_with_5pct_divergence` - 10kb with 5% divergence
- `bench_10kb_with_10pct_divergence` - 10kb with 10% divergence
- `bench_varying_lengths` - 100bp to 10kb
- `bench_varying_divergence` - 0% to 20% divergence
- `bench_varying_seed_positions` - Different seed positions
- `bench_short_sequences` - 100bp baseline
- `bench_empty_sequences` - Empty sequences

**Additional Unit Test**:
- `test_10kb_alignment_completes` - Ensures 10kb alignment completes without panic

## Additional Edge Cases Tested

**Test File**: `tests/test_wfa_scalar_acceptance.rs`

- `test_single_base_sequences` - Single nucleotide sequences
- `test_single_base_mismatch` - Single base mismatch
- `test_very_long_indel` - 50bp indel
- `test_sequences_with_ambiguous_bases` - IUPAC codes handling

## Test Execution

### Run All Unit Tests
```bash
cargo test --package phraya-align --test test_wfa_scalar_acceptance
```

**Expected**: 56 tests fail with compilation errors (RED phase)

### Run Benchmarks
```bash
cargo bench --package phraya-align --bench bench_wfa_scalar
```

**Expected**: Benchmarks fail to compile (RED phase)

## Summary

**Total Coverage**:
- 56 unit tests covering all acceptance criteria
- 9 criterion benchmarks for performance measurement
- All edge cases covered (empty, single base, long indels, etc.)

**Test Status**: RED PHASE ✓
- All tests fail as expected
- Tests define the correct API from issue #70 acceptance criteria
- Implementation phase will make tests GREEN

## Next Steps for Implementation Agent

1. Update `Alignment` struct in `phraya-align/src/lib.rs`:
   - Change `score: i32` to `edit_distance: usize`
   - Add `query_start: usize`, `query_end: usize`
   - Add `target_start: usize`, `target_end: usize`

2. Update `wfa_extend` function signature:
   - Change `seed: SeedAnchor` parameter to `seed_pos: usize`
   - Change return type from `WfaResult` to `Alignment` (direct return, not Result)

3. Implement scalar WFA algorithm:
   - Use simple dynamic programming approach
   - No SIMD optimizations (baseline correctness)
   - Generate CIGAR string with M/I/D operations
   - Calculate edit distance from alignment

4. Handle edge cases:
   - Empty sequences (return appropriate alignment)
   - Single base sequences
   - Very long indels

5. Verify all 56 tests pass (GREEN phase)

6. Run benchmarks to establish baseline performance for 10kb alignments

## Files Created

1. `tests/test_wfa_scalar_acceptance.rs` - 56 acceptance tests
2. `benches/bench_wfa_scalar.rs` - 9 performance benchmarks
3. `TEST_PLAN_ISSUE_70.md` - This test plan (documentation)

## Notes

This test plan follows Test-Driven Development (TDD) principles:
- **RED phase**: Tests written first, all fail (current state)
- **GREEN phase**: Implementation makes tests pass (next agent)
- **REFACTOR phase**: Optimize without changing behavior (future)

The tests are immutable contracts - implementation agent cannot modify them.
