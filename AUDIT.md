# Phraya Codebase Audit

**Date**: 2026-06-08  
**Method**: Code inspection + test execution. Spec document: CLAUDE.md Phase 1 MVP claims (shipped 2026-06-06).  
All conclusions derived from source only; documentation ignored.

## Executive Summary

Phase 1 MVP is largely implemented and functional. Cases 2, 3, 4 work end-to-end, filter presets operational, local coverage genuine, tandem repeats wired. **Critical finding**: SIMD claim overstated — WFA diagonal fill uses `wide` crate portable SIMD, not raw SSE4.2/NEON intrinsics as spec implies. Multiple stale RED/implementation comments remain. One timing test fails. Overall: genuine implementation with documentation debt.

## 1. What Is Genuinely Implemented

### 1.1 End-to-End Pipeline (Cases 2, 3, 4)
- **Case 2** (N reads + ref → N alignments): `integration_test_e2e_case2.rs` — 8/8 tests pass (SNP calling, VCF output, coverage filters)
- **Case 3** (M contigs + N reads, no ref → auto-centroid): `integration_test_e2e_case3.rs` — `issue_89_e2e_case3_pipeline_succeeds` passes (133s wall time)
- **Case 4** (M contigs ± ref): `integration_test_plan.rs::issue_68_plan_case4_contigs_only` passes

### 1.2 Real WFA Alignment
- WFA O(s·n) algorithm implemented at `phraya-align/src/wfa_simd.rs:95` (`fill_wfa_fitting`).
- Fitting alignment mode (query fully consumed, target end free) correctly used for read-vs-reference.
- 84/85 WFA tests pass. Timing test `wfa_is_faster_than_on2_for_sparse_edits_windowed` fails: 12.5ms vs <10ms target (150bp×300bp).

### 1.3 Local Coverage (±50bp Window, Not Stubbed)
- Real implementation at `phraya-align/src/executor.rs:401` (`compute_raw_coverage`): increments per-position from alignment start..end.
- `local_coverage_tests.rs` confirms NOT placeholder `vec![1]` (issue #130 regression guard).
- Coverage quantized to nearest 5, RLE-compressed in `.phraya` format.

### 1.4 Variation Hotspot Detection
- Implemented at `phraya-core/src/types.rs:1278` (`detect_hotspot_intervals`).
- All 8 hotspot tests pass (`issue_145_*`).
- **Stale comment** at `phraya-core/src/hotspot_tests.rs:1`: "All tests call `detect_hotspot_intervals` which is `unimplemented!()` in production code" — FALSE, impl exists.

### 1.5 Tandem Repeat Detection
- Real impl at `phraya-core/src/lib.rs:79` (`detect_tandem_repeats`): period 2–4, min 3 repeats.
- Wired end-to-end: `phraya-align/src/executor.rs:168` calls it, variants annotated with `in_tandem_repeat` flag.
- Test `executor::tests::tandem_repeat_variants_are_annotated` passes.

### 1.6 Filter Presets & Library API
- **Conservative** (min_cov=10, mapq=30, allele_freq≥10%, exclude_tandem_repeats=true) at `phraya-filter/src/lib.rs:27`.
- **Sensitive** (min_cov=3, mapq=20, allele_freq≥2%) at `lib.rs:32`.
- CLI integration: `phraya-cli/src/main.rs:542` (`--preset conservative|sensitive`).
- 3/3 preset tests pass (`conservative_*`).
- Library API (`FilterBuilder`, `ThresholdFilter`) exposed in `phraya-filter` crate.

### 1.7 K-mer Sketching & Reuse
- `.phrayaplan` v2 (MessagePack+zstd) stores `kmer_index: HashMap<String, MinimizerSketch>` (per-sequence ID).
- Version check at `phraya-io/src/plan.rs:118` rejects v1 files.
- Alignment reuses via `plan.get_sketch(query.id())` at `phraya-align/src/executor.rs:100`.
- Sketching via `simd-minimizers` crate (k=21, w=11).

### 1.8 BAM/CRAM Input (Pure Rust, No htslib)
- `noodles_bam` / `noodles_cram` at `phraya-io/src/bam_cram.rs:37,104`.
- mapq extracted at `bam_cram.rs:68`, quality scores at `:60`.
- 25/25 `test_issue_61_*` tests pass (BAM/CRAM parsing, indexed files, mapq, quality scores).
- **Stale RED comments** at `bam_cram.rs:347,362,374`: "Feature not yet implemented" — FALSE, tests pass.

### 1.9 Confidence / MAPQ / Avg Base Quality
- Confidence = `1 - edit_dist/query_len` (alignment-derived), computed at `executor.rs:159`.
- MAPQ / avg_base_quality extracted from BAM records (`bam_cram.rs:68,60`) or FASTQ Phred scores.
- Propagated into `VariantObservation` at `executor.rs:266`.

### 1.10 VCF/TSV Output
- VCF writer at `phraya-filter/src/vcf.rs`, TSV at `phraya-filter/src/tsv.rs`.
- E2E tests confirm VCF contains correct REF/ALT alleles (`issue_88_vcf_correct_ref_alt_alleles`).

## 2. Critical Gaps

### 2.1 SIMD Intrinsics Claim Overstated

**Claim** (CLAUDE.md):  
> "Real WFA O(s·n) alignment with SIMD diagonal fill (SSE4.2 / NEON)"  
> "Platform-optimized SIMD (SSE4.2/AVX2/NEON via simd-minimizers)"  
> "WFA diagonal fill uses SSE4.2/NEON intrinsics"

**Reality**:  
WFA diagonal fill at `phraya-align/src/wfa_simd.rs:674` (`fill_simd`) uses `wide::i32x8` (portable SIMD), NOT raw `_mm_*` / `vmovq` intrinsic calls.

```rust
// phraya-align/src/wfa_simd.rs:674
fn fill_simd(q: &[u8], t: &[u8]) -> DiagMatrix {
    use wide::i32x8;  // ← portable SIMD, not intrinsics
    // ...
    let m = (up + one).min(left + one).min(diag + cost);
    // ↑ wide's min() may lower to SSE4.2 pminsd IF compiled with -C target-cpu=native,
    // but code itself is arch-neutral
}
```

- Zero occurrences of `_mm_`, `vld1`, `vmovq`, or `#[cfg(target_feature)]` unsafe blocks in `wfa_simd.rs`.
- Test at `wfa_simd.rs:1064` confirms: "portable-SIMD kernel currently uses no `unsafe` at all".
- `wide` crate (`safe_arch` backend) CAN lower to SSE4.2 / NEON, but only if build flags enable it — the code is portable, not intrinsic-based.

**Consequence**:  
Spec's repeated claims of "SSE4.2 intrinsics" mislead. The implementation is `wide`-based portable SIMD with *potential* platform lowering, not a hand-written intrinsic kernel. On a CPU without SSE4.2 or when not compiled with `-C target-cpu=native`, `wide` falls back to scalar emulation. The "SIMD" is abstract, not guaranteed.

### 2.2 Stale "RED" and "unimplemented!()" Comments

**Finding**: Multiple comments claim features are stubs or RED (test-driven development initial failure state), but tests pass and implementations exist.

| File | Line | Comment | Reality |
|------|------|---------|---------|
| `phraya-core/src/hotspot_tests.rs` | 1 | "All tests call `detect_hotspot_intervals` which is `unimplemented!()` in production code" | Impl exists at `types.rs:1278`, 8/8 tests pass |
| `phraya-io/src/bam_cram.rs` | 347 | "RED: Feature not yet implemented" (auto-detect .bam) | Test passes, feature works |
| `phraya-io/src/bam_cram.rs` | 362 | "RED: Feature not yet implemented" (auto-detect .cram) | Test passes, feature works |
| `phraya-io/src/bam_cram.rs` | 374 | "RED: Feature not yet implemented" (.bam.gz) | Test passes, feature works |
| `phraya-cli/tests/issue_58_mvp_acceptance_tests.rs` | 5 | "RED phase: All tests fail initially (feature not yet implemented)" | 28/28 tests pass |

**Consequence**:  
A code reader sees "unimplemented!" comments and assumes the feature is missing. In fact, the implementation is present and passing tests. Comments are stale artifacts from TDD red-green-refactor cycle; they were never updated after GREEN phase.

## 3. Missing Features

All Phase 1 MVP claims are met. The following are correctly marked as deferred:

| Feature | Spec Status | Code Status |
|---------|-------------|-------------|
| Expression-based filters (`--expr`) | Phase 2+ | Library impl exists (`phraya-filter/src/lib.rs:616` `ExprFilter`), CLI flag missing |
| CRAM with external reference | Phase 2+ | Not implemented (comment at `bam_cram.rs:90`) |
| Variation hotspot estimation in plan | Phase 2+ | Detection impl exists, plan-time estimation not wired |
| Two-tier evidence (k-mer → alignment refinement) | Phase 2+ | Not implemented |
| Python/R bindings | Phase 2+ | Not implemented |

## 4. Test Suite Integrity Summary

**Total**: 84 passed (phraya-align lib), 8 passed (case 2 E2E), 1 passed (case 3 E2E), 2 passed (case 4), 8 passed (hotspots), 3 passed (filter presets), 25 passed (BAM/CRAM), 28 passed (MVP acceptance).  
**1 failure**: `wfa_is_faster_than_on2_for_sparse_edits_windowed` (timing regression: 12.5ms vs <10ms target).  
**0 ignored** (except `wfa_perf.rs` microbenchmarks, correctly ignored in debug build).

### What Tests Prove
- WFA correctness: 84 tests compare SIMD vs naive DP; edit distances and CIGARs match reference.
- End-to-end pipelines: real `phraya plan → align → merge → filter` workflows succeed.
- Local coverage: NOT placeholder `vec![1]`; real ±50bp window from alignment.
- Tandem repeats: variants correctly annotated, excluded by conservative preset.
- Filter presets: threshold combinations (conservative/sensitive) tested.

### What Tests Don't Prove
- **SIMD performance**: only one timing test, and it FAILS. No benchmark comparing `wide` SIMD vs scalar across sequence sizes.
- **Platform dispatch correctness**: no test confirms SSE4.2 codepath actually runs on x86_64 or NEON on aarch64; tests only check dispatch selection string matches arch.
- **Multi-mapping score threshold justification**: 0.95 score_ratio threshold hard-coded; no test validates it's the right value for bacterial genomics.

---

## Appendix A: Bad Tests

### A.1 Timing Test Failure

**Test**: `phraya-align/src/wfa_simd.rs:2456` (`wfa_is_faster_than_on2_for_sparse_edits_windowed`)  
**Status**: FAILS  
**Why**: `fill_wfa` takes 12.5ms for 150bp×300bp; test requires <10ms.  
**Root cause**: Either (1) WFA not fast enough, or (2) timing threshold too strict for debug build, or (3) portable SIMD (`wide`) not lowering to real intrinsics on this CPU.

### A.2 Vacuous Tests

None detected. All passing tests check real behavior against expected values or non-trivial assertions.

### A.3 Stale "RED" Comments

See §2.2. Tests marked "RED: Feature not yet implemented" pass, implying feature IS implemented. Comment maintenance debt, not test failure.

---

## Appendix B: Dead Code

**Functions never called in prod**:
- `phraya-align/src/wfa_simd.rs:533` `fill_scalar`: O(n×m) DP reference. Used only in tests for correctness comparison.
- `wfa_simd.rs:598` `traceback`: Row-major DP traceback. Used only in tests via `fill_scalar`.

**Verdict**: Test helpers, not prod dead code. Compiler warns "never used" because no prod path calls them. Normal for test-only utilities not marked `#[cfg(test)]`. Not a bug.

---

## Appendix C: Implementation Notes for Maintainers

### C.1 SIMD Reality Check
If performance matters, verify `wide` is actually using SIMD:
```bash
RUSTFLAGS="-C target-cpu=native" cargo build --release
objdump -d target/release/phraya-align | grep -i pminsd  # SSE4.2
objdump -d target/release/phraya-align | grep -i smin    # NEON
```
If no SIMD instructions appear, `wide` is emulating. Consider raw intrinsics or `std::simd` (nightly).

### C.2 Stale Comment Cleanup
Grep for `RED|unimplemented!|not yet|Phase 2` in comments and verify each against current code state. Update or remove.

### C.3 Timing Test Fix
Either:
1. Relax threshold to `<15ms` (matches observed 12.5ms in debug).
2. Mark `#[ignore]` for debug builds, run only in release (`#[cfg(not(debug_assertions))]`).
3. Optimize WFA fitting mode for short reads (current impl may do unnecessary work).

### C.4 Expression Filter CLI Wiring
`ExprFilter` library exists. Add CLI flag:
```rust
// phraya-cli/src/main.rs, Filter subcommand:
#[arg(long)]
expr: Option<String>,
```
Then call `ExprFilter::new(expr.as_ref().unwrap())?` when `--expr` provided.
