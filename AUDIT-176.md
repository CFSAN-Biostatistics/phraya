# Phraya Codebase Audit

**Date**: 2026-06-29
**Method**: Code inspection + test execution. Spec document: `CLAUDE.md` (project instructions).
All conclusions derived from source only; documentation ignored unless corroborated by code.
**Auditor stance**: assume implementation is incomplete until code proves otherwise.

---

## Executive Summary

Phraya's Phase 1 MVP claims are genuinely and correctly implemented. The core algorithms (WFA, Myers, SIMD match extension, RLE coverage track, strategy ladder) are real — not stubs — and backed by thorough differential testing. The single most important finding is a coverage gap in the paired-end filter pipeline: the BAM mate-info extraction path has no non-vacuous integration test. The three "passing" tests in `test_bam_mate_extraction.rs` never call `BamCramParser::from_bam_path()` at all; the only test that would (`bam_parser_extracts_mate_info`) is explicitly ignored as a placeholder. The entire insert-size / discordant-pair filter stack depends on this extraction being correct under real BAM input. Expression-based filters (`--expr`) are fully implemented in the library but absent from the CLI — consistent with the Phase 2+ label in the spec, but the gap between library and binary deserves explicit tracking.

---

## 1. What Is Genuinely Implemented

### 1.1 WFA O(s·n) Alignment

`wfa_simd.rs:220` — `fill_wfa` is a real wavefront algorithm. It maintains per-diagonal frontiers, advances via `count_matching_prefix` SIMD extension, and backtracks through wavefront history to emit CIGAR. Performance guard at `wfa_simd.rs:2572` asserts 150bp vs. 300bp query completes in <10ms, ruling out O(n×m) behavior. **Genuine.**

### 1.2 Myers Bit-Parallel O(nm/w)

`wfa_simd.rs:2621-2880` — `myers_forward` builds per-column bitvectors (Peq, Pv, Mv), `myers_backtrace` traces CIGAR, `myers_fitting_impl` enforces fitting semantics. Differential suite at `issue_148_strategy_tests.rs` proves Myers ≡ WFA edit distance on 100 random pairs. **Genuine.**

### 1.3 SIMD Match Extension (`count_matching_prefix`)

`wfa_simd.rs:52-89` (x86_64 SSE2), `wfa_simd.rs:96-117` (aarch64 NEON). Both branches contain real intrinsic calls: `_mm_loadu_si128`, `_mm_cmpeq_epi8`, `_mm_movemask_epi8` (x86_64); `vld1q_u8`, `vst1q_u8` (aarch64). u64-XOR fallback for other platforms. **Genuine SIMD — not a wide-crate emulation wrapper.**

### 1.4 Strategy Ladder (exact / balanced / fast)

`executor.rs:10-107` — `Strategy` enum, `AlignConfig::new` sets coverage-window radius per strategy (25/50/150bp), `FAST_MAX_DIVERGENCE = 0.20` at line 107. `build_anchors` subsamples seeds for `Fast` (best-voted target-start) vs. all anchors for `Balanced`/`Exact`. `Exact` forces WFA; `Balanced`/`Fast` choose Myers ≤500bp, WFA fallback. **Genuine.**

### 1.5 Coverage Track (quantized + RLE)

`types.rs:433-521` — `CoverageTrack::new` quantizes to nearest 5 (explicit `quantize` function), then RLE-compresses into `Vec<(u8, u32)>`. 7 unit tests including `compression_ratio` prove the track compresses 10,000-element uniform arrays to 1 run. **Genuine.**

### 1.6 Filter Presets (conservative / sensitive)

`phraya-filter/src/lib.rs:18-36` — `FilterPreset::Conservative` and `FilterPreset::Sensitive` return seeded `FilterBuilder` instances with distinct thresholds. Integration tests `conservative_preset_rejects_low_quality` and `sensitive_preset_passes_low_coverage_variant` verify threshold semantics with real observations. **Genuine.**

### 1.7 VCF and TSV Output

`phraya-filter/src/vcf.rs` and `phraya-filter/src/tsv.rs` each contain real formatters with 7+ unit tests covering position encoding, multi-allelic sites, header format, and 1-indexed output. **Genuine.**

### 1.8 Multi-mapping (score_ratio ≥ 0.95)

`phraya-align/src/lib.rs:166-174` — alternatives filtered by normalized edit-distance ratio against the primary. Three unit tests at lines 217, 242, 274 verify inclusion/exclusion at the boundary. **Genuine.**

### 1.9 Paired-End Filters (merge-stable)

`phraya-core/src/types.rs` — `insert_size_sum`, `insert_size_count`, `unmapped_mate_count` aggregate fields stamped per-variant at alignment time, summed during merge. Four filter modes (`require_proper_pairs`, `exclude_discordant_pairs`, `min/max_insert_size`, `require_both_mates_mapped`) operate on aggregates. 12 unit tests in `test_paired_end_filters.rs`. **Genuine — with coverage caveat noted below.**

### 1.10 Variation Hotspot Estimation

`detect_hotspot_intervals` (`phraya-core/src/lib.rs:1368`) called from `run_plan` at `main.rs:712` and stored in `plan.hotspot_intervals`. 5 integration tests in `issue_145_hotspot_intervals_test.rs` covering interval detection, merging, and empty-input edge cases. **Genuine.**

### 1.11 Expression-Based Filter (Library Only)

`phraya-filter/src/lib.rs:848-1499` — `ExpressionFilter` parses boolean expressions into an AST (`Expr::Or/And/Not`), evaluates against all `VariantObservation` fields. 20 unit tests under `issue_150_*` cover operators, parentheses, unknown fields, and complex compositions. **Library is genuinely implemented** — see §2.1 for the CLI gap.

---

## 2. Critical Gaps

### 2.1 `--expr` Absent from CLI

**Claim (spec)**: `--expr "coverage >= 10 && mapq > 30"` [Phase 2+]

**Reality**: `phraya-cli/src/main.rs` — the `Filter { ... }` struct (lines 112–173) has no `expr` field. Running `phraya filter --help` produces no `--expr` flag. `run_filter` at line 1119 does not instantiate `ExpressionFilter`. The library is complete; the binary surface is not wired up.

**Consequence**: Users cannot use expression filters from the CLI. The spec explicitly marks this Phase 2+, so the gap is documented but it means the library is ahead of the binary by one integration point.

### 2.2 BAM Mate-Info Extraction Has No Integration Test

**Claim (spec)**: "BAM/CRAM input via `noodles` (pure Rust, no htslib)" + paired-end filters that depend on TLEN and flag extraction from BAM records.

**Reality**: `test_bam_mate_extraction.rs` has 4 tests:
- `bam_parser_extracts_mate_info` — **ignored** ("placeholder for future BAM file creation")
- `parsed_reads_structure` — constructs a `ParsedReads` struct by hand, never calls `BamCramParser::from_bam_path()`
- `mate_id_toggling` — duplicates ID-toggling logic in-test, never calls the parser
- `sam_flags_interpretation` — tests SAM bitfield arithmetic inline, never calls noodles

The actual extraction code in `bam_cram.rs:51-130` reads `flags.is_properly_segmented()`, `record.template_length()`, and builds `MateInfo` — but this path is never exercised by a non-trivial test. All four paired-end filters (`require_proper_pairs`, `exclude_discordant_pairs`, `min/max_insert_size`, `require_both_mates_mapped`) depend on this extraction to produce non-zero aggregate counts at alignment time. If extraction silently fails or produces wrong values (e.g., `TLEN=0` for all records), the filters would silently pass everything.

**Consequence**: The paired-end filter pipeline has an untested seam between BAM parsing and filter application. A regression in `bam_cram.rs` flag/TLEN extraction would not be caught by the test suite.

---

## 3. Missing Features (Phase 2+ Deferred — Expected Absent)

| Feature | Spec Section | Status | Notes |
|---------|-------------|--------|-------|
| Expression filter CLI (`--expr`) | Filter Operations | Library done, CLI not wired | `ExpressionFilter` exists with 20 tests; no `--arg` in `Filter` subcommand |
| CRAM + external reference | Phase 2 deferred | Not implemented | `bam_cram.rs:137`: "reference-compressed mapped reads require external reference (not yet supported)" |
| Two-tier evidence (k-mer → alignment refinement) | Phase 2 deferred | Not implemented | Plan produces k-mer evidence; no refinement loop |
| Python bindings / R integration | Phase 2 deferred | Not implemented | Expected absent |
| Per-base coverage breakdown | Phase 2 deferred | Not implemented | Only ±N bp window; no per-base |

All items above are explicitly labeled Phase 2+ in the spec. Their absence is correct.

---

## 4. Test Suite Integrity Summary

**Total**: 625 passing, 0 failing, 3 ignored, wall time ~80s (non-trivial; rules out all-stub behavior).

| Category | Count |
|----------|-------|
| Passing | 625 |
| Failing | 0 |
| Ignored | 3 |
| Vacuous (see §A.2) | 5 |

The 80-second wall time includes end-to-end integration tests (plan → align → filter pipeline) in `issue_58_mvp_acceptance_tests.rs` and `integration_test_e2e_case2/3.rs`. These tests are genuine: they write real FASTA/BAM files, invoke the full pipeline, and assert on `.phraya` file contents.

---

## Appendix A: Bad Tests

### A.1 Unconditionally Ignored Tests

| Test | File | Reason | Does the feature work? |
|------|------|--------|----------------------|
| `bam_parser_extracts_mate_info` | `phraya-io/tests/test_bam_mate_extraction.rs:8` | "Placeholder for future BAM file creation, which is complex" | Code exists in `bam_cram.rs`, but **untested end-to-end** |
| `detect_tandem_repeats` (doctest) | `phraya-core/src/lib.rs:72` | Doctest marked `ignore` (likely scoping issue) | Function is real and tested by `phraya-cli` integration tests |
| `wfa_perf_myers_vs_wfa_speed_comparison` | `phraya-align/tests/wfa_perf.rs:37` | "Release-only microbenchmark" — meaningless in debug | Correctly documented; not a coverage gap |

**The `bam_parser_extracts_mate_info` ignored test is the only one that represents a genuine coverage gap.**

### A.2 Vacuous Tests

| Test | File | Why Vacuous |
|------|------|-------------|
| `every_unsafe_in_impl_has_safety_comment` | `wfa_simd.rs:1202` | The comment in the test itself says it "passes vacuously" when there are no `unsafe` blocks. It's a guard for the future, not evidence of current correctness. |
| `test_sse42_feature_detection` | `wfa_simd_dispatch.rs` | Asserts `is_sse42_available() == is_sse42_available()` — determinism check, not a capability test. |
| `parsed_reads_structure` | `test_bam_mate_extraction.rs` | Constructs `ParsedReads` by hand and checks `.len()`. Never invokes `BamCramParser`. |
| `sam_flags_interpretation` | `test_bam_mate_extraction.rs` | Tests SAM bitfield arithmetic copied inline; never calls `bam_cram.rs`. |
| `mate_id_toggling` | `test_bam_mate_extraction.rs` | Duplicates the ID-toggling logic in-test; never calls `bam_cram.rs`. |

### A.3 Tests Named for Feature X That Test Feature Y

| Test | Claimed Feature | Actually Tests |
|------|----------------|---------------|
| `sam_flags_interpretation` | BAM flag parsing (implied by file name `test_bam_mate_extraction`) | SAM bitfield arithmetic in isolation — no noodles calls |
| `mate_id_toggling` | BAM mate ID extraction | String manipulation logic inline — no noodles calls |

### A.4 Stale "Currently Fails" Comments

None found. No test has a "CURRENTLY FAILS" comment that passes.

---

## Appendix B: Non-Implemented Features

| Feature | Spec Location | Status | Notes |
|---------|-------------|--------|-------|
| `--expr` CLI flag | "Expression-based: `--expr` [Phase 2+]" | Library done; CLI surface missing | Consistent with Phase 2+ label |
| CRAM + external reference | Phase 1 "BAM/CRAM input" note | Partial | Reference-free CRAM works; external reference not supported |
| Expression filter CLI wire-up | Filter Operations §2 | Not done | `run_filter` in `main.rs` does not call `ExpressionFilter` |
