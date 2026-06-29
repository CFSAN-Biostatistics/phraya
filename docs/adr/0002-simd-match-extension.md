# 2. SIMD-accelerated match extension primitive

- **Status**: Accepted
- **Date**: 2026-06-19

## Context

The hot path of both the WFA and Myers aligners is the *match extension* step: advancing
along a diagonal while query and target bytes are equal. Profiling and the algorithm's
structure both point here — it runs `O(s · match_run_length)` times per alignment, one
byte comparison per iteration in the original scalar loop.

An earlier experiment ("SIMD diagonal fill") implemented a full `O(n·m)` anti-diagonal DP
with the `wide` crate. It was ~28× slower than `O(s·n)` WFA and is retained only as a
test oracle. The lesson: SIMD belongs in the *inner primitive*, not in replacing the
sub-quadratic algorithm.

Two complications shaped the design:

1. **Endianness** — reading bytes as a word and using `trailing_zeros` to locate the
   first mismatch is only correct if byte index 0 is the least-significant byte.
2. **Portability vs. determinism** — `multiversion`-style runtime feature dispatch adds
   per-call overhead and non-determinism in which code path runs; but SSE2 is mandatory
   on every x86-64 CPU and NEON is mandatory on every AArch64 CPU.

## Decision

Introduce `count_matching_prefix(a, b) -> usize`, the length of the longest common prefix,
as the single match-extension primitive, used by both WFA `extend` closures and Myers.

Three tiers, all required to return bit-identical results:

1. `count_matching_prefix_scalar` — byte-by-byte reference; the semantic ground truth.
2. `count_matching_prefix_u64` — 8 bytes/step via **little-endian** word XOR
   (`from_le_bytes`), so `trailing_zeros()/8` is the first mismatch index on any host.
   No `unsafe`. Portable fallback for non-x86-64/AArch64 targets.
3. Architecture SIMD — 16 bytes/step using SSE2 (`_mm_cmpeq_epi8` + `_mm_movemask_epi8`)
   on x86-64 and NEON (`vceqq_u8` + `vminvq_u8`) on AArch64.

Tier 3 is selected at **compile time** via `cfg(target_arch)`, not runtime dispatch,
because SSE2/NEON are architecture baselines and are therefore always available.

Equivalence is enforced by a differential test that runs all tiers (and the dispatched
entry point) over a battery of inputs straddling the 8- and 16-byte boundaries, a
mismatch at every position, and empty/unequal-length slices. The whole-aligner WFA test
suite is the behavioral backstop.

## Consequences

- Measured **5.4×–8.9×** speedup over scalar on long match runs (32 B → 8.2→1.5 ns,
  256 B → 61.8→7.0 ns, 4096 B → 906→114 ns), well above the 2× target.
- Zero behavior change: all pre-existing WFA/Myers tests pass unchanged.
- No runtime SIMD dispatch and no `unsafe` in the common (u64) path; the architecture
  tiers contain small, bounded, differential-tested `unsafe` blocks (justified by SSE2/
  NEON being guaranteed on their targets).
- AVX2/AVX-512 (32/64 bytes/step) are intentionally **not** added here. They require
  runtime detection and yield diminishing returns over a 16-byte primitive; deferred
  until profiling on real workloads justifies the added dispatch complexity.
- The misleadingly named `wfa_extend_simd`/`wfa_extend_neon` entry points (historical,
  from the diagonal-DP era) are now genuinely SIMD-accelerated transitively, since the
  acceleration lives inside the shared `extend` primitive.

## Alternatives considered

- **`multiversion` runtime dispatch**: rejected for the inner loop — adds per-call
  overhead and nondeterminism for no benefit when the target baselines already guarantee
  SSE2/NEON.
- **Rely on LLVM autovectorization of the scalar loop**: rejected — the data-dependent
  early exit defeats reliable vectorization; performance would vary with compiler version.
- **Promote the `O(n·m)` SIMD diagonal fill to production**: rejected — quadratic, ~28×
  slower than WFA for the low-divergence alignments Phraya targets.
