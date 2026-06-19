# 4. Fast strategy: seed subsampling + divergence cutoff

- **Status**: Accepted
- **Date**: 2026-06-19

## Context

The strategy ladder (ADR-0003) needs a genuinely *faster, lower-sensitivity* tier for
survey/screening workloads (classification, presence/absence, QC triage) where missing a
few reads or alternative mappings is acceptable in exchange for throughput.

The known per-read bottleneck is not the aligner itself but **anchor count**: a minimizer
landing in a tandem repeat or low-complexity region produces thousands of shared seeds,
hence thousands of distinct candidate target starts, hence thousands of WFA/Myers calls.
The default (`balanced`/`exact`) deliberately keeps every distinct seed position to
preserve multi-mapping sensitivity.

## Decision

`fast` differs from `balanced` not in extension engine but in two cheap, well-precedented
levers:

1. **Seed-vote subsampling.** Group shared minimizers by implied target start
   (`target_pos − query_pos`), count votes per start, and keep only the single
   best-supported anchor (ties broken toward the earliest position). Falls back to
   `(0,0)` when no seeds are shared. This collapses the per-read anchor count to `O(1)`
   regardless of repetitiveness — the same chaining intuition minimap2/bwa use.

2. **Divergence cutoff.** After scoring, drop the read entirely if its best alignment
   exceeds `FAST_MAX_DIVERGENCE` (0.20 — mismatches+indels per base). Hard hopeless reads
   are abandoned instead of reported.

Both thresholds are hard-coded opinions, consistent with Phraya's existing
score-ratio-0.95 precedent. They are not user-tunable today.

## Consequences

- Per-read work is bounded even in pathological repeat regions: one anchor, one extension.
- **Sensitivity is sacrificed deliberately and measurably** (pinned by tests):
  - Multi-mapping is under-reported — against a tandem-duplicated target, `exact` records
    both equally-good positions while `fast` keeps only the best-voted one.
  - Divergent reads (>20%) are dropped, so they contribute no variants and no coverage.
- For a uniquely-mapping read within the cutoff, `fast` produces the same variant calls as
  `balanced` — the speed comes from doing less bookkeeping, not from a worse answer on
  easy reads.
- `fast` is therefore appropriate for screening/triage, **not** for outbreak analysis or
  any application that depends on complete, multi-mapping-aware variant evidence.

## Alternatives considered

- **Top-K anchors (K>1) instead of top-1**: deferred — top-1 is the simplest defensible
  point and the largest speed win; K could become a tunable later if screening accuracy
  demands it.
- **X-drop / banded WFA early termination as the primary lever**: useful but secondary —
  Myers (the default engine) computes the full column DP with no early-exit benefit, and
  anchor explosion, not extension cost, is the measured bottleneck. The divergence cutoff
  captures most of the "abandon hopeless reads" value without an engine change.
- **Read subsampling (align every Nth read)**: rejected here — it is a pipeline-level
  concern, not a per-alignment strategy, and would belong upstream of the executor.
