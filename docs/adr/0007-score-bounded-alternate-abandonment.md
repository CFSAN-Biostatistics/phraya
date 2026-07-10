# 7. Score-bounded early abandonment of alternate anchors

- **Status**: Superseded by ADR-0012
- **Date**: 2026-07-03

## Context

The measured per-read bottleneck is anchor count, not extension depth (ADR-0003,
ADR-0004). `balanced`/`exact` enumerate every distinct seed-derived target start and run
a full extension on each; in a repeat that sprays thousands of seeds this is thousands of
extensions per read. The Phase-3 benchmark quantifies it: on the chicken target,
`fast` (1 anchor) ran in 142 s while `balanced` (all anchors) took 1737 s at *identical*
plan cost — a ~12× gap that is almost entirely spurious-anchor extension work.

Two facts about the extension engines constrain any fix:

1. **WFA is boundable; Myers is not.** `fill_wfa_fitting` is an `O(s·n)` wavefront with a
   `for s in 1..=max_s` loop that already exits at the first fitting-end (wfa_simd.rs).
   Capping `max_s` makes it abandon a hopeless alignment mid-search. Myers
   (`myers_fitting_impl`) is a bit-parallel *column sweep* whose cost is
   `O(n·⌈m/64⌉)` **independent of edit distance** — the score only emerges after the full
   sweep, so there is no cheap early-exit.

2. **The 0.95 reporting filter runs *after* extension.** `score_alignments` computes the
   score ratio and drops sub-0.95 alternatives only once every anchor is fully aligned
   (lib.rs). So an alternate destined for the trash still costs its entire extension.
   "Discard alternates by score before extending" is impossible as literally stated — the
   score *is* the extension.

To skip work losslessly we need a discard signal available *before* an anchor's extension
completes. The 0.95 threshold supplies exactly that bound.

## Decision

Extend alternates through a **score-capped WFA** whose cap is derived from the 0.95
reporting threshold, so an anchor abandons the instant it provably cannot be reported.

**Mechanism (C): split engines by role.**

- The **primary** anchor (highest seed vote) is extended with **Myers** — fast, fixed
  cost, no bound needed. This yields the incumbent edit distance `d_best`.
- Each remaining **alternate** anchor is extended with **WFA**, its `s`-loop capped at

  ```
  max_s = floor(0.05·L + 0.95·d_best)        L = query length
  ```

  the exact edit distance at which the score ratio `(1 − d_alt/L)/(1 − d_best/L)` falls to
  0.95. An anchor that reaches no fitting-end within `max_s` is abandoned.

Because Myers and WFA are proven to produce identical edit distances (ADR-0003), an
alternate produced by WFA is result-identical to what Myers would have produced — the
split costs no accuracy.

**Dynamic re-tightening (branch-and-bound, ii).** Anchors are processed in vote order,
but `d_best` is updated whenever an extension returns a smaller distance, and `max_s` is
recomputed from the new `d_best`. The bound is monotonically non-increasing, so later
anchors get the tightest possible cap. This closes the "top vote was wrong" hole: even if
a repeat copy is extended first, the moment the true low-divergence locus is reached its
small `d` tightens the bound and the remaining junk anchors die fast.

**The cap is dual-purpose and provably safe.** Since
`0.05·L + 0.95·d_best − d_best = 0.05·(L − d_best) ≥ 0`, the cap is always `≥ d_best`.
Therefore any anchor that could *beat* the incumbent (a better primary) is strictly inside
the cap and always runs to completion — **the bound can never prune a potential new
primary.** One number serves both jobs: "don't miss a better primary" and "don't waste
time on a filter-doomed alternate."

**Prerequisite — abandonment sentinel.** `fill_wfa_fitting`'s loop-exhaustion fallback
currently returns `(String::new(), 0, t.len())` — edit distance **0**. A capped WFA that
abandons must not report a *perfect* alignment. The fallback must return an explicit
*abandoned* sentinel that `score_alignments` discards outright (never a candidate, never a
primary). We mirror WFA2-lib's `max_alignment_score` "unaligned" semantics rather than
inventing our own.

This is sequenced **before** the anchor-cap `K` work (ADR-0008): it is
correctness-preserving (differential tests pass unchanged — only work whose result was
destined for the filter is skipped) and it de-risks `K` by decoupling "how many placements
we report" from "how much time we spend."

## Consequences

- Spurious anchors on a well-placed read (small `d_best` → tiny `max_s`, e.g. `max_s = 7`
  for a perfect 150 bp read) abandon after ~7 wavefront steps instead of a full extension.
  Uniquely-mapping reads — the common case — get their repeat-noise annihilated for free.
- **The `sensitive`/exhaustive strategy is the largest beneficiary.** It keeps every anchor
  by design and thus has no anchor-cap speedup available; score-bounded abandonment is its
  *only* free lunch. It converts that mode from *impossibly slow* (full extension per
  anchor) to *expensively complete* (cheap abandon per anchor) — still `O(cap·n)` per
  spurious anchor, not free, but tractable.
- **The hard-coded 0.95 is now load-bearing for *speed*, not just for *what is reported*.**
  Loosening it toward 0.90 makes every read do more work; tightening toward 0.99 prunes
  more aggressively. This coupling is accepted and deliberate — "report more" genuinely
  *should* cost more — but the constant must be commented in code as having two masters, so
  a future maintainer does not retune it for output reasons without realising they are also
  changing the performance profile.
- Two engines now participate in a single read's alignment (Myers primary, WFA alternates).
  Acceptable because they are result-identical; the split is invisible downstream.

## Alternatives considered

- **Banded Myers (Ukkonen).** Restrict the bit-parallel sweep to a diagonal band of width
  ~`2·max_s`, giving Myers an edit-distance bound too. Strictly better asymptotically, but
  it is a genuinely new, error-prone algorithm (band/block alignment, edge correctness) —
  exactly the implementation risk we sequenced this decision first to avoid. **Deferred**:
  revisit only if profiling shows the *primary* Myers sweep (not the anchor multiplier) is
  the bottleneck.
- **Cheap pre-filter before extending** (estimate a lower bound on each anchor's edit
  distance, skip if it exceeds `max_s`). The purest "do nothing," but a *sound* lower bound
  is hard — seed-vote count is not one — so it risks being lossy, violating the
  correctness-preserving contract. Rejected.
- **Static bound (extend primary once, fix `max_s`, never re-tighten).** Simpler, but loose
  whenever the top-voted anchor is not the true locus — precisely the repetitive reads that
  motivate the change. Rejected in favour of dynamic re-tightening.
- **`⌊d_best / 0.95⌋` as the cap.** The formula considered first and discarded: it ignores
  query length and would wrongly prune *every* alternate of a perfect read (`d_best = 0 →
  cap 0`). The correct, length-aware cap is `floor(0.05·L + 0.95·d_best)`.
