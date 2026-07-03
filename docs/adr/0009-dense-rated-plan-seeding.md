# 9. Dense-rated minimizer plan; strategies subselect by density

- **Status**: Proposed
- **Date**: 2026-07-03

## Context

The `sensitive` strategy (ADR-0008) wants maximal recall. Anchor cap `K = ∞` recovers true
loci that were *out-voted* by repeat copies, plus the complete paralog placement set — but
it cannot recover a read whose true locus has **zero surviving shared minimizers** because
the variants we are hunting knocked them all out. That read simply has no anchor to
enumerate. Recovering it needs **denser seeding** (smaller `w`, more minimizers).

But denser seeding collides with a headline design point. Every strategy's seeds come from
`find_seeds_indexed` over the query's minimizer sketch, and that sketch is computed *once
at plan time* (k=21, w=11) and stored in `.phrayaplan` for reuse across `align`. Denser
seeding at align time would re-sketch every read, discarding the sketch-reuse optimization
and paying plan-overhead-per-read again.

The resolution: compute dense seeds **in the plan** and let each strategy subselect the
density it needs. The align-time cost of this is essentially zero — the dominant align cost
is *extension*, and extension is untouched for `fast`/`balanced` (see below). The cost is
paid once, at plan time, in plan compute and plan file size. For Phraya's target domain
(bacterial genomes, 2–6 Mbp) that is KB-scale — negligible. For out-of-domain large genomes
a `--sparse` escape hatch avoids the plan bloat.

## Decision

**`phraya plan` computes a dense sketch by default and tags which minimizers belong to the
canonical w=11 set.** The plan carries one dense minimizer set per sequence with a
per-minimizer flag (or equivalent membership set) marking the w=11 subset.

- `fast` and `balanced` use the **tagged w=11 subset** — provably the same seed set as
  today.
- `sensitive` uses the **full dense set** — recovering divergent-region reads that the
  sparse set would miss, *while sketch reuse is preserved* (nothing re-sketches at align
  time; strategies subselect from one stored plan).

**Byte-equivalence for `fast`/`balanced` is guaranteed by construction, not by theorem.**
Standard minimizers nest (the w=11 set ⊆ a denser w=5 set), but canonical minimizers in
`simd-minimizers` have canonicalization and tie-breaking whose nesting at window edges we do
not rely on. Therefore we **compute the w=11 sketch exactly as today AND the dense sketch,
and explicitly tag** which dense minimizers are the w=11 ones. `fast`/`balanced` read the
tagged subset — identical to today by construction. We do **not** threshold the dense set
down by density (that shortcut loses the window-coverage guarantee — every `w+k−1` window
contains a seed — producing seeding *gaps* and *lower* recall; it is strictly worse, not a
valid relaxation).

**Density ceiling.** Canonicality requires `l = w + k − 1` odd (CLAUDE.md). With k=21 the
dense `w` must be odd: `w ∈ {5, 7, 9}`, roughly up to 2× density at w=5. The dense level is
chosen within this bound.

**Default dense, `--sparse` opt-out.**

- Default dense: in-domain the plan tax is KB-scale noise; it honours plan-once/align-many
  (one canonical plan supports *any* strategy without re-planning — the re-alignment
  workflow that motivated this); and it removes a footgun (opt-in would make "sensitive on a
  sparse plan" either error or silently under-seed — silent degradation of the
  maximal-recall mode is the worst failure).
- `--sparse` opt-out serves the out-of-domain large-genome case (human/wheat) where plan
  bloat could matter — and those already OOM for unrelated reasons.

## Consequences

- One plan supports every strategy. A workflow can re-align the same plan under `fast`,
  `balanced`, and `sensitive` with no re-planning — the payoff that justifies the dense
  default.
- `fast`/`balanced` results are unchanged (byte-identical seed set) — no re-benchmark of
  existing modes required.
- Align-time tax for `fast`/`balanced` is two `O(minimizers)` integer-filter passes (one per
  target, one per read) — noise against a single extension. The real costs are plan-time
  compute (~2× a cheap simd-minimizers pass, once) and plan file size (2–3× of already-tiny
  sketches; the plan is transmitted to every worker, so this is the number that matters at
  scale — KB-scale in-domain).
- `sensitive` gets *both* recalls in one mode: complete paralog enumeration (from `K=∞`) and
  divergent-region rescue (from dense seeds), without breaking sketch reuse.

## Alternatives considered

- **Re-sketch reads at align time for `sensitive`.** Rejected — breaks the sketch-reuse
  design point and pays plan overhead per read. The dense-plan approach preserves reuse.
- **Threshold the dense set down to w=11 density (skip tagging).** Rejected — an arbitrary
  density-thresholded subset loses the window-coverage guarantee, causing seeding gaps and
  *lower* `fast`/`balanced` recall. Strictly worse; the tag-and-subset construction avoids
  it.
- **Opt-in dense (`plan --dense`), sparse by default.** Rejected — reintroduces the footgun
  of running `sensitive` on a sparse plan, and the in-domain plan tax is too small to
  justify defaulting to the degraded plan.
- **Dense-informed anchor *ranking* for `balanced`** (rank the top-`K` on the finer
  dense-vote signal rather than coarse w=11 votes; same extension count, better-chosen `K`).
  This breaks byte-equivalence for `balanced` and could raise placement accuracy — a real
  lever. **Deferred as a measured experiment**, gated on actual PA gains against the
  T3–T7b harness: ship `fast`/`balanced` byte-identical, relax byte-equivalence for
  `balanced` *only if* the dense-ranking experiment demonstrably improves PA. The dense
  signal is already in the plan by default, so the experiment costs no architecture to
  defer. `fast` stays byte-identical permanently — its value is the confident point
  estimate, for which finer ranking is irrelevant at `K=1`.
