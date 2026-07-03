# 8. Strategy ladder as a single recall axis; `exact` → `sensitive`

- **Status**: Proposed
- **Date**: 2026-07-03
- **Supersedes**: ADR-0003 (redefines the ladder), extends ADR-0004

## Context

ADR-0003 defined `exact`/`balanced`/`fast` as presets bundling several independent knobs
(anchor enumeration, extension engine, divergence cutoff, coverage-window radius). The
Phase-3 benchmark and a review of `executor.rs` exposed two problems:

1. **`exact` is not distinguishable from `balanced` in results.** `build_anchors` gives
   `Balanced` and `Exact` the *identical* anchor set (a single shared match arm), and both
   run `score_alignments` with the same 0.95 filter. Since Myers ≡ WFA on edit distance
   (ADR-0003), the two produce byte-identical primary and alternative alignments. The
   benchmark confirms it: per-target CAS for `exact` equals `balanced` on every non-OOM
   target. `exact` is currently "balanced, but slower (always-WFA), with a narrower
   coverage window." It surfaces **zero** additional alternate placements. The name is
   downstream of a mechanism that does not exist.

2. **"Exact" is the wrong word regardless.** Exactness is *saturated* — all three
   strategies compute the optimal edit distance. The distinction that matters for Phraya's
   deferred-filtering thesis is not precision but **how much placement/variant ambiguity
   each mode preserves versus discards.**

A design principle pins down the missing mode: a downstream filter can only ever
*subtract* from the alignment superposition. You cannot filter *in* a read you never
placed or an alternate locus you never enumerated. **Recall must be bought at alignment
time or it is gone forever.** So the exhaustive end of the ladder is definitionally an
alignment-strategy concern, not a filter preset over `balanced` output.

## Decision

**The strategy ladder is sorted on one axis: recall / ambiguity preservation.** Each rung
differs by exactly one primary lever — the **anchor cap `K`** (how many candidate loci are
extended and reported) — plus, at the sensitive end, removal of the divergence cutoff.
Coverage-window radius ceases to be a strategy-defining trait (it was already orthogonally
overridable) and extension-engine choice is an implementation detail chosen for
speed-at-equal-results, not a strategy trait.

| Strategy | Anchor cap `K` | Divergence cutoff | Role |
|----------|----------------|-------------------|------|
| `fast` | 1 (best-voted, collapse) | 0.20 (drop hard reads) | allelic typing — precise point estimate, ambiguity-*intolerant* |
| `balanced` *(default)* | 5 (minimap2 `-N` parity) | none (BWA-faithful) | general BWA/minimap2 drop-in substitute |
| `sensitive` *(was `exact`)* | ∞ (all anchors) | none | maximal-recall variant/SNP discovery — enumerate every placement, drop nothing |

**`exact` is renamed `sensitive`.** The mode's identity is now concrete: `K = ∞`, dense
seeds (ADR-0009), no divergence cutoff, 0.95 retention — it enumerates every placement it
can seed and discards nothing. Its defining verb is *enumerate exhaustively*; its use case
is variant discovery via the full evidence superposition. This is the genuine antipode of
`fast` — not "more exact," but "maximally undecided": `fast` throws ambiguity away,
`balanced` picks a winner and notes the contenders, `sensitive` refuses to decide and hands
the whole superposition to the filter layer. "Sensitive" is the domain-native term for the
recall axis (minimap2/BLAST usage) and now fits the mode exactly.

**`balanced`'s `K = 5` mirrors minimap2's `-N` default.** Post-ADR-0007, `K` no longer
governs speed (junk anchors self-abandon cheaply regardless), so `K` is purely a
*reporting* knob: "report up to 5 placements." This also makes `balanced` a *more* faithful
BWA/minimap2 substitute — those tools cap reported alternates (BWA `XA`, minimap2 `-N`),
whereas the old uncapped enumeration was *more* exhaustive than the tools it emulates.

**Multi-mapping signal is preserved at `K = 5`.** Flagging a read as multi-mapping for
downstream filtering ("exclude variants where >50 % of supporting reads multi-map") needs
only *>1* placement, not all copies. Exhaustive per-copy enumeration becomes the exclusive
job of `sensitive`.

**The divergence cutoff stays `fast`-only and is unchanged.** It is not a speed lever — it
runs *after* extension (executor.rs) and saves no compute; it only drops reads. A read
carrying several nearby SNPs *is* a divergent read — the exact signal `sensitive` must
retain — so `sensitive` (and `balanced`, for BWA fidelity) never applies it.

**The 0.95 retention threshold is fixed across all three strategies.** The recall axis is
`K`-alone; `sensitive` retains more by *looking in more places* (K=∞) and *keeping divergent
reads*, not by lowering the reporting bar. Paralog-suspicion signal for variant discovery is
already computed at a better layer (`kmer_uniqueness`, `hotspot_intervals` in executor.rs),
so the filter can flag paralog-suspect SNPs without the aligner loosening retention.

## Consequences

- The ladder is legible: one axis, three named points, each a single-sentence promise.
  `fast` = typing, `balanced` = general substitute, `sensitive` = variant discovery.
- `sensitive` is the reason to reach for Phraya over BWA at all: not a better single
  alignment, but a complete, trustworthy evidence superposition that a BAM-then-filter
  pipeline cannot reconstruct after the fact.
- Renaming `exact` → `sensitive` is a **breaking CLI/API change** (`--strategy exact`, the
  `Strategy::Exact` variant, `AlignConfig::exact()`). Hard rename, no deprecation alias:
  Phraya is pre-release, shipped Phase-1 2026-06-06, and no operational code depends on it.
- ADR-0004's "top-K deferred" alternative is now resolved: `K` is the ladder's organizing
  parameter (`fast`=1, `balanced`=5, `sensitive`=∞).
- `sensitive` is only tractable because of ADR-0007 (score-bounded abandonment); without it
  `K=∞` would be the old `exact`'s 4600 s+ runtimes.

## Alternatives considered

- **Keep strategies as multi-dimensional presets (a speed/accuracy/precision cube).**
  Rejected as the *organizing* model — bundling several knobs per rung is what made `exact`
  incoherent. Presets remain multi-dimensional only to the extent their *utility* is legible
  (a named use case each); the underlying sort key is the single recall axis.
- **`sensitive` loosens the retention threshold (e.g. keep alternates to 0.90).** Rejected:
  keeps 0.95 a single opinionated constant, keeps one score-cap formula (ADR-0007) across
  all strategies, and the paralog signal lives better in the filter layer.
- **Name it `exhaustive` or `thorough`.** `exhaustive` is literally accurate but `sensitive`
  is the domain-native recall term; `thorough` is vaguer about *what*. Chose `sensitive`,
  resolving the resulting collision with the filter preset by renaming the preset
  (ADR-0010), because "sensitive" fits the *alignment recall* axis better than it fits a
  filter-stringency knob.
- **Adaptive `K` (keep anchors within X % of top vote).** Rejected for `balanced`: a second
  opaque parameter that still blows up on perfect tandem repeats where many loci tie. A
  fixed `K=5` is legible and bounds worst-case reporting unconditionally.
