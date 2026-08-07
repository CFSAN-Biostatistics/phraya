# 14. Gap-affine scoring, defaulted in `sensitive`

- **Status**: Accepted
- **Date**: 2026-08-07

## Context

Phraya's current WFA/Myers scoring is uniform (linear-cost) edit distance: every
mismatch costs 1, every inserted/deleted base costs 1 regardless of run length — a 5bp
deletion costs the same as 5 isolated substitutions. This matches Phraya's haploid,
low-divergence bacterial-surveillance target and is what lets Myers and WFA stay
drop-in-equivalent (ADR-0003, ADR-0008: "Myers and WFA are both exact and produce
identical edit distances").

Real indel biology — and amino-acid alignment (ADR-0013) — is better modeled by
gap-affine cost: a k-bp indel is usually one mutational event, not k independent ones, so
`gap_open + k * gap_extend` scores it more faithfully than `k * mismatch_cost`. This is
the standard model (Smith-Waterman-affine; WFA2-lib's gap-affine mode) used by
general-purpose aligners (minimap2, BWA, the BLAST family) precisely because linear-gap
cost systematically misprices indel-heavy regions and multi-residue protein indels.

Myers' bit-parallel algorithm is fundamentally a linear-cost algorithm — there is no
affine-gap Myers; the recurrence has no notion of a differential open/extend cost. WFA,
by contrast, is the general algorithm family: gap-affine WFA (separate M/I/D wavefronts,
open cost `o` + extend cost `e` per indel) is a well-documented, standard extension
(WFA2-lib ships it as a mode, not a different algorithm).

The strategy ladder (ADR-0008) already sorts on one axis — recall/ambiguity preservation
via anchor cap `K` — and explicitly treats extension-engine choice as "an implementation
detail chosen for speed-at-equal-results, not a strategy trait." `Strategy::Sensitive`
(K=∞, no divergence cutoff) already runs WFA unconditionally — "the reference path."
`Balanced`/`Fast` are Myers-primary, with a WFA fallback only for queries exceeding
`MYERS_MAX_QUERY_LEN` (500bp), chosen purely because it was cheaper at equivalent results,
never because the ladder wanted a different cost model at that tier.

Team direction (this conversation, superseding the "engine choice never changes results"
framing insofar as it implied a hard invariant): it is acceptable for engine/cost-model
choice to change results between strategies — this is already implicit in the strategy
ladder's design, and won't astonish users, since choosing a strategy is already choosing
a different algorithm. Users who reach for `sensitive` already know they're paying a time
cost for a better/more-complete answer, so gap-affine becoming `sensitive`'s default (not
merely an opt-in) doesn't violate any existing expectation.

## Decision

- Implement a gap-affine WFA extension mode (separate M/I/D wavefronts, `gap_open` +
  `gap_extend` costs) alongside the existing linear-cost WFA. Both emit the same
  `Alignment`/`WfaResult`/CIGAR shape, so downstream (variant extraction, `.phraya`
  output, MAPQ) stays engine/cost-model agnostic — consistent with the existing
  convention that "all variants carry the same CIGAR convention regardless of engine"
  (ADR-0003).
- **`Strategy::Sensitive` defaults to gap-affine WFA.**
- **`Strategy::Balanced` and `Strategy::Fast` are unaffected and stay linear-cost,
  unconditionally.** Myers remains their primary engine (no affine Myers exists), and
  their WFA fallback for queries >500bp stays on the existing **linear**-cost WFA — that
  fallback exists purely to mirror Myers' behavior at lower cost for long queries, not to
  introduce a different cost model at that tier. Linear-cost scoring stays in the
  Myers-using strategies because Myers has no other mode to offer, full stop — not
  because of any residual belief that engine choice must preserve identical results.
- Add `--gap-model {linear|affine}` as an override, valid only with `--strategy
  sensitive` (clap `requires = "strategy=sensitive"`), defaulting to `affine` under
  `sensitive`. Passing `--gap-model affine` with `balanced`/`fast` is a hard CLI error at
  parse time, not a silent no-op or engine swap — consistent with Phraya's fail-fast
  ("sealed") philosophy for invocation-vs-plan mismatches.
- Default `gap_open`/`gap_extend` constants are an implementation decision made and
  documented in code (mirroring how `FAST_MAX_DIVERGENCE` and `SCORE_REPORT_THRESHOLD`
  are hard-coded, non-configurable Phraya opinions today), not fixed by this ADR. They
  should be chosen so an isolated SNP is still cheaper than opening a gap (gap_open >
  mismatch cost) and a long indel is charged close to linear once opened (gap_extend ≈
  mismatch cost) — appropriate to Phraya's close-strain-comparison domain rather than
  BLAST-style remote-homology defaults.
- `score_alignments`'s 0.95 retention threshold and MAPQ, currently defined over
  normalized (linear) edit distance, need an equivalent normalized-score definition under
  the affine cost model before `sensitive`'s default can ship — an implementation task
  this ADR calls out but does not resolve.

## Consequences

- `sensitive` becomes strictly more expensive per anchor (M/I/D wavefronts vs one
  wavefront). Accepted per team direction: `sensitive` users already pay for K=∞ recall;
  this doesn't touch `balanced`/`fast`, where throughput-sensitive workloads
  (screening/triage, general BWA/minimap2 substitute) actually live.
- `sensitive`'s results are no longer edit-distance-equivalent to `balanced`/`fast` on
  indel-bearing reads — a cost-model difference, not just a `K` difference. Explicitly
  endorsed by this ADR, extending ADR-0008's existing tolerance for engine choice as an
  implementation detail into cost-model choice for `sensitive` specifically.
- The Myers≡WFA differential suite is unchanged for `balanced`/`fast` and for
  `sensitive`'s linear path (kept reachable via `--gap-model linear` for anyone who wants
  the old exact-equivalence guarantee). A new differential suite is needed for gap-affine
  WFA — Myers cannot serve as the affine oracle, so correctness needs its own reference
  (e.g. exhaustive affine DP on small cases).
- Enables ADR-0013's protein-space alignment to eventually pair with a cost model
  appropriate to amino-acid indels, without forcing that pairing now — the two features
  are orthogonal and composable later, never coupled by this decision.
- `extract_variants_from_cigar` and downstream evidence extraction will see longer,
  consolidated indel runs more often under affine (it actively prefers merging adjacent
  edits into one gap over scattering them) — CIGAR parsing is already op-based so this
  should be a non-issue, but flagged for verification during implementation.

## Alternatives considered

- **Make gap-affine available under any strategy, with `balanced`/`fast` silently falling
  back to WFA whenever affine is requested (even for sub-500bp queries)**: rejected — this
  turns "Myers-using strategies" into "Myers-using strategies, except when a scoring flag
  silently swaps the whole engine," a larger behavioral surprise than a hard CLI error,
  and defeats `balanced`'s/`fast`'s entire cost premise for short reads.
- **`sensitive` supports both cost models with no default preference (`--gap-model`
  required, no default)**: rejected per explicit team direction — `sensitive` should
  default to affine since it is the more biologically faithful model and its users are
  already opting into the slowest, most-complete tier.
- **Approximate affine via post-hoc CIGAR compaction on the linear Myers/WFA result**
  (heuristic gap-merging instead of a true affine DP): rejected — produces a different
  edit distance than an actual affine optimum would find, breaking the "Myers and WFA are
  both exact" invariant for any path that used it, for uncertain benefit over just
  implementing gap-affine WFA properly.
