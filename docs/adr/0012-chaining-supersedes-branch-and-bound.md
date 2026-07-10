# 12. Seed chaining supersedes ADR-0007's branch-and-bound

- **Status**: Accepted
- **Date**: 2026-07-10
- **Supersedes**: ADR-0007 (score-bounded early abandonment of alternate anchors)

## Context

ADR-0007 solved a real problem: `sensitive`'s K=∞ anchor cap enumerated every distinct
seed-derived target start, so a read landing in a repeat-dense region could spray hundreds
of spurious alternate anchors, each paying for a full WFA extension only to be discarded by
`score_alignments`'s 0.95 multi-mapping filter. ADR-0007's fix — extend the top-voted anchor
fully to establish an incumbent edit distance `d_best`, then cap every subsequent
alternate's WFA wavefront at the exact edit distance where the 0.95 score ratio is lost —
made `sensitive` tractable without changing its output. The fix was correct *by
construction*: an abandoned alternate provably could never have passed the existing filter,
so nothing reportable was lost.

Benchmarking in mid-2026 found phraya's `balanced` strategy running ~500x slower than a
comparable aligner (minibwa) on the same workload. Two allocation/scan bugs accounted for
most of that gap (see the O(n) repeat-region-scan and coverage-window-union fixes, issues
#219/#220) and closed it to ~5x. Reading minibwa's source for the remaining gap found the
actual cause: minibwa collapses raw seed hits into co-linear **chains** before running any
DP extension at all, so only ~2-4 extensions happen per read total. Phraya had no chaining
step — `build_anchors` voted each seed for an independent target-start and extended up to K
of those raw positions from scratch, each a full DP call.

The user explicitly authorized dropping `sensitive`'s algebraic byte-equivalence guarantee
(the exact property ADR-0007 was designed to preserve) in favor of empirical validation via
the simulated-read PA/CBS benchmark, on the stated principle that `sensitive` needs to be
"comprehensive at the expense of speed, to the extent it's useful" rather than provably
identical to a prior baseline. This unlocked replacing the three separate
Fast/Balanced/Sensitive raw-vote anchor-selection code paths with one seed→chain→extend
pipeline, parameterized by how many top chains get extended (`chain_cap`), rather than
scoping chaining to `balanced`/`fast` only and leaving `sensitive`'s branch-and-bound intact.

## Decision

**Seed chaining (`phraya-align/src/chaining.rs`, `chain_seeds`) replaces raw anchor voting
as the anchor-selection mechanism for all three strategies.** Chaining collapses co-linear
seeds into candidate loci *before* extension, so a repeat family that would have sprayed
hundreds of raw votes now produces one chain per genuine copy. `sensitive` additionally
takes a finite chain cap (`chain_cap(Strategy::Sensitive) = 50`), replacing its former
`K = ∞` raw-vote enumeration — a deliberate move away from ADR-0008's "enumerate literally
everything" framing, trading unbounded worst-case enumeration for predictable cost.

**`sensitive`'s extension loop in `align_oriented` now extends every chain-derived anchor
plainly, identically to `balanced`/`fast`.** The Sensitive-specific branch that ordered
extension by vote count and capped alternates via `wfa_extend_capped`/`score_bound_max_s`
is removed from `align_oriented`. Chaining already bounds the candidate count structurally
(each repeat copy → one chain, not one anchor per seed within it), and the K=50 cap bounds
it further, so the anchor set entering extension is small before ADR-0007's cheap-abandon
optimization would even have anything meaningful to prune.

**ADR-0007's functions (`score_bound_max_s`, `wfa_extend_capped`, `extend_alternates_bounded`)
remain defined and unit-tested**, but are no longer called from `align_oriented`. They are
slated for deletion once a hidden debug-only legacy anchor-selection toggle
(`AlignConfig::use_legacy_anchors`, added to A/B the pre-chaining and chained paths during
HPC validation) is removed — see the chaining redesign's implementation plan.

## Consequences

- `sensitive` is no longer provably byte-identical to its pre-chaining baseline. Its
  correctness contract is now empirical: the simulated-read placement-accuracy (PA) and
  call-based-sensitivity (CBS) benchmark (`scripts/benchmark/slurm/utils/sam_accuracy.py`,
  `phraya_accuracy.py`), not an algebraic proof. This matches how every other aligner in the
  comparison set is actually validated.
- `align_oriented` loses its per-strategy branch entirely: all three strategies now share
  one code path (seed → chain → cap to K → extend every survivor → score), differing only
  in `chain_cap(strategy)` and `extend_anchor`'s engine choice (Myers vs WFA), both already
  established as orthogonal implementation details per ADR-0008.
- The branch-and-bound machinery's removal simplifies `align_oriented` substantially (the
  Sensitive-specific block, including vote-ordering and incumbent-tracking, is deleted) with
  no loss of test coverage: `issue_183_score_bounded_branch_and_bound.rs` tests the ADR-0007
  functions directly (not via `align_oriented`), so it continues to validate them in
  isolation even though production code no longer calls them through that path.
- `sensitive`'s worst-case cost is now bounded by `chain_cap = 50` chains rather than by how
  cheaply a hopeless alternate can be abandoned — a stronger and simpler guarantee for a
  genuinely pathological repeat-saturated genome, at the cost of a hard (if generous)
  ceiling on enumerated placements that ADR-0008's original "enumerate everything" framing
  did not have.

## Alternatives considered

- **Keep ADR-0007's branch-and-bound as a fallback inside the chained path, gated on chain
  count exceeding some threshold.** Rejected as unnecessary complexity: chaining's own
  structural collapse plus the K=50 cap already bound the candidate count well below where
  branch-and-bound would meaningfully help; reintroducing it would resurrect the byte-
  equivalence bookkeeping the redesign explicitly moved away from, for a case that empirical
  benchmarking has not shown to occur.
- **Scope chaining to `balanced`/`fast` only, leave `sensitive`'s raw-vote + branch-and-bound
  untouched.** This was the plan's original framing before the user's explicit authorization
  to drop `sensitive`'s byte-equivalence guarantee. Rejected once that authorization was
  given: maintaining two anchor-selection mechanisms (chained for balanced/fast, raw-vote +
  branch-and-bound for sensitive) would have been more code for less benefit than
  collapsing all three onto one pipeline.
