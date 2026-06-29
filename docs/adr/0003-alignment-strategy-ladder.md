# 3. Alignment strategy ladder and Myers fitting default

- **Status**: Accepted
- **Date**: 2026-06-19

## Context

Phraya exposes a `--strategy` flag (`exact`/`balanced`/`fast`) so users pick an
accuracy/speed tradeoff with a flag instead of switching tools. Before this decision the
flag only changed the coverage-annotation window radius; every strategy ran the same
seeded WFA. We now have a second exact engine — Myers bit-parallel edit distance — and a
genuine need for a speed/sensitivity continuum.

Three facts drove the design:

1. **Myers and WFA are both exact** and produce identical edit distances (verified by a
   320-case randomized SNP+indel differential sweep, and by identical variant calls
   across the entire integration suite). Myers is `O(nm/w)`; WFA is `O(s·n)`.
2. **Reads are windowed against a longer reference** (~2× read length around a seed).
   The correct mode is *fitting* alignment — query fully consumed, target end free — so
   the unconsumed reference tail is not charged as edits. Both engines must do this.
3. The original spec's `exact` ("single `(0,0)` anchor, full target, no seeds") is
   **broken**: WFA fitting fixes the target *start* at the anchor, so `(0,0)` only finds
   reads mapping at reference position 0. A read mapping deep in a megabase reference
   would be charged thousands of edits. Truly seedless search needs free-*start*
   semi-global alignment, which we do not yet have.

## Decision

`Strategy` is a convenience preset that selects an **algorithm** and a **default coverage
window radius**. The window radius is independently overridable
(`--coverage-window N` / `AlignConfig::with_coverage_window_radius`), so choosing an
algorithm never silently changes annotation width.

| Strategy | Algorithm | Anchors | Window | Role |
|----------|-----------|---------|--------|------|
| `exact` | seeded WFA | all distinct seed target-starts + `(0,0)` | ±25 bp | canonical reference path |
| `balanced` *(default)* | Myers fitting ≤500 bp, WFA fallback above | same as exact | ±50 bp | exact results, faster engine |
| `fast` | Myers/WFA + cutoff | single best-voted target-start | ±150 bp | low-sensitivity survey (see ADR-0004) |

Both engines run in **fitting** mode above the `n ≤ m + m/2 + 10` length threshold and
fall back to global below it (the threshold `fill_wfa_fitting` already used), so Myers and
WFA stay mutually consistent. The global `myers_edit_distance` primitive is left untouched
(its contract — e.g. empty query vs. target = N insertions — is intentionally global);
fitting is a separate `myers_extend` path.

**`exact` is redefined** as the most-sensitive *seeded* WFA path (all anchors, no cutoff)
— the trusted reference engine — not a from-scratch seedless search. Seedless/free-start
alignment is deferred to its own design (see Consequences).

The Myers cutoff (`MYERS_MAX_QUERY_LEN = 500`) is a hard-coded opinion: above it WFA's
`O(s·n)` beats Myers' length-quadratic term for the low-divergence reads Phraya targets.

## Consequences

- The default engine becomes Myers (via `balanced`). Safe because it is provably
  equivalent to WFA on results; faster on short and higher-divergence reads.
- `exact` remains available as the WFA reference for users who want the canonical
  algorithm with no approximation — valuable for validation and regulatory contexts even
  though `balanced` matches it.
- All variants carry the same CIGAR convention regardless of engine (M/X consume both,
  `I` = target-only, `D` = query-only), so downstream variant extraction is
  engine-agnostic.
- **Deferred**: true seedless / free-start semi-global alignment (find a read anywhere in
  the reference without seeds). This is inherently `~O(reference length)` per read and
  needs its own ADR and validation; it is *not* what `exact` provides today.
- The PRD's original strategy mapping (FM-index seeding; `exact` = full seedless WFA) is
  superseded by this minimizer-seeded ladder; the PRD has been annotated accordingly.

## Alternatives considered

- **Keep WFA as the default, Myers as opt-in**: rejected — Myers is equivalent and
  faster; making it the default delivers the win transparently, and `exact` preserves the
  WFA reference for anyone who wants it.
- **One enum carrying both algorithm and window radius with no override**: rejected —
  conflates two unrelated concerns; a user choosing an algorithm should not have their
  coverage annotation silently resized. Hence the orthogonal override.
- **Implement seedless `exact` now (free-start semi-global WFA)**: rejected for this
  change — correct but substantial, with `~O(reference)` per-read cost and its own
  validation burden. Shipping the broken `(0,0)` version would have been worse than
  deferring.
