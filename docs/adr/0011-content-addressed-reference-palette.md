# 11. Content-addressed reference palette and cross-space superposition

- **Status**: Accepted
- **Date**: 2026-07-22

## Context

Phraya's original pipeline aligns a read set against a *single* coordinate space (one
reference, or a Case-3 centroid). A growing class of use cases — host/pathogen depletion,
best-space in-silico typing, competitive mapping against a panel of candidate references —
needs each read considered against *several* reference spaces at once, with the results
comparable across spaces.

Two things made this awkward:

1. **Identity by path, not content.** A plan referred to its reference by file path. The same
   reference presented from a different path (or a benchmark re-run) could not be recognized
   as "already planned," and there was no way to carry more than one reference space in a
   plan.
2. **Per-space score is relative.** The `.phraya.queries` sidecar stored a *primary-relative*
   score ratio, which cannot distinguish a primary of 0.99 from one of 0.75 and is meaningless
   across spaces (each space anchors its own primary). Cross-space reasoning was therefore
   impossible from the artifacts.

This ADR was implemented incrementally across issues #196 (plan-side palette type), #197/#222
(content-hash reference space in the plan, repeatable `plan --reference`), #200/#227 (content-
hash-keyed read sketches), #226 (align-side `--reference` wiring, hit/miss resolution,
composability), #198 (cross-space sidecar with absolute identity), and #199 (`--sealed`).

## Decision

### Content-addressed reference spaces

A reference is identified by a **strong content hash** of its bytes, not its path. A plan
carries a *palette* — `PhrayaPlan.reference_space: Vec<ReferenceSpace>` — where each
`ReferenceSpace` is `{ content_hash, name: Option<String>, sketches }`. `plan --reference` is
repeatable and folds each presented reference into the palette; `get_reference_space(hash)` is
the lookup.

### Align-side resolution: tolerant by default, sealed on request

`phraya align --reference <ref>...` is a dedicated mode (mutually exclusive with the
traditional `QUERY_ID`/`TARGET_ID`, `--worker`, and `--ensure` modes; `--output` becomes a
directory). Each presented reference resolves by content hash:

- **hit** → reuse the planned sketch (no recompute; `TargetContext::build_with_sketch`);
- **miss** → **tolerant default**: warn and sketch on the fly, then align anyway;
- **miss under `--sealed`** → hard error naming the offending reference.

"Sealed" is the canonical fail-fast term for "this invocation diverged from what was planned,"
deliberately distinct from the filter-stringency term "strict" (ADR-0010). Tolerant is the
default so ad-hoc runs Just Work; sealed is for production pipelines that want an unplanned
reference to fail loudly once rather than warn quietly across many workers.

### Per-space artifacts + one cross-space sidecar

Align writes one `<output>/<label>.phraya` per space (`label` = palette name, else a content-
hash prefix) — each a clean single-coordinate artifact — plus one
`<output>/cross_space.phraya.queries` spanning the palette.

- Per-space `.phraya`/`.phraya.queries` store **absolute normalized identity** `(1 - edit/len)`
  per placement, for primary and alternates alike. Absolute identity is invocation-independent
  (a fact about `(read, space)`), so it survives per-space anchoring and is comparable across
  reads and spaces. This replaces the old primary-relative ratio.
- The cross-space sidecar is `CrossSpaceQueryIndex = query -> [(space, pos, identity)]`. A read
  placed in N spaces appears once per space with directly comparable identities, so a
  downstream filter can read a cross-space margin (e.g. host 0.96 vs target 0.98) straight off
  the sidecar. All cross-space coupling lives here.

### Composability

The read pool excludes every sequence whose content hash is in the plan's palette (not merely
the subset presented this invocation). This makes the read pool invocation-independent, so:

```
align({A, B}) = align({A}) ∪ align({B})
```

with each per-space `.phraya` byte-identical whether produced alone or alongside another space.
(Necessary because `phraya plan` folds its first reference into `input_files` for backward-
compatible task generation, so a palette reference can reappear among the inputs; without the
palette-wide exclusion it would be aligned *as a read* against the other spaces.)

### Non-goals

The align stage enacts **no** decision over the superposition: no `--best`, no margin cutoff,
no per-read host/target classification. Binning is deferred to `phraya filter` (issue #202).
Align produces the competitive superposition; filter decides what it means.

## Consequences

- Benchmark/regression harnesses can use a content hash as an equality gate for a reference,
  and multi-reference runs are decomposable for parallelism.
- The stale "primary-relative ratio" description of the sidecar (CLAUDE.md) is corrected to
  absolute identity.
- `.phraya` output remains non-byte-deterministic run-to-run within a single space when a read
  multi-maps (HashMap iteration order leaks into serialization); the composability guarantee is
  same-invocation determinism plus alone-vs-combined byte-identity, both of which hold.
- Filter-side cross-reference operations (depletion, best-space typing) build directly on the
  cross-space sidecar (issue #202, still open / ready-for-human).
