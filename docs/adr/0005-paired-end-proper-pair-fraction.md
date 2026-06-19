# 5. Paired-end proper-pair fraction filtering

- **Status**: Accepted
- **Date**: 2026-06-19

## Context

Phraya supports filtering variants by paired-end structural signals (proper pairs, insert
size, discordant pairs). The first implementation exposed `--require-proper-pairs` as a
boolean and read the flag from a single read's `MateInfo` attached to each
`VariantObservation`.

Two problems surfaced:

1. **Granularity.** "Require proper pairs" is not a per-read yes/no for a *position* — a
   site is covered by many reads, some properly paired and some not. A boolean cannot
   express "at least 90% of covering reads are proper pairs."
2. **Merge loses mate info.** `merge_phraya_files` rebuilds observations via
   `VariantObservation::new()`, which sets `mate_info: None`. Any filter reading
   per-read `mate_info` therefore silently stops working after a multi-sample merge —
   exactly when cohort-level filtering matters most.

A third issue appeared while fixing the above: adding aggregate count fields after the
existing `mate_info: Option<MateInfo>` (which carried
`skip_serializing_if = "Option::is_none"`) corrupted the `.phraya` MessagePack files.
rmp-serde encodes structs as positional arrays; skipping a *middle* field shifts every
later field, so the new counts deserialized into the wrong slots
(`invalid type: integer 0, expected struct MateInfo`).

## Decision

1. **`--require-proper-pairs` takes a numeric threshold**, not a boolean: a fraction
   (`0 < N ≤ 1`) or a percentage (`1 < N ≤ 100`). It is the minimum fraction of a
   position's *paired covering reads* that must be properly paired. Internally
   `require_proper_pairs: Option<f64>`.

2. **Aggregate counts live on the observation and survive merge.**
   `VariantObservation` gains `total_paired_count` and `proper_pair_count` (u32). The
   executor stamps `(1, proper?1:0)` per contributing read; `merge_phraya_files` *sums*
   them per position. `proper_pair_fraction()` derives the ratio. These counts are the
   merge-stable source of truth, replacing per-read `mate_info` for the proper-pair
   filter. Per-read `mate_info` is retained only for the pre-merge insert-size/discordant
   filters.

3. **`mate_info` is the final struct field**, with `skip_serializing_if` restored. As the
   trailing field, omitting it when `None` (the common, unpaired case) only shortens the
   end of the MessagePack array; `#[serde(default)]` restores it to `None` on read. This
   keeps the array compact *and* keeps JSON output free of a null `mate_info`, while the
   always-present count fields sit at fixed positions ahead of it.

## Consequences

- Proper-pair filtering works correctly after merge (the common cohort case), and is
  expressible as a confidence threshold rather than an all-or-nothing flag.
- `.phraya` files round-trip correctly with the new fields; field *order* is now
  load-bearing for the MessagePack encoding — new persisted fields must be appended ahead
  of `mate_info`, or `mate_info` kept last.
- Slightly more state per observation (two u32s); negligible given existing fields.

## Alternatives considered

- **Keep the boolean filter**: rejected — cannot express a fraction, and the per-read
  `mate_info` it depended on is destroyed by merge.
- **Preserve `mate_info` through merge**: rejected — merge is position-centric and
  order-independent by design; carrying per-read mate structures through it would bloat
  the format and re-introduce per-read coupling. Aggregate counts are the right
  granularity for a position-level filter.
- **Serialize `mate_info` always (as null) instead of reordering**: workable but adds a
  nil per unpaired observation and a noisy null in JSON; reordering keeps both encodings
  clean at the cost of making field order significant.
