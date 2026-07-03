# 10. Filter presets renamed `conservative`/`sensitive` → `strict`/`tolerant`

- **Status**: Proposed
- **Date**: 2026-07-03

## Context

`phraya-filter` ships two named post-alignment variant-call presets (lib.rs,
`FilterPreset`): bundles of thresholds that decide which `VariantObservation`s survive into
VCF/TSV output.

| Threshold | `conservative` | `sensitive` |
|-----------|----------------|-------------|
| min coverage | 10× | 3× |
| min MAPQ | 30 | 20 |
| min allele frequency | 10 % | 2 % |
| exclude tandem repeats | yes | no |

- `conservative` = high-confidence calls only (deep coverage, well-placed reads, substantial
  allele fraction, repeats excluded). Few false positives, misses low-frequency real
  variants.
- `sensitive` = catch low-frequency variants at the cost of noise (thin coverage, lower
  MAPQ, alleles to 2 %, repeats kept).

ADR-0008 renames the `exact` **alignment strategy** to `sensitive`. That would put
"sensitive" at two layers meaning related-but-different things — a legibility footgun. The
word fits the *alignment recall* axis better than a *filter-stringency* knob (the filter is
really about stringency, not sensitivity per se), so the filter pair is what should move.

## Decision

Rename the filter presets:

- `conservative` → **`strict`**
- `sensitive` → **`tolerant`**

Rationale:

- Frees "sensitive" for the alignment strategy layer (ADR-0008); zero name overload.
  "Tolerant" appears in the alignment layer only as descriptive prose ("sensitive mode is
  divergence-tolerant"), never as a mode name.
- `strict`/`tolerant` is a clean single-axis antonym pairing (better than `strict`/`lenient`,
  which reads lopsided).
- `strict` is 6 letters vs `conservative`'s 12 — fixes the "anti-human length" ergonomics.
- **`tolerant` over `lenient` deliberately, for spelling.** "Lenient" fails the
  type-it-under-pressure test (the *e/ie* trap); "tolerant" is phonetic and unambiguous.
- The faint thematic resonance (both `sensitive` alignment and `tolerant` filtering
  "tolerate more") signals the same recall-maximizing philosophy across layers *without*
  reusing a token — a feature, not a collision.

**Hard rename, no deprecation alias.** Update the `FilterPreset` variants
(`Conservative`→`Strict`, `Sensitive`→`Tolerant`), the CLI parser (main.rs — `--preset`
matches and the error message), doc-strings, and tests. Phraya is pre-release; no
operational code depends on the old spellings. A lingering alias would mean two spellings in
docs/tests forever.

## Consequences

- Non-overloaded vocabulary across both layers:
  - **Alignment strategies:** `fast` / `balanced` / `sensitive`
  - **Filter presets:** `strict` / `tolerant`
- Breaking change: `--preset conservative` and `--preset sensitive` stop working. Documented
  in the changelog.
- The threshold *values* are unchanged — this is a pure rename of the two bundles.

## Alternatives considered

- **`strict` / `permissive`** — `permissive` is marginally more precise than `tolerant` but
  longer and no easier to spell. Rejected on ergonomics.
- **`strict` / `lenient`** — the initial proposal; rejected because `lenient` is a common
  misspelling risk.
- **Rename the alignment strategy instead (`exhaustive`/`thorough`), keep the filter
  `sensitive`.** Rejected: "sensitive" is the more accurate word for the *alignment recall*
  axis, and the filter preset (older, domain-standard, paired with `conservative`) is really
  a stringency knob better named on a strict/tolerant axis anyway.
- **Keep both, add a deprecation alias.** Rejected — pre-release with no external callers;
  a clean break avoids permanent dual spellings.
