# Architecture Decision Records

This directory records the **why** behind non-obvious engineering decisions in Phraya —
the choices a future maintainer (or a regulatory reviewer) would otherwise have to
reverse-engineer from the code.

Each ADR uses the [Michael Nygard format](https://cognitect.com/blog/2011/11/15/documenting-architecture-decisions):

- **Status** — Proposed | Accepted | Superseded by ADR-NNNN
- **Context** — the forces at play; what made this a decision rather than a default
- **Decision** — what we chose, stated plainly
- **Consequences** — what becomes easier and harder, including the costs we accept
- **Alternatives considered** — what we rejected and why

ADRs are immutable once Accepted. To change a decision, write a new ADR that supersedes
the old one (and update the old one's Status). Do not edit history.

## Index

| ADR | Title | Status |
|-----|-------|--------|
| [0001](0001-record-architecture-decisions.md) | Record architecture decisions | Accepted |
| [0002](0002-simd-match-extension.md) | SIMD-accelerated match extension primitive | Accepted |
| [0003](0003-alignment-strategy-ladder.md) | Alignment strategy ladder and Myers fitting default | Superseded by ADR-0008 |
| [0004](0004-fast-strategy-low-sensitivity.md) | Fast strategy: seed subsampling + divergence cutoff | Accepted (extended by ADR-0008) |
| [0005](0005-paired-end-proper-pair-fraction.md) | Paired-end proper-pair fraction filtering | Accepted |
| [0006](0006-simd-build-requirements-and-packaging.md) | SIMD build requirements and packaging (`target-cpu` / `ensure_simd`) | Proposed |
| [0007](0007-score-bounded-alternate-abandonment.md) | Score-bounded early abandonment of alternate anchors | Superseded by ADR-0012 |
| [0008](0008-strategy-ladder-recall-axis.md) | Strategy ladder as a single recall axis; `exact` → `sensitive` | Proposed |
| [0009](0009-dense-rated-plan-seeding.md) | Dense-rated minimizer plan; strategies subselect by density | Proposed |
| [0010](0010-filter-preset-rename-strict-tolerant.md) | Filter presets renamed `conservative`/`sensitive` → `strict`/`tolerant` | Proposed |
| [0012](0012-chaining-supersedes-branch-and-bound.md) | Seed chaining supersedes ADR-0007's branch-and-bound | Accepted |
