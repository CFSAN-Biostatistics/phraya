# 1. Record architecture decisions

- **Status**: Accepted
- **Date**: 2026-06-19

## Context

Phraya is a SNP caller intended for use in scientific research and public-health
genomics. Decisions about alignment algorithms, scoring, and filtering have to be
*defensible* — in code review, in a methods section, and potentially to a regulator.
Several such decisions had already been made (hard-coded score-ratio threshold, k-mer
parameters, fitting alignment) and lived only as comments scattered through the code.

As the alignment engine grows (SIMD, multiple algorithms, speed/sensitivity strategies),
the rationale behind each choice matters as much as the choice itself. Code shows *what*;
it rarely shows *why this and not the obvious alternative*.

## Decision

We keep Architecture Decision Records in `docs/adr/`, one file per significant decision,
in Michael Nygard's format. A decision is "significant" if reversing it would be
expensive, if it trades off a quality attribute (speed vs. sensitivity, simplicity vs.
performance), or if a reasonable engineer would otherwise ask "why did they do it that
way?"

ADRs are immutable once Accepted; a decision is changed by superseding it with a new ADR.

## Consequences

- New significant decisions get an ADR alongside the implementation PR.
- Reviewers and downstream users have a single place to find rationale and rejected
  alternatives, instead of inferring intent from diffs.
- Minor or easily-reversible decisions stay out of ADRs to avoid noise; code comments
  remain the right tool for local detail.
