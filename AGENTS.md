# Phraya — Domain Glossary

Canonical terms for the Phraya sequence aligner. Definitions only — no
implementation details, specs, or scratch notes. When code and this glossary
disagree, one of them is wrong; fix it.

## Terms

### Reference space

A **named, content-hashed coordinate space** that reads can be aligned
against. A plan may hold several, each heterogeneous in purpose (a host genome
to deplete against, an allele database to type against, a genome to call
variants on). Each reference space has:

- a **content hash** — the identity. A reference space *is* its content.
  Lookup and the align-time robustness check are both keyed on this, never on
  filesystem path or filename.
- a **name** — an optional human-facing label, assigned at plan time and
  stored *in the plan*. Convenience for selection and error messages only;
  never trusted as identity, because workflow engines stage, symlink, and
  rename files freely.

Identity is **content-dependent, not location-dependent**: the same bytes are
the same reference space no matter what the file is called or where it lives.

Supersedes the old positional convention where "the reference" was whatever
sequence sat at index 0 of the pool and every read was "everything after
index 0." Multiple reference spaces per plan is the new default; a single
reference is the degenerate case.

Reference spaces are **purely mechanical and role-free**. A space does not know
or record whether it exists to deplete against, type against, or call variants
on — that intent lives in the pipeline author's head and the surrounding
workflow, never in the plan.

### Sealed

**Fail-fast when an invocation diverges from what was planned.** The canonical
term for any check that turns "you're doing something the plan didn't
anticipate" from a tolerated warning into a hard error. First use:
`align --sealed` hard-errors on a presented reference whose content hash is not
in the plan's palette (vs. the tolerant default — warn and sketch on the fly).
Reach for "sealed" for every future plan-vs-invocation divergence guard;
deliberately distinct from "strict", which is a filter-stringency term.

### Cross-reference superposition

When reads are aligned against several reference spaces at once, alignment does
**not** pick a winning reference. It preserves the *superposition*: each read
retains its placements in every space it hits well enough, exactly as
multi-mapping already retains within-reference alternatives. Choosing a
reference — "this read is host, deplete it"; "this read types as allele X" — is
a **filter-step decision**, made post-hoc on the superposition, never enacted
in align. This is the multi-reference extension of Phraya's deferred-filtering
commitment: align surfaces possibilities, filter decides.

Consequences: "unassigned" (a read with no qualifying placement in any space)
and "ambiguous" (a read placed near-equally across spaces) are things a filter
*reads off* the superposition — not bins align computes. There is no `--best`,
no margin cutoff, and no per-read classification produced by align.

**Retention is per-space, not globally anchored.** A read's placements in space
X are retained or dropped by comparison to that read's best hit *within X*
(the existing single-reference `score_ratio ≥ 0.95` rule), independent of what
it does in any other space. This makes competitive alignment **composable**:
`align({A, B}) = align({A}) ∪ align({B})`. A read's presence and scores in a
space are a stable function of (read, space), never contingent on which other
references shared the invocation. A global anchor was rejected because it makes
align encode a decision (which space to drop) and, worse, discards exactly the
ambiguous reads that depletion must preserve.

**The queries sidecar stores absolute normalized identity `(1 − edit/len)` per
placement**, not the primary-relative ratio. Absolute identity is
invocation-independent, so cross-space margins survive per-space anchoring: a
filter reads host-0.96 vs target-0.98 straight off the sidecar. It is also
strictly more informative than the relative ratio, which cannot distinguish a
primary of 0.99 from one of 0.75.

### Use case (explanatory only — not a code construct)

The "Case 2 / 3 / 4" labels (reads-with-ref, contigs-with-reads, contigs-only)
are **teaching language** from the PRD for describing input shapes. They are
*not* a type, field, or control-flow construct. The old plan-level `use_case`
field was a leak of pedagogy into the type system and has been removed. The
mechanical shape of an alignment is derived at align time from what is actually
presented (reads vs contigs, which reference space) — never frozen into the
plan. Do not reintroduce a `UseCase` enum into stored formats.
