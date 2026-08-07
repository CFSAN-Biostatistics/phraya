# 13. Protein-space alphabet support

- **Status**: Accepted
- **Date**: 2026-08-07

## Context

Phraya's alignment core is already largely alphabet-agnostic, confirmed by reading the
kernels rather than assuming:

- `Sequence` stores raw `Vec<u8>` — never 2-bit packed.
- WFA's SIMD match-extension (`count_matching_prefix`) is a pure byte-compare intrinsic.
- Myers' bit-parallel DP already builds a `[[0u64; 256]; num_blocks]` Peq table
  (`wfa_simd.rs`) — full-byte-alphabet, not a 4-symbol table. Protein costs nothing extra
  here.
- Chaining (`chaining.rs`) operates on abstract `(query_pos, target_pos)` seed
  coordinates with no base semantics.

Three places do assume nucleotide DNA, all in the seeding/orientation layer:

1. `sketch()` (`phraya-core/src/types.rs`) calls `simd_minimizers::canonical_minimizers`,
   which picks `min(kmer, revcomp(kmer))`. Meaningless without a defined complement
   alphabet.
2. `align_read` (`phraya-align/src/executor.rs`) unconditionally tries both the forward
   sequence and `reverse_complement()` (types.rs), keeping whichever orientation scores
   better — a DNA-strand concept with no protein analog (no reverse-frame search is in
   scope here).
3. `DEFAULT_K = 21, DEFAULT_W = 11` (types.rs) are tuned for 4-symbol entropy (~2
   bits/base). At 20-symbol amino-acid entropy (~4.3 bits/aa), k=21 is >90 bits of
   specificity — homologous but divergent proteins would share essentially no seeds.

Research into the already-vendored `simd-minimizers` 2.3.1 crate confirms a non-canonical
entry point exists as a first-class peer of `canonical_minimizers`:
`simd_minimizers::minimizers(k, w)` returns the same `Builder` type (`CANONICAL = false`),
exposes the same `.run()` API, and is generic over any `Seq` implementation including
`packed_seq::AsciiSeq` — i.e. arbitrary-byte SIMD minimizers, not a scalar fallback. This
is a call-site swap, not new integration work.

`phraya filter`'s VCF emission is spec'd to nucleotide IUPAC REF/ALT alleles (VCF 4.2) —
redefining that for amino acids is out of scope here.

## Decision

- Add an `Alphabet` type (`Dna | Protein`) detected automatically at `phraya plan` time
  from sequence byte content, using the same mechanism as the existing
  `detect_use_case` FASTA/FASTQ/contig-vs-read content sniffing. Store it in
  `PhrayaPlan` so `phraya align` recovers it without re-sniffing files.
- Add `--alphabet {auto|dna|protein}` on `plan`, default `auto`, as an override for the
  one genuine detection ambiguity (below).
- For `Alphabet::Protein`: sketch with `simd_minimizers::minimizers()` (non-canonical) at
  protein-tuned `k`/`w` defaults (in the range BLASTP/DIAMOND use for amino-acid seeding,
  k≈5–7 — exact constants are an implementation decision, not fixed here); skip
  `reverse_complement`/dual-strand search entirely. `Alphabet::Dna` keeps the existing
  canonical/dual-strand path byte-for-byte unchanged.
- Extension engines (WFA, Myers) are untouched — both are already alphabet-blind.
  Protein alignment reuses Phraya's existing uniform edit-distance model, not a
  substitution matrix — same rationale as DNA today, appropriate for close-strain
  comparison, and keeps this a plumbing change rather than a new algorithm class.
  (Gap-affine scoring, ADR-0014, is an orthogonal axis that composes with this freely
  once both exist.)
- `phraya filter`'s VCF output stays DNA-only. Protein alignment output is
  `.phraya`/TSV only (CIGAR + edit distance) — no VCF REF/ALT emission for amino acids.
  Protein *variant calling* (an amino-acid-substitution evidence model) is explicitly out
  of scope for this ADR.

## Consequences

- Zero changes to WFA/Myers/chaining — the perf-critical path is untouched, so DNA
  performance is provably unaffected (identical code path, identical branches).
- Protein alignment does strictly *less* work per query than DNA (no second-orientation
  seed/chain/extend pass), so it carries no risk of being slower than the DNA baseline.
- Composes cleanly with ADR-0014: alphabet only affects seeding/orientation, gap-model
  only affects extension cost — no shared state, no interaction to design around.
- Known detection ambiguity: a short peptide spelled entirely in A/C/G/T-coding residues
  (Ala/Cys/Gly/Thr) is indistinguishable from DNA by content alone. `--alphabet` override
  exists for this; `auto` should require either a minimum sequence length or the presence
  of a non-ACGT amino-acid-only byte before concluding DNA, and this limitation should be
  documented, not silently wrong.
- No CLI subcommand or flag changes for existing DNA workflows — `--alphabet auto` is the
  default and behaves identically to today for DNA input.

## Alternatives considered

- **`--protein` boolean flag instead of content auto-detection**: rejected — breaks the
  existing convention (`detect_use_case`) of inferring input shape from content, and asks
  users to declare something the tool can determine itself in the overwhelming majority
  of cases.
- **Bundle BLOSUM62/PAM scoring into this same change**: rejected — conflates a small,
  safe plumbing change with a large, algorithmic scoring change (see ADR-0014's Myers vs
  WFA discussion for why substitution-matrix scoring is a materially bigger project).
  Tracked separately; not needed for alignment to work correctly in protein space.
- **Packed sub-byte amino-acid encoding, mirroring DNA's 2-bit-packing precedent**:
  rejected — there is no such precedent to preserve. `Sequence` already stores raw bytes
  for DNA; byte storage is uniformly fine at both alphabet sizes.
