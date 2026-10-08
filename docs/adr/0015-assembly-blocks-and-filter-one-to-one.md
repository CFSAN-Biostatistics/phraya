# 15. Assembly alignment blocks, with one-to-one selection in the filter

- **Status**: Accepted
- **Date**: 2026-10-08

## Context

CSP2 wants to replace MUMmer `dnadiff` with `phraya plan → align → filter --format snpdiffs`.
On a real pair of draft assemblies (`phraya-csp2-snp-discrepancy.md`) Phraya recovered 2.7% of
dnadiff's SNPs, aligned 6.6% of the reference against dnadiff's 91%, and emitted noise from
alignments under 60% identity. Three causes sat in align and one in snpdiffs:

- **A1.** One whole-contig WFA per (query, space) is capped at
  `512 MiB / (15·(q+t+1))` edits, so a homologous contig longer than about 50–60 kb at
  0.5–0.8% divergence was abandoned.
- **A2.** `select_anchors` projects every chain to `query_pos 0` and extends the *whole*
  query; only the primary alignment produces variants. A query contig that overlaps a reference
  contig only partially, or spans several (the normal case for two draft assemblies), cannot be
  represented.
- **A3.** `fill_wfa_fitting_impl` falls back to global alignment when the target window is at
  most `1.5·|q| + 10`, which the end of a reference contig always causes, so a contained contig is
  charged an edit for every reference base it does not cover.
- **B.** `.phraya` stored observations and a coverage track but no alignment records, so
  snpdiffs could not build per-alignment BED rows, query coverage, block identity or one-to-one
  selection. It emitted one aggregate row per contig and divided reference coverage by one
  query contig's length (`Query_Percent_Aligned` 4090%).

dnadiff's model is: exact matches, clustered, extended into several local alignment blocks per
contig pair, then `delta-filter -1` keeps each reference and query base in at most one block.

## Decision

1. **`Strategy::Assembly` (`--strategy assembly`)** aligns contigs as local blocks. Every chain,
   in either orientation, becomes one block (or several, if a gap cannot be filled):
   - its seeds are verified byte-for-byte (canonical minimizers also hit the opposite strand)
     and grown into maximal exact segments;
   - consecutive segments are joined by a global WFA of the gap between them;
   - blocks under 100 bp or 80% identity are dropped.

   Chaining already breaks chains at gaps over ~121 bp and at diagonal drift over 100, so no
   separate splitting step is needed. There is no anchor selection, no orientation winner and no
   `(0,0)` fallback. Variants come from every kept block. Align still picks no winner across
   blocks or spaces (glossary: cross-reference superposition), so
   `align({A,B}) = align({A}) ∪ align({B})` holds.

   It is a strategy rather than something inferred from a `UseCase` or from query length: the
   caller states the alignment shape (like `long-read`). The plan does not freeze it, and other
   strategies' results do not change.

2. **`.phraya` v2 stores the blocks** (`PhrayaFile::alignments: Vec<AlignmentBlock>`):
   - each block has the query ID, its length, the query interval in forward coordinates, the
     reference interval, strand, edit distance and CIGAR;
   - the field is trailing and `#[serde(default)]`, so v1 files still read;
   - blocks are not attached to `VariantObservation`: its trailing `skip_serializing_if` field
     (`mate_info`) makes any new optional field shift MessagePack slots. Instead, an observation
     belongs to the block with the same query and strand whose reference and strand-local query
     intervals contain it (`AlignmentBlock::contains_observation`);
   - only `Strategy::Assembly` emits blocks, so read-alignment outputs do not grow.

3. **The filter decides** (`phraya_filter::blocks::select_one_to_one`):
   - blocks from all input spaces are taken greedily by `aligned columns × identity`, with a
     deterministic tie-break;
   - a block overlapping already-claimed reference or query bases is cut into its maximal
     unclaimed runs, and each run is kept if it still passes the floors.

   `--format snpdiffs` is computed from the kept blocks:
   - one BED row per block, with real query contig and interval;
   - `.` rows for uncovered regions;
   - variant rows only for observations inside a kept block, carrying that block's
     coordinates and identity, with query positions in forward coordinates;
   - header coverage as kept bases over assembly bases;
   - medians over kept blocks.

   Inputs with observations but no blocks are rejected with a message to re-run with
   `--strategy assembly`; the previous aggregate output is not kept as a fallback.

## Consequences

- On synthetic assembly pairs (identical breakpoints, shifted breakpoints,
  reverse-complemented query contigs, decoys, a duplicated region) callable SNP recall is 1.00
  with zero false positives and zero duplicated reference positions. Under `balanced` it was
  0.44 and 0.00.
- Each block loses about k+w bp at its ends, because there is no end extension yet: a SNP that
  close to a contig end or homology break is not reported. dnadiff trims terminal mismatches
  too. If the real-data comparison shows edge loss above 0.5% of SNPs, the next step is an
  X-drop extension from the outermost segments.
- Greedy selection approximates `delta-filter -1` (the intersection of its `-r` and `-q`
  weighted-LIS filters). Where the two disagree, it is in repeated regions such as rRNA
  operons, which copy gets a query region. If the P4 gate against dnadiff fails there, port
  MUMmer's LIS.
- `phraya merge` unions blocks, but snpdiffs is defined over per-space files from palette mode
  (`phraya align --reference`), whose `sample_id` is the reference contig.
- Structural counts (breakpoints, relocations, ...) and k-mer fields stay `0`/`NA`.

## Alternatives considered

- **Switch every strategy to block extension for long queries.** No CLI change for CSP2, but it
  silently changes long-read and contig results everywhere. Rejected in favour of an explicit
  strategy.
- **Key the mode on `UseCase::ContigsOnly`.** Freezes input shape into the plan, against the
  glossary.
- **Pick blocks one-to-one in align.** Makes align encode a decision and breaks composability
  across spaces. Selection belongs in the filter.
- **A block index on `VariantObservation`.** Rejected because of the MessagePack field-order
  hazard described above.
