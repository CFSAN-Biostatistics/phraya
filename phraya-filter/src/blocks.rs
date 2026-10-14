//! One-to-one selection of `--strategy assembly` alignment blocks (ADR-0015).
//!
//! Align surfaces every homologous block of a query contig in every reference space
//! (cross-reference superposition). For a dnadiff-style comparison each reference base and
//! each query base must be explained by at most one alignment — MUMmer's `delta-filter -1`.
//! [`select_one_to_one`] does this greedily: blocks are taken by descending weight (aligned
//! columns × identity); a block that overlaps bases already claimed on either axis is cut
//! into its maximal runs of unclaimed alignment columns, and each run is kept while it still
//! spans [`MIN_BLOCK_LEN`] reference bases at [`MIN_BLOCK_IDENTITY`].

use phraya_core::cigar::CigarStats;
use phraya_core::types::{AlignmentBlock, Strand, VariantObservation};
use std::collections::HashMap;
use std::sync::Arc;

/// Minimum reference span of a kept block (bp); the same floor `Strategy::Assembly` applies.
pub const MIN_BLOCK_LEN: u32 = 100;

/// Minimum identity of a kept block; the same floor `Strategy::Assembly` applies.
pub const MIN_BLOCK_IDENTITY: f64 = 0.80;

/// The blocks of one reference space (one `.phraya` file).
pub struct SpaceBlocks<'a> {
    pub ref_contig: &'a str,
    pub ref_len: u32,
    pub blocks: &'a [AlignmentBlock],
}

/// A block kept by one-to-one selection, possibly trimmed to its unclaimed part.
#[derive(Debug, Clone, PartialEq)]
pub struct KeptBlock {
    /// Index of the reference space in the [`select_one_to_one`] input.
    pub space: usize,
    pub block: AlignmentBlock,
    /// `block.identity()`, cached.
    pub identity: f64,
}

/// The kept blocks and, per reference space, which kept block owns each reference base.
pub struct OneToOne {
    /// In input-space order, then by reference start.
    pub kept: Vec<KeptBlock>,
    ref_owner: Vec<Vec<u32>>,
    query_claimed: HashMap<Arc<str>, Vec<bool>>,
}

const UNOWNED: u32 = u32::MAX;

impl OneToOne {
    /// The kept block `obs` (from reference space `space`) belongs to, if its reference base
    /// survived selection and the block it was extracted from is the one that owns it.
    pub fn block_for(&self, space: usize, obs: &VariantObservation) -> Option<&KeptBlock> {
        let owner = *self.ref_owner.get(space)?.get(obs.position() as usize)?;
        let kept = self.kept.get(owner as usize)?;
        kept.block.contains_observation(obs).then_some(kept)
    }

    /// Reference bases owned by a kept block, over all spaces.
    pub fn reference_bases_covered(&self) -> u64 {
        self.ref_owner.iter().flatten().filter(|&&o| o != UNOWNED).count() as u64
    }

    /// Query bases covered by a kept block, over all query contigs.
    pub fn query_bases_covered(&self) -> u64 {
        self.query_claimed.values().flatten().filter(|&&c| c).count() as u64
    }

    /// Maximal uncovered `[start, end)` runs of reference space `space`.
    pub fn uncovered_reference(&self, space: usize) -> Vec<(u32, u32)> {
        runs(self.ref_owner.get(space).map(|v| v.iter().map(|&o| o == UNOWNED).collect()).unwrap_or_default())
    }

    /// Maximal uncovered `[start, end)` runs of a query contig of length `query_len`.
    pub fn uncovered_query(&self, query: &str, query_len: u32) -> Vec<(u32, u32)> {
        match self.query_claimed.get(query) {
            Some(claimed) => runs(claimed.iter().map(|&c| !c).collect()),
            None if query_len > 0 => vec![(0, query_len)],
            None => Vec::new(),
        }
    }
}

fn runs(flags: Vec<bool>) -> Vec<(u32, u32)> {
    let mut out = Vec::new();
    let mut start = None;
    for (i, &f) in flags.iter().enumerate() {
        match (f, start) {
            (true, None) => start = Some(i as u32),
            (false, Some(s)) => {
                out.push((s, i as u32));
                start = None;
            }
            _ => {}
        }
    }
    if let Some(s) = start {
        out.push((s, flags.len() as u32));
    }
    out
}

/// One alignment column: the reference base and strand-local query base it consumes.
#[derive(Clone, Copy)]
struct Column {
    t: Option<u32>,
    q: Option<u32>,
    op: char,
}

fn columns(block: &AlignmentBlock) -> Vec<Column> {
    let mut out = Vec::new();
    let (mut t, mut q) = (block.t_start, block.strand_local_query_range().start);
    for (n, op) in phraya_core::cigar::parse_ops(&block.cigar) {
        for _ in 0..n {
            let col = match op {
                'M' | 'X' => {
                    let c = Column { t: Some(t), q: Some(q), op };
                    t += 1;
                    q += 1;
                    c
                }
                // This codebase's alphabet: 'I' = target-extra, 'D' = query-extra.
                'I' => {
                    let c = Column { t: Some(t), q: None, op };
                    t += 1;
                    c
                }
                'D' => {
                    let c = Column { t: None, q: Some(q), op };
                    q += 1;
                    c
                }
                _ => continue,
            };
            out.push(col);
        }
    }
    out
}

/// Rebuild a block from a run of its columns (which starts and ends on an M/X column).
fn sub_block(parent: &AlignmentBlock, run: &[Column]) -> AlignmentBlock {
    let mut cigar = String::new();
    let mut edit_distance = 0u32;
    let mut i = 0;
    while i < run.len() {
        let op = run[i].op;
        let n = run[i..].iter().take_while(|c| c.op == op).count();
        cigar.push_str(&format!("{n}{op}"));
        if op != 'M' {
            edit_distance += n as u32;
        }
        i += n;
    }
    let t_start = run.iter().find_map(|c| c.t).unwrap();
    let t_end = run.iter().rev().find_map(|c| c.t).unwrap() + 1;
    let local_start = run.iter().find_map(|c| c.q).unwrap();
    let local_end = run.iter().rev().find_map(|c| c.q).unwrap() + 1;
    let (q_start, q_end) = match parent.strand {
        Strand::Forward => (local_start, local_end),
        Strand::Reverse => (parent.query_len - local_end, parent.query_len - local_start),
    };
    AlignmentBlock {
        q_start,
        q_end,
        t_start,
        t_end,
        edit_distance,
        cigar: cigar.into(),
        ..parent.clone()
    }
}

fn passes_floors(block: &AlignmentBlock) -> bool {
    block.t_end - block.t_start >= MIN_BLOCK_LEN && block.identity() >= MIN_BLOCK_IDENTITY
}

/// Assign every reference base and every query base to at most one block.
pub fn select_one_to_one(spaces: &[SpaceBlocks<'_>]) -> OneToOne {
    let mut candidates: Vec<(usize, &AlignmentBlock, f64)> = spaces
        .iter()
        .enumerate()
        .flat_map(|(i, s)| s.blocks.iter().map(move |b| (i, b)))
        .filter(|(_, b)| passes_floors(b))
        .map(|(i, b)| {
            let weight = CigarStats::parse(&b.cigar).alignment_columns() as f64 * b.identity();
            (i, b, weight)
        })
        .collect();
    // Heaviest first; ties broken by position so the result is input-order independent.
    candidates.sort_by(|(ia, a, wa), (ib, b, wb)| {
        wb.total_cmp(wa)
            .then_with(|| spaces[*ia].ref_contig.cmp(spaces[*ib].ref_contig))
            .then_with(|| a.t_start.cmp(&b.t_start))
            .then_with(|| a.query_id.cmp(&b.query_id))
            .then_with(|| a.q_start.cmp(&b.q_start))
            .then_with(|| (a.strand == Strand::Reverse).cmp(&(b.strand == Strand::Reverse)))
    });

    let mut ref_owner: Vec<Vec<u32>> = spaces.iter().map(|s| vec![UNOWNED; s.ref_len as usize]).collect();
    let mut query_claimed: HashMap<Arc<str>, Vec<bool>> = HashMap::new();
    let mut kept: Vec<KeptBlock> = Vec::new();

    for (space, block, _) in candidates {
        let qlen = block.query_len as usize;
        let cols = columns(block);
        // A block's own columns never repeat a base, so freedom is decided once up front.
        let free: Vec<bool> = {
            let owner = &ref_owner[space];
            let claimed = query_claimed.get(&block.query_id);
            cols.iter()
                .map(|c| {
                    c.t.is_none_or(|t| owner.get(t as usize) == Some(&UNOWNED))
                        && c.q.is_none_or(|q| {
                            let fwd = block.forward_query_position(q) as usize;
                            claimed.is_none_or(|v| !v[fwd])
                        })
                })
                .collect()
        };
        let mut i = 0;
        while i < cols.len() {
            if !free[i] {
                i += 1;
                continue;
            }
            let start = i;
            while i < cols.len() && free[i] {
                i += 1;
            }
            // A kept run begins and ends on an aligned (M/X) column.
            let run = &cols[start..i];
            let Some(first) = run.iter().position(|c| c.t.is_some() && c.q.is_some()) else { continue };
            let last = run.iter().rposition(|c| c.t.is_some() && c.q.is_some()).unwrap();
            let run = &run[first..=last];
            let piece = sub_block(block, run);
            if !passes_floors(&piece) {
                continue;
            }
            let id = kept.len() as u32;
            let claimed = query_claimed.entry(block.query_id.clone()).or_insert_with(|| vec![false; qlen]);
            for c in run {
                if let Some(t) = c.t {
                    ref_owner[space][t as usize] = id;
                }
                if let Some(q) = c.q {
                    claimed[block.forward_query_position(q) as usize] = true;
                }
            }
            let identity = piece.identity();
            kept.push(KeptBlock { space, block: piece, identity });
        }
    }

    // Report in input-space order, then by reference start; renumber the owner map to match.
    let mut order: Vec<usize> = (0..kept.len()).collect();
    order.sort_by_key(|&i| (kept[i].space, kept[i].block.t_start, kept[i].block.q_start));
    let mut renumber = vec![0u32; kept.len()];
    for (new, &old) in order.iter().enumerate() {
        renumber[old] = new as u32;
    }
    for owners in &mut ref_owner {
        for o in owners.iter_mut().filter(|o| **o != UNOWNED) {
            *o = renumber[*o as usize];
        }
    }
    let mut slots: Vec<Option<KeptBlock>> = kept.into_iter().map(Some).collect();
    let kept = order.iter().map(|&i| slots[i].take().unwrap()).collect();
    OneToOne { kept, ref_owner, query_claimed }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    /// A gap-free block of `len` columns with mismatches at the given column offsets.
    fn block(query: &str, query_len: u32, q_start: u32, t_start: u32, len: u32, strand: Strand, mismatches: &[u32]) -> AlignmentBlock {
        let mut cigar = String::new();
        let mut prev = 0;
        for &m in mismatches {
            if m > prev {
                cigar.push_str(&format!("{}M", m - prev));
            }
            cigar.push_str("1X");
            prev = m + 1;
        }
        if len > prev {
            cigar.push_str(&format!("{}M", len - prev));
        }
        AlignmentBlock {
            query_id: query.into(),
            query_len,
            q_start,
            q_end: q_start + len,
            t_start,
            t_end: t_start + len,
            strand,
            edit_distance: mismatches.len() as u32,
            cigar: cigar.into(),
        }
    }

    fn space<'a>(ref_contig: &'a str, ref_len: u32, blocks: &'a [AlignmentBlock]) -> SpaceBlocks<'a> {
        SpaceBlocks { ref_contig, ref_len, blocks }
    }

    fn spans(r: &OneToOne) -> Vec<(usize, &str, u32, u32, u32, u32)> {
        r.kept
            .iter()
            .map(|k| (k.space, &*k.block.query_id, k.block.t_start, k.block.t_end, k.block.q_start, k.block.q_end))
            .collect()
    }

    #[test]
    fn disjoint_blocks_are_kept_whole() {
        let a = [block("q", 10_000, 0, 0, 1_000, Strand::Forward, &[])];
        let b = [block("q", 10_000, 1_000, 500, 2_000, Strand::Forward, &[10])];
        let r = select_one_to_one(&[space("r0", 1_000, &a), space("r1", 3_000, &b)]);
        assert_eq!(spans(&r), vec![(0, "q", 0, 1_000, 0, 1_000), (1, "q", 500, 2_500, 1_000, 3_000)]);
        assert_eq!(&*r.kept[1].block.cigar, "10M1X1989M");
        assert_eq!((r.reference_bases_covered(), r.query_bases_covered()), (3_000, 3_000));
        assert_eq!(r.uncovered_reference(1), vec![(0, 500), (2_500, 3_000)]);
        assert_eq!(r.uncovered_query("q", 10_000), vec![(3_000, 10_000)]);
        assert_eq!(r.uncovered_query("unaligned", 700), vec![(0, 700)]);
    }

    #[test]
    fn reference_overlap_trims_the_lighter_block() {
        let blocks = [
            block("long", 5_000, 0, 0, 2_000, Strand::Forward, &[]),
            block("short", 5_000, 0, 1_500, 1_000, Strand::Forward, &[]),
        ];
        let r = select_one_to_one(&[space("r", 3_000, &blocks)]);
        assert_eq!(spans(&r), vec![(0, "long", 0, 2_000, 0, 2_000), (0, "short", 2_000, 2_500, 500, 1_000)]);
        assert_eq!(r.reference_bases_covered(), 2_500, "every reference base counted once");
    }

    #[test]
    fn query_bases_aligned_to_two_spaces_are_claimed_once() {
        // One query region matching two reference copies keeps only the better copy.
        let copy1 = [block("q", 3_000, 0, 0, 1_000, Strand::Forward, &[])];
        let copy2 = [block("q", 3_000, 0, 0, 1_000, Strand::Forward, &[5, 50, 500])];
        let r = select_one_to_one(&[space("copy1", 1_000, &copy1), space("copy2", 1_000, &copy2)]);
        assert_eq!(spans(&r), vec![(0, "q", 0, 1_000, 0, 1_000)]);
        assert_eq!(r.uncovered_reference(1), vec![(0, 1_000)]);
    }

    #[test]
    fn reverse_strand_trim_maps_query_coordinates() {
        let blocks = [
            block("a", 5_000, 0, 0, 2_000, Strand::Forward, &[]),
            // Forward query 1000..2000 on the reverse strand: strand-local 3000..4000.
            block("b", 5_000, 1_000, 1_500, 1_000, Strand::Reverse, &[]),
        ];
        let r = select_one_to_one(&[space("r", 3_000, &blocks)]);
        // b loses reference 1500..2000 (its first 500 columns, strand-local 3000..3500), so it
        // keeps strand-local 3500..4000, which is forward 1000..1500.
        assert_eq!(spans(&r), vec![(0, "a", 0, 2_000, 0, 2_000), (0, "b", 2_000, 2_500, 1_000, 1_500)]);
        assert_eq!(r.kept[1].block.strand_local_query_range(), 3_500..4_000);
    }

    #[test]
    fn remnants_below_the_length_floor_are_dropped() {
        let blocks = [
            block("long", 5_000, 0, 0, 2_000, Strand::Forward, &[]),
            block("short", 5_000, 0, 1_000, 1_050, Strand::Forward, &[]),
        ];
        let r = select_one_to_one(&[space("r", 3_000, &blocks)]);
        assert_eq!(spans(&r), vec![(0, "long", 0, 2_000, 0, 2_000)], "a 50 bp remnant is not a block");
    }

    #[test]
    fn blocks_below_the_floors_are_ignored() {
        let low_identity = block("low", 5_000, 0, 0, 1_000, Strand::Forward, &(0..300).map(|i| i * 3).collect::<Vec<_>>());
        let too_short = block("short", 5_000, 0, 2_000, 80, Strand::Forward, &[]);
        let r = select_one_to_one(&[space("r", 3_000, &[low_identity, too_short])]);
        assert!(r.kept.is_empty());
        assert_eq!(r.reference_bases_covered(), 0);
    }

    #[test]
    fn equal_weight_ties_do_not_depend_on_input_order() {
        let q1 = block("q1", 1_000, 0, 0, 1_000, Strand::Forward, &[]);
        let q2 = block("q2", 1_000, 0, 0, 1_000, Strand::Forward, &[]);
        let forward = select_one_to_one(&[space("r", 1_000, &[q1.clone(), q2.clone()])]);
        let backward = select_one_to_one(&[space("r", 1_000, &[q2, q1])]);
        assert_eq!(spans(&forward), vec![(0, "q1", 0, 1_000, 0, 1_000)]);
        assert_eq!(spans(&forward), spans(&backward));
    }

    #[test]
    fn observations_map_to_the_block_that_owns_their_base() {
        let blocks = [block("a", 2_000, 0, 0, 1_000, Strand::Forward, &[500])];
        let r = select_one_to_one(&[space("r", 1_000, &blocks)]);
        let obs = |provenance: &str, query_position: u32| {
            VariantObservation::new(500, b'A', HashMap::from([(b'C', 1)]), 1.0, "500M1X499M".into(), 60, 1, vec![1], 60.0, provenance.into())
                .with_query_position(query_position)
        };
        assert!(r.block_for(0, &obs("a", 500)).is_some());
        assert!(r.block_for(0, &obs("b", 500)).is_none(), "another query's variant at the same base");
        assert!(r.block_for(0, &obs("a", 1_500)).is_none(), "the same query aligned elsewhere");
        assert!(r.block_for(0, &obs("a", 500).with_strand(Strand::Reverse)).is_none(), "other strand");
    }
}
