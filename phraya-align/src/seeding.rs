use phraya_core::types::MinimizerSketch;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// A minimizer index: value -> the target positions where it occurs.
///
/// Positions live in a single flat `Vec<u32>`, grouped by value, with `slots` mapping a
/// value to its `(start, len)` span. This replaces a `HashMap<u64, Vec<u32>>` that
/// allocated one `Vec` per distinct minimizer — ~18M allocations for one vertebrate
/// chromosome, measured at 0.146-0.213 s per Mb of reference and the dominant cost of
/// `align`'s per-space setup.
#[derive(Debug, Default)]
pub struct MinimizerIndex {
    slots: HashMap<u64, (u32, u32)>,
    positions: Vec<u32>,
}

impl MinimizerIndex {
    /// Group `target`'s minimizers by value.
    ///
    /// Within a value, positions keep the order they appear in the sketch — *not* ascending
    /// order. That is load-bearing: it is the order the old per-value `Vec` produced, so
    /// seeds, chains and the resulting `.phraya` bytes are unchanged. Sketch positions are
    /// not ascending to begin with (`simd-minimizers` emits them interleaved across its 8
    /// parallel streams, and duplicate positions occur), so this cannot be reconstructed by
    /// sorting.
    pub fn build(target: &MinimizerSketch) -> Self {
        // One entry per distinct value; sized for the all-distinct case to avoid rehashing.
        let mut slots: HashMap<u64, (u32, u32)> =
            HashMap::with_capacity(target.minimizers.len());

        // Pass 1: occurrence count per value.
        for &(val, _) in &target.minimizers {
            slots.entry(val).or_insert((0, 0)).1 += 1;
        }

        // Prefix-sum: give each value a contiguous span of `positions`.
        let mut cursor = 0u32;
        for (start, len) in slots.values_mut() {
            *start = cursor;
            cursor += *len;
        }

        // Pass 2: place each position at its value's next free slot, walking the sketch in
        // its original order so per-value order is preserved. `start` doubles as the write
        // cursor here and is rewound below.
        let mut positions = vec![0u32; cursor as usize];
        for &(val, pos) in &target.minimizers {
            let slot = slots.get_mut(&val).expect("counted in pass 1");
            positions[slot.0 as usize] = pos;
            slot.0 += 1;
        }
        for (start, len) in slots.values_mut() {
            *start -= *len;
        }

        MinimizerIndex { slots, positions }
    }

    /// Target positions for `value`, or an empty slice when it does not occur.
    pub fn positions(&self, value: u64) -> &[u32] {
        match self.slots.get(&value) {
            Some(&(start, len)) => &self.positions[start as usize..(start + len) as usize],
            None => &[],
        }
    }

    /// Occurrence count of each distinct value, in unspecified order.
    pub fn occurrence_counts(&self) -> impl Iterator<Item = usize> + '_ {
        self.slots.values().map(|&(_, len)| len as usize)
    }

    /// Number of distinct minimizer values.
    pub fn len(&self) -> usize {
        self.slots.len()
    }

    pub fn is_empty(&self) -> bool {
        self.slots.is_empty()
    }
}

/// A seed: a shared minimizer between query and target that anchors WFA extension.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct Seed {
    pub query_pos: u32,
    pub target_pos: u32,
    pub minimizer: u64,
}

/// Find shared minimizer seeds between two sketches, sorted by query position.
///
/// Convenience for one-off pairs. When aligning many queries against a single
/// target, build a [`MinimizerIndex`] from the target once with
/// [`MinimizerIndex::build`] and call [`find_seeds_indexed`] per query instead —
/// this function rebuilds the target-side index on every call.
pub fn find_seeds(query: &MinimizerSketch, target: &MinimizerSketch) -> Vec<Seed> {
    find_seeds_indexed(query, &MinimizerIndex::build(target))
}

/// Find shared minimizer seeds against a prebuilt target [`MinimizerIndex`],
/// sorted by query position. Equivalent to [`find_seeds`] but without rebuilding
/// the target-side map, so O(query minimizers) per call instead of O(target).
pub fn find_seeds_indexed(query: &MinimizerSketch, index: &MinimizerIndex) -> Vec<Seed> {
    find_seeds_indexed_capped(query, index, usize::MAX)
}

/// Like [`find_seeds_indexed`], but skips any minimizer whose target occurrence count
/// exceeds `max_occ` (repeat masking).
///
/// A minimizer that occurs thousands of times in the target (a homopolymer or
/// microsatellite k-mer in an AT-rich genome) contributes thousands of near-useless
/// seeds — an O(occurrences) blow-up in seed generation and diagonal voting that can
/// stall alignment (issue #194). Such hyper-frequent minimizers also carry almost no
/// positional information: the read's true locus is still anchored by its rarer
/// minimizers. Dropping them bounds the work with negligible recall cost — only a read
/// lying *entirely* within a hyper-repeat (no rarer minimizer to anchor on) is lost, and
/// that read is genuinely unmappable to a unique locus. `max_occ = usize::MAX` disables
/// masking (identical to [`find_seeds_indexed`]).
pub fn find_seeds_indexed_capped(
    query: &MinimizerSketch,
    index: &MinimizerIndex,
    max_occ: usize,
) -> Vec<Seed> {
    let mut seeds = Vec::new();
    for &(val, qpos) in &query.minimizers {
        let tposs = index.positions(val);
        if tposs.len() > max_occ {
            continue; // hyper-frequent minimizer: mask it
        }
        for &tpos in tposs {
            seeds.push(Seed {
                query_pos: qpos,
                target_pos: tpos,
                minimizer: val,
            });
        }
    }
    seeds.sort_by_key(|s| s.query_pos);
    seeds
}

/// Multiplier on the median minimizer occurrence used to set the repeat-masking cap.
const SEED_CAP_MEDIAN_MULT: usize = 8;

/// Choose a repeat-masking occurrence cap from a target [`MinimizerIndex`]'s own
/// occurrence-count distribution (a plan-phase redundancy signal).
///
/// Self-normalizing: returns `max(floor, SEED_CAP_MEDIAN_MULT × median occurrence)`. The
/// median is robust to the heavy tail (unlike a high percentile, which on few-distinct-value
/// indices lands *on* the outlier we mean to mask). On a clean genome the median is 1, so the
/// `floor` governs and nothing is masked; on a genome so repetitive that its *typical*
/// minimizer already occurs many times, the cap lifts so masking stays proportionate rather
/// than stripping real signal. `floor` is the no-op guard for clean/small genomes.
pub fn seed_occurrence_cap(index: &MinimizerIndex, floor: usize) -> usize {
    if index.is_empty() {
        return floor;
    }
    let mut counts: Vec<usize> = index.occurrence_counts().collect();
    counts.sort_unstable();
    let median = counts[counts.len() / 2];
    floor.max(median.saturating_mul(SEED_CAP_MEDIAN_MULT))
}

#[cfg(test)]
mod tests {
    use super::*;
    use phraya_core::types::MinimizerSketch;

    fn sketch(mins: &[(u64, u32)]) -> MinimizerSketch {
        MinimizerSketch { minimizers: mins.to_vec(), k: 21, w: 11 }
    }

    #[test]
    fn capped_masks_hyperfrequent_minimizer() {
        // Minimizer value 7 occurs 5× in the target; value 9 occurs once.
        let target = sketch(&[(7, 10), (7, 20), (7, 30), (7, 40), (7, 50), (9, 100)]);
        let index = MinimizerIndex::build(&target);
        let query = sketch(&[(7, 0), (9, 5)]);

        // Uncapped: both contribute (5 + 1 = 6 seeds).
        assert_eq!(find_seeds_indexed_capped(&query, &index, usize::MAX).len(), 6);
        // Cap at 4: the 5×-occurring value 7 is masked, only value 9's single seed remains.
        let capped = find_seeds_indexed_capped(&query, &index, 4);
        assert_eq!(capped.len(), 1);
        assert_eq!(capped[0].minimizer, 9);
    }

    #[test]
    fn cap_is_a_noop_on_a_clean_index() {
        // Every minimizer unique → percentile is 1 → cap == floor → nothing maskable.
        let mins: Vec<(u64, u32)> = (0..1000u64).map(|v| (v, v as u32)).collect();
        let index = MinimizerIndex::build(&sketch(&mins));
        assert_eq!(seed_occurrence_cap(&index, 256), 256);
    }

    #[test]
    fn floor_governs_when_median_is_low_and_hyperrepeat_is_maskable() {
        // 999 unique values + one value occurring 10_000×. Median occurrence is 1, so the
        // floor (256) governs and the hyper-repeat (10_000 > 256) is maskable — the outlier
        // does not pull the cap up to shelter itself.
        let mut mins: Vec<(u64, u32)> = (0..999u64).map(|v| (v, v as u32)).collect();
        mins.extend(std::iter::repeat((9999u64, 0u32)).take(10_000));
        let index = MinimizerIndex::build(&sketch(&mins));
        let cap = seed_occurrence_cap(&index, 256);
        assert_eq!(cap, 256);
        assert!(10_000 > cap, "the hyper-repeat must exceed the cap and be maskable");
    }

    #[test]
    fn cap_lifts_when_typical_minimizer_is_repetitive() {
        // Every minimizer occurs 100× (a uniformly repetitive genome). Median is 100, so the
        // cap lifts to 8×100 = 800, above the floor — masking stays proportionate instead of
        // stripping the genome's normal signal.
        let mut mins: Vec<(u64, u32)> = Vec::new();
        for v in 0..500u64 {
            mins.extend(std::iter::repeat((v, 0u32)).take(100));
        }
        let index = MinimizerIndex::build(&sketch(&mins));
        assert_eq!(seed_occurrence_cap(&index, 256), 800);
    }

    #[test]
    fn empty_index_returns_floor() {
        let index = MinimizerIndex::build(&sketch(&[]));
        assert_eq!(seed_occurrence_cap(&index, 256), 256);
    }

    /// Positions within a value keep sketch order, not ascending order. Seeds are emitted in
    /// this order and then stably sorted by query position, so reordering here would change
    /// `.phraya` bytes.
    #[test]
    fn build_preserves_sketch_order_within_a_value() {
        let index = MinimizerIndex::build(&sketch(&[(7, 30), (7, 10), (9, 5), (7, 20)]));
        assert_eq!(index.positions(7), &[30, 10, 20]);
        assert_eq!(index.positions(9), &[5]);
        assert!(index.positions(11).is_empty());
    }
}
