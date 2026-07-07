use phraya_core::types::MinimizerSketch;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// A minimizer index: value → target positions where it occurs. Built once from a
/// target sketch and reused across many queries (see [`MinimizerIndex`]).
pub type MinimizerIndex = HashMap<u64, Vec<u32>>;

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
/// [`build_minimizer_index`] and call [`find_seeds_indexed`] per query instead —
/// this function rebuilds the target-side hash map on every call.
pub fn find_seeds(query: &MinimizerSketch, target: &MinimizerSketch) -> Vec<Seed> {
    find_seeds_indexed(query, &build_minimizer_index(target))
}

/// Build a reusable minimizer index from a target sketch: value → target positions.
///
/// This is the per-target work that [`find_seeds`] otherwise repeats on every call.
/// Hoist it out of a per-query loop and pass the result to [`find_seeds_indexed`].
pub fn build_minimizer_index(target: &MinimizerSketch) -> MinimizerIndex {
    let mut index: MinimizerIndex = HashMap::new();
    for &(val, pos) in &target.minimizers {
        index.entry(val).or_default().push(pos);
    }
    index
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
        if let Some(tposs) = index.get(&val) {
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
    let mut counts: Vec<usize> = index.values().map(|v| v.len()).collect();
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
        let mut index: MinimizerIndex = HashMap::new();
        index.insert(7, vec![10, 20, 30, 40, 50]);
        index.insert(9, vec![100]);
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
        let mut index: MinimizerIndex = HashMap::new();
        for v in 0..1000u64 {
            index.insert(v, vec![v as u32]);
        }
        assert_eq!(seed_occurrence_cap(&index, 256), 256);
    }

    #[test]
    fn floor_governs_when_median_is_low_and_hyperrepeat_is_maskable() {
        // 999 unique values + one value occurring 10_000×. Median occurrence is 1, so the
        // floor (256) governs and the hyper-repeat (10_000 > 256) is maskable — the outlier
        // does not pull the cap up to shelter itself.
        let mut index: MinimizerIndex = HashMap::new();
        for v in 0..999u64 {
            index.insert(v, vec![v as u32]);
        }
        index.insert(9999, vec![0u32; 10_000]);
        let cap = seed_occurrence_cap(&index, 256);
        assert_eq!(cap, 256);
        assert!(10_000 > cap, "the hyper-repeat must exceed the cap and be maskable");
    }

    #[test]
    fn cap_lifts_when_typical_minimizer_is_repetitive() {
        // Every minimizer occurs 100× (a uniformly repetitive genome). Median is 100, so the
        // cap lifts to 8×100 = 800, above the floor — masking stays proportionate instead of
        // stripping the genome's normal signal.
        let mut index: MinimizerIndex = HashMap::new();
        for v in 0..500u64 {
            index.insert(v, vec![0u32; 100]);
        }
        assert_eq!(seed_occurrence_cap(&index, 256), 800);
    }

    #[test]
    fn empty_index_returns_floor() {
        let index: MinimizerIndex = HashMap::new();
        assert_eq!(seed_occurrence_cap(&index, 256), 256);
    }
}
