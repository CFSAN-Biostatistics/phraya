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
    let mut seeds = Vec::new();
    for &(val, qpos) in &query.minimizers {
        if let Some(tposs) = index.get(&val) {
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
