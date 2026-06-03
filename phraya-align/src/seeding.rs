use phraya_core::types::MinimizerSketch;
use serde::{Deserialize, Serialize};

/// A seed: a shared minimizer between query and target that anchors WFA extension.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct Seed {
    pub query_pos: u32,
    pub target_pos: u32,
    pub minimizer: u64,
}

/// Find shared minimizer seeds between two sketches, sorted by query position.
pub fn find_seeds(query: &MinimizerSketch, target: &MinimizerSketch) -> Vec<Seed> {
    query
        .find_shared_minimizers(target)
        .into_iter()
        .map(|(m, qp, tp)| Seed {
            query_pos: qp,
            target_pos: tp,
            minimizer: m,
        })
        .collect()
}
