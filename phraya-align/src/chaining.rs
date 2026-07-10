//! Seed chaining: collapse raw minimizer seeds into co-linear candidate regions
//! ("chains") before DP extension, so extension runs once per genuine locus instead
//! of once per raw anchor vote.
//!
//! Ported from minimap2/minibwa's chaining DP (`mb_lchain_dp` in minibwa's `lchain.c`),
//! adapted to Phraya's single-strand `Seed` list (Phraya seeds each orientation
//! separately — see `align_read` in `executor.rs` — so unlike minimap2/minibwa there is
//! no `sid`/strand field to check here; every seed passed to [`chain_seeds`] is already
//! known to share an orientation).
//!
//! This module is purely additive: it introduces `Chain`/`ChainParams`/`chain_seeds` but
//! nothing in `executor.rs` calls it yet.

use crate::seeding::Seed;

/// A co-linear run of seeds collapsed into one candidate alignment region.
///
/// Seeds are sorted by `query_pos`, strictly increasing (a valid chain never revisits a
/// query or target position). `score` is a chaining-DP score (sum of matched seed
/// lengths minus gap penalties) — a proxy for how well-anchored the candidate locus is,
/// not a DP alignment score; the actual edit-distance/CIGAR is computed later by WFA or
/// Myers extension against the target window this chain implies.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Chain {
    pub seeds: Vec<Seed>,
    pub score: i64,
}

impl Chain {
    /// Implied target start for full-query (`query_pos = 0`) extension: the first seed's
    /// target position, projected back along the chain's dominant diagonal.
    pub fn target_start(&self) -> usize {
        let first = &self.seeds[0];
        (first.target_pos as i64 - first.query_pos as i64).max(0) as usize
    }
}

/// Tunables for [`chain_seeds`], analogous to minimap2/minibwa's chaining options
/// (`mm_mapopt_t`/`mb_opt_t`). Initial values are ported from minibwa's short-read
/// preset (`options.c`); Phraya's minimizer density differs from minibwa's own seed
/// regime (canonical minimizers vs SMEMs), so these are a starting point to tune
/// empirically against the local benchmark harness, not a final answer.
#[derive(Debug, Clone, Copy)]
pub struct ChainParams {
    /// Seed length used for scoring and the diagonal-consistency check — the sketch's
    /// `k` (minimizer length), since Phraya's `Seed` doesn't carry an explicit length.
    pub seed_len: i64,
    /// Max diagonal drift `|dr - dq|` allowed between a seed and its chained
    /// predecessor (minibwa's `bw`).
    pub max_band: i64,
    /// Max target-position gap before giving up scanning further back for a
    /// predecessor (minibwa's `max_dist_x`; bounds the backward scan so it isn't
    /// O(n^2) on repeat-dense seed sets).
    pub max_dist_x: i64,
    /// Max number of backward predecessor candidates scanned per seed before the
    /// early-exit heuristic kicks in (minibwa's `max_chain_iter`).
    pub max_iter: usize,
    /// Once a seed has been dominated by `max_skip` better-scoring predecessors
    /// without itself improving the running best, stop scanning further back
    /// (minibwa's `max_chain_skip`).
    pub max_skip: usize,
    /// Per-base gap penalty (minibwa's `chn_pen_gap`, derived from `chain_gap_scale`).
    pub gap_penalty: f64,
    /// Minimum chain score to keep a chain at all (minibwa's `min_chain_score`; drops
    /// noise chains formed from a single spurious seed).
    pub min_chain_score: i64,
}

impl Default for ChainParams {
    /// minibwa short-read preset (`options.c`: `bw=100`, `max_gap=100`,
    /// `max_chain_skip=25`, `max_chain_iter=5000`, `min_chain_score=25`), with
    /// `seed_len` set to Phraya's default minimizer `k` (21) and `gap_penalty` set to
    /// minibwa's `chain_gap_scale * .01 * min_len` with `chain_gap_scale=1.0,
    /// min_len≈seed_len` collapsing to a small per-base value.
    fn default() -> Self {
        ChainParams {
            seed_len: 21,
            max_band: 100,
            max_dist_x: 100,
            max_iter: 5000,
            max_skip: 25,
            gap_penalty: 0.21,
            min_chain_score: 25,
        }
    }
}

/// Chaining score between a seed `i` and a candidate predecessor `j` (both already
/// known to be strictly increasing in query_pos, `i` after `j`).
///
/// Mirrors minibwa's `comput_sc` (`lchain.c`): requires the predecessor to also advance
/// in target position, requires the diagonal drift to stay within `max_band`, and
/// charges a linear+log gap penalty proportional to the smaller of the query/target
/// gaps. Returns `None` if `j` cannot chain before `i` at all (non-increasing target
/// position, or drift beyond the band).
fn chain_score(seed_i: &Seed, seed_j: &Seed, params: &ChainParams) -> Option<i64> {
    let dq = seed_i.query_pos as i64 - seed_j.query_pos as i64;
    if dq <= 0 {
        return None;
    }
    let dr = seed_i.target_pos as i64 - seed_j.target_pos as i64;
    if dr <= 0 || dr > params.max_dist_x + params.seed_len {
        return None;
    }
    let dd = (dr - dq).abs();
    if dd > params.max_band {
        return None;
    }
    let dg = dr.min(dq);
    let mut sc = params.seed_len.min(dg);
    if dd > 0 {
        let lin_pen = params.gap_penalty * dd as f64;
        let log_pen = 0.5 * (dd as f64 + 1.0).log2();
        sc -= (lin_pen + log_pen).round() as i64;
    }
    Some(sc)
}

/// Collapse raw seeds into co-linear chains, sorted by score descending.
///
/// `seeds` need not be pre-sorted; this sorts by `target_pos` internally (chaining
/// requires scanning in target order, unlike [`crate::seeding::find_seeds_indexed_capped`]'s
/// query-position order). Every seed in `seeds` is assumed to share one orientation —
/// Phraya seeds forward and reverse-complement bytes as two separate calls (see
/// `align_read` in `executor.rs`), so there is no strand field to check here, unlike
/// minimap2/minibwa's `sid` (target-id<<1|strand).
///
/// Ported from minibwa's `mb_lchain_dp` (`lchain.c`): for each seed in target-position
/// order, scan backward through prior seeds (bounded by `max_iter` and `max_dist_x`) for
/// the best-scoring co-linear predecessor, with a `max_skip` early-exit once a seed has
/// been dominated by that many better-scoring alternatives. Backtrack from the
/// highest-scoring endpoints to reconstruct chains; a seed used in a higher-scoring
/// chain is never reused as part of a lesser one.
pub fn chain_seeds(seeds: &[Seed], params: &ChainParams) -> Vec<Chain> {
    if seeds.is_empty() {
        return Vec::new();
    }

    let mut sorted: Vec<Seed> = seeds.to_vec();
    sorted.sort_by_key(|s| (s.target_pos, s.query_pos));
    let n = sorted.len();

    // f[i] = best chain score ending at seed i. pred[i] = predecessor seed index, or
    // None if i starts its own chain.
    let mut f: Vec<i64> = vec![0; n];
    let mut pred: Vec<Option<usize>> = vec![None; n];
    // used_as_pred[j] tracks how many times j has been dominated by a better-scoring
    // successor without itself becoming part of the running best (minibwa's `t[]`/
    // `n_skip` bookkeeping, simplified to a per-seed skip counter reset per outer i).
    let mut dominated_count: Vec<usize> = vec![0; n];

    for i in 0..n {
        let mut best_score = params.seed_len; // a chain of just seed i alone
        let mut best_pred: Option<usize> = None;
        let mut n_skip = 0usize;

        let scan_floor = i.saturating_sub(params.max_iter);
        for j in (scan_floor..i).rev() {
            if sorted[i].target_pos as i64 - sorted[j].target_pos as i64
                >= params.max_dist_x + params.seed_len
            {
                break; // matches minibwa's early positional cutoff
            }
            let Some(sc) = chain_score(&sorted[i], &sorted[j], params) else {
                continue;
            };
            let total = sc + f[j];
            if total > best_score {
                best_score = total;
                best_pred = Some(j);
                if n_skip > 0 {
                    n_skip -= 1;
                }
            } else {
                dominated_count[j] += 1;
                if dominated_count[j] > params.max_skip {
                    n_skip += 1;
                    if n_skip > params.max_skip {
                        break;
                    }
                }
            }
        }

        f[i] = best_score;
        pred[i] = best_pred;
    }

    // Backtrack from the highest-scoring endpoints. A seed already claimed by a
    // higher-scoring chain (`used[seed]`) is skipped as a starting point — this is
    // minibwa's `compact_a`/backtrack dedup, simplified since we don't need the
    // radix-sort machinery C uses for its oversized scratch buffers.
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by_key(|&i| std::cmp::Reverse(f[i]));

    let mut used = vec![false; n];
    let mut chains = Vec::new();

    for &end in &order {
        if used[end] || f[end] < params.min_chain_score {
            continue;
        }
        let mut chain_seeds_rev = Vec::new();
        let mut cur = Some(end);
        while let Some(idx) = cur {
            if used[idx] {
                break;
            }
            used[idx] = true;
            chain_seeds_rev.push(sorted[idx]);
            cur = pred[idx];
        }
        if chain_seeds_rev.is_empty() {
            continue;
        }
        chain_seeds_rev.reverse();
        chain_seeds_rev.sort_by_key(|s| s.query_pos);
        chains.push(Chain {
            seeds: chain_seeds_rev,
            score: f[end],
        });
    }

    chains.sort_by_key(|c| std::cmp::Reverse(c.score));
    chains
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Deterministic filler seeds along a clean diagonal: `target_pos = target_start +
    /// query_pos` for each of `count` seeds spaced `spacing` apart, starting at
    /// `query_start`.
    fn diagonal_seeds(target_start: usize, query_start: u32, count: u32, spacing: u32) -> Vec<Seed> {
        (0..count)
            .map(|i| {
                let qpos = query_start + i * spacing;
                Seed {
                    query_pos: qpos,
                    target_pos: target_start as u32 + qpos,
                    minimizer: 0x1000 + i as u64,
                }
            })
            .collect()
    }

    #[test]
    fn empty_seeds_produce_no_chains() {
        assert!(chain_seeds(&[], &ChainParams::default()).is_empty());
    }

    #[test]
    fn single_seed_below_min_score_produces_no_chain() {
        // seed_len=21 alone must clear min_chain_score=25 by default — a lone seed
        // scores exactly seed_len, so it's dropped unless min_chain_score <= seed_len.
        let seeds = vec![Seed { query_pos: 0, target_pos: 1000, minimizer: 1 }];
        let params = ChainParams { min_chain_score: 25, seed_len: 21, ..Default::default() };
        assert!(chain_seeds(&seeds, &params).is_empty());
    }

    #[test]
    fn single_seed_produces_single_chain_when_above_min_score() {
        let seeds = vec![Seed { query_pos: 0, target_pos: 1000, minimizer: 1 }];
        let params = ChainParams { min_chain_score: 10, seed_len: 21, ..Default::default() };
        let chains = chain_seeds(&seeds, &params);
        assert_eq!(chains.len(), 1);
        assert_eq!(chains[0].seeds.len(), 1);
        assert_eq!(chains[0].target_start(), 1000);
    }

    #[test]
    fn clean_diagonal_run_collapses_into_one_chain() {
        // 5 seeds on a clean diagonal (no indel, no drift) must chain into ONE chain
        // covering all 5 seeds, not 5 separate chains.
        let seeds = diagonal_seeds(1000, 0, 5, 40);
        let chains = chain_seeds(&seeds, &ChainParams::default());
        assert_eq!(chains.len(), 1, "clean co-linear seeds must collapse into one chain");
        assert_eq!(chains[0].seeds.len(), 5);
        assert_eq!(chains[0].target_start(), 1000);
    }

    #[test]
    fn seeds_spanning_a_short_indel_still_chain_into_one() {
        // Seeds before a 5bp deletion (target_pos = 1000 + query_pos), then seeds after
        // the deletion shifted by +5 in target (target_pos = 1005 + query_pos). This is
        // the core minibwa-parity behavior: a real short indel must not fragment a read
        // into two separate chains.
        let mut seeds = diagonal_seeds(1000, 0, 3, 40); // query 0,40,80; target 1000,1040,1080
        let mut after = diagonal_seeds(1005, 160, 3, 40); // query 160,200,240; target 1165,1205,1245
        seeds.append(&mut after);

        let chains = chain_seeds(&seeds, &ChainParams::default());
        assert_eq!(
            chains.len(),
            1,
            "seeds spanning a short indel must chain into one chain, not fragment"
        );
        assert_eq!(chains[0].seeds.len(), 6);
    }

    #[test]
    fn seeds_from_unrelated_repeat_copies_do_not_merge() {
        // Two internally-consistent diagonal runs at target positions far enough apart
        // (5000bp) that they must NOT be treated as one chain — this is the
        // repeat-copy-confusion case chaining must avoid.
        let mut seeds = diagonal_seeds(1000, 0, 4, 40);
        let mut copy2 = diagonal_seeds(6000, 0, 4, 40);
        seeds.append(&mut copy2);

        let chains = chain_seeds(&seeds, &ChainParams::default());
        assert_eq!(
            chains.len(),
            2,
            "seeds from unrelated repeat copies must produce two separate chains"
        );
        let starts: Vec<usize> = chains.iter().map(|c| c.target_start()).collect();
        assert!(starts.contains(&1000));
        assert!(starts.contains(&6000));
    }

    #[test]
    fn longer_chain_outscores_isolated_spurious_seed() {
        let mut seeds = diagonal_seeds(1000, 0, 6, 40); // a strong 6-seed chain
        seeds.push(Seed { query_pos: 500, target_pos: 9000, minimizer: 0xdead }); // isolated noise

        let chains = chain_seeds(&seeds, &ChainParams::default());
        assert!(chains.len() >= 1);
        // Highest-scoring chain must be the 6-seed run, not the isolated seed.
        assert_eq!(chains[0].seeds.len(), 6);
        assert_eq!(chains[0].target_start(), 1000);
    }

    #[test]
    fn out_of_band_drift_prevents_chaining() {
        // Second seed's diagonal drifts by 200 (target_pos jumps far more than query_pos
        // advances) — beyond max_band=100 default, so it must not chain with the first.
        let seeds = vec![
            Seed { query_pos: 0, target_pos: 1000, minimizer: 1 },
            Seed { query_pos: 40, target_pos: 1240, minimizer: 2 }, // drift = |240-40|=200
        ];
        let params = ChainParams { min_chain_score: 10, ..Default::default() };
        let chains = chain_seeds(&seeds, &params);
        assert_eq!(chains.len(), 2, "seeds beyond the band must not chain together");
    }
}
