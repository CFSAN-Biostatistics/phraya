/// Regression test for CSP2 spike finding B4: `~100kb` contig-scale extension must be
/// memory-bounded, not explore an unbounded wavefront search that OOMs on divergent or
/// unrelated pairs.
///
/// Root cause (traced in `phraya-align/src/wfa_simd.rs`): Case 4 (contigs-only)
/// compares two *comparable-length* sequences, which both the linear
/// (`fill_wfa_fitting_impl`) and affine (`fill_wfa_affine_generic`) engines route
/// through their "global" alignment mode (`tn <= qn + qn/2 + 10` always holds for
/// equal-length inputs). The linear global path (`fill_wfa`) had no cap parameter at
/// all — it explored up to `s = qn + tn` wavefront phases, each allocating an
/// `O(qn + tn)` array retained for traceback, before any cap was checked. Production
/// entry points (`wfa_extend`, `wfa_extend_affine`) also never supplied a cap in the
/// first place. Combined: any divergent/unrelated ~100kb pair could allocate tens of
/// GB. Fixed by making `fill_wfa` itself enforce a cap during its loop, and giving
/// every production entry point a derived default cap (`default_max_s_cap`).
///
/// This test proves the fix via *bounded time*, not by trying to provoke an actual OOM
/// (which would be a hazardous, non-portable thing to assert on in CI): before the fix,
/// the unrelated-pair case below would explore the full uncapped wavefront search;
/// after the fix, it must abandon (`Err`) well within a couple of seconds.
use phraya_align::{wfa_extend, wfa_extend_affine, wfa_simd::AffineCosts, SeedAnchor};
use std::time::{Duration, Instant};

/// Deterministic pseudo-random DNA (LCG). Independent seeds are ~75% divergent from
/// each other (matches the documented convention in
/// `phraya-cli/tests/issue_210_fixture_tests.rs`).
fn diverse_dna(len: usize, seed: u64) -> Vec<u8> {
    let mut state = seed;
    let bases = [b'A', b'C', b'G', b'T'];
    (0..len)
        .map(|_| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            bases[((state >> 32) as usize) % 4]
        })
        .collect()
}

/// Mutated copy of `base` with `mutation_rate` fraction of positions substituted.
fn mutate_sequence(base: &[u8], mutation_rate: f64, seed: u64) -> Vec<u8> {
    let mut x = seed;
    let bases = [b'A', b'C', b'G', b'T'];
    base.iter()
        .map(|&b| {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            let rand_01 = ((x >> 32) as u32) as f64 / (u32::MAX as f64 + 1.0);
            if rand_01 < mutation_rate {
                // Filter out the current base *before* indexing, matching
                // `issue_210_fixture_tests.rs`'s convention — indexing a fixed-size
                // "keep retrying" loop with an already-consumed random draw can pick
                // the same excluded value forever (found via B4's own regression test
                // hanging on a 40kb+ fixture).
                let candidates: Vec<u8> = bases.iter().copied().filter(|&c| c != b).collect();
                let idx = ((x >> 16) as usize) % candidates.len();
                candidates[idx]
            } else {
                b
            }
        })
        .collect()
}

const CONTIG_LEN: usize = 100_000;
const BOUND: Duration = Duration::from_secs(5);
const ORIGIN: SeedAnchor = SeedAnchor { query_pos: 0, target_pos: 0 };

/// A ~0.05%-diverged pair (closely related isolate contigs, well within realistic
/// bacterial-outbreak divergence) aligns successfully through the production linear
/// entry point, bounded in time.
#[test]
fn findings_b4_related_100kb_pair_aligns_via_linear_engine() {
    let contig_a = diverse_dna(CONTIG_LEN, 1);
    let contig_b = mutate_sequence(&contig_a, 0.0005, 2);

    let start = Instant::now();
    let result = wfa_extend(&contig_a, &contig_b, ORIGIN);
    let elapsed = start.elapsed();

    assert!(
        result.is_ok(),
        "a closely-related 100kb pair should align successfully, got: {:?}",
        result
    );
    assert!(
        elapsed < BOUND,
        "related-pair alignment took {:?}, expected well under {:?}",
        elapsed,
        BOUND
    );
}

/// Two independent (~75%-diverged, unrelated) ~100kb sequences must abandon quickly via
/// the production linear entry point, not explore an unbounded wavefront search.
#[test]
fn findings_b4_unrelated_100kb_pair_abandons_via_linear_engine() {
    let contig_a = diverse_dna(CONTIG_LEN, 10);
    let contig_b = diverse_dna(CONTIG_LEN, 20); // independent seed, no real homology

    let start = Instant::now();
    let result = wfa_extend(&contig_a, &contig_b, ORIGIN);
    let elapsed = start.elapsed();

    assert!(
        result.is_err(),
        "an unrelated 100kb pair should abandon (Err), not succeed with a fabricated \
         near-total-edit-distance alignment"
    );
    assert!(
        elapsed < BOUND,
        "unrelated-pair abandonment took {:?}, expected well under {:?} — an uncapped \
         wavefront search would take vastly longer (and allocate far more memory) \
         before exhausting s = qn + tn phases",
        elapsed,
        BOUND
    );
}

/// The same bound applies to the affine engine (`Strategy::Sensitive`'s default path).
#[test]
fn findings_b4_unrelated_100kb_pair_abandons_via_affine_engine() {
    let contig_a = diverse_dna(CONTIG_LEN, 30);
    let contig_b = diverse_dna(CONTIG_LEN, 40);

    let start = Instant::now();
    let result = wfa_extend_affine(&contig_a, &contig_b, ORIGIN, AffineCosts::default());
    let elapsed = start.elapsed();

    assert!(
        result.is_err(),
        "an unrelated 100kb pair should abandon (Err) via the affine engine too"
    );
    assert!(
        elapsed < BOUND,
        "unrelated-pair abandonment (affine) took {:?}, expected well under {:?}",
        elapsed,
        BOUND
    );
}

/// A ~0.05%-diverged pair also aligns successfully through the affine engine.
#[test]
fn findings_b4_related_100kb_pair_aligns_via_affine_engine() {
    let contig_a = diverse_dna(CONTIG_LEN, 50);
    let contig_b = mutate_sequence(&contig_a, 0.0005, 60);

    let start = Instant::now();
    let result = wfa_extend_affine(&contig_a, &contig_b, ORIGIN, AffineCosts::default());
    let elapsed = start.elapsed();

    assert!(
        result.is_ok(),
        "a closely-related 100kb pair should align successfully via the affine engine, \
         got: {:?}",
        result
    );
    assert!(
        elapsed < BOUND,
        "related-pair alignment (affine) took {:?}, expected well under {:?}",
        elapsed,
        BOUND
    );
}
