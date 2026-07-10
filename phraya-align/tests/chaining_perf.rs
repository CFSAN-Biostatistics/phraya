//! Performance comparison: legacy raw-vote anchor selection vs seed-chaining, on a
//! realistic-scale synthetic reference and read batch.
//!
//! `#[ignore]` by default (microbenchmarks are meaningless in the debug build
//! `cargo test --all` uses). Run in release with native SIMD enabled:
//!   RUSTFLAGS="-C target-cpu=native" cargo test --release --test chaining_perf -- --ignored --nocapture

use phraya_align::executor::{align_read, AlignConfig, Strategy, TargetContext};
use phraya_core::types::Sequence;
use phraya_io::plan::{PhrayaPlan, UseCase};
use std::collections::HashMap;
use std::time::Instant;

fn make_plan() -> PhrayaPlan {
    PhrayaPlan::new(
        UseCase::ReadsWithRef,
        vec!["test".to_string()],
        "2026-07-10T00:00:00Z".to_string(),
        HashMap::new(),
        HashMap::new(),
        vec![],
    )
}

fn diverse_dna(len: usize, seed: u64) -> Vec<u8> {
    let mut x = seed;
    (0..len)
        .map(|_| {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            b"ACGT"[((x >> 33) & 3) as usize]
        })
        .collect()
}

/// Simulate a read from `target` at `pos`, with roughly `divergence` fraction of bases
/// mutated (a crude per-base SNP injector, not a full simulator — enough for a wall-time
/// comparison, not for accuracy scoring).
fn simulate_read(target: &[u8], pos: usize, len: usize, divergence: f64, seed: u64) -> Vec<u8> {
    let mut read = target[pos..pos + len].to_vec();
    let mut x = seed;
    for base in read.iter_mut() {
        x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        if (x >> 40) as f64 / (1u64 << 24) as f64 % 1.0 < divergence {
            *base = b"ACGT"[((x >> 33) & 3) as usize];
        }
    }
    read
}

#[test]
#[ignore = "release-only microbenchmark; run with --ignored in release"]
fn chained_anchors_faster_than_legacy_on_realistic_batch() {
    const GENOME_LEN: usize = 200_000;
    const NUM_READS: usize = 2_000;
    const READ_LEN: usize = 150;
    const DIVERGENCE: f64 = 0.01;

    let target_bases = diverse_dna(GENOME_LEN, 42);
    let target = Sequence::new(target_bases.clone(), None, "ref".to_string(), None);
    let plan = make_plan();

    let reads: Vec<Sequence> = (0..NUM_READS)
        .map(|i| {
            let pos = (i * 37) % (GENOME_LEN - READ_LEN);
            let bases = simulate_read(&target_bases, pos, READ_LEN, DIVERGENCE, i as u64 + 1);
            Sequence::new(bases, None, format!("read{i}"), None)
        })
        .collect();

    let legacy_config = AlignConfig::new(Strategy::Balanced).with_legacy_anchors(true);
    let chained_config = AlignConfig::new(Strategy::Balanced).with_legacy_anchors(false);

    // Build TargetContext once per config, matching real usage (phraya-cli builds it
    // once per worker and reuses it across every read in the batch — see
    // run_align_worker_with_plan in phraya-cli/src/main.rs). Rebuilding it per-read
    // (as align_task_with_config's one-off convenience wrapper does) would re-run
    // detect_tandem_repeats' O(target) scan every iteration, an identical fixed cost
    // for both configs that would swamp the actual per-read anchor-selection signal
    // this test exists to measure.
    let legacy_ctx = TargetContext::build(&target, &plan, legacy_config.strategy);
    let chained_ctx = TargetContext::build(&target, &plan, chained_config.strategy);

    let t0 = Instant::now();
    let mut legacy_placed = 0;
    for read in &reads {
        if align_read(&legacy_ctx, read, &plan, &legacy_config, None).is_some() {
            legacy_placed += 1;
        }
    }
    let legacy_elapsed = t0.elapsed();

    let t1 = Instant::now();
    let mut chained_placed = 0;
    for read in &reads {
        if align_read(&chained_ctx, read, &plan, &chained_config, None).is_some() {
            chained_placed += 1;
        }
    }
    let chained_elapsed = t1.elapsed();

    println!(
        "legacy: {:?} ({} placed) | chained: {:?} ({} placed) | speedup: {:.2}x",
        legacy_elapsed,
        legacy_placed,
        chained_elapsed,
        chained_placed,
        legacy_elapsed.as_secs_f64() / chained_elapsed.as_secs_f64().max(1e-9)
    );

    // Placement counts should be close (chaining shouldn't lose real placements); a
    // generous tolerance since exact equality isn't required (see the differential
    // test module's documented divergence for multi-mapping cases).
    let diff = (legacy_placed as i64 - chained_placed as i64).abs();
    assert!(
        diff <= (NUM_READS / 20) as i64,
        "chained placement count should be close to legacy's: legacy={} chained={}",
        legacy_placed,
        chained_placed
    );

    assert!(
        chained_elapsed < legacy_elapsed,
        "chaining should be at least as fast as legacy raw-vote anchor selection on a \
         Balanced-strategy batch (legacy extends up to 6 independent anchors per \
         orientation; chaining should extend far fewer per genuine locus)"
    );
}
