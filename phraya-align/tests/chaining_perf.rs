//! Performance gate for the chaining redesign (ADR-0012): `balanced`-strategy alignment
//! throughput on a realistic-scale synthetic reference and read batch stays under a
//! generous wall-time ceiling. Not a comparison against the removed legacy raw-vote path
//! (validated once on the full HPC benchmark ladder during the redesign, see
//! `~/idiot-binfie/benchmarking-5.md`) — this is a plain throughput floor to catch a
//! future regression in the anchor-selection hot path.
//!
//! `#[ignore]` by default (microbenchmarks are meaningless in the debug build
//! `cargo test --all` uses). Run in release with native SIMD enabled:
//!   RUSTFLAGS="-C target-cpu=native" cargo test --release --test chaining_perf -- --ignored --nocapture

use phraya_align::executor::{align_read, AlignConfig, Strategy, TargetContext};
use phraya_core::types::Sequence;
use phraya_io::plan::{PhrayaPlan, UseCase};
use std::collections::HashMap;
use std::time::{Duration, Instant};

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
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"ACGT"[((x >> 33) & 3) as usize]
        })
        .collect()
}

/// Simulate a read from `target` at `pos`, with roughly `divergence` fraction of bases
/// mutated (a crude per-base SNP injector, not a full simulator — enough for a wall-time
/// floor test, not for accuracy scoring).
fn simulate_read(target: &[u8], pos: usize, len: usize, divergence: f64, seed: u64) -> Vec<u8> {
    let mut read = target[pos..pos + len].to_vec();
    let mut x = seed;
    for base in read.iter_mut() {
        x = x
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        if (x >> 40) as f64 / (1u64 << 24) as f64 % 1.0 < divergence {
            *base = b"ACGT"[((x >> 33) & 3) as usize];
        }
    }
    read
}

#[test]
#[ignore = "release-only microbenchmark; run with --ignored in release"]
fn balanced_strategy_throughput_stays_under_ceiling() {
    const GENOME_LEN: usize = 200_000;
    const NUM_READS: usize = 2_000;
    const READ_LEN: usize = 150;
    const DIVERGENCE: f64 = 0.01;
    // Generous ceiling: on the reference machine this batch ran in ~3.2s post-chaining
    // (down from ~4.4s pre-chaining on the same synthetic scale). 10s leaves ample
    // headroom for slower CI/dev hardware while still catching a real regression.
    const WALL_CEILING: Duration = Duration::from_secs(10);

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

    let config = AlignConfig::new(Strategy::Balanced);

    // Build TargetContext once and reuse across every read, matching real usage
    // (phraya-cli builds it once per worker — see run_align_worker_with_plan in
    // phraya-cli/src/main.rs). Rebuilding it per-read would re-run
    // detect_tandem_repeats' O(target) scan every iteration, swamping the per-read
    // anchor-selection signal this test exists to measure.
    let ctx = TargetContext::build(&target, &plan, config.strategy);

    let t0 = Instant::now();
    let mut placed = 0;
    for read in &reads {
        if align_read(&ctx, read, &plan, &config, None).is_some() {
            placed += 1;
        }
    }
    let elapsed = t0.elapsed();

    println!(
        "balanced: {:?} ({} / {} placed)",
        elapsed, placed, NUM_READS
    );

    assert!(
        placed as f64 / NUM_READS as f64 > 0.9,
        "at least 90% of reads should place on this low-divergence synthetic batch, got {}/{}",
        placed,
        NUM_READS
    );
    assert!(
        elapsed < WALL_CEILING,
        "balanced-strategy throughput regressed: {:?} exceeds the {:?} ceiling",
        elapsed,
        WALL_CEILING
    );
}
