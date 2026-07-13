//! Performance gate for issue #133 criterion: the SIMD diagonal fill must be
//! measurably faster than the scalar naive path on sequences >= 500bp.
//!
//! `#[ignore]` by default (microbenchmarks are meaningless in the debug build
//! `cargo test --all` uses). CI runs it in release with native SIMD enabled:
//!   RUSTFLAGS="-C target-cpu=native" cargo test --release --test wfa_perf -- --ignored
use phraya_align::{wfa_extend, wfa_extend_naive, SeedAnchor};
use std::time::{Duration, Instant};

fn make_pair(len: usize, div_pct: usize) -> (Vec<u8>, Vec<u8>) {
    const BASES: &[u8; 4] = b"ACGT";
    let q: Vec<u8> = (0..len).map(|i| BASES[(i * 7 + 3) % 4]).collect();
    let mut t = q.clone();
    let step = if div_pct == 0 {
        usize::MAX
    } else {
        100 / div_pct
    };
    if step != usize::MAX {
        let mut i = step;
        while i < t.len() {
            t[i] = BASES[(t[i] as usize + 1) % 4];
            i += step;
        }
    }
    (q, t)
}

/// Best-of-`runs` wall time for `f` (min reduces scheduler/cache noise).
fn best_of(runs: usize, mut f: impl FnMut()) -> Duration {
    let mut best = Duration::MAX;
    for _ in 0..runs {
        let start = Instant::now();
        f();
        best = best.min(start.elapsed());
    }
    best
}

#[test]
#[ignore = "release-only microbenchmark; run with --ignored in release"]
fn simd_faster_than_naive_ge_500bp() {
    let seed = SeedAnchor {
        query_pos: 0,
        target_pos: 0,
    };
    // Sizes >= 500bp are the issue #133 criterion. Measured ~1.2-1.5x across
    // the board on x86 (AVX2) and ARM (NEON); assert a conservative margin to
    // stay robust against CI noise while still catching a real regression.
    for &len in &[500usize, 1000, 2000, 5000] {
        let (q, t) = make_pair(len, 5);
        let naive = best_of(15, || {
            std::hint::black_box(wfa_extend_naive(&q, &t, seed).unwrap());
        });
        let simd = best_of(15, || {
            std::hint::black_box(wfa_extend(&q, &t, seed).unwrap());
        });
        let speedup = naive.as_secs_f64() / simd.as_secs_f64();
        eprintln!("len={len:>5}: naive={naive:>10.2?} simd={simd:>10.2?} speedup={speedup:.2}x");
        assert!(
            speedup >= 1.1,
            "SIMD must be measurably faster than naive at {len}bp: \
             naive={naive:?} simd={simd:?} ({speedup:.2}x, want >= 1.10x)"
        );
    }
}
