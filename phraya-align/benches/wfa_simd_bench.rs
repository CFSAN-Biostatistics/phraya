use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use phraya_align::wfa_simd::{count_matching_prefix, count_matching_prefix_scalar};
use phraya_align::{wfa_extend, wfa_extend_naive, SeedAnchor};

/// Micro-benchmark the WFA inner-loop primitive: scalar byte-by-byte vs the
/// architecture-dispatched SIMD longest-common-prefix. Run with
/// `RUSTFLAGS="-C target-feature=+sse2" cargo bench -p phraya-align count_matching_prefix`.
fn bench_count_matching_prefix(c: &mut Criterion) {
    let mut group = c.benchmark_group("count_matching_prefix");
    // Long matching runs with a single trailing mismatch — the common extend case.
    for run_len in [32usize, 256, 4096] {
        let a: Vec<u8> = (0..run_len).map(|i| b"ACGT"[i % 4]).collect();
        let mut b = a.clone();
        *b.last_mut().unwrap() ^= 0x01; // force the mismatch at the final byte

        group.bench_with_input(
            BenchmarkId::new("scalar", run_len),
            &(&a, &b),
            |bn, (a, b)| bn.iter(|| count_matching_prefix_scalar(black_box(a), black_box(b))),
        );
        group.bench_with_input(
            BenchmarkId::new("simd", run_len),
            &(&a, &b),
            |bn, (a, b)| bn.iter(|| count_matching_prefix(black_box(a), black_box(b))),
        );
    }
    group.finish();
}

// Test will fail: functions do not exist yet
fn bench_10kb_alignment_naive(c: &mut Criterion) {
    // Generate a 10kb query sequence
    let query: Vec<u8> = (0..10_000)
        .map(|i| match i % 4 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        })
        .collect();

    // Target with ~95% identity (5% divergence)
    let mut target = query.clone();
    for i in (0..target.len()).step_by(20) {
        if i < target.len() {
            target[i] = match target[i] {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                _ => b'A',
            };
        }
    }

    let seed = SeedAnchor {
        query_pos: 0,
        target_pos: 0,
    };

    c.bench_function("wfa_10kb_naive", |b| {
        b.iter(|| {
            wfa_extend_naive(
                black_box(&query),
                black_box(&target),
                black_box(seed.clone()),
            )
        });
    });
}

// Test will fail: functions do not exist yet
fn bench_10kb_alignment_simd(c: &mut Criterion) {
    // Generate a 10kb query sequence
    let query: Vec<u8> = (0..10_000)
        .map(|i| match i % 4 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        })
        .collect();

    // Target with ~95% identity (5% divergence)
    let mut target = query.clone();
    for i in (0..target.len()).step_by(20) {
        if i < target.len() {
            target[i] = match target[i] {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                _ => b'A',
            };
        }
    }

    let seed = SeedAnchor {
        query_pos: 0,
        target_pos: 0,
    };

    c.bench_function("wfa_10kb_simd", |b| {
        b.iter(|| {
            wfa_extend(
                black_box(&query),
                black_box(&target),
                black_box(seed.clone()),
            )
        });
    });
}

// Test will fail: functions do not exist yet
fn bench_comparison(c: &mut Criterion) {
    let mut group = c.benchmark_group("wfa_10kb_comparison");

    // Generate a 10kb query sequence
    let query: Vec<u8> = (0..10_000)
        .map(|i| match i % 4 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        })
        .collect();

    // Target with ~95% identity (5% divergence)
    let mut target = query.clone();
    for i in (0..target.len()).step_by(20) {
        if i < target.len() {
            target[i] = match target[i] {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                _ => b'A',
            };
        }
    }

    let seed = SeedAnchor {
        query_pos: 0,
        target_pos: 0,
    };

    group.bench_with_input(
        BenchmarkId::new("naive", "10kb"),
        &(&query, &target, &seed),
        |b, (q, t, s)| {
            b.iter(|| wfa_extend_naive(black_box(q), black_box(t), black_box((*s).clone())))
        },
    );

    group.bench_with_input(
        BenchmarkId::new("simd", "10kb"),
        &(&query, &target, &seed),
        |b, (q, t, s)| b.iter(|| wfa_extend(black_box(q), black_box(t), black_box((*s).clone()))),
    );

    group.finish();
}

// Test will fail: functions do not exist yet
fn bench_varying_sizes(c: &mut Criterion) {
    let mut group = c.benchmark_group("wfa_varying_sizes");

    for size in [100, 500, 1_000, 5_000, 10_000].iter() {
        let query: Vec<u8> = (0..*size)
            .map(|i| match i % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            })
            .collect();

        let mut target = query.clone();
        for i in (0..target.len()).step_by(20) {
            if i < target.len() {
                target[i] = match target[i] {
                    b'A' => b'T',
                    b'C' => b'G',
                    b'G' => b'C',
                    _ => b'A',
                };
            }
        }

        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        group.bench_with_input(
            BenchmarkId::new("naive", size),
            &(&query, &target, &seed),
            |b, (q, t, s)| {
                b.iter(|| wfa_extend_naive(black_box(q), black_box(t), black_box((*s).clone())))
            },
        );

        group.bench_with_input(
            BenchmarkId::new("simd", size),
            &(&query, &target, &seed),
            |b, (q, t, s)| {
                b.iter(|| wfa_extend(black_box(q), black_box(t), black_box((*s).clone())))
            },
        );
    }

    group.finish();
}

// Test will fail: functions do not exist yet
fn bench_varying_divergence(c: &mut Criterion) {
    let mut group = c.benchmark_group("wfa_varying_divergence");

    let query: Vec<u8> = (0..10_000)
        .map(|i| match i % 4 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        })
        .collect();

    // Test different divergence levels: 1%, 5%, 10%, 20%
    for divergence in [1, 5, 10, 20].iter() {
        let mut target = query.clone();
        let step = 100 / divergence;
        for i in (0..target.len()).step_by(step) {
            if i < target.len() {
                target[i] = match target[i] {
                    b'A' => b'T',
                    b'C' => b'G',
                    b'G' => b'C',
                    _ => b'A',
                };
            }
        }

        let seed = SeedAnchor {
            query_pos: 0,
            target_pos: 0,
        };

        group.bench_with_input(
            BenchmarkId::new("naive", format!("{}%", divergence)),
            &(&query, &target, &seed),
            |b, (q, t, s)| {
                b.iter(|| wfa_extend_naive(black_box(q), black_box(t), black_box((*s).clone())))
            },
        );

        group.bench_with_input(
            BenchmarkId::new("simd", format!("{}%", divergence)),
            &(&query, &target, &seed),
            |b, (q, t, s)| {
                b.iter(|| wfa_extend(black_box(q), black_box(t), black_box((*s).clone())))
            },
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_count_matching_prefix,
    bench_10kb_alignment_naive,
    bench_10kb_alignment_simd,
    bench_comparison,
    bench_varying_sizes,
    bench_varying_divergence
);
criterion_main!(benches);
