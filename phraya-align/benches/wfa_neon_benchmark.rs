//! Benchmark for NEON SIMD WFA diagonal fill
//!
//! Acceptance criteria: ≥1.5× speedup on ARM64 vs naive for 10kb alignment

use criterion::{BenchmarkId, Criterion, black_box, criterion_group, criterion_main};
use phraya_align::wfa_neon::wfa_diagonal_fill_neon;

/// Benchmark NEON implementation on 10kb sequences
fn bench_neon_10kb(c: &mut Criterion) {
    let query = vec![b'A'; 10_000];
    let target = vec![b'A'; 10_000];
    let prev_wavefront = vec![0; 10_000];

    c.bench_function("neon_10kb", |b| {
        b.iter(|| {
            let mut diagonal = vec![0; 10_000];
            wfa_diagonal_fill_neon(
                black_box(&mut diagonal),
                black_box(&prev_wavefront),
                black_box(&query),
                black_box(&target),
            );
        });
    });
}

/// Benchmark NEON implementation with various sequence sizes
fn bench_neon_sizes(c: &mut Criterion) {
    let mut group = c.benchmark_group("neon_by_size");

    for size in [100, 1_000, 5_000, 10_000, 50_000].iter() {
        let query = vec![b'A'; *size];
        let target = vec![b'A'; *size];
        let prev_wavefront = vec![0; *size];

        group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, &_size| {
            b.iter(|| {
                let mut diagonal = vec![0; *size];
                wfa_diagonal_fill_neon(
                    black_box(&mut diagonal),
                    black_box(&prev_wavefront),
                    black_box(&query),
                    black_box(&target),
                );
            });
        });
    }
    group.finish();
}

/// Benchmark NEON vs naive (requires naive implementation from #5)
fn bench_neon_vs_naive(c: &mut Criterion) {
    let query = vec![b'A'; 10_000];
    let target = vec![b'A'; 10_000];
    let prev_wavefront = vec![0; 10_000];

    let mut group = c.benchmark_group("neon_vs_naive");

    group.bench_function("neon", |b| {
        b.iter(|| {
            let mut diagonal = vec![0; 10_000];
            wfa_diagonal_fill_neon(
                black_box(&mut diagonal),
                black_box(&prev_wavefront),
                black_box(&query),
                black_box(&target),
            );
        });
    });

    // Naive implementation benchmark - will fail until #5 is implemented
    // group.bench_function("naive", |b| {
    //     b.iter(|| {
    //         let mut diagonal = vec![0; 10_000];
    //         wfa_diagonal_fill_naive(
    //             black_box(&mut diagonal),
    //             black_box(&prev_wavefront),
    //             black_box(&query),
    //             black_box(&target),
    //         );
    //     });
    // });

    group.finish();
}

/// Benchmark with mismatches (more realistic)
fn bench_neon_with_mismatches(c: &mut Criterion) {
    let query = vec![b'A'; 10_000];
    let mut target = vec![b'A'; 10_000];
    // Introduce 1% mismatches
    for i in (0..10_000).step_by(100) {
        target[i] = b'T';
    }
    let prev_wavefront = vec![0; 10_000];

    c.bench_function("neon_10kb_with_mismatches", |b| {
        b.iter(|| {
            let mut diagonal = vec![0; 10_000];
            wfa_diagonal_fill_neon(
                black_box(&mut diagonal),
                black_box(&prev_wavefront),
                black_box(&query),
                black_box(&target),
            );
        });
    });
}

criterion_group!(
    benches,
    bench_neon_10kb,
    bench_neon_sizes,
    bench_neon_vs_naive,
    bench_neon_with_mismatches
);
criterion_main!(benches);
