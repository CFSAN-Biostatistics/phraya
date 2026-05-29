use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use phraya_core::CoverageTrack;

fn benchmark_compression_ratio(c: &mut Criterion) {
    let mut group = c.benchmark_group("coverage_track_compression");

    // Test 1: Uniform coverage (best case - single RLE run)
    let uniform_coverage = vec![30; 4_600_000]; // E. coli genome size at 30x
    group.bench_with_input(
        BenchmarkId::new("uniform", "4.6Mbp_30x"),
        &uniform_coverage,
        |b, cov| {
            b.iter(|| {
                let track = CoverageTrack::from_coverage(black_box(cov.clone()));
                black_box(track)
            })
        },
    );

    // Test 2: Realistic coverage with variation (multiple regions)
    let mut realistic_coverage = Vec::new();
    realistic_coverage.extend(vec![30; 1_000_000]); // 1Mbp at 30x
    realistic_coverage.extend(vec![10; 500_000]); // 500kbp at 10x
    realistic_coverage.extend(vec![35; 2_000_000]); // 2Mbp at 35x
    realistic_coverage.extend(vec![0; 100_000]); // 100kbp no coverage
    realistic_coverage.extend(vec![30; 1_000_000]); // 1Mbp at 30x

    group.bench_with_input(
        BenchmarkId::new("realistic", "4.6Mbp_mixed"),
        &realistic_coverage,
        |b, cov| {
            b.iter(|| {
                let track = CoverageTrack::from_coverage(black_box(cov.clone()));
                black_box(track)
            })
        },
    );

    // Test 3: High variation coverage (worst case - many runs)
    let mut high_variation = Vec::new();
    for i in 0..100_000 {
        high_variation.push(if i % 2 == 0 { 10 } else { 20 });
    }

    group.bench_with_input(
        BenchmarkId::new("high_variation", "100k_alternating"),
        &high_variation,
        |b, cov| {
            b.iter(|| {
                let track = CoverageTrack::from_coverage(black_box(cov.clone()));
                black_box(track)
            })
        },
    );

    // Test 4: Random coverage (realistic noise)
    fn simple_rand(seed: u64) -> u64 {
        seed.wrapping_mul(6364136223846793005).wrapping_add(1)
    }

    let mut random_coverage = Vec::new();
    let mut seed = 42u64;
    for _ in 0..1_000_000 {
        seed = simple_rand(seed);
        let coverage = 20 + (seed % 20) as usize; // 20-40x with variation
        random_coverage.push(coverage);
    }

    group.bench_with_input(
        BenchmarkId::new("random", "1Mbp_20-40x"),
        &random_coverage,
        |b, cov| {
            b.iter(|| {
                let track = CoverageTrack::from_coverage(black_box(cov.clone()));
                black_box(track)
            })
        },
    );

    group.finish();
}

fn benchmark_random_access(c: &mut Criterion) {
    let mut group = c.benchmark_group("coverage_track_random_access");

    // Build a realistic track
    let mut coverage = Vec::new();
    coverage.extend(vec![30; 1_000_000]);
    coverage.extend(vec![10; 500_000]);
    coverage.extend(vec![35; 2_000_000]);
    coverage.extend(vec![0; 100_000]);
    coverage.extend(vec![30; 1_000_000]);

    let track = CoverageTrack::from_coverage(coverage.clone());

    // Benchmark random access via binary search
    group.bench_function("coverage_at_position", |b| {
        b.iter(|| {
            let pos = black_box(1_234_567);
            let cov = track.coverage_at(pos);
            black_box(cov)
        })
    });

    // Benchmark sequential access via iterator
    group.bench_function("iterate_all_positions", |b| {
        b.iter(|| {
            let sum: usize = track.iter().map(|(_, cov)| cov).sum();
            black_box(sum)
        })
    });

    group.finish();
}

fn benchmark_decompression(c: &mut Criterion) {
    let mut group = c.benchmark_group("coverage_track_decompression");

    // Build a realistic track
    let mut coverage = Vec::new();
    coverage.extend(vec![30; 1_000_000]);
    coverage.extend(vec![10; 500_000]);
    coverage.extend(vec![35; 2_000_000]);

    let track = CoverageTrack::from_coverage(coverage.clone());

    // Benchmark full decompression
    group.bench_function("to_vec_full_decompression", |b| {
        b.iter(|| {
            let decompressed = track.to_vec();
            black_box(decompressed)
        })
    });

    group.finish();
}

criterion_group!(
    benches,
    benchmark_compression_ratio,
    benchmark_random_access,
    benchmark_decompression
);
criterion_main!(benches);
