/// Benchmark for Issue #70: WFA extension (scalar)
///
/// Acceptance criterion: Benchmark 10kb sequences (measure time)
///
/// This benchmark measures the performance of the scalar WFA implementation
/// on 10kb sequence alignments to establish a baseline before SIMD optimization.
use criterion::{BenchmarkId, Criterion, black_box, criterion_group, criterion_main};
use phraya_align::wfa_extend;

/// Generate a DNA sequence of specified length with a repeating ACGT pattern
fn generate_sequence(length: usize) -> Vec<u8> {
    (0..length).map(|i| b"ACGT"[i % 4]).collect()
}

/// Generate a DNA sequence with specified divergence from a reference
fn generate_diverged_sequence(reference: &[u8], divergence_percent: usize) -> Vec<u8> {
    let mut result = reference.to_vec();
    let num_changes = (reference.len() * divergence_percent) / 100;

    for i in 0..num_changes {
        let pos = (i * reference.len()) / num_changes;
        if pos < result.len() {
            // Introduce a mismatch
            result[pos] = match result[pos] {
                b'A' => b'T',
                b'T' => b'A',
                b'C' => b'G',
                b'G' => b'C',
                _ => result[pos],
            };
        }
    }

    result
}

// ============================================================================
// ACCEPTANCE CRITERION: Benchmark 10kb alignment
// ============================================================================

fn bench_10kb_exact_match(c: &mut Criterion) {
    let query = generate_sequence(10000);
    let target = generate_sequence(10000);

    c.bench_function("wfa_scalar_10kb_exact_match", |b| {
        b.iter(|| wfa_extend(black_box(&query), black_box(&target), black_box(0)))
    });
}

fn bench_10kb_with_1pct_divergence(c: &mut Criterion) {
    let query = generate_sequence(10000);
    let target = generate_diverged_sequence(&query, 1);

    c.bench_function("wfa_scalar_10kb_1pct_divergence", |b| {
        b.iter(|| wfa_extend(black_box(&query), black_box(&target), black_box(0)))
    });
}

fn bench_10kb_with_5pct_divergence(c: &mut Criterion) {
    let query = generate_sequence(10000);
    let target = generate_diverged_sequence(&query, 5);

    c.bench_function("wfa_scalar_10kb_5pct_divergence", |b| {
        b.iter(|| wfa_extend(black_box(&query), black_box(&target), black_box(0)))
    });
}

fn bench_10kb_with_10pct_divergence(c: &mut Criterion) {
    let query = generate_sequence(10000);
    let target = generate_diverged_sequence(&query, 10);

    c.bench_function("wfa_scalar_10kb_10pct_divergence", |b| {
        b.iter(|| wfa_extend(black_box(&query), black_box(&target), black_box(0)))
    });
}

// ============================================================================
// Varying sequence lengths
// ============================================================================

fn bench_varying_lengths(c: &mut Criterion) {
    let mut group = c.benchmark_group("wfa_scalar_varying_lengths");

    for length in [100, 500, 1000, 5000, 10000].iter() {
        let query = generate_sequence(*length);
        let target = generate_diverged_sequence(&query, 1);

        group.bench_with_input(
            BenchmarkId::from_parameter(length),
            length,
            |b, &_length| {
                b.iter(|| wfa_extend(black_box(&query), black_box(&target), black_box(0)))
            },
        );
    }

    group.finish();
}

// ============================================================================
// Varying divergence levels
// ============================================================================

fn bench_varying_divergence(c: &mut Criterion) {
    let mut group = c.benchmark_group("wfa_scalar_varying_divergence");

    let query = generate_sequence(10000);

    for divergence in [0, 1, 5, 10, 20].iter() {
        let target = generate_diverged_sequence(&query, *divergence);

        group.bench_with_input(
            BenchmarkId::new("divergence", divergence),
            divergence,
            |b, &_div| b.iter(|| wfa_extend(black_box(&query), black_box(&target), black_box(0))),
        );
    }

    group.finish();
}

// ============================================================================
// Seed position benchmarks
// ============================================================================

fn bench_varying_seed_positions(c: &mut Criterion) {
    let mut group = c.benchmark_group("wfa_scalar_varying_seed_pos");

    let query = generate_sequence(10000);
    let target = generate_diverged_sequence(&query, 5);

    for seed_pos in [0, 1000, 5000, 8000].iter() {
        group.bench_with_input(
            BenchmarkId::from_parameter(seed_pos),
            seed_pos,
            |b, &pos| b.iter(|| wfa_extend(black_box(&query), black_box(&target), black_box(pos))),
        );
    }

    group.finish();
}

// ============================================================================
// Edge cases
// ============================================================================

fn bench_short_sequences(c: &mut Criterion) {
    let query = generate_sequence(100);
    let target = generate_diverged_sequence(&query, 5);

    c.bench_function("wfa_scalar_100bp", |b| {
        b.iter(|| wfa_extend(black_box(&query), black_box(&target), black_box(0)))
    });
}

fn bench_empty_sequences(c: &mut Criterion) {
    let query = b"";
    let target = b"";

    c.bench_function("wfa_scalar_empty", |b| {
        b.iter(|| wfa_extend(black_box(query), black_box(target), black_box(0)))
    });
}

// ============================================================================
// Criterion setup
// ============================================================================

criterion_group!(
    benches,
    bench_10kb_exact_match,
    bench_10kb_with_1pct_divergence,
    bench_10kb_with_5pct_divergence,
    bench_10kb_with_10pct_divergence,
    bench_varying_lengths,
    bench_varying_divergence,
    bench_varying_seed_positions,
    bench_short_sequences,
    bench_empty_sequences,
);

criterion_main!(benches);
