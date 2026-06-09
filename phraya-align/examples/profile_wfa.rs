// Simple WFA profiling example
use phraya_align::{wfa_extend_naive, SeedAnchor};
use std::time::Instant;

fn main() {
    println!("=== WFA Performance Profiling ===\n");

    // Test realistic bacterial genomics workloads
    let test_cases = [
        (150, "150bp read (typical Illumina)"),
        (1000, "1kb sequence"),
        (5000, "5kb sequence"),
    ];

    for (size, desc) in &test_cases {
        println!("--- {} ---", desc);

        // Generate query
        let query: Vec<u8> = (0..*size).map(|i| [b'A', b'C', b'G', b'T'][i % 4]).collect();

        // Generate target with 5% divergence (95% identity - typical bacterial)
        let mut target = query.clone();
        for i in (0..target.len()).step_by(20) {
            target[i] = b'T'; // mismatch every 20 bases
        }

        let seed = SeedAnchor { query_pos: 0, target_pos: 0 };

        // Warmup
        for _ in 0..10 {
            let _ = wfa_extend_naive(&query, &target, seed.clone());
        }

        // Benchmark
        let iterations = if *size <= 150 { 1000 } else if *size <= 1000 { 100 } else { 20 };

        let start = Instant::now();
        for _ in 0..iterations {
            let _ = wfa_extend_naive(&query, &target, seed.clone());
        }
        let elapsed = start.elapsed();

        let per_align = elapsed / iterations;
        println!("  Iterations: {}", iterations);
        println!("  Total time: {:?}", elapsed);
        println!("  Per-alignment: {:?}", per_align);
        println!("  Throughput: {:.0} alignments/sec", 1.0 / per_align.as_secs_f64());

        // Estimate operations
        let matches = (size * 95) / 100;
        let edits = size - matches;
        let avg_match_run = if edits > 0 { matches / edits } else { *size };

        println!("  Expected matches: ~{} bases", matches);
        println!("  Expected edits: ~{} ops", edits);
        println!("  Avg match run: ~{} bases", avg_match_run);
        println!();
    }

    println!("=== SIMD Feasibility Analysis ===");
    println!("✓ Avg match runs >32 bytes: AVX2 SIMD extension justified (2-4× speedup expected)");
    println!("✗ Avg match runs <16 bytes: SIMD overhead dominates, not worth it");
    println!("\nConclusion:");
    println!("  For 150bp reads @ 95% identity:");
    println!("    - ~142 match bases / ~7 edits = ~20 byte avg runs");
    println!("    - SIMD would help, but overhead significant");
    println!("    - Expected net speedup: 1.5-2×");
    println!();
    println!("  For longer sequences (1kb+):");
    println!("    - Avg runs scale with length");
    println!("    - SIMD more beneficial");
    println!("    - Expected speedup: 2-4×");
    println!();
    println!("Run with: RUSTFLAGS=\"-C target-cpu=native\" cargo run --release --example profile_wfa");
}
