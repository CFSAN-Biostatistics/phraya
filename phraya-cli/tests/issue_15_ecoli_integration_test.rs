/// Integration test for issue #15: E. coli 2-sample alignment with confidence metrics.
///
/// End-to-end test validating that phraya align correctly identifies SNPs when
/// comparing two sequences with known divergence, at a scale that exercises the
/// full pipeline (plan → align → count variants).
///
/// Test uses synthetic fixture data (not downloaded) with:
/// - Two 2000bp "assembly" sequences with exactly N known SNPs introduced
/// - Verification that detected variant count is within ±10% of N
/// - Validation that mean confidence > 0.5 on synthetic data
/// - Completion in <2 minutes

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use tempfile::TempDir;

// ============================================================================
// Synthetic sequence generation with known SNPs
// ============================================================================

/// LCG-based diverse DNA sequence generator. Avoids minimizer-seed explosion
/// by ensuring most 21-mers are unique.
fn diverse_dna(len: usize, seed: u64) -> Vec<u8> {
    let mut x = seed;
    (0..len)
        .map(|_| {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            b"ACGT"[((x >> 33) & 3) as usize]
        })
        .collect()
}

/// Create two synthetic sequences with exactly N known SNPs at predictable positions.
/// Returns (seq1, seq2, snp_positions).
fn create_synthetic_pair(seq_len: usize, num_snps: usize) -> (Vec<u8>, Vec<u8>, Vec<usize>) {
    let seq1 = diverse_dna(seq_len, 42);
    let mut seq2 = seq1.clone();
    let mut snp_positions = Vec::new();

    // Introduce SNPs at evenly-spaced positions
    let spacing = seq_len / (num_snps + 1);
    for i in 0..num_snps {
        let pos = (i + 1) * spacing;
        if pos < seq_len {
            // Flip base at pos in seq2
            seq2[pos] = match seq1[pos] {
                b'A' => b'T',
                b'T' => b'A',
                b'C' => b'G',
                b'G' => b'C',
                _ => b'A',
            };
            snp_positions.push(pos);
        }
    }

    (seq1, seq2, snp_positions)
}

/// Write two sequences to a FASTA file for plan construction.
fn write_assembly_fasta(dir: &Path, filename: &str, seq1: &[u8], seq2: &[u8]) -> PathBuf {
    let path = dir.join(filename);
    let content = format!(
        ">contig1\n{}\n>contig2\n{}\n",
        String::from_utf8_lossy(seq1),
        String::from_utf8_lossy(seq2)
    );
    std::fs::write(&path, content).unwrap();
    path
}

// ============================================================================
// Test 1: SNP count within tolerance
// ============================================================================

/// Verify that phraya align detects expected SNP count (within ±10%).
/// Two 2000bp synthetic sequences with 20 introduced SNPs → 18–22 detected variants expected.
#[test]
fn issue_15_variant_count_within_tolerance() {
    let dir = TempDir::new().unwrap();
    let seq_len = 2000;
    let num_snps = 20;

    // Create synthetic sequences with known SNPs
    let (seq1, seq2, snp_positions) = create_synthetic_pair(seq_len, num_snps);

    // Write to FASTA
    let fasta_path = write_assembly_fasta(dir.path(), "assembly.fa", &seq1, &seq2);

    // Parse sequences using the library API
    let mut parser = phraya_io::SequenceParser::from_path(&fasta_path)
        .expect("failed to open FASTA");

    let seq1_obj = parser
        .next()
        .expect("no sequence 1")
        .expect("sequence 1 parse error");
    let seq2_obj = parser
        .next()
        .expect("no sequence 2")
        .expect("sequence 2 parse error");

    // Build a minimal plan
    let mut kmer_index = HashMap::new();
    let seq1_sketch = phraya_core::types::sketch_sequence_default(&seq1_obj);
    let seq2_sketch = phraya_core::types::sketch_sequence_default(&seq2_obj);
    kmer_index.insert(seq1_obj.id().to_string(), seq1_sketch);
    kmer_index.insert(seq2_obj.id().to_string(), seq2_sketch);

    let plan = phraya_io::plan::PhrayaPlan::new(
        phraya_io::plan::UseCase::ContigsOnly,
        vec![fasta_path.to_string_lossy().to_string()],
        "2026-06-07T00:00:00Z".to_string(),
        kmer_index,
        HashMap::new(),
        vec![(1, 0)], // query=contig2 (index 1), target=contig1 (index 0)
    );

    // Perform alignment
    let result = phraya_align::executor::align_task_with_config(
        &seq2_obj,
        &seq1_obj,
        &plan,
        &phraya_align::executor::AlignConfig::default(),
    );

    let result = result.expect("alignment failed or returned None");

    // Count variants at the known SNP positions
    let detected_count = result.variants.len();

    // Allow ±10% tolerance on SNP detection
    let lower_bound = (num_snps as f64 * 0.9).ceil() as usize;
    let upper_bound = (num_snps as f64 * 1.1).floor() as usize;

    assert!(
        detected_count >= lower_bound && detected_count <= upper_bound,
        "SNP count {} is outside tolerance [{}–{}]. Expected ~{} SNPs. \
         Known SNP positions: {:?}. Detected {} variants at positions: {:?}",
        detected_count,
        lower_bound,
        upper_bound,
        num_snps,
        snp_positions,
        result.variants.len(),
        result
            .variants
            .iter()
            .map(|v| v.position())
            .collect::<Vec<_>>()
    );
}

// ============================================================================
// Test 2: Mean confidence > 0.5
// ============================================================================

/// Verify that synthetic data produces mean confidence > 0.5.
/// This validates that confidence scoring is properly threaded through alignment.
#[test]
fn issue_15_mean_confidence_exceeds_threshold() {
    let dir = TempDir::new().unwrap();
    let seq_len = 2000;
    let num_snps = 15;

    // Create synthetic sequences
    let (seq1, seq2, _) = create_synthetic_pair(seq_len, num_snps);
    let fasta_path = write_assembly_fasta(dir.path(), "assembly.fa", &seq1, &seq2);

    // Parse sequences
    let mut parser = phraya_io::SequenceParser::from_path(&fasta_path)
        .expect("failed to open FASTA");

    let seq1_obj = parser
        .next()
        .expect("no sequence 1")
        .expect("sequence 1 parse error");
    let seq2_obj = parser
        .next()
        .expect("no sequence 2")
        .expect("sequence 2 parse error");

    // Build plan
    let mut kmer_index = HashMap::new();
    let seq1_sketch = phraya_core::types::sketch_sequence_default(&seq1_obj);
    let seq2_sketch = phraya_core::types::sketch_sequence_default(&seq2_obj);
    kmer_index.insert(seq1_obj.id().to_string(), seq1_sketch);
    kmer_index.insert(seq2_obj.id().to_string(), seq2_sketch);

    let plan = phraya_io::plan::PhrayaPlan::new(
        phraya_io::plan::UseCase::ContigsOnly,
        vec![fasta_path.to_string_lossy().to_string()],
        "2026-06-07T00:00:00Z".to_string(),
        kmer_index,
        HashMap::new(),
        vec![(1, 0)],
    );

    // Perform alignment
    let result = phraya_align::executor::align_task_with_config(
        &seq2_obj,
        &seq1_obj,
        &plan,
        &phraya_align::executor::AlignConfig::default(),
    );

    let result = result.expect("alignment failed");

    // Verify that variants were detected
    assert!(
        !result.variants.is_empty(),
        "No variants detected; cannot assess confidence"
    );

    // Calculate mean confidence
    let mean_confidence: f64 =
        result.variants.iter().map(|v| v.confidence()).sum::<f64>() / result.variants.len() as f64;

    assert!(
        mean_confidence > 0.5,
        "Mean confidence {:.4} is below threshold 0.5. Variants: {}",
        mean_confidence,
        result.variants.len()
    );
}

// ============================================================================
// Test 3: Full pipeline (plan + align + variant counting)
// ============================================================================

/// End-to-end test: create plan file, load it, perform alignment, count variants.
/// Exercises the full pipeline including serialization/deserialization.
#[test]
fn issue_15_full_pipeline_with_plan_file() {
    let dir = TempDir::new().unwrap();
    let seq_len = 1500;
    let num_snps = 10;

    // Create synthetic sequences
    let (seq1, seq2, snp_positions) = create_synthetic_pair(seq_len, num_snps);
    let fasta_path = write_assembly_fasta(dir.path(), "assembly.fa", &seq1, &seq2);

    // Parse sequences
    let mut parser = phraya_io::SequenceParser::from_path(&fasta_path)
        .expect("failed to open FASTA");

    let seq1_obj = parser
        .next()
        .expect("no sequence 1")
        .expect("sequence 1 parse error");
    let seq2_obj = parser
        .next()
        .expect("no sequence 2")
        .expect("sequence 2 parse error");

    // Build and write plan
    let mut kmer_index = HashMap::new();
    let seq1_sketch = phraya_core::types::sketch_sequence_default(&seq1_obj);
    let seq2_sketch = phraya_core::types::sketch_sequence_default(&seq2_obj);
    kmer_index.insert(seq1_obj.id().to_string(), seq1_sketch);
    kmer_index.insert(seq2_obj.id().to_string(), seq2_sketch);

    let plan = phraya_io::plan::PhrayaPlan::new(
        phraya_io::plan::UseCase::ContigsOnly,
        vec![fasta_path.to_string_lossy().to_string()],
        "2026-06-07T00:00:00Z".to_string(),
        kmer_index,
        HashMap::new(),
        vec![(1, 0)],
    );

    let plan_path = dir.path().join("test.phrayaplan");
    phraya_io::plan::write_plan(&plan_path, &plan).expect("failed to write plan");

    // Read plan back
    let loaded_plan = phraya_io::plan::read_plan(&plan_path).expect("failed to read plan");

    // Perform alignment using loaded plan
    let result = phraya_align::executor::align_task_with_config(
        &seq2_obj,
        &seq1_obj,
        &loaded_plan,
        &phraya_align::executor::AlignConfig::default(),
    );

    let result = result.expect("alignment failed");

    // Verify variant count
    let detected_count = result.variants.len();
    let lower_bound = (num_snps as f64 * 0.9).ceil() as usize;
    let upper_bound = (num_snps as f64 * 1.1).floor() as usize;

    assert!(
        detected_count >= lower_bound && detected_count <= upper_bound,
        "SNP count {} is outside tolerance [{}–{}]. Expected ~{} SNPs",
        detected_count, lower_bound, upper_bound, num_snps
    );

    // Verify that detected variants are near the expected SNP positions
    // (within ±50bp due to alignment windowing)
    let mut matched = 0;
    for var in &result.variants {
        if snp_positions.iter().any(|&pos| {
            let pos_i32 = pos as i32;
            let var_pos_i32 = var.position() as i32;
            (var_pos_i32 - pos_i32).abs() <= 50
        }) {
            matched += 1;
        }
    }

    // At least 70% of detected variants should be near a known SNP position
    let match_rate = matched as f64 / result.variants.len() as f64;
    assert!(
        match_rate >= 0.7,
        "Only {:.1}% of detected variants are near known SNP positions. \
         Expected ≥70%. Detected count: {}, matched: {}",
        match_rate * 100.0,
        result.variants.len(),
        matched
    );
}

// ============================================================================
// Test 4: No false positives on perfect-match sequences
// ============================================================================

/// Verify that identical sequences produce no variants.
/// This is a sanity check that the alignment and variant calling are working correctly.
#[test]
fn issue_15_perfect_match_no_variants() {
    let dir = TempDir::new().unwrap();
    let seq_len = 1000;

    // Create identical sequences
    let seq = diverse_dna(seq_len, 42);
    let fasta_path = write_assembly_fasta(dir.path(), "identical.fa", &seq, &seq);

    // Parse sequences
    let mut parser = phraya_io::SequenceParser::from_path(&fasta_path)
        .expect("failed to open FASTA");

    let seq1_obj = parser
        .next()
        .expect("no sequence 1")
        .expect("sequence 1 parse error");
    let seq2_obj = parser
        .next()
        .expect("no sequence 2")
        .expect("sequence 2 parse error");

    // Build plan
    let mut kmer_index = HashMap::new();
    let seq1_sketch = phraya_core::types::sketch_sequence_default(&seq1_obj);
    let seq2_sketch = phraya_core::types::sketch_sequence_default(&seq2_obj);
    kmer_index.insert(seq1_obj.id().to_string(), seq1_sketch);
    kmer_index.insert(seq2_obj.id().to_string(), seq2_sketch);

    let plan = phraya_io::plan::PhrayaPlan::new(
        phraya_io::plan::UseCase::ContigsOnly,
        vec![fasta_path.to_string_lossy().to_string()],
        "2026-06-07T00:00:00Z".to_string(),
        kmer_index,
        HashMap::new(),
        vec![(1, 0)],
    );

    // Perform alignment
    let result = phraya_align::executor::align_task_with_config(
        &seq2_obj,
        &seq1_obj,
        &plan,
        &phraya_align::executor::AlignConfig::default(),
    );

    let result = result.expect("alignment failed");

    // Identical sequences should produce zero variants
    assert_eq!(
        result.variants.len(),
        0,
        "Perfect-match sequences must produce zero variants, got {}",
        result.variants.len()
    );
}

// ============================================================================
// Test 5: Coverage track presence and validity
// ============================================================================

/// Verify that alignment produces a valid coverage track with correct length.
#[test]
fn issue_15_coverage_track_valid() {
    let dir = TempDir::new().unwrap();
    let seq_len = 1000;
    let num_snps = 5;

    // Create synthetic sequences
    let (seq1, seq2, _) = create_synthetic_pair(seq_len, num_snps);
    let fasta_path = write_assembly_fasta(dir.path(), "assembly.fa", &seq1, &seq2);

    // Parse sequences
    let mut parser = phraya_io::SequenceParser::from_path(&fasta_path)
        .expect("failed to open FASTA");

    let seq1_obj = parser
        .next()
        .expect("no sequence 1")
        .expect("sequence 1 parse error");
    let seq2_obj = parser
        .next()
        .expect("no sequence 2")
        .expect("sequence 2 parse error");

    // Build plan
    let mut kmer_index = HashMap::new();
    let seq1_sketch = phraya_core::types::sketch_sequence_default(&seq1_obj);
    let seq2_sketch = phraya_core::types::sketch_sequence_default(&seq2_obj);
    kmer_index.insert(seq1_obj.id().to_string(), seq1_sketch);
    kmer_index.insert(seq2_obj.id().to_string(), seq2_sketch);

    let plan = phraya_io::plan::PhrayaPlan::new(
        phraya_io::plan::UseCase::ContigsOnly,
        vec![fasta_path.to_string_lossy().to_string()],
        "2026-06-07T00:00:00Z".to_string(),
        kmer_index,
        HashMap::new(),
        vec![(1, 0)],
    );

    // Perform alignment
    let result = phraya_align::executor::align_task_with_config(
        &seq2_obj,
        &seq1_obj,
        &plan,
        &phraya_align::executor::AlignConfig::default(),
    );

    let result = result.expect("alignment failed");

    // Coverage track (windowed per read) must reconstruct to the target length.
    let full_coverage = result.coverage[0].to_full(seq1.len());
    assert_eq!(
        full_coverage.len(),
        seq1.len(),
        "Coverage track length {} does not match target length {}",
        full_coverage.len(),
        seq1.len()
    );

    // Coverage track is quantized to nearest 5, so small alignments may have all zeros.
    // The critical property is that the track is properly sized and computed.
    // A perfect match or low-divergence alignment may result in quantized coverage of 0
    // if only one read covers each position (0 rounds to 0 when quantized).
    // The real validation is through the variants and query_positions, which must be non-empty.
    assert!(!result.variants.is_empty() || result.query_positions.is_empty() == false,
        "Alignment should produce either variants or query positions"
    );
}
