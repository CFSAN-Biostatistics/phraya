//! Integration tests for `phraya align` command
//!
//! These tests verify the CLI accepts inputs, dispatches alignments, and writes .phraya files.
//! All tests should FAIL until implementation is complete (RED phase of TDD).

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use tempfile::TempDir;

/// Helper to get the phraya binary path
fn phraya_bin() -> PathBuf {
    let mut path = std::env::current_exe().unwrap();
    path.pop(); // remove test binary name
    path.pop(); // remove deps/
    path.push("phraya");
    path
}

/// Helper to create test FASTA file
fn create_test_fasta(path: &Path, name: &str, seq: &str) {
    fs::write(path, format!(">{}\n{}\n", name, seq)).unwrap();
}

/// Helper to create test FASTQ file
fn create_test_fastq(path: &Path, name: &str, seq: &str, qual: &str) {
    fs::write(path, format!("@{}\n{}\n+\n{}\n", name, seq, qual)).unwrap();
}

#[test]
fn test_align_accepts_required_arguments() {
    // Acceptance: CLI accepts: <input_files>..., --reference <path> (optional),
    //             --strategy <balanced>, --outdir <path>

    let tmpdir = TempDir::new().unwrap();
    let ref_path = tmpdir.path().join("reference.fasta");
    let sample1 = tmpdir.path().join("sample1.fasta");
    let outdir = tmpdir.path().join("results");

    create_test_fasta(&ref_path, "ref", "ATCGATCG");
    create_test_fasta(&sample1, "s1", "ATCGATCG");
    fs::create_dir(&outdir).unwrap();

    let output = Command::new(phraya_bin())
        .arg("align")
        .arg(&sample1)
        .arg("--reference")
        .arg(&ref_path)
        .arg("--strategy")
        .arg("balanced")
        .arg("--outdir")
        .arg(&outdir)
        .output()
        .expect("Failed to execute phraya align");

    assert!(
        output.status.success(),
        "phraya align should accept required arguments. stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn test_align_rejects_invalid_strategy() {
    // Acceptance: Phase 1 only supports --strategy balanced

    let tmpdir = TempDir::new().unwrap();
    let sample1 = tmpdir.path().join("sample1.fasta");
    let outdir = tmpdir.path().join("results");

    create_test_fasta(&sample1, "s1", "ATCGATCG");
    fs::create_dir(&outdir).unwrap();

    let output = Command::new(phraya_bin())
        .arg("align")
        .arg(&sample1)
        .arg("--strategy")
        .arg("exact") // Not implemented in Phase 1
        .arg("--outdir")
        .arg(&outdir)
        .output()
        .expect("Failed to execute phraya align");

    assert!(
        !output.status.success(),
        "phraya align should reject non-balanced strategy in Phase 1"
    );

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("balanced") || stderr.contains("strategy"),
        "Error message should mention strategy requirement. stderr: {}",
        stderr
    );
}

#[test]
fn test_reference_mode_creates_n_phraya_files() {
    // Acceptance: Reference mode: N alignments (each input → reference) → N .phraya files

    let tmpdir = TempDir::new().unwrap();
    let ref_path = tmpdir.path().join("reference.fasta");
    let sample1 = tmpdir.path().join("sample1.fasta");
    let sample2 = tmpdir.path().join("sample2.fasta");
    let sample3 = tmpdir.path().join("sample3.fasta");
    let outdir = tmpdir.path().join("results");

    create_test_fasta(&ref_path, "ref", "ATCGATCGATCGATCG");
    create_test_fasta(&sample1, "s1", "ATCGATCGATCGATCG");
    create_test_fasta(&sample2, "s2", "ATCGATCGTTCGATCG"); // 1 SNP
    create_test_fasta(&sample3, "s3", "ATCGATCGATCGATCG");
    fs::create_dir(&outdir).unwrap();

    let output = Command::new(phraya_bin())
        .arg("align")
        .arg(&sample1)
        .arg(&sample2)
        .arg(&sample3)
        .arg("--reference")
        .arg(&ref_path)
        .arg("--strategy")
        .arg("balanced")
        .arg("--outdir")
        .arg(&outdir)
        .output()
        .expect("Failed to execute phraya align");

    assert!(
        output.status.success(),
        "phraya align reference mode should succeed. stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Verify 3 .phraya files exist
    let phraya_files: Vec<_> = fs::read_dir(&outdir)
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| e.path().extension().is_some_and(|ext| ext == "phraya"))
        .collect();

    assert_eq!(
        phraya_files.len(),
        3,
        "Reference mode with 3 inputs should create 3 .phraya files. Found: {:?}",
        phraya_files.iter().map(|e| e.path()).collect::<Vec<_>>()
    );

    // Verify files are named after input samples
    let expected_names = ["sample1.phraya", "sample2.phraya", "sample3.phraya"];
    for name in &expected_names {
        let path = outdir.join(name);
        assert!(
            path.exists(),
            "Expected .phraya file {} to exist",
            path.display()
        );

        // Verify file is not empty
        let metadata = fs::metadata(&path).unwrap();
        assert!(
            metadata.len() > 0,
            ".phraya file {} should not be empty",
            path.display()
        );
    }
}

#[test]
fn test_msa_mode_creates_n_phraya_files_all_vs_all() {
    // Acceptance: MSA mode: N×(N-1)/2 alignments → N .phraya files
    //             (each aggregates all pairs involving that sample)

    let tmpdir = TempDir::new().unwrap();
    let sample1 = tmpdir.path().join("sample1.fasta");
    let sample2 = tmpdir.path().join("sample2.fasta");
    let sample3 = tmpdir.path().join("sample3.fasta");
    let outdir = tmpdir.path().join("results");

    create_test_fasta(&sample1, "s1", "ATCGATCGATCGATCG");
    create_test_fasta(&sample2, "s2", "ATCGATCGTTCGATCG"); // 1 SNP vs s1
    create_test_fasta(&sample3, "s3", "ATCGAGCGATCGATCG"); // 1 SNP vs s1, different from s2
    fs::create_dir(&outdir).unwrap();

    let output = Command::new(phraya_bin())
        .arg("align")
        .arg(&sample1)
        .arg(&sample2)
        .arg(&sample3)
        .arg("--strategy")
        .arg("balanced")
        .arg("--outdir")
        .arg(&outdir)
        .output()
        .expect("Failed to execute phraya align");

    assert!(
        output.status.success(),
        "phraya align MSA mode should succeed. stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Verify 3 .phraya files exist (one per sample)
    let phraya_files: Vec<_> = fs::read_dir(&outdir)
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| e.path().extension().is_some_and(|ext| ext == "phraya"))
        .collect();

    assert_eq!(
        phraya_files.len(),
        3,
        "MSA mode with 3 inputs should create 3 .phraya files (N×(N-1)/2 = 3 comparisons). Found: {:?}",
        phraya_files.iter().map(|e| e.path()).collect::<Vec<_>>()
    );

    // Verify each file contains observations from all pairwise comparisons involving that sample
    // sample1 should have observations from: s1-vs-s2, s1-vs-s3
    // sample2 should have observations from: s1-vs-s2, s2-vs-s3
    // sample3 should have observations from: s1-vs-s3, s2-vs-s3

    for name in &["sample1.phraya", "sample2.phraya", "sample3.phraya"] {
        let path = outdir.join(name);
        assert!(
            path.exists(),
            "Expected .phraya file {} to exist in MSA mode",
            path.display()
        );

        let metadata = fs::metadata(&path).unwrap();
        assert!(
            metadata.len() > 0,
            ".phraya file {} should not be empty in MSA mode",
            path.display()
        );
    }
}

#[test]
fn test_auto_detects_fasta_format() {
    // Acceptance: Auto-detects input formats via phraya-io (FASTA/FASTQ/BAM/CRAM, gz/bz2/raw)

    let tmpdir = TempDir::new().unwrap();
    let ref_path = tmpdir.path().join("reference.fasta");
    let sample_fasta = tmpdir.path().join("sample.fasta");
    let outdir = tmpdir.path().join("results");

    create_test_fasta(&ref_path, "ref", "ATCGATCG");
    create_test_fasta(&sample_fasta, "s1", "ATCGATCG");
    fs::create_dir(&outdir).unwrap();

    let output = Command::new(phraya_bin())
        .arg("align")
        .arg(&sample_fasta)
        .arg("--reference")
        .arg(&ref_path)
        .arg("--strategy")
        .arg("balanced")
        .arg("--outdir")
        .arg(&outdir)
        .output()
        .expect("Failed to execute phraya align");

    assert!(
        output.status.success(),
        "phraya align should auto-detect FASTA format. stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    assert!(
        outdir.join("sample.phraya").exists(),
        "Should create .phraya file for FASTA input"
    );
}

#[test]
fn test_auto_detects_fastq_format() {
    // Acceptance: Auto-detects FASTQ format

    let tmpdir = TempDir::new().unwrap();
    let ref_path = tmpdir.path().join("reference.fasta");
    let sample_fastq = tmpdir.path().join("sample.fastq");
    let outdir = tmpdir.path().join("results");

    create_test_fasta(&ref_path, "ref", "ATCGATCG");
    create_test_fastq(&sample_fastq, "read1", "ATCGATCG", "IIIIIIII");
    fs::create_dir(&outdir).unwrap();

    let output = Command::new(phraya_bin())
        .arg("align")
        .arg(&sample_fastq)
        .arg("--reference")
        .arg(&ref_path)
        .arg("--strategy")
        .arg("balanced")
        .arg("--outdir")
        .arg(&outdir)
        .output()
        .expect("Failed to execute phraya align");

    assert!(
        output.status.success(),
        "phraya align should auto-detect FASTQ format. stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    assert!(
        outdir.join("sample.phraya").exists(),
        "Should create .phraya file for FASTQ input"
    );
}

#[test]
fn test_auto_detects_gzipped_input() {
    // Acceptance: Auto-detects gz compressed input

    use flate2::Compression;
    use flate2::write::GzEncoder;
    use std::io::Write;

    let tmpdir = TempDir::new().unwrap();
    let ref_path = tmpdir.path().join("reference.fasta");
    let sample_gz = tmpdir.path().join("sample.fasta.gz");
    let outdir = tmpdir.path().join("results");

    create_test_fasta(&ref_path, "ref", "ATCGATCG");

    // Create gzipped FASTA
    let fasta_content = b">s1\nATCGATCG\n";
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(fasta_content).unwrap();
    let compressed = encoder.finish().unwrap();
    fs::write(&sample_gz, compressed).unwrap();

    fs::create_dir(&outdir).unwrap();

    let output = Command::new(phraya_bin())
        .arg("align")
        .arg(&sample_gz)
        .arg("--reference")
        .arg(&ref_path)
        .arg("--strategy")
        .arg("balanced")
        .arg("--outdir")
        .arg(&outdir)
        .output()
        .expect("Failed to execute phraya align");

    assert!(
        output.status.success(),
        "phraya align should auto-detect gzipped FASTA. stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    assert!(
        outdir.join("sample.phraya").exists(),
        "Should create .phraya file for gzipped FASTA input"
    );
}

#[test]
fn test_mixed_format_inputs() {
    // Acceptance: Can handle mixed FASTA/FASTQ inputs in same invocation

    let tmpdir = TempDir::new().unwrap();
    let ref_path = tmpdir.path().join("reference.fasta");
    let sample_fasta = tmpdir.path().join("sample1.fasta");
    let sample_fastq = tmpdir.path().join("sample2.fastq");
    let outdir = tmpdir.path().join("results");

    create_test_fasta(&ref_path, "ref", "ATCGATCG");
    create_test_fasta(&sample_fasta, "s1", "ATCGATCG");
    create_test_fastq(&sample_fastq, "read1", "ATCGATCG", "IIIIIIII");
    fs::create_dir(&outdir).unwrap();

    let output = Command::new(phraya_bin())
        .arg("align")
        .arg(&sample_fasta)
        .arg(&sample_fastq)
        .arg("--reference")
        .arg(&ref_path)
        .arg("--strategy")
        .arg("balanced")
        .arg("--outdir")
        .arg(&outdir)
        .output()
        .expect("Failed to execute phraya align");

    assert!(
        output.status.success(),
        "phraya align should handle mixed format inputs. stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    assert!(
        outdir.join("sample1.phraya").exists(),
        "Should create .phraya file for FASTA input"
    );
    assert!(
        outdir.join("sample2.phraya").exists(),
        "Should create .phraya file for FASTQ input"
    );
}

#[test]
fn test_progress_reporting() {
    // Acceptance: Progress reporting (log which pair is being processed)

    let tmpdir = TempDir::new().unwrap();
    let ref_path = tmpdir.path().join("reference.fasta");
    let sample1 = tmpdir.path().join("sample1.fasta");
    let sample2 = tmpdir.path().join("sample2.fasta");
    let outdir = tmpdir.path().join("results");

    create_test_fasta(&ref_path, "ref", "ATCGATCGATCGATCG");
    create_test_fasta(&sample1, "s1", "ATCGATCGATCGATCG");
    create_test_fasta(&sample2, "s2", "ATCGATCGTTCGATCG");
    fs::create_dir(&outdir).unwrap();

    let output = Command::new(phraya_bin())
        .arg("align")
        .arg(&sample1)
        .arg(&sample2)
        .arg("--reference")
        .arg(&ref_path)
        .arg("--strategy")
        .arg("balanced")
        .arg("--outdir")
        .arg(&outdir)
        .output()
        .expect("Failed to execute phraya align");

    assert!(
        output.status.success(),
        "phraya align should succeed. stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let stderr = String::from_utf8_lossy(&output.stderr);

    // Check for progress indicators mentioning alignment pairs
    assert!(
        stderr.contains("sample1")
            || stderr.contains("sample2")
            || stderr.contains("Aligning")
            || stderr.contains("Processing"),
        "Progress output should mention which samples are being processed. stderr: {}",
        stderr
    );
}

#[test]
fn test_continues_on_individual_alignment_failure() {
    // Acceptance: Error handling: continues if one pair fails

    let tmpdir = TempDir::new().unwrap();
    let ref_path = tmpdir.path().join("reference.fasta");
    let sample1 = tmpdir.path().join("sample1.fasta");
    let sample2 = tmpdir.path().join("sample2.fasta"); // Will be invalid/missing
    let sample3 = tmpdir.path().join("sample3.fasta");
    let outdir = tmpdir.path().join("results");

    create_test_fasta(&ref_path, "ref", "ATCGATCG");
    create_test_fasta(&sample1, "s1", "ATCGATCG");
    // sample2 deliberately not created to cause failure
    create_test_fasta(&sample3, "s3", "ATCGATCG");
    fs::create_dir(&outdir).unwrap();

    let output = Command::new(phraya_bin())
        .arg("align")
        .arg(&sample1)
        .arg(&sample2) // Missing file
        .arg(&sample3)
        .arg("--reference")
        .arg(&ref_path)
        .arg("--strategy")
        .arg("balanced")
        .arg("--outdir")
        .arg(&outdir)
        .output()
        .expect("Failed to execute phraya align");

    // Should continue and process sample1 and sample3 even though sample2 fails
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("sample2") && (stderr.contains("error") || stderr.contains("failed")),
        "Should report error for missing sample2. stderr: {}",
        stderr
    );

    // Verify that sample1 and sample3 were still processed
    assert!(
        outdir.join("sample1.phraya").exists() || outdir.join("sample3.phraya").exists(),
        "Should process other samples even when one fails"
    );
}

#[test]
fn test_parallel_dispatch_via_rayon() {
    // Acceptance: Parallel dispatch via Rayon for pairwise alignments
    // This test verifies that multiple alignments can happen concurrently

    use std::time::Instant;

    let tmpdir = TempDir::new().unwrap();
    let ref_path = tmpdir.path().join("reference.fasta");
    let outdir = tmpdir.path().join("results");

    // Create reference with longer sequence to ensure some processing time
    let long_seq = "ATCGATCG".repeat(100); // 800bp
    create_test_fasta(&ref_path, "ref", &long_seq);

    // Create multiple samples
    let mut samples = Vec::new();
    for i in 1..=4 {
        let sample = tmpdir.path().join(format!("sample{}.fasta", i));
        create_test_fasta(&sample, &format!("s{}", i), &long_seq);
        samples.push(sample);
    }

    fs::create_dir(&outdir).unwrap();

    let start = Instant::now();

    let mut cmd = Command::new(phraya_bin());
    cmd.arg("align");
    for sample in &samples {
        cmd.arg(sample);
    }
    cmd.arg("--reference")
        .arg(&ref_path)
        .arg("--strategy")
        .arg("balanced")
        .arg("--outdir")
        .arg(&outdir);

    let output = cmd.output().expect("Failed to execute phraya align");

    let duration = start.elapsed();

    assert!(
        output.status.success(),
        "phraya align should succeed. stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Verify all output files exist
    for i in 1..=4 {
        let phraya_file = outdir.join(format!("sample{}.phraya", i));
        assert!(
            phraya_file.exists(),
            "Expected .phraya file for sample{} to exist",
            i
        );
    }

    // Note: We can't directly test Rayon usage, but the implementation should use it.
    // This test at least verifies that multiple inputs are processed successfully.
    println!("Processed 4 alignments in {:?}", duration);
}

#[test]
fn test_requires_at_least_one_input() {
    // Edge case: No input files provided

    let tmpdir = TempDir::new().unwrap();
    let ref_path = tmpdir.path().join("reference.fasta");
    let outdir = tmpdir.path().join("results");

    create_test_fasta(&ref_path, "ref", "ATCGATCG");
    fs::create_dir(&outdir).unwrap();

    let output = Command::new(phraya_bin())
        .arg("align")
        .arg("--reference")
        .arg(&ref_path)
        .arg("--strategy")
        .arg("balanced")
        .arg("--outdir")
        .arg(&outdir)
        .output()
        .expect("Failed to execute phraya align");

    assert!(
        !output.status.success(),
        "phraya align should fail when no input files provided"
    );

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("input") || stderr.contains("required"),
        "Error message should indicate missing input files. stderr: {}",
        stderr
    );
}

#[test]
fn test_requires_outdir() {
    // Edge case: --outdir is required

    let tmpdir = TempDir::new().unwrap();
    let ref_path = tmpdir.path().join("reference.fasta");
    let sample1 = tmpdir.path().join("sample1.fasta");

    create_test_fasta(&ref_path, "ref", "ATCGATCG");
    create_test_fasta(&sample1, "s1", "ATCGATCG");

    let output = Command::new(phraya_bin())
        .arg("align")
        .arg(&sample1)
        .arg("--reference")
        .arg(&ref_path)
        .arg("--strategy")
        .arg("balanced")
        .output()
        .expect("Failed to execute phraya align");

    assert!(
        !output.status.success(),
        "phraya align should fail when --outdir not provided"
    );

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("outdir") || stderr.contains("required"),
        "Error message should indicate missing --outdir. stderr: {}",
        stderr
    );
}

#[test]
fn test_msa_mode_requires_multiple_samples() {
    // Edge case: MSA mode with only 1 sample should error or handle gracefully

    let tmpdir = TempDir::new().unwrap();
    let sample1 = tmpdir.path().join("sample1.fasta");
    let outdir = tmpdir.path().join("results");

    create_test_fasta(&sample1, "s1", "ATCGATCG");
    fs::create_dir(&outdir).unwrap();

    let output = Command::new(phraya_bin())
        .arg("align")
        .arg(&sample1)
        .arg("--strategy")
        .arg("balanced")
        .arg("--outdir")
        .arg(&outdir)
        .output()
        .expect("Failed to execute phraya align");

    assert!(
        !output.status.success(),
        "MSA mode with single sample should fail or warn"
    );

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("multiple") || stderr.contains("at least"),
        "Should indicate MSA mode requires multiple samples. stderr: {}",
        stderr
    );
}

#[test]
fn test_creates_outdir_if_missing() {
    // Acceptance: Should create output directory if it doesn't exist

    let tmpdir = TempDir::new().unwrap();
    let ref_path = tmpdir.path().join("reference.fasta");
    let sample1 = tmpdir.path().join("sample1.fasta");
    let outdir = tmpdir.path().join("results/nested/path");

    create_test_fasta(&ref_path, "ref", "ATCGATCG");
    create_test_fasta(&sample1, "s1", "ATCGATCG");
    // Don't create outdir

    let output = Command::new(phraya_bin())
        .arg("align")
        .arg(&sample1)
        .arg("--reference")
        .arg(&ref_path)
        .arg("--strategy")
        .arg("balanced")
        .arg("--outdir")
        .arg(&outdir)
        .output()
        .expect("Failed to execute phraya align");

    assert!(
        output.status.success(),
        "phraya align should create outdir if missing. stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    assert!(outdir.exists(), "Output directory should be created");

    assert!(
        outdir.join("sample1.phraya").exists(),
        "Should create .phraya file in new directory"
    );
}

#[test]
fn test_reference_mode_explicit_with_two_samples() {
    // Integration test from acceptance criteria: 3 samples in reference mode → verify 3 .phraya files
    // This is a more thorough version testing the full workflow

    let tmpdir = TempDir::new().unwrap();
    let ref_path = tmpdir.path().join("reference.fasta");
    let sample1 = tmpdir.path().join("sample1.fasta");
    let sample2 = tmpdir.path().join("sample2.fasta");
    let sample3 = tmpdir.path().join("sample3.fasta");
    let outdir = tmpdir.path().join("results");

    // Create reference
    let ref_seq = "ATCGATCGATCGATCGATCGATCGATCGATCG"; // 32bp
    create_test_fasta(&ref_path, "reference", ref_seq);

    // Create samples with known variants
    create_test_fasta(&sample1, "sample1", "ATCGATCGATCGATCGATCGATCGATCGATCG"); // identical
    create_test_fasta(&sample2, "sample2", "ATCGATCGATCGATGGATCGATCGATCGATCG"); // C→G at pos 15
    create_test_fasta(&sample3, "sample3", "ATCGATCGATCGATCGATCGATCTATCGATCG"); // G→T at pos 22

    fs::create_dir(&outdir).unwrap();

    let output = Command::new(phraya_bin())
        .arg("align")
        .arg(&sample1)
        .arg(&sample2)
        .arg(&sample3)
        .arg("--reference")
        .arg(&ref_path)
        .arg("--strategy")
        .arg("balanced")
        .arg("--outdir")
        .arg(&outdir)
        .output()
        .expect("Failed to execute phraya align");

    assert!(
        output.status.success(),
        "phraya align reference mode integration test should succeed. stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Verify exactly 3 .phraya files
    let phraya_files: Vec<_> = fs::read_dir(&outdir)
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| e.path().extension().is_some_and(|ext| ext == "phraya"))
        .collect();

    assert_eq!(
        phraya_files.len(),
        3,
        "Should create exactly 3 .phraya files. Found: {:?}",
        phraya_files.iter().map(|e| e.path()).collect::<Vec<_>>()
    );

    // Verify each file has content
    for entry in phraya_files {
        let metadata = fs::metadata(entry.path()).unwrap();
        assert!(
            metadata.len() > 0,
            ".phraya file {:?} should have content",
            entry.path()
        );
    }
}

#[test]
fn test_msa_mode_integration_three_samples() {
    // Integration test from acceptance criteria: 3 samples in MSA mode → verify 3 .phraya files
    // with all-vs-all observations

    let tmpdir = TempDir::new().unwrap();
    let sample1 = tmpdir.path().join("sample1.fasta");
    let sample2 = tmpdir.path().join("sample2.fasta");
    let sample3 = tmpdir.path().join("sample3.fasta");
    let outdir = tmpdir.path().join("results");

    // Create samples with pairwise differences to detect in MSA mode
    let base_seq = "ATCGATCGATCGATCGATCGATCGATCGATCG";
    create_test_fasta(&sample1, "sample1", base_seq); // reference
    create_test_fasta(&sample2, "sample2", "ATCGATCGATCGATGGATCGATCGATCGATCG"); // differs from s1
    create_test_fasta(&sample3, "sample3", "ATCGATCGATCGATCGATCGATCTATCGATCG"); // differs from s1 and s2

    fs::create_dir(&outdir).unwrap();

    let output = Command::new(phraya_bin())
        .arg("align")
        .arg(&sample1)
        .arg(&sample2)
        .arg(&sample3)
        .arg("--strategy")
        .arg("balanced")
        .arg("--outdir")
        .arg(&outdir)
        .output()
        .expect("Failed to execute phraya align");

    assert!(
        output.status.success(),
        "phraya align MSA mode integration test should succeed. stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Verify exactly 3 .phraya files (one per sample, aggregating all comparisons)
    let phraya_files: Vec<_> = fs::read_dir(&outdir)
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| e.path().extension().is_some_and(|ext| ext == "phraya"))
        .collect();

    assert_eq!(
        phraya_files.len(),
        3,
        "MSA mode should create exactly 3 .phraya files. Found: {:?}",
        phraya_files.iter().map(|e| e.path()).collect::<Vec<_>>()
    );

    // Verify each file has content (should contain observations from pairwise alignments)
    for entry in phraya_files {
        let metadata = fs::metadata(entry.path()).unwrap();
        assert!(
            metadata.len() > 0,
            ".phraya file {:?} should have content from all-vs-all comparisons",
            entry.path()
        );
    }

    // Log progress output for verification
    let stderr = String::from_utf8_lossy(&output.stderr);
    println!("MSA mode progress output:\n{}", stderr);
}
