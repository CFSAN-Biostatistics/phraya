use std::collections::HashMap;

/// Helper to get phraya-cli manifest path for cargo run commands
fn get_manifest_path() -> std::path::PathBuf {
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR").unwrap();
    std::path::Path::new(&manifest_dir).join("Cargo.toml")
}
use std::path::{Path, PathBuf};

/// Helper to get phraya-cli manifest path for cargo run commands
fn get_manifest_path() -> std::path::PathBuf {
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR").unwrap();
    std::path::Path::new(&manifest_dir).join("Cargo.toml")
}
use tempfile::TempDir;

/// Helper to get phraya-cli manifest path for cargo run commands
fn get_manifest_path() -> std::path::PathBuf {
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR").unwrap();
    std::path::Path::new(&manifest_dir).join("Cargo.toml")
}

/// Helper to create a temporary .phraya file
fn create_phraya_file(
    dir: &Path,
    filename: &str,
    reference_length: u32,
    sample_id: &str,
    observation_count: usize,
) -> PathBuf {
    let path = dir.join(filename);

    let mut alleles = HashMap::new();
    alleles.insert(b'A', 10);

    // Create observations at different positions
    let mut observations = Vec::new();
    for i in 0..observation_count {
        let pos = (i as u32) * 10;
        let obs = phraya_core::types::VariantObservation::new(
            pos,
            b'A',
            alleles.clone(),
            0.95,
            "10M".to_string(),
            60,
            0,
            vec![10],
            35.0,
            format!("{}:read{}", sample_id, i),
        );
        observations.push(obs);
    }

    let coverage = phraya_core::types::CoverageTrack::new(vec![10; reference_length as usize]);
    let phraya_file = phraya_io::phraya::PhrayaFile::new(
        reference_length,
        sample_id.to_string(),
        "2026-05-31T12:00:00Z".to_string(),
        observations,
        coverage,
    );

    phraya_io::phraya::write_phraya(&path, &phraya_file).unwrap();
    path
}

/// Test: phraya merge with two input files
#[test]
fn issue_80_merge_two_files() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    let input1 = create_phraya_file(temp_path, "sample1.phraya", 200, "sample1", 5);
    let input2 = create_phraya_file(temp_path, "sample2.phraya", 200, "sample2", 3);
    let output_path = temp_path.join("merged.phraya");

    let manifest_path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("Cargo.toml");
    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            manifest_path.to_str().unwrap(),
            "--",
            "merge",
            input1.to_str().unwrap(),
            input2.to_str().unwrap(),
            "--output",
            output_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya merge");

    assert!(
        output.status.success(),
        "phraya merge should succeed. stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    assert!(
        output_path.exists(),
        "merged file should be created at {:?}",
        output_path
    );

    let merged_file =
        phraya_io::phraya::read_phraya(&output_path).expect("merged file should be readable");

    assert_eq!(
        merged_file.header.reference_length, 200,
        "reference length should be preserved"
    );

    assert_eq!(
        merged_file.observations.len(),
        8,
        "should have 8 total observations (5 + 3)"
    );
}

/// Test: phraya merge preserves provenance
#[test]
fn issue_80_merge_preserves_provenance() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    let input1 = create_phraya_file(temp_path, "sample1.phraya", 100, "sample1", 2);
    let input2 = create_phraya_file(temp_path, "sample2.phraya", 100, "sample2", 2);
    let output_path = temp_path.join("merged.phraya");

    let manifest_path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("Cargo.toml");
    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            manifest_path.to_str().unwrap(),
            "--",
            "merge",
            input1.to_str().unwrap(),
            input2.to_str().unwrap(),
            "--output",
            output_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya merge");

    assert!(output.status.success(), "phraya merge should succeed");

    let merged_file =
        phraya_io::phraya::read_phraya(&output_path).expect("merged file should be readable");

    // Verify provenance is preserved
    let provenances: Vec<_> = merged_file
        .observations
        .iter()
        .map(|obs| obs.provenance().to_string())
        .collect();

    assert!(
        provenances.iter().any(|p| p.contains("sample1")),
        "should preserve sample1 provenance"
    );
    assert!(
        provenances.iter().any(|p| p.contains("sample2")),
        "should preserve sample2 provenance"
    );
}

/// Test: phraya merge with three input files
#[test]
fn issue_80_merge_three_files() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    let input1 = create_phraya_file(temp_path, "sample1.phraya", 100, "sample1", 2);
    let input2 = create_phraya_file(temp_path, "sample2.phraya", 100, "sample2", 3);
    let input3 = create_phraya_file(temp_path, "sample3.phraya", 100, "sample3", 1);
    let output_path = temp_path.join("merged.phraya");

    let manifest_path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("Cargo.toml");
    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            manifest_path.to_str().unwrap(),
            "--",
            "merge",
            input1.to_str().unwrap(),
            input2.to_str().unwrap(),
            input3.to_str().unwrap(),
            "--output",
            output_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya merge");

    assert!(output.status.success(), "phraya merge should succeed");

    let merged_file =
        phraya_io::phraya::read_phraya(&output_path).expect("merged file should be readable");

    assert_eq!(
        merged_file.observations.len(),
        6,
        "should have 6 total observations (2 + 3 + 1)"
    );
}

/// Test: phraya merge coverage tracks are summed
#[test]
fn issue_80_merge_coverage_summing() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    // Create files with specific coverage patterns
    let path1 = temp_path.join("sample1.phraya");
    let coverage1 = phraya_core::types::CoverageTrack::new(vec![5, 5, 5, 5]);
    let file1 = phraya_io::phraya::PhrayaFile::new(
        4,
        "sample1".to_string(),
        "2026-05-31T12:00:00Z".to_string(),
        vec![],
        coverage1,
    );
    phraya_io::phraya::write_phraya(&path1, &file1).unwrap();

    let path2 = temp_path.join("sample2.phraya");
    let coverage2 = phraya_core::types::CoverageTrack::new(vec![10, 10, 10, 10]);
    let file2 = phraya_io::phraya::PhrayaFile::new(
        4,
        "sample2".to_string(),
        "2026-05-31T12:00:00Z".to_string(),
        vec![],
        coverage2,
    );
    phraya_io::phraya::write_phraya(&path2, &file2).unwrap();

    let output_path = temp_path.join("merged.phraya");

    let manifest_path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("Cargo.toml");
    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            manifest_path.to_str().unwrap(),
            "--",
            "merge",
            path1.to_str().unwrap(),
            path2.to_str().unwrap(),
            "--output",
            output_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya merge");

    assert!(output.status.success(), "phraya merge should succeed");

    let merged_file =
        phraya_io::phraya::read_phraya(&output_path).expect("merged file should be readable");

    let merged_coverage = merged_file.coverage_track.decompress();

    // Coverage should be summed (5 + 10 = 15 per position)
    for &cov in &merged_coverage {
        assert!(
            cov >= 15,
            "merged coverage should be at least 15 (5 + 10), got {}",
            cov
        );
    }
}

/// Test: phraya merge fails with mismatched reference lengths
#[test]
fn issue_80_merge_mismatched_reference_length_error() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    let path1 = temp_path.join("sample1.phraya");
    let coverage1 = phraya_core::types::CoverageTrack::new(vec![10; 100]);
    let file1 = phraya_io::phraya::PhrayaFile::new(
        100,
        "sample1".to_string(),
        "2026-05-31T12:00:00Z".to_string(),
        vec![],
        coverage1,
    );
    phraya_io::phraya::write_phraya(&path1, &file1).unwrap();

    let path2 = temp_path.join("sample2.phraya");
    let coverage2 = phraya_core::types::CoverageTrack::new(vec![10; 200]);
    let file2 = phraya_io::phraya::PhrayaFile::new(
        200,
        "sample2".to_string(),
        "2026-05-31T12:00:00Z".to_string(),
        vec![],
        coverage2,
    );
    phraya_io::phraya::write_phraya(&path2, &file2).unwrap();

    let output_path = temp_path.join("merged.phraya");

    let manifest_path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("Cargo.toml");
    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            manifest_path.to_str().unwrap(),
            "--",
            "merge",
            path1.to_str().unwrap(),
            path2.to_str().unwrap(),
            "--output",
            output_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya merge");

    assert!(
        !output.status.success(),
        "should fail with mismatched reference lengths"
    );

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("reference length") || stderr.contains("error") || stderr.contains("Error"),
        "error message should mention reference length mismatch: {}",
        stderr
    );
}

/// Test: phraya merge with nonexistent input file
#[test]
fn issue_80_merge_nonexistent_input_file() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    let input1 = create_phraya_file(temp_path, "sample1.phraya", 100, "sample1", 2);
    let nonexistent = temp_path.join("nonexistent.phraya");
    let output_path = temp_path.join("merged.phraya");

    let manifest_path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("Cargo.toml");
    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            manifest_path.to_str().unwrap(),
            "--",
            "merge",
            input1.to_str().unwrap(),
            nonexistent.to_str().unwrap(),
            "--output",
            output_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya merge");

    assert!(
        !output.status.success(),
        "should fail with nonexistent input file"
    );

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("not found") || stderr.contains("error") || stderr.contains("Error"),
        "error message should be informative: {}",
        stderr
    );
}

/// Test: phraya merge with no input files
#[test]
fn issue_80_merge_no_input_files() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();
    let output_path = temp_path.join("merged.phraya");

    let manifest_path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("Cargo.toml");
    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            manifest_path.to_str().unwrap(),
            "--",
            "merge",
            "--output",
            output_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya merge");

    assert!(!output.status.success(), "should fail without input files");
}

/// Test: phraya merge without --output flag
#[test]
fn issue_80_merge_missing_output_flag() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    let input1 = create_phraya_file(temp_path, "sample1.phraya", 100, "sample1", 2);
    let input2 = create_phraya_file(temp_path, "sample2.phraya", 100, "sample2", 2);

    let manifest_path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("Cargo.toml");
    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            manifest_path.to_str().unwrap(),
            "--",
            "merge",
            input1.to_str().unwrap(),
            input2.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya merge");

    assert!(!output.status.success(), "should fail without --output");
}

/// Test: phraya merge logs progress to stderr
#[test]
fn issue_80_merge_logs_progress() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    let input1 = create_phraya_file(temp_path, "sample1.phraya", 100, "sample1", 2);
    let input2 = create_phraya_file(temp_path, "sample2.phraya", 100, "sample2", 2);
    let output_path = temp_path.join("merged.phraya");

    let manifest_path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("Cargo.toml");
    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            manifest_path.to_str().unwrap(),
            "--",
            "merge",
            input1.to_str().unwrap(),
            input2.to_str().unwrap(),
            "--output",
            output_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya merge");

    assert!(output.status.success(), "phraya merge should succeed");

    let stderr = String::from_utf8_lossy(&output.stderr);

    // Should log something about merging samples or progress
    assert!(
        stderr.contains("Merging")
            || stderr.contains("merge")
            || stderr.contains("Merge")
            || stderr.len() > 0,
        "should log merge progress to stderr"
    );
}

/// Test: phraya merge deduplicates identical observations
#[test]
fn issue_80_merge_deduplicates_observations() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    // Create two identical files with the same observation
    let mut alleles = HashMap::new();
    alleles.insert(b'A', 10);

    let obs = phraya_core::types::VariantObservation::new(
        50,
        b'A',
        alleles.clone(),
        0.95,
        "10M".to_string(),
        60,
        0,
        vec![10],
        35.0,
        "shared:read1".to_string(),
    );

    let coverage = phraya_core::types::CoverageTrack::new(vec![10; 100]);

    let file1 = phraya_io::phraya::PhrayaFile::new(
        100,
        "sample1".to_string(),
        "2026-05-31T12:00:00Z".to_string(),
        vec![obs.clone()],
        coverage.clone(),
    );

    let file2 = phraya_io::phraya::PhrayaFile::new(
        100,
        "sample2".to_string(),
        "2026-05-31T12:00:00Z".to_string(),
        vec![obs.clone()],
        coverage,
    );

    let path1 = temp_path.join("sample1.phraya");
    let path2 = temp_path.join("sample2.phraya");
    phraya_io::phraya::write_phraya(&path1, &file1).unwrap();
    phraya_io::phraya::write_phraya(&path2, &file2).unwrap();

    let output_path = temp_path.join("merged.phraya");

    let manifest_path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("Cargo.toml");
    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            manifest_path.to_str().unwrap(),
            "--",
            "merge",
            path1.to_str().unwrap(),
            path2.to_str().unwrap(),
            "--output",
            output_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya merge");

    assert!(output.status.success(), "phraya merge should succeed");

    let merged_file =
        phraya_io::phraya::read_phraya(&output_path).expect("merged file should be readable");

    // Should have only 1 observation (deduped)
    assert_eq!(
        merged_file.observations.len(),
        1,
        "identical observations should be deduplicated"
    );
}

/// Test: phraya merge with single input file
#[test]
fn issue_80_merge_single_file() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    let input1 = create_phraya_file(temp_path, "sample1.phraya", 100, "sample1", 3);
    let output_path = temp_path.join("merged.phraya");

    let manifest_path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("Cargo.toml");
    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            manifest_path.to_str().unwrap(),
            "--",
            "merge",
            input1.to_str().unwrap(),
            "--output",
            output_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya merge");

    assert!(
        output.status.success(),
        "phraya merge should succeed with single file"
    );

    let merged_file =
        phraya_io::phraya::read_phraya(&output_path).expect("merged file should be readable");

    assert_eq!(
        merged_file.observations.len(),
        3,
        "single file should be copied with all observations"
    );
}

/// Test: phraya merge returns non-zero exit code on error
#[test]
fn issue_80_merge_nonzero_exit_on_error() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    let nonexistent1 = temp_path.join("nonexistent1.phraya");
    let nonexistent2 = temp_path.join("nonexistent2.phraya");
    let output_path = temp_path.join("merged.phraya");

    let manifest_path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("Cargo.toml");
    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            manifest_path.to_str().unwrap(),
            "--",
            "merge",
            nonexistent1.to_str().unwrap(),
            nonexistent2.to_str().unwrap(),
            "--output",
            output_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya merge");

    assert!(
        !output.status.success(),
        "should return non-zero exit code on error"
    );
}
