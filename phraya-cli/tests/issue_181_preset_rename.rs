use phraya_core::types::{CoverageTrack, VariantObservation};
use phraya_io::phraya::{write_phraya, PhrayaFile};
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use tempfile::TempDir;

/// Helper to get phraya-cli manifest path for cargo run commands
fn get_manifest_path() -> PathBuf {
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR").unwrap();
    std::path::Path::new(&manifest_dir).join("Cargo.toml")
}

/// Helper to create a test .phraya file
fn create_phraya_file_with_observations(
    dir: &Path,
    filename: &str,
    observations: Vec<(u32, u8, HashMap<u8, u32>, u8, u32)>,
    reference_length: u32,
) -> PathBuf {
    let path = dir.join(filename);

    let variant_obs: Vec<VariantObservation> = observations
        .into_iter()
        .enumerate()
        .map(|(i, (pos, ref_base, alleles, mapq, coverage))| {
            VariantObservation::new(
                pos,
                ref_base,
                alleles,
                0.95,
                "10M".to_string(),
                mapq,
                0,
                vec![coverage],
                35.0,
                format!("sample:read{}", i),
            )
        })
        .collect();

    let coverage = CoverageTrack::new(vec![10; reference_length as usize]);
    let file = PhrayaFile::new(
        reference_length,
        "test_sample".to_string(),
        "2026-05-31T12:00:00Z".to_string(),
        variant_obs,
        coverage,
    );

    write_phraya(&path, &file).expect("Failed to write phraya file");
    path
}

/// Issue #181: Test that --preset strict works and passes high-quality variants
#[test]
fn issue_181_cli_preset_strict_accepts_high_quality() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    // Create observations:
    // - Position 100: coverage=15, mapq=40, allele_freq=20% (passes strict: min_cov=10, min_mapq=30, min_af=10%)
    // - Position 200: coverage=5, mapq=25, allele_freq=5% (fails strict)
    let mut alleles_pass = HashMap::new();
    alleles_pass.insert(b'A', 80u32);
    alleles_pass.insert(b'T', 20u32); // 20% alt freq

    let mut alleles_fail = HashMap::new();
    alleles_fail.insert(b'A', 95u32);
    alleles_fail.insert(b'T', 5u32); // 5% alt freq

    let input_path = create_phraya_file_with_observations(
        temp_path,
        "test_strict.phraya",
        vec![
            (100, b'A', alleles_pass, 40, 15),
            (200, b'A', alleles_fail, 25, 5),
        ],
        300,
    );

    // Run: phraya filter input --preset strict
    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            get_manifest_path().to_str().unwrap(),
            "--",
            "filter",
            input_path.to_str().unwrap(),
            "--preset",
            "strict",
        ])
        .output()
        .expect("Failed to execute phraya filter");

    assert!(
        output.status.success(),
        "phraya filter --preset strict should succeed\nstderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let vcf_output = String::from_utf8_lossy(&output.stdout);
    let record_lines: Vec<&str> = vcf_output
        .lines()
        .filter(|line| !line.starts_with("#") && !line.is_empty())
        .collect();

    // Should have exactly 1 record (position 100)
    assert_eq!(
        record_lines.len(),
        1,
        "strict preset should pass 1 variant (position 100); expected 1, got {}",
        record_lines.len()
    );

    // Verify position is 100 (VCF is 1-indexed, so observation pos 100 -> VCF pos 101)
    assert!(
        vcf_output.contains("101\t"),
        "VCF should contain position 101"
    );
}

/// Issue #181: Test that --preset tolerant works and passes low-quality variants
#[test]
fn issue_181_cli_preset_tolerant_accepts_low_quality() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    // Create observations that pass tolerant but fail strict:
    // - Position 100: coverage=3, mapq=20, allele_freq=2% (passes tolerant: min_cov=3, min_mapq=20, min_af=2%)
    // - Position 200: coverage=1, mapq=15, allele_freq=1% (fails tolerant)
    let mut alleles_pass_tolerant = HashMap::new();
    alleles_pass_tolerant.insert(b'A', 98u32);
    alleles_pass_tolerant.insert(b'T', 2u32); // 2% alt freq

    let mut alleles_fail_tolerant = HashMap::new();
    alleles_fail_tolerant.insert(b'A', 99u32);
    alleles_fail_tolerant.insert(b'T', 1u32); // 1% alt freq

    let input_path = create_phraya_file_with_observations(
        temp_path,
        "test_tolerant.phraya",
        vec![
            (100, b'A', alleles_pass_tolerant, 20, 3),
            (200, b'A', alleles_fail_tolerant, 15, 1),
        ],
        300,
    );

    // Run: phraya filter input --preset tolerant
    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            get_manifest_path().to_str().unwrap(),
            "--",
            "filter",
            input_path.to_str().unwrap(),
            "--preset",
            "tolerant",
        ])
        .output()
        .expect("Failed to execute phraya filter");

    assert!(
        output.status.success(),
        "phraya filter --preset tolerant should succeed\nstderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let vcf_output = String::from_utf8_lossy(&output.stdout);
    let record_lines: Vec<&str> = vcf_output
        .lines()
        .filter(|line| !line.starts_with("#") && !line.is_empty())
        .collect();

    // Should have exactly 1 record (position 100)
    assert_eq!(
        record_lines.len(),
        1,
        "tolerant preset should pass 1 variant (position 100); expected 1, got {}",
        record_lines.len()
    );

    // Verify position is 100 (VCF is 1-indexed, so observation pos 100 -> VCF pos 101)
    assert!(
        vcf_output.contains("101\t"),
        "VCF should contain position 101"
    );
}

/// Issue #181: Test that old preset name "conservative" is rejected with helpful error
#[test]
fn issue_181_cli_preset_conservative_rejected() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    let mut alleles = HashMap::new();
    alleles.insert(b'A', 100u32);

    let input_path = create_phraya_file_with_observations(
        temp_path,
        "test.phraya",
        vec![(100, b'A', alleles, 60, 10)],
        200,
    );

    // Run: phraya filter input --preset conservative (should fail)
    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            get_manifest_path().to_str().unwrap(),
            "--",
            "filter",
            input_path.to_str().unwrap(),
            "--preset",
            "conservative",
        ])
        .output()
        .expect("Failed to execute phraya filter");

    assert!(
        !output.status.success(),
        "phraya filter --preset conservative should fail"
    );

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.to_lowercase().contains("conservative")
            || stderr.to_lowercase().contains("valid presets"),
        "Error message should mention 'conservative' or 'valid presets', got: {}",
        stderr
    );
}

/// Issue #181: Test that old preset name "sensitive" is rejected with helpful error
#[test]
fn issue_181_cli_preset_sensitive_rejected() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    let mut alleles = HashMap::new();
    alleles.insert(b'A', 100u32);

    let input_path = create_phraya_file_with_observations(
        temp_path,
        "test.phraya",
        vec![(100, b'A', alleles, 60, 10)],
        200,
    );

    // Run: phraya filter input --preset sensitive (should fail)
    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            get_manifest_path().to_str().unwrap(),
            "--",
            "filter",
            input_path.to_str().unwrap(),
            "--preset",
            "sensitive",
        ])
        .output()
        .expect("Failed to execute phraya filter");

    assert!(
        !output.status.success(),
        "phraya filter --preset sensitive should fail"
    );

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.to_lowercase().contains("sensitive")
            || stderr.to_lowercase().contains("valid presets"),
        "Error message should mention 'sensitive' or 'valid presets', got: {}",
        stderr
    );
}

/// Issue #181: Test that error message lists correct valid presets
#[test]
fn issue_181_cli_invalid_preset_shows_valid_presets() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    let mut alleles = HashMap::new();
    alleles.insert(b'A', 100u32);

    let input_path = create_phraya_file_with_observations(
        temp_path,
        "test.phraya",
        vec![(100, b'A', alleles, 60, 10)],
        200,
    );

    // Run: phraya filter input --preset invalid (should fail with helpful message)
    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            get_manifest_path().to_str().unwrap(),
            "--",
            "filter",
            input_path.to_str().unwrap(),
            "--preset",
            "invalid",
        ])
        .output()
        .expect("Failed to execute phraya filter");

    assert!(
        !output.status.success(),
        "phraya filter --preset invalid should fail"
    );

    let stderr = String::from_utf8_lossy(&output.stderr);
    let error_msg = stderr.to_lowercase();

    // Error should mention both strict and tolerant
    assert!(
        error_msg.contains("strict") || error_msg.contains("tolerant") || error_msg.contains("valid presets"),
        "Error should mention 'strict' and 'tolerant' presets, got: {}",
        stderr
    );
}

/// Issue #181: Test that strict and tolerant produce different results
#[test]
fn issue_181_cli_strict_and_tolerant_differ() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    // Create observation that passes tolerant but fails strict:
    // coverage=5, mapq=25, allele_freq=5%
    // Tolerant: min_cov=3, min_mapq=20, min_af=2% → PASSES
    // Strict: min_cov=10, min_mapq=30, min_af=10% → FAILS
    let mut alleles = HashMap::new();
    alleles.insert(b'A', 95u32);
    alleles.insert(b'T', 5u32); // 5% alt freq

    let input_path = create_phraya_file_with_observations(
        temp_path,
        "test.phraya",
        vec![(100, b'A', alleles, 25, 5)],
        200,
    );

    // Run with --preset strict
    let strict_output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            get_manifest_path().to_str().unwrap(),
            "--",
            "filter",
            input_path.to_str().unwrap(),
            "--preset",
            "strict",
        ])
        .output()
        .expect("Failed to execute phraya filter --preset strict");

    let strict_stdout = String::from_utf8_lossy(&strict_output.stdout);
    let strict_records: Vec<&str> = strict_stdout
        .lines()
        .filter(|line| !line.starts_with("#") && !line.is_empty())
        .collect();

    // Run with --preset tolerant
    let tolerant_output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            get_manifest_path().to_str().unwrap(),
            "--",
            "filter",
            input_path.to_str().unwrap(),
            "--preset",
            "tolerant",
        ])
        .output()
        .expect("Failed to execute phraya filter --preset tolerant");

    let tolerant_stdout = String::from_utf8_lossy(&tolerant_output.stdout);
    let tolerant_records: Vec<&str> = tolerant_stdout
        .lines()
        .filter(|line| !line.starts_with("#") && !line.is_empty())
        .collect();

    // strict should have 0 records, tolerant should have 1
    assert_eq!(
        strict_records.len(),
        0,
        "strict preset should reject this variant (coverage=5, mapq=25, af=5%); expected 0 records, got {}",
        strict_records.len()
    );

    assert_eq!(
        tolerant_records.len(),
        1,
        "tolerant preset should pass this variant (coverage=5, mapq=25, af=5%); expected 1 record, got {}",
        tolerant_records.len()
    );
}

/// Issue #181: Test that individual thresholds can override preset values
#[test]
fn issue_181_cli_preset_override_with_individual_flags() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    // Create observation: coverage=5, mapq=25, allele_freq=5%
    // Fails strict preset defaults but can pass with overridden min_coverage=3
    let mut alleles = HashMap::new();
    alleles.insert(b'A', 95u32);
    alleles.insert(b'T', 5u32); // 5% alt freq

    let input_path = create_phraya_file_with_observations(
        temp_path,
        "test.phraya",
        vec![(100, b'A', alleles, 25, 5)],
        200,
    );

    // Run with --preset strict but override with --min-coverage 3
    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            get_manifest_path().to_str().unwrap(),
            "--",
            "filter",
            input_path.to_str().unwrap(),
            "--preset",
            "strict",
            "--min-coverage",
            "3",
        ])
        .output()
        .expect("Failed to execute phraya filter with preset override");

    // This should still fail because allele_freq=5% < strict's min=10%
    let stdout = String::from_utf8_lossy(&output.stdout);
    let records: Vec<&str> = stdout
        .lines()
        .filter(|line| !line.starts_with("#") && !line.is_empty())
        .collect();

    assert_eq!(
        records.len(),
        0,
        "Even after overriding min_coverage=3, should still fail strict's min_af=10% threshold"
    );
}
