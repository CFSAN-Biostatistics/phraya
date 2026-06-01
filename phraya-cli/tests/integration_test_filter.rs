use std::path::{Path, PathBuf};
use std::collections::HashMap;
use tempfile::TempDir;
use phraya_filter::FilterBuilder;
use phraya_core::types::VariantObservation;

#[test]
fn unit_test_filter_directly() {
    // Test that the filter works directly without file I/O
    let mut alleles = HashMap::new();
    alleles.insert(b'A', 5u32);
    let obs1 = VariantObservation::new(50, b'A', alleles.clone(), 0.95, "10M".to_string(), 60, 0, vec![5], 35.0, "sample:read0".to_string());

    let mut alleles = HashMap::new();
    alleles.insert(b'A', 15u32);
    let obs2 = VariantObservation::new(100, b'A', alleles.clone(), 0.95, "10M".to_string(), 60, 0, vec![15], 35.0, "sample:read1".to_string());

    let mut alleles = HashMap::new();
    alleles.insert(b'A', 10u32);
    let obs3 = VariantObservation::new(150, b'A', alleles.clone(), 0.95, "10M".to_string(), 60, 0, vec![10], 35.0, "sample:read2".to_string());

    let observations = vec![obs1, obs2, obs3];

    // Apply filter with min-coverage 10
    let filter = FilterBuilder::new().min_coverage(10).build();
    let filtered: Vec<_> = filter.filter(&observations).cloned().collect();

    // Should have 2 observations (positions 100 and 150, both have coverage >= 10)
    assert_eq!(filtered.len(), 2, "Filter should pass obs with coverage >= 10, got {} instead of 2", filtered.len());
    assert_eq!(filtered[0].position(), 100);
    assert_eq!(filtered[1].position(), 150);
}

#[test]
fn unit_test_observations_through_file_io() {
    // Test that observations preserve coverage values through file I/O
    use phraya_io::phraya::{PhrayaFile, write_phraya, read_phraya};
    use phraya_core::types::CoverageTrack;

    let temp_dir = TempDir::new().unwrap();
    let test_file = temp_dir.path().join("test.phraya");

    // Create observations with specific coverage values
    let mut alleles = HashMap::new();
    alleles.insert(b'A', 5u32);
    let obs1 = VariantObservation::new(50, b'A', alleles.clone(), 0.95, "10M".to_string(), 60, 0, vec![5], 35.0, "sample:read0".to_string());

    let mut alleles = HashMap::new();
    alleles.insert(b'A', 15u32);
    let obs2 = VariantObservation::new(100, b'A', alleles.clone(), 0.95, "10M".to_string(), 60, 0, vec![15], 35.0, "sample:read1".to_string());

    let observations = vec![obs1, obs2];

    // Write to file
    let coverage_track = CoverageTrack::new(vec![10; 200]);
    let phraya_file = PhrayaFile::new(200, "test_sample".to_string(), "2026-06-01T00:00:00Z".to_string(), observations, coverage_track);
    write_phraya(&test_file, &phraya_file).expect("Failed to write");

    // Read back from file
    let read_file = read_phraya(&test_file).expect("Failed to read");

    // Check that coverage values are preserved
    assert_eq!(read_file.observations.len(), 2);
    assert_eq!(read_file.observations[0].local_coverage()[0], 5, "First observation coverage should be 5, got {}", read_file.observations[0].local_coverage()[0]);
    assert_eq!(read_file.observations[1].local_coverage()[0], 15, "Second observation coverage should be 15, got {}", read_file.observations[1].local_coverage()[0]);
}

/// Helper to create a temporary .phraya file with observations
fn create_phraya_file(
    dir: &Path,
    filename: &str,
    observations: Vec<(u32, u8, HashMap<u8, u32>, u8, u32)>, // position, ref_base, alleles, mapq, coverage
    reference_length: u32,
) -> PathBuf {
    use phraya_core::types::{VariantObservation, CoverageTrack};
    use phraya_io::phraya::{PhrayaFile, write_phraya};

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

/// Test: phraya filter reads .phraya file and outputs VCF (default format)
#[test]
fn issue_85_filter_basic_vcf_output() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    // Create a test .phraya file with some observations
    let mut alleles = HashMap::new();
    alleles.insert(b'A', 10);
    alleles.insert(b'T', 5);

    let input_path = create_phraya_file(
        temp_path,
        "test.phraya",
        vec![
            (50, b'A', alleles.clone(), 60, 10),
            (100, b'C', {
                let mut a = HashMap::new();
                a.insert(b'C', 8);
                a.insert(b'G', 2);
                a
            }, 50, 10),
        ],
        200,
    );

    // Command: phraya filter test.phraya (no format specified, defaults to VCF)
    let manifest_path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("Cargo.toml");
    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            manifest_path.to_str().unwrap(),
            "--",
            "filter",
            input_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya filter");

    assert!(
        output.status.success(),
        "phraya filter should succeed. stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let vcf_output = String::from_utf8_lossy(&output.stdout);

    // Verify VCF header is present
    assert!(vcf_output.contains("##fileformat=VCFv4.2"), "VCF header should be present");
    assert!(vcf_output.contains("#CHROM"), "VCF column header should be present");

    // Verify at least one VCF record exists
    let record_lines: Vec<&str> = vcf_output
        .lines()
        .filter(|line| !line.starts_with("#") && !line.is_empty())
        .collect();
    assert!(record_lines.len() >= 2, "Should have at least 2 VCF records");
}

/// Test: phraya filter applies min-coverage threshold
#[test]
fn issue_85_filter_min_coverage_threshold() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    // Create observations with different coverages
    let input_path = create_phraya_file(
        temp_path,
        "test_coverage.phraya",
        vec![
            (50, b'A', {
                let mut a = HashMap::new();
                a.insert(b'A', 5);
                a
            }, 60, 5), // low coverage
            (100, b'A', {
                let mut a = HashMap::new();
                a.insert(b'A', 15);
                a
            }, 60, 15), // high coverage
            (150, b'A', {
                let mut a = HashMap::new();
                a.insert(b'A', 10);
                a
            }, 60, 10), // medium coverage
        ],
        200,
    );

    // Filter with min-coverage 10
    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "/home/crash/phraya/phraya-cli/Cargo.toml",
            "--",
            "filter",
            input_path.to_str().unwrap(),
            "--min-coverage",
            "10",
        ])
        .output()
        .expect("Failed to execute phraya filter");

    assert!(
        output.status.success(),
        "phraya filter should succeed with --min-coverage"
    );

    let vcf_output = String::from_utf8_lossy(&output.stdout);
    let record_lines: Vec<&str> = vcf_output
        .lines()
        .filter(|line| !line.starts_with("#") && !line.is_empty())
        .collect();

    // Should have 2 records (positions 100 and 150, both >= coverage 10)
    // Position 50 has coverage 5, should be filtered out
    assert_eq!(
        record_lines.len(),
        2,
        "Should have exactly 2 records after min-coverage=10 filter"
    );
}

/// Test: phraya filter applies min-mapq threshold
#[test]
fn issue_85_filter_min_mapq_threshold() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    let input_path = create_phraya_file(
        temp_path,
        "test_mapq.phraya",
        vec![
            (50, b'A', {
                let mut a = HashMap::new();
                a.insert(b'A', 10);
                a
            }, 20, 10), // low mapq
            (100, b'A', {
                let mut a = HashMap::new();
                a.insert(b'A', 10);
                a
            }, 50, 10), // good mapq
            (150, b'A', {
                let mut a = HashMap::new();
                a.insert(b'A', 10);
                a
            }, 60, 10), // high mapq
        ],
        200,
    );

    // Filter with min-mapq 40
    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "/home/crash/phraya/phraya-cli/Cargo.toml",
            "--",
            "filter",
            input_path.to_str().unwrap(),
            "--min-mapq",
            "40",
        ])
        .output()
        .expect("Failed to execute phraya filter");

    assert!(output.status.success(), "phraya filter should succeed with --min-mapq");

    let vcf_output = String::from_utf8_lossy(&output.stdout);
    let record_lines: Vec<&str> = vcf_output
        .lines()
        .filter(|line| !line.starts_with("#") && !line.is_empty())
        .collect();

    // Should have 2 records (positions 100 and 150, both have mapq >= 40)
    assert_eq!(
        record_lines.len(),
        2,
        "Should have exactly 2 records after min-mapq=40 filter"
    );
}

/// Test: phraya filter supports --format vcf (explicit)
#[test]
fn issue_85_filter_format_vcf_explicit() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    let input_path = create_phraya_file(
        temp_path,
        "test_format.phraya",
        vec![(50, b'A', {
            let mut a = HashMap::new();
            a.insert(b'A', 10);
            a
        }, 60, 10)],
        100,
    );

    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "/home/crash/phraya/phraya-cli/Cargo.toml",
            "--",
            "filter",
            input_path.to_str().unwrap(),
            "--format",
            "vcf",
        ])
        .output()
        .expect("Failed to execute phraya filter");

    assert!(output.status.success(), "phraya filter should support --format vcf");

    let output_str = String::from_utf8_lossy(&output.stdout);
    assert!(output_str.contains("##fileformat=VCFv4.2"), "Output should be VCF format");
}

/// Test: phraya filter supports --format tsv
#[test]
fn issue_85_filter_format_tsv() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    let input_path = create_phraya_file(
        temp_path,
        "test_tsv.phraya",
        vec![(50, b'A', {
            let mut a = HashMap::new();
            a.insert(b'A', 10);
            a
        }, 60, 10)],
        100,
    );

    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "/home/crash/phraya/phraya-cli/Cargo.toml",
            "--",
            "filter",
            input_path.to_str().unwrap(),
            "--format",
            "tsv",
        ])
        .output()
        .expect("Failed to execute phraya filter");

    assert!(
        output.status.success(),
        "phraya filter should support --format tsv"
    );

    let output_str = String::from_utf8_lossy(&output.stdout);

    // TSV should have header with columns
    let lines: Vec<&str> = output_str.lines().collect();
    assert!(!lines.is_empty(), "TSV output should not be empty");

    // First line should be header with tabs
    if !lines.is_empty() {
        assert!(
            lines[0].contains('\t'),
            "TSV output should have tab-separated columns"
        );
    }
}

/// Test: phraya filter --format phraya writes valid .phraya file
#[test]
fn issue_85_filter_format_phraya() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    let input_path = create_phraya_file(
        temp_path,
        "test_phraya_format.phraya",
        vec![(50, b'A', {
            let mut a = HashMap::new();
            a.insert(b'A', 10);
            a
        }, 60, 10)],
        100,
    );

    let output_path = temp_path.join("filtered.phraya");

    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "/home/crash/phraya/phraya-cli/Cargo.toml",
            "--",
            "filter",
            input_path.to_str().unwrap(),
            "--format",
            "phraya",
            "--output",
            output_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya filter");

    assert!(
        output.status.success(),
        "phraya filter should support --format phraya"
    );

    // Verify output file was created
    assert!(
        output_path.exists(),
        "filtered .phraya file should be created"
    );

    // Verify it's a valid .phraya file by reading it back
    use phraya_io::phraya::read_phraya;
    let result = read_phraya(&output_path);
    assert!(
        result.is_ok(),
        "filtered .phraya file should be readable as valid .phraya format"
    );

    let phraya_file = result.unwrap();
    assert!(
        phraya_file.observations.len() > 0,
        "filtered .phraya file should contain observations"
    );
}

/// Test: phraya filter combines multiple thresholds (AND logic)
#[test]
fn issue_85_filter_multiple_thresholds() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    let input_path = create_phraya_file(
        temp_path,
        "test_multi.phraya",
        vec![
            (50, b'A', {
                let mut a = HashMap::new();
                a.insert(b'A', 15);
                a
            }, 20, 15), // good coverage, low mapq
            (100, b'A', {
                let mut a = HashMap::new();
                a.insert(b'A', 5);
                a
            }, 50, 5), // low coverage, good mapq
            (150, b'A', {
                let mut a = HashMap::new();
                a.insert(b'A', 15);
                a
            }, 50, 15), // good coverage, good mapq
        ],
        200,
    );

    // Filter with both min-coverage 10 AND min-mapq 40
    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "/home/crash/phraya/phraya-cli/Cargo.toml",
            "--",
            "filter",
            input_path.to_str().unwrap(),
            "--min-coverage",
            "10",
            "--min-mapq",
            "40",
        ])
        .output()
        .expect("Failed to execute phraya filter");

    assert!(output.status.success(), "phraya filter should combine thresholds");

    let vcf_output = String::from_utf8_lossy(&output.stdout);
    let record_lines: Vec<&str> = vcf_output
        .lines()
        .filter(|line| !line.starts_with("#") && !line.is_empty())
        .collect();

    // Only position 150 should pass both filters
    assert_eq!(
        record_lines.len(),
        1,
        "Should have exactly 1 record passing both min-coverage=10 AND min-mapq=40"
    );
}

/// Test: phraya filter logs statistics to stderr
#[test]
fn issue_85_filter_logs_statistics() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    let input_path = create_phraya_file(
        temp_path,
        "test_stats.phraya",
        vec![
            (50, b'A', {
                let mut a = HashMap::new();
                a.insert(b'A', 10);
                a
            }, 60, 10),
            (100, b'A', {
                let mut a = HashMap::new();
                a.insert(b'A', 5);
                a
            }, 60, 5), // will be filtered
            (150, b'A', {
                let mut a = HashMap::new();
                a.insert(b'A', 10);
                a
            }, 60, 10),
        ],
        200,
    );

    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "/home/crash/phraya/phraya-cli/Cargo.toml",
            "--",
            "filter",
            input_path.to_str().unwrap(),
            "--min-coverage",
            "10",
        ])
        .output()
        .expect("Failed to execute phraya filter");

    assert!(output.status.success(), "phraya filter should succeed");

    let stderr = String::from_utf8_lossy(&output.stderr);

    // Should log "Filtered X → Y observations" pattern
    assert!(
        stderr.contains("Filtered") || stderr.contains("filtered"),
        "Should log filter statistics to stderr. stderr: {}",
        stderr
    );
}

/// Test: phraya filter returns non-zero exit code on error
#[test]
fn issue_85_filter_error_handling_nonexistent_file() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    let nonexistent_path = temp_path.join("nonexistent.phraya");

    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "/home/crash/phraya/phraya-cli/Cargo.toml",
            "--",
            "filter",
            nonexistent_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya filter");

    assert!(
        !output.status.success(),
        "phraya filter should fail with nonexistent input file"
    );

    // Should have error message
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        !stderr.is_empty(),
        "should print error message to stderr"
    );
}

/// Test: phraya filter supports chaining (filtered .phraya → filter)
#[test]
fn issue_85_filter_chaining_support() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    // Create initial .phraya file with 10 observations
    let input_path = create_phraya_file(
        temp_path,
        "test_chain_input.phraya",
        (0..10)
            .map(|i| {
                (
                    i * 10,
                    b'A',
                    {
                        let mut a = HashMap::new();
                        a.insert(b'A', 10 + i as u32);
                        a
                    },
                    30 + i as u8,
                    10 + i as u32,
                )
            })
            .collect(),
        150,
    );

    let filtered1_path = temp_path.join("filtered1.phraya");

    // First filter: min-coverage 15
    let output1 = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "/home/crash/phraya/phraya-cli/Cargo.toml",
            "--",
            "filter",
            input_path.to_str().unwrap(),
            "--min-coverage",
            "15",
            "--format",
            "phraya",
            "--output",
            filtered1_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute first phraya filter");

    assert!(output1.status.success(), "First filter should succeed");
    assert!(filtered1_path.exists(), "First filter should create output file");

    // Second filter: min-mapq 45, input is the output of first filter
    let filtered2_path = temp_path.join("filtered2.phraya");
    let output2 = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "/home/crash/phraya/phraya-cli/Cargo.toml",
            "--",
            "filter",
            filtered1_path.to_str().unwrap(),
            "--min-mapq",
            "45",
            "--format",
            "phraya",
            "--output",
            filtered2_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute second phraya filter");

    assert!(output2.status.success(), "Second filter should succeed");
    assert!(filtered2_path.exists(), "Second filter should create output file");

    // Verify result is valid .phraya
    use phraya_io::phraya::read_phraya;
    let result = read_phraya(&filtered2_path);
    assert!(result.is_ok(), "Chained filter output should be valid .phraya");
}

/// Test: phraya filter with no observations produces valid output
#[test]
fn issue_85_filter_empty_result() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    let input_path = create_phraya_file(
        temp_path,
        "test_empty.phraya",
        vec![(50, b'A', {
            let mut a = HashMap::new();
            a.insert(b'A', 5);
            a
        }, 20, 5)],
        100,
    );

    // Use very strict filter that filters out everything
    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "/home/crash/phraya/phraya-cli/Cargo.toml",
            "--",
            "filter",
            input_path.to_str().unwrap(),
            "--min-coverage",
            "100", // Very high, will filter out all observations
        ])
        .output()
        .expect("Failed to execute phraya filter");

    assert!(output.status.success(), "phraya filter should succeed even with no results");

    let vcf_output = String::from_utf8_lossy(&output.stdout);

    // Should have VCF header but no data records
    assert!(
        vcf_output.contains("##fileformat=VCFv4.2"),
        "Should have VCF header even with no observations"
    );

    let record_lines: Vec<&str> = vcf_output
        .lines()
        .filter(|line| !line.starts_with("#") && !line.is_empty())
        .collect();
    assert_eq!(
        record_lines.len(),
        0,
        "Should have 0 records when filter excludes all observations"
    );
}

/// Test: phraya filter CLI argument validation (missing required args)
#[test]
fn issue_85_filter_argument_validation() {
    // Test: missing input file argument
    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "/home/crash/phraya/phraya-cli/Cargo.toml",
            "--",
            "filter",
        ])
        .output()
        .expect("Failed to execute phraya filter");

    assert!(
        !output.status.success(),
        "phraya filter should fail without input file"
    );

    // Test: invalid format argument
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    let input_path = create_phraya_file(
        temp_path,
        "test_arg.phraya",
        vec![(50, b'A', {
            let mut a = HashMap::new();
            a.insert(b'A', 10);
            a
        }, 60, 10)],
        100,
    );

    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "/home/crash/phraya/phraya-cli/Cargo.toml",
            "--",
            "filter",
            input_path.to_str().unwrap(),
            "--format",
            "invalid_format",
        ])
        .output()
        .expect("Failed to execute phraya filter");

    assert!(
        !output.status.success(),
        "phraya filter should fail with invalid format"
    );
}

/// Test: phraya filter preserves observation order in output
#[test]
fn issue_85_filter_preserves_position_order() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    let input_path = create_phraya_file(
        temp_path,
        "test_order.phraya",
        vec![
            (10, b'A', {
                let mut a = HashMap::new();
                a.insert(b'A', 10);
                a
            }, 60, 10),
            (50, b'A', {
                let mut a = HashMap::new();
                a.insert(b'A', 10);
                a
            }, 60, 10),
            (100, b'A', {
                let mut a = HashMap::new();
                a.insert(b'A', 10);
                a
            }, 60, 10),
        ],
        150,
    );

    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "/home/crash/phraya/phraya-cli/Cargo.toml",
            "--",
            "filter",
            input_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya filter");

    assert!(output.status.success(), "Filter should succeed");

    let vcf_output = String::from_utf8_lossy(&output.stdout);
    let record_lines: Vec<&str> = vcf_output
        .lines()
        .filter(|line| !line.starts_with("#") && !line.is_empty())
        .collect();

    // Extract positions from VCF records (POS column is 1-indexed in VCF, 0-indexed in internal)
    let positions: Vec<u32> = record_lines
        .iter()
        .filter_map(|line| {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 2 {
                parts[1].parse::<u32>().ok()
            } else {
                None
            }
        })
        .collect();

    // Positions should be in ascending order: 11, 51, 101 (VCF 1-indexed)
    assert!(
        positions.windows(2).all(|w| w[0] <= w[1]),
        "Positions should be in ascending order"
    );
}
