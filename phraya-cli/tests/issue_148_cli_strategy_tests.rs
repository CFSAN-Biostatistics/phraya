/// Issue #148: CLI tests for --strategy flag
///
/// These RED acceptance tests verify that the CLI accepts and correctly propagates
/// the --strategy flag to the align command, and that different strategies produce
/// different local_coverage window sizes in the output.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use tempfile::TempDir;

fn get_manifest_path() -> PathBuf {
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR").unwrap();
    Path::new(&manifest_dir).join("Cargo.toml")
}

fn create_fasta(dir: &Path, filename: &str, sequences: &[(&str, &str)]) -> PathBuf {
    let path = dir.join(filename);
    let mut content = String::new();
    for (id, seq) in sequences {
        content.push_str(&format!(">{id}\n{seq}\n"));
    }
    std::fs::write(&path, content).unwrap();
    path
}

fn write_test_plan(plan_path: &Path, fasta_path: &Path) {
    use phraya_io::plan::{write_plan, PhrayaPlan, UseCase};
    let plan = PhrayaPlan::new(
        UseCase::ReadsWithRef,
        vec![fasta_path.to_string_lossy().to_string()],
        "2026-06-02T00:00:00Z".to_string(),
        HashMap::new(),
        HashMap::new(),
        vec![(1, 0)],
    );
    write_plan(plan_path, &plan).unwrap();
}

/// Test that phraya align accepts --strategy flag with value "fast".
/// The command must parse and recognize the flag without error.
#[test]
fn issue_148_align_accepts_strategy_fast() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();

    let fasta = create_fasta(
        p,
        "seqs.fa",
        &[
            ("ref", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"),
            ("read1", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"),
        ],
    );

    let plan_path = p.join("test.phrayaplan");
    write_test_plan(&plan_path, &fasta);

    let output_path = p.join("out.phraya");

    // TODO: Once --strategy flag is added to align subcommand:
    // let status = std::process::Command::new("cargo")
    //     .args([
    //         "run",
    //         "--manifest-path",
    //         get_manifest_path().to_str().unwrap(),
    //         "--",
    //         "align",
    //         plan_path.to_str().unwrap(),
    //         "read1",
    //         "ref",
    //         "--output",
    //         output_path.to_str().unwrap(),
    //         "--strategy",
    //         "fast",
    //     ])
    //     .output()
    //     .expect("cargo run failed");
    //
    // assert!(
    //     status.status.success(),
    //     "phraya align --strategy fast should succeed.\nstderr: {}",
    //     String::from_utf8_lossy(&status.stderr)
    // );

    assert!(
        false,
        "CLI align subcommand must accept --strategy flag with values: fast, balanced, exact"
    );
}

#[test]
fn issue_148_align_accepts_strategy_balanced() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();

    let fasta = create_fasta(
        p,
        "seqs.fa",
        &[
            ("ref", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"),
            ("read1", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"),
        ],
    );

    let plan_path = p.join("test.phrayaplan");
    write_test_plan(&plan_path, &fasta);

    let output_path = p.join("out.phraya");

    // TODO: Once --strategy flag is added:
    // let status = std::process::Command::new("cargo")
    //     .args([
    //         "run",
    //         "--manifest-path",
    //         get_manifest_path().to_str().unwrap(),
    //         "--",
    //         "align",
    //         plan_path.to_str().unwrap(),
    //         "read1",
    //         "ref",
    //         "--output",
    //         output_path.to_str().unwrap(),
    //         "--strategy",
    //         "balanced",
    //     ])
    //     .output()
    //     .expect("cargo run failed");
    //
    // assert!(
    //     status.status.success(),
    //     "phraya align --strategy balanced should succeed.\nstderr: {}",
    //     String::from_utf8_lossy(&status.stderr)
    // );

    assert!(
        false,
        "CLI align subcommand must accept --strategy balanced"
    );
}

#[test]
fn issue_148_align_accepts_strategy_exact() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();

    let fasta = create_fasta(
        p,
        "seqs.fa",
        &[
            ("ref", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"),
            ("read1", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"),
        ],
    );

    let plan_path = p.join("test.phrayaplan");
    write_test_plan(&plan_path, &fasta);

    let output_path = p.join("out.phraya");

    // TODO: Once --strategy flag is added:
    // let status = std::process::Command::new("cargo")
    //     .args([
    //         "run",
    //         "--manifest-path",
    //         get_manifest_path().to_str().unwrap(),
    //         "--",
    //         "align",
    //         plan_path.to_str().unwrap(),
    //         "read1",
    //         "ref",
    //         "--output",
    //         output_path.to_str().unwrap(),
    //         "--strategy",
    //         "exact",
    //     ])
    //     .output()
    //     .expect("cargo run failed");
    //
    // assert!(
    //     status.status.success(),
    //     "phraya align --strategy exact should succeed.\nstderr: {}",
    //     String::from_utf8_lossy(&status.stderr)
    // );

    assert!(
        false,
        "CLI align subcommand must accept --strategy exact"
    );
}

/// Test that phraya align rejects invalid strategy values.
/// Only "fast", "balanced", "exact" are valid; others must be rejected with a clear error.
#[test]
fn issue_148_align_rejects_invalid_strategy() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();

    let fasta = create_fasta(
        p,
        "seqs.fa",
        &[
            ("ref", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"),
            ("read1", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"),
        ],
    );

    let plan_path = p.join("test.phrayaplan");
    write_test_plan(&plan_path, &fasta);

    let output_path = p.join("out.phraya");

    // TODO: Once --strategy flag is added:
    // let status = std::process::Command::new("cargo")
    //     .args([
    //         "run",
    //         "--manifest-path",
    //         get_manifest_path().to_str().unwrap(),
    //         "--",
    //         "align",
    //         plan_path.to_str().unwrap(),
    //         "read1",
    //         "ref",
    //         "--output",
    //         output_path.to_str().unwrap(),
    //         "--strategy",
    //         "invalid_strategy",
    //     ])
    //     .output()
    //     .expect("cargo run failed");
    //
    // assert!(
    //     !status.status.success(),
    //     "phraya align should reject invalid strategy value"
    // );
    //
    // let stderr = String::from_utf8_lossy(&status.stderr);
    // assert!(
    //     stderr.contains("invalid") || stderr.contains("strategy") || stderr.contains("one of"),
    //     "error message should mention invalid strategy: {stderr}"
    // );

    assert!(
        false,
        "CLI must validate --strategy values and reject invalid ones"
    );
}

/// Test that omitting --strategy defaults to "balanced".
/// The current behavior (±50bp window) must be preserved when the flag is not provided.
#[test]
fn issue_148_align_uses_default_strategy_without_flag() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();

    let fasta = create_fasta(
        p,
        "seqs.fa",
        &[
            ("ref", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"),
            ("read1", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"),
        ],
    );

    let plan_path = p.join("test.phrayaplan");
    write_test_plan(&plan_path, &fasta);

    let output_path = p.join("out.phraya");

    // TODO: Once --strategy flag is added:
    // let status = std::process::Command::new("cargo")
    //     .args([
    //         "run",
    //         "--manifest-path",
    //         get_manifest_path().to_str().unwrap(),
    //         "--",
    //         "align",
    //         plan_path.to_str().unwrap(),
    //         "read1",
    //         "ref",
    //         "--output",
    //         output_path.to_str().unwrap(),
    //         // No --strategy flag: should default to balanced
    //     ])
    //     .output()
    //     .expect("cargo run failed");
    //
    // assert!(
    //     status.status.success(),
    //     "phraya align without --strategy should succeed with default.\nstderr: {}",
    //     String::from_utf8_lossy(&status.stderr)
    // );
    //
    // assert!(
    //     output_path.exists(),
    //     ".phraya output file should be created with default strategy"
    // );

    assert!(
        false,
        "Omitting --strategy must default to balanced (±50bp) for backward compatibility"
    );
}

/// Integration test: align with different strategies and verify window sizes differ.
/// Run alignment with --strategy fast and --strategy balanced, then compare the
/// resulting .phraya files' local_coverage window sizes.
#[test]
fn issue_148_cli_different_strategies_produce_different_windows() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();

    // Create a test sequence with a SNP at position 50 (in a 100bp read vs 200bp ref)
    let ref_seq = "A".repeat(200);
    let mut read_seq = "A".repeat(100);
    // Introduce a SNP at position 50 (ref: A -> C, read: T)
    unsafe {
        // SAFETY: we're modifying a string we just created
        read_seq.as_bytes_mut()[50] = b'T';
    }

    let fasta = create_fasta(
        p,
        "seqs.fa",
        &[
            ("ref", &ref_seq),
            ("read1", &read_seq),
        ],
    );

    let plan_path = p.join("test.phrayaplan");
    write_test_plan(&plan_path, &fasta);

    let output_fast = p.join("out_fast.phraya");
    let output_balanced = p.join("out_balanced.phraya");

    // TODO: Once --strategy flag is fully integrated and working:
    // Run alignment with fast strategy
    // let status_fast = std::process::Command::new("cargo")
    //     .args([
    //         "run",
    //         "--manifest-path",
    //         get_manifest_path().to_str().unwrap(),
    //         "--",
    //         "align",
    //         plan_path.to_str().unwrap(),
    //         "read1",
    //         "ref",
    //         "--output",
    //         output_fast.to_str().unwrap(),
    //         "--strategy",
    //         "fast",
    //     ])
    //     .output()
    //     .expect("cargo run failed");
    //
    // assert!(status_fast.status.success(), "fast alignment should succeed");
    //
    // Run alignment with balanced strategy
    // let status_balanced = std::process::Command::new("cargo")
    //     .args([
    //         "run",
    //         "--manifest-path",
    //         get_manifest_path().to_str().unwrap(),
    //         "--",
    //         "align",
    //         plan_path.to_str().unwrap(),
    //         "read1",
    //         "ref",
    //         "--output",
    //         output_balanced.to_str().unwrap(),
    //         "--strategy",
    //         "balanced",
    //     ])
    //     .output()
    //     .expect("cargo run failed");
    //
    // assert!(status_balanced.status.success(), "balanced alignment should succeed");
    //
    // Read both .phraya files and compare local_coverage window sizes
    // let phraya_fast = phraya_io::phraya::read_phraya(&output_fast)
    //     .expect("should read fast .phraya");
    // let phraya_balanced = phraya_io::phraya::read_phraya(&output_balanced)
    //     .expect("should read balanced .phraya");
    //
    // assert!(!phraya_fast.observations.is_empty(), "fast alignment should produce variants");
    // assert!(!phraya_balanced.observations.is_empty(), "balanced alignment should produce variants");
    //
    // Find the variant at position 50 in both files
    // let var_fast = phraya_fast.observations
    //     .iter()
    //     .find(|v| v.position() == 50)
    //     .expect("variant at position 50 in fast alignment");
    // let var_balanced = phraya_balanced.observations
    //     .iter()
    //     .find(|v| v.position() == 50)
    //     .expect("variant at position 50 in balanced alignment");
    //
    // Compare window sizes: fast should be wider (more positions) than balanced
    // let window_fast = var_fast.local_coverage().len();
    // let window_balanced = var_balanced.local_coverage().len();
    //
    // assert!(
    //     window_fast > window_balanced,
    //     "fast strategy window ({}) should be larger than balanced window ({})",
    //     window_fast,
    //     window_balanced
    // );
    // Fast: ±150bp → window length ≈ 201 (or 200 if clamped)
    // Balanced: ±50bp → window length = 101
    // assert_eq!(
    //     window_fast, 200,
    //     "fast strategy at pos 50 should produce window length 200 (±150bp with clamping)"
    // );
    // assert_eq!(
    //     window_balanced, 101,
    //     "balanced strategy at pos 50 should produce window length 101 (±50bp)"
    // );

    assert!(
        false,
        "Different strategies must produce different local_coverage window sizes in output"
    );
}
