/// Integration tests for `phraya evidence` subcommand
///
/// These tests verify the two-phase workflow:
/// 1. Extract evidence with `phraya evidence`
/// 2. Use evidence with `phraya align --evidence`
///
/// All tests should FAIL until the feature is implemented.
use std::fs;
use std::path::PathBuf;
use std::process::Command;
use tempfile::TempDir;

/// Helper to get the phraya binary path
fn phraya_bin() -> PathBuf {
    let mut path = std::env::current_exe().unwrap();
    path.pop(); // Remove test binary name
    path.pop(); // Remove 'deps'
    path.push("phraya");
    path
}

/// Helper to create a temporary test directory
fn setup_temp_dir() -> TempDir {
    TempDir::new().expect("Failed to create temp dir")
}

/// Helper to create a minimal FASTA reference file
fn create_reference_fasta(dir: &TempDir) -> PathBuf {
    let path = dir.path().join("reference.fa");
    fs::write(
        &path,
        ">ref_seq1\nACGTACGTACGTACGTACGT\n>ref_seq2\nTTGGCCAAGGCCTTAA\n",
    )
    .expect("Failed to write reference");
    path
}

/// Helper to create a minimal assembly FASTA input
fn create_assembly_input(dir: &TempDir) -> PathBuf {
    let path = dir.path().join("assembly.fa");
    fs::write(
        &path,
        ">contig1\nACGTACGTACGTACGTACGT\n>contig2\nTTGGCCAAGGCCTTAA\n",
    )
    .expect("Failed to write assembly");
    path
}

/// Helper to create a minimal FASTQ read input
fn create_fastq_input(dir: &TempDir) -> PathBuf {
    let path = dir.path().join("reads.fq");
    fs::write(
        &path,
        "@read1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n@read2\nTTGGCCAAGGCC\n+\nIIIIIIIIIIII\n",
    )
    .expect("Failed to write FASTQ");
    path
}

#[test]
fn test_evidence_subcommand_exists() {
    // Test that `phraya evidence` is recognized as a valid subcommand
    let output = Command::new(phraya_bin())
        .arg("evidence")
        .arg("--help")
        .output()
        .expect("Failed to execute phraya");

    assert!(
        output.status.success(),
        "phraya evidence --help should succeed"
    );

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        stdout.contains("evidence"),
        "Help text should mention evidence subcommand"
    );
}

#[test]
fn test_evidence_requires_reference() {
    // Test that `phraya evidence` requires --reference flag
    let temp_dir = setup_temp_dir();
    let assembly = create_assembly_input(&temp_dir);
    let outdir = temp_dir.path().join("evidence_out");

    let output = Command::new(phraya_bin())
        .arg("evidence")
        .arg(assembly)
        .arg("--outdir")
        .arg(&outdir)
        .output()
        .expect("Failed to execute phraya");

    assert!(
        !output.status.success(),
        "Should fail without --reference flag"
    );

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("reference") || stderr.contains("required"),
        "Error should mention missing reference"
    );
}

#[test]
fn test_evidence_requires_outdir() {
    // Test that `phraya evidence` requires --outdir flag
    let temp_dir = setup_temp_dir();
    let reference = create_reference_fasta(&temp_dir);
    let assembly = create_assembly_input(&temp_dir);

    let output = Command::new(phraya_bin())
        .arg("evidence")
        .arg(assembly)
        .arg("--reference")
        .arg(reference)
        .output()
        .expect("Failed to execute phraya");

    assert!(
        !output.status.success(),
        "Should fail without --outdir flag"
    );

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("outdir") || stderr.contains("required"),
        "Error should mention missing outdir"
    );
}

#[test]
fn test_evidence_creates_evidence_json() {
    // Test that `phraya evidence` creates evidence.json in outdir
    let temp_dir = setup_temp_dir();
    let reference = create_reference_fasta(&temp_dir);
    let assembly = create_assembly_input(&temp_dir);
    let outdir = temp_dir.path().join("evidence_out");

    let output = Command::new(phraya_bin())
        .arg("evidence")
        .arg(assembly)
        .arg("--reference")
        .arg(&reference)
        .arg("--outdir")
        .arg(&outdir)
        .output()
        .expect("Failed to execute phraya");

    assert!(
        output.status.success(),
        "phraya evidence should succeed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let evidence_json = outdir.join("evidence.json");
    assert!(
        evidence_json.exists(),
        "evidence.json should be created in outdir"
    );

    // Verify it's valid JSON
    let content = fs::read_to_string(&evidence_json).expect("Failed to read evidence.json");
    assert!(
        serde_json::from_str::<serde_json::Value>(&content).is_ok(),
        "evidence.json should contain valid JSON"
    );
}

#[test]
fn test_evidence_accepts_multiple_inputs() {
    // Test that `phraya evidence` accepts multiple input files
    let temp_dir = setup_temp_dir();
    let reference = create_reference_fasta(&temp_dir);
    let assembly = create_assembly_input(&temp_dir);
    let reads = create_fastq_input(&temp_dir);
    let outdir = temp_dir.path().join("evidence_out");

    let output = Command::new(phraya_bin())
        .arg("evidence")
        .arg(&assembly)
        .arg(&reads)
        .arg("--reference")
        .arg(&reference)
        .arg("--outdir")
        .arg(&outdir)
        .output()
        .expect("Failed to execute phraya");

    assert!(
        output.status.success(),
        "phraya evidence should accept multiple inputs: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let evidence_json = outdir.join("evidence.json");
    assert!(
        evidence_json.exists(),
        "evidence.json should be created for multiple inputs"
    );
}

#[test]
fn test_align_accepts_evidence_flag() {
    // Test that `phraya align` accepts --evidence flag
    let temp_dir = setup_temp_dir();
    let reference = create_reference_fasta(&temp_dir);
    let assembly = create_assembly_input(&temp_dir);
    let evidence_dir = temp_dir.path().join("evidence_out");
    let align_outdir = temp_dir.path().join("align_out");

    // First, extract evidence
    let evidence_output = Command::new(phraya_bin())
        .arg("evidence")
        .arg(&assembly)
        .arg("--reference")
        .arg(&reference)
        .arg("--outdir")
        .arg(&evidence_dir)
        .output()
        .expect("Failed to extract evidence");

    assert!(
        evidence_output.status.success(),
        "Evidence extraction should succeed"
    );

    let evidence_json = evidence_dir.join("evidence.json");

    // Then, use evidence with align
    let output = Command::new(phraya_bin())
        .arg("align")
        .arg(&assembly)
        .arg("--evidence")
        .arg(&evidence_json)
        .arg("--outdir")
        .arg(&align_outdir)
        .output()
        .expect("Failed to execute phraya align");

    assert!(
        output.status.success(),
        "phraya align with --evidence should succeed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Verify output was created
    let result_file = align_outdir.join("results.phraya");
    assert!(
        result_file.exists(),
        "phraya align with --evidence should produce results.phraya"
    );
}

#[test]
fn test_align_works_without_evidence_flag() {
    // Test that `phraya align` still works without --evidence (single-phase mode)
    let temp_dir = setup_temp_dir();
    let reference = create_reference_fasta(&temp_dir);
    let assembly = create_assembly_input(&temp_dir);
    let outdir = temp_dir.path().join("align_out");

    let output = Command::new(phraya_bin())
        .arg("align")
        .arg(&assembly)
        .arg("--reference")
        .arg(&reference)
        .arg("--outdir")
        .arg(&outdir)
        .output()
        .expect("Failed to execute phraya align");

    assert!(
        output.status.success(),
        "phraya align without --evidence should succeed (auto-extract): {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Verify output was created
    let result_file = outdir.join("results.phraya");
    assert!(
        result_file.exists(),
        "phraya align should produce results.phraya"
    );
}

#[test]
fn test_two_phase_matches_single_phase() {
    // Test that two-phase workflow produces same results as single-phase
    let temp_dir = setup_temp_dir();
    let reference = create_reference_fasta(&temp_dir);
    let assembly = create_assembly_input(&temp_dir);

    // Single-phase workflow
    let single_phase_outdir = temp_dir.path().join("single_phase");
    let single_output = Command::new(phraya_bin())
        .arg("align")
        .arg(&assembly)
        .arg("--reference")
        .arg(&reference)
        .arg("--outdir")
        .arg(&single_phase_outdir)
        .output()
        .expect("Failed single-phase align");

    assert!(
        single_output.status.success(),
        "Single-phase should succeed"
    );

    // Two-phase workflow
    let evidence_dir = temp_dir.path().join("evidence");
    let two_phase_outdir = temp_dir.path().join("two_phase");

    // Phase 1: Extract evidence
    let evidence_output = Command::new(phraya_bin())
        .arg("evidence")
        .arg(&assembly)
        .arg("--reference")
        .arg(&reference)
        .arg("--outdir")
        .arg(&evidence_dir)
        .output()
        .expect("Failed evidence extraction");

    assert!(evidence_output.status.success(), "Evidence should succeed");

    // Phase 2: Align with evidence
    let two_phase_output = Command::new(phraya_bin())
        .arg("align")
        .arg(&assembly)
        .arg("--evidence")
        .arg(evidence_dir.join("evidence.json"))
        .arg("--outdir")
        .arg(&two_phase_outdir)
        .output()
        .expect("Failed two-phase align");

    assert!(
        two_phase_output.status.success(),
        "Two-phase should succeed"
    );

    // Compare results - both should produce .phraya output files
    let single_result = single_phase_outdir.join("results.phraya");
    let two_phase_result = two_phase_outdir.join("results.phraya");

    assert!(
        single_result.exists(),
        "Single-phase should produce results.phraya"
    );
    assert!(
        two_phase_result.exists(),
        "Two-phase should produce results.phraya"
    );

    // Results should be identical (bit-for-bit or logically equivalent)
    let single_content = fs::read(&single_result).expect("Failed to read single-phase result");
    let two_phase_content = fs::read(&two_phase_result).expect("Failed to read two-phase result");

    assert_eq!(
        single_content, two_phase_content,
        "Single-phase and two-phase results should be identical"
    );
}

#[test]
fn test_evidence_reused_across_multiple_align_calls() {
    // Test that evidence extracted once can be used for multiple align invocations
    let temp_dir = setup_temp_dir();
    let reference = create_reference_fasta(&temp_dir);
    let assembly1 = create_assembly_input(&temp_dir);
    let assembly2 = {
        let path = temp_dir.path().join("assembly2.fa");
        fs::write(&path, ">contig3\nGGCCTTAAGGCC\n").expect("Failed to write assembly2");
        path
    };

    let evidence_dir = temp_dir.path().join("evidence");

    // Extract evidence once with multiple inputs
    let evidence_output = Command::new(phraya_bin())
        .arg("evidence")
        .arg(&assembly1)
        .arg(&assembly2)
        .arg("--reference")
        .arg(&reference)
        .arg("--outdir")
        .arg(&evidence_dir)
        .output()
        .expect("Failed evidence extraction");

    assert!(evidence_output.status.success(), "Evidence should succeed");

    let evidence_json = evidence_dir.join("evidence.json");
    assert!(evidence_json.exists(), "evidence.json should exist");

    // Use the same evidence for multiple separate align calls
    let outdir1 = temp_dir.path().join("align_out1");
    let align1_output = Command::new(phraya_bin())
        .arg("align")
        .arg(&assembly1)
        .arg("--evidence")
        .arg(&evidence_json)
        .arg("--outdir")
        .arg(&outdir1)
        .output()
        .expect("Failed first align");

    assert!(
        align1_output.status.success(),
        "First align with evidence should succeed"
    );

    let outdir2 = temp_dir.path().join("align_out2");
    let align2_output = Command::new(phraya_bin())
        .arg("align")
        .arg(&assembly2)
        .arg("--evidence")
        .arg(&evidence_json)
        .arg("--outdir")
        .arg(&outdir2)
        .output()
        .expect("Failed second align");

    assert!(
        align2_output.status.success(),
        "Second align with same evidence should succeed"
    );

    // Both should produce output
    assert!(
        outdir1.join("results.phraya").exists(),
        "First align should produce output"
    );
    assert!(
        outdir2.join("results.phraya").exists(),
        "Second align should produce output"
    );
}

#[test]
fn test_align_fails_with_nonexistent_evidence_file() {
    // Test that `phraya align` fails gracefully with nonexistent evidence file
    let temp_dir = setup_temp_dir();
    let assembly = create_assembly_input(&temp_dir);
    let outdir = temp_dir.path().join("align_out");
    let nonexistent_evidence = temp_dir.path().join("nonexistent/evidence.json");

    let output = Command::new(phraya_bin())
        .arg("align")
        .arg(&assembly)
        .arg("--evidence")
        .arg(&nonexistent_evidence)
        .arg("--outdir")
        .arg(&outdir)
        .output()
        .expect("Failed to execute phraya align");

    assert!(
        !output.status.success(),
        "Should fail with nonexistent evidence file"
    );

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("evidence")
            || stderr.contains("not found")
            || stderr.contains("No such file"),
        "Error should mention evidence file problem"
    );
}

#[test]
fn test_align_fails_with_invalid_evidence_json() {
    // Test that `phraya align` fails gracefully with invalid JSON in evidence file
    let temp_dir = setup_temp_dir();
    let assembly = create_assembly_input(&temp_dir);
    let evidence_dir = temp_dir.path().join("evidence");
    fs::create_dir_all(&evidence_dir).expect("Failed to create evidence dir");

    let evidence_json = evidence_dir.join("evidence.json");
    fs::write(&evidence_json, "{ invalid json }").expect("Failed to write invalid JSON");

    let outdir = temp_dir.path().join("align_out");

    let output = Command::new(phraya_bin())
        .arg("align")
        .arg(&assembly)
        .arg("--evidence")
        .arg(&evidence_json)
        .arg("--outdir")
        .arg(&outdir)
        .output()
        .expect("Failed to execute phraya align");

    assert!(
        !output.status.success(),
        "Should fail with invalid evidence JSON"
    );

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("parse") || stderr.contains("invalid") || stderr.contains("JSON"),
        "Error should mention JSON parsing problem"
    );
}

#[test]
fn test_evidence_creates_outdir_if_not_exists() {
    // Test that `phraya evidence` creates outdir if it doesn't exist
    let temp_dir = setup_temp_dir();
    let reference = create_reference_fasta(&temp_dir);
    let assembly = create_assembly_input(&temp_dir);
    let outdir = temp_dir.path().join("nested/deep/evidence_out");

    assert!(!outdir.exists(), "Outdir should not exist initially");

    let output = Command::new(phraya_bin())
        .arg("evidence")
        .arg(assembly)
        .arg("--reference")
        .arg(&reference)
        .arg("--outdir")
        .arg(&outdir)
        .output()
        .expect("Failed to execute phraya");

    assert!(
        output.status.success(),
        "Should succeed and create nested outdir: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    assert!(outdir.exists(), "Outdir should be created");
    assert!(
        outdir.join("evidence.json").exists(),
        "evidence.json should exist in created outdir"
    );
}

#[test]
fn test_evidence_json_contains_expected_structure() {
    // Test that evidence.json has expected high-level structure
    let temp_dir = setup_temp_dir();
    let reference = create_reference_fasta(&temp_dir);
    let assembly = create_assembly_input(&temp_dir);
    let outdir = temp_dir.path().join("evidence_out");

    let output = Command::new(phraya_bin())
        .arg("evidence")
        .arg(assembly)
        .arg("--reference")
        .arg(&reference)
        .arg("--outdir")
        .arg(&outdir)
        .output()
        .expect("Failed to execute phraya");

    assert!(
        output.status.success(),
        "Evidence extraction should succeed"
    );

    let evidence_json = outdir.join("evidence.json");
    let content = fs::read_to_string(&evidence_json).expect("Failed to read evidence.json");
    let json: serde_json::Value =
        serde_json::from_str(&content).expect("Failed to parse evidence.json");

    // Evidence should have some expected top-level keys
    // (exact structure depends on implementation, but should have metadata)
    assert!(
        json.is_object(),
        "evidence.json should be a JSON object at top level"
    );

    // Should contain reference to the reference genome
    let json_str = json.to_string();
    assert!(
        json_str.contains("reference") || json_str.contains("ref"),
        "evidence.json should reference the reference genome"
    );
}

#[test]
fn test_align_with_evidence_skips_extraction() {
    // Test that when --evidence is provided, align skips evidence extraction
    // This is a behavioral test - implementation should log or indicate skipping
    let temp_dir = setup_temp_dir();
    let reference = create_reference_fasta(&temp_dir);
    let assembly = create_assembly_input(&temp_dir);
    let evidence_dir = temp_dir.path().join("evidence");
    let align_outdir = temp_dir.path().join("align_out");

    // Extract evidence
    let evidence_output = Command::new(phraya_bin())
        .arg("evidence")
        .arg(&assembly)
        .arg("--reference")
        .arg(&reference)
        .arg("--outdir")
        .arg(&evidence_dir)
        .output()
        .expect("Failed evidence extraction");

    assert!(evidence_output.status.success());

    let evidence_json = evidence_dir.join("evidence.json");

    // Align with evidence - should be faster/skip extraction
    // (We can't directly measure "skipping", but it should succeed)
    let align_output = Command::new(phraya_bin())
        .arg("align")
        .arg(&assembly)
        .arg("--evidence")
        .arg(&evidence_json)
        .arg("--outdir")
        .arg(&align_outdir)
        .output()
        .expect("Failed align with evidence");

    assert!(
        align_output.status.success(),
        "Align with evidence should succeed"
    );

    // Verify output was created
    let result_file = align_outdir.join("results.phraya");
    assert!(
        result_file.exists(),
        "Align with evidence should produce results.phraya"
    );

    // If verbose logging is available, could check for "using pre-computed evidence"
    let stderr = String::from_utf8_lossy(&align_output.stderr);
    // This is optional - depends on implementation verbosity
    // Just verify it succeeds without errors
    assert!(
        !stderr.contains("error") && !stderr.contains("failed"),
        "Should not have errors when using evidence"
    );
}
