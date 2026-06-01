use std::fs;
use std::path::{Path, PathBuf};
use tempfile::TempDir;

/// Helper to create a temporary FASTA file
fn create_fasta_file(dir: &Path, filename: &str, sequences: &[(&str, &str)]) -> PathBuf {
    let path = dir.join(filename);
    let mut content = String::new();
    for (id, seq) in sequences {
        content.push('>');
        content.push_str(id);
        content.push('\n');
        content.push_str(seq);
        content.push('\n');
    }
    fs::write(&path, content).unwrap();
    path
}

/// Helper to create a temporary FASTQ file
fn create_fastq_file(dir: &Path, filename: &str, sequences: &[(&str, &str, &str)]) -> PathBuf {
    let path = dir.join(filename);
    let mut content = String::new();
    for (id, seq, qual) in sequences {
        content.push('@');
        content.push_str(id);
        content.push('\n');
        content.push_str(seq);
        content.push_str("\n+\n");
        content.push_str(qual);
        content.push('\n');
    }
    fs::write(&path, content).unwrap();
    path
}

/// Test: phraya plan with FASTA reference and FASTQ reads (Case 2)
#[test]
#[ignore = "test: implement phraya plan CLI"]
fn issue_68_plan_case2_reads_with_reference() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    // Create a reference FASTA
    let ref_path = create_fasta_file(
        temp_path,
        "reference.fa",
        &[(
            "ref",
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
        )],
    );

    // Create reads FASTQ
    let reads_path = create_fastq_file(
        temp_path,
        "reads.fq",
        &[
            (
                "read1",
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
                "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
            ),
            (
                "read2",
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
                "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
            ),
        ],
    );

    let output_path = temp_path.join("plan.phrayaplan");

    // Command: phraya plan --inputs reads.fq --reference reference.fa --output plan.phrayaplan
    // This test verifies the CLI parses arguments correctly and invokes the plan logic.
    // Expected: non-zero exit code on success (success == file written)
    // Expected: task list contains N=(num_reads) tasks since this is Case 2 (reads + ref)
    // Expected: use_case should be ReadsWithRef
    // Expected: output file exists and is valid .phrayaplan format

    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "/home/crash/phraya/phraya-cli/Cargo.toml",
            "--",
            "plan",
            "--inputs",
            reads_path.to_str().unwrap(),
            "--reference",
            ref_path.to_str().unwrap(),
            "--output",
            output_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya plan");

    // Verify successful execution
    assert!(
        output.status.success(),
        "phraya plan should succeed. stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Verify plan file was created
    assert!(
        output_path.exists(),
        "plan file should be created at {:?}",
        output_path
    );

    // Verify file is non-empty
    let plan_size = fs::metadata(&output_path).unwrap().len();
    assert!(plan_size > 0, "plan file should not be empty");

    // Parse the plan and verify use_case
    let plan = phraya_io::plan::read_plan(&output_path).expect("plan file should be readable");
    assert_eq!(
        plan.use_case,
        phraya_io::plan::UseCase::ReadsWithRef,
        "should detect Case 2: reads with reference"
    );

    // Verify task list has 2 tasks (one per read)
    assert_eq!(
        plan.task_list.len(),
        2,
        "task list should contain 2 tasks (one per read)"
    );

    // Verify input files are recorded
    assert!(
        plan.input_files
            .contains(&reads_path.to_string_lossy().to_string()),
        "input files should include reads path"
    );
    assert!(
        plan.input_files
            .contains(&ref_path.to_string_lossy().to_string()),
        "input files should include reference path"
    );
}

/// Test: phraya plan with contigs and reads, no reference (Case 3)
#[test]
#[ignore = "test: implement phraya plan CLI"]
fn issue_68_plan_case3_contigs_with_reads() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    // Create contigs FASTA
    let contigs_path = create_fasta_file(
        temp_path,
        "contigs.fa",
        &[
            (
                "contig1",
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            ),
            (
                "contig2",
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGA",
            ),
        ],
    );

    // Create reads FASTQ
    let reads_path = create_fastq_file(
        temp_path,
        "reads.fq",
        &[
            (
                "read1",
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
                "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
            ),
            (
                "read2",
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGA",
                "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
            ),
        ],
    );

    let output_path = temp_path.join("plan.phrayaplan");

    // Command: phraya plan --inputs contigs.fa reads.fq --output plan.phrayaplan
    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "/home/crash/phraya/phraya-cli/Cargo.toml",
            "--",
            "plan",
            "--inputs",
            contigs_path.to_str().unwrap(),
            "--inputs",
            reads_path.to_str().unwrap(),
            "--output",
            output_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya plan");

    // Verify successful execution
    assert!(
        output.status.success(),
        "phraya plan should succeed. stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Verify plan file was created
    assert!(output_path.exists(), "plan file should be created");

    // Parse the plan and verify use_case
    let plan = phraya_io::plan::read_plan(&output_path).expect("plan file should be readable");
    assert_eq!(
        plan.use_case,
        phraya_io::plan::UseCase::ContigsWithReads,
        "should detect Case 3: contigs with reads"
    );

    // Verify task list is populated (should have tasks for alignment)
    assert!(
        plan.task_list.len() > 0,
        "task list should not be empty for Case 3"
    );

    // Verify we have 4 sequences total (2 contigs + 2 reads)
    assert_eq!(
        plan.kmer_index.len(),
        4,
        "kmer_index should contain 4 sketches (2 contigs + 2 reads)"
    );
}

/// Test: phraya plan with contigs only (Case 4)
#[test]
#[ignore = "test: implement phraya plan CLI"]
fn issue_68_plan_case4_contigs_only() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    // Create contigs FASTA
    let contigs_path = create_fasta_file(
        temp_path,
        "contigs.fa",
        &[
            (
                "contig1",
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            ),
            (
                "contig2",
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGA",
            ),
            (
                "contig3",
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGC",
            ),
        ],
    );

    let output_path = temp_path.join("plan.phrayaplan");

    // Command: phraya plan --inputs contigs.fa --output plan.phrayaplan
    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "/home/crash/phraya/phraya-cli/Cargo.toml",
            "--",
            "plan",
            "--inputs",
            contigs_path.to_str().unwrap(),
            "--output",
            output_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya plan");

    // Verify successful execution
    assert!(
        output.status.success(),
        "phraya plan should succeed. stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Parse the plan and verify use_case
    let plan = phraya_io::plan::read_plan(&output_path).expect("plan file should be readable");
    assert_eq!(
        plan.use_case,
        phraya_io::plan::UseCase::ContigsOnly,
        "should detect Case 4: contigs only"
    );

    // Verify all contigs are sketched
    assert_eq!(
        plan.kmer_index.len(),
        3,
        "kmer_index should contain 3 sketches"
    );
}

/// Test: phraya plan with invalid input file
#[test]
#[ignore = "test: implement phraya plan CLI"]
fn issue_68_plan_invalid_input_file() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();
    let nonexistent_path = temp_path.join("nonexistent.fa");
    let output_path = temp_path.join("plan.phrayaplan");

    // Command: phraya plan --inputs nonexistent.fa --output plan.phrayaplan
    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "/home/crash/phraya/phraya-cli/Cargo.toml",
            "--",
            "plan",
            "--inputs",
            nonexistent_path.to_str().unwrap(),
            "--output",
            output_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya plan");

    // Verify it fails
    assert!(
        !output.status.success(),
        "phraya plan should fail with invalid input file"
    );

    // Verify error message is informative
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("invalid") || stderr.contains("not found") || stderr.contains("error"),
        "error message should be clear: {}",
        stderr
    );

    // Verify no plan file was created
    assert!(
        !output_path.exists(),
        "plan file should not be created on error"
    );
}

/// Test: phraya plan parses command-line arguments correctly
#[test]
#[ignore = "test: implement phraya plan CLI"]
fn issue_68_plan_cli_argument_parsing() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    // Create a minimal FASTA
    let reads_path = create_fasta_file(
        temp_path,
        "reads.fa",
        &[(
            "read1",
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
        )],
    );

    let output_path = temp_path.join("plan.phrayaplan");

    // Test: Missing required --output argument
    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "/home/crash/phraya/phraya-cli/Cargo.toml",
            "--",
            "plan",
            "--inputs",
            reads_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya plan");

    assert!(!output.status.success(), "should fail without --output");

    // Test: Missing required --inputs argument
    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "/home/crash/phraya/phraya-cli/Cargo.toml",
            "--",
            "plan",
            "--output",
            output_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya plan");

    assert!(!output.status.success(), "should fail without --inputs");
}

/// Test: phraya plan computes k-mer uniqueness
#[test]
#[ignore = "test: implement phraya plan CLI"]
fn issue_68_plan_kmer_uniqueness_computed() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    // Create sequences with varying k-mer uniqueness
    // Sequence 1 and 2 share common k-mers
    // Sequence 3 is unique
    let reads_path = create_fasta_file(
        temp_path,
        "reads.fa",
        &[
            (
                "read1",
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            ),
            (
                "read2",
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            ),
            (
                "read3",
                "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
            ),
        ],
    );

    let output_path = temp_path.join("plan.phrayaplan");

    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "/home/crash/phraya/phraya-cli/Cargo.toml",
            "--",
            "plan",
            "--inputs",
            reads_path.to_str().unwrap(),
            "--output",
            output_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya plan");

    assert!(output.status.success(), "should succeed");

    let plan = phraya_io::plan::read_plan(&output_path).expect("plan file should be readable");

    // Verify k-mer uniqueness is computed and populated
    assert!(
        !plan.kmer_uniqueness.is_empty(),
        "kmer_uniqueness should be computed and not empty"
    );

    // All uniqueness values should be in [0.0, 1.0]
    for (_, uniqueness) in plan.kmer_uniqueness.iter() {
        assert!(
            *uniqueness >= 0.0 && *uniqueness <= 1.0,
            "uniqueness score should be in [0.0, 1.0], got {}",
            uniqueness
        );
    }
}

/// Test: phraya plan logs detected use case to stderr
#[test]
#[ignore = "test: implement phraya plan CLI"]
fn issue_68_plan_logs_use_case() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    let ref_path = create_fasta_file(
        temp_path,
        "reference.fa",
        &[(
            "ref",
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
        )],
    );

    let reads_path = create_fastq_file(
        temp_path,
        "reads.fq",
        &[(
            "read1",
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        )],
    );

    let output_path = temp_path.join("plan.phrayaplan");

    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "/home/crash/phraya/phraya-cli/Cargo.toml",
            "--",
            "plan",
            "--inputs",
            reads_path.to_str().unwrap(),
            "--reference",
            ref_path.to_str().unwrap(),
            "--output",
            output_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya plan");

    assert!(output.status.success(), "should succeed");

    let stderr = String::from_utf8_lossy(&output.stderr);

    // Should log which case was detected
    assert!(
        stderr.contains("Case") || stderr.contains("case") || stderr.contains("Detected"),
        "should log detected use case to stderr: {}",
        stderr
    );
}

/// Test: phraya plan with multiple input files
#[test]
#[ignore = "test: implement phraya plan CLI"]
fn issue_68_plan_multiple_input_files() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    let ref_path = create_fasta_file(
        temp_path,
        "reference.fa",
        &[(
            "ref",
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
        )],
    );

    let reads1_path = create_fastq_file(
        temp_path,
        "reads1.fq",
        &[(
            "read1",
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        )],
    );

    let reads2_path = create_fastq_file(
        temp_path,
        "reads2.fq",
        &[(
            "read2",
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        )],
    );

    let output_path = temp_path.join("plan.phrayaplan");

    // Command with multiple --inputs
    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "/home/crash/phraya/phraya-cli/Cargo.toml",
            "--",
            "plan",
            "--inputs",
            reads1_path.to_str().unwrap(),
            "--inputs",
            reads2_path.to_str().unwrap(),
            "--reference",
            ref_path.to_str().unwrap(),
            "--output",
            output_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya plan");

    assert!(output.status.success(), "should succeed");

    let plan = phraya_io::plan::read_plan(&output_path).expect("plan file should be readable");

    // Should have 2 reads + 1 reference = 3 sketches
    assert_eq!(
        plan.kmer_index.len(),
        3,
        "should have sketches for both reads files and reference"
    );

    // Should have 2 tasks (one per read)
    assert_eq!(
        plan.task_list.len(),
        2,
        "should have task for each read file"
    );
}

/// Test: phraya plan generates valid task list
#[test]
#[ignore = "test: implement phraya plan CLI"]
fn issue_68_plan_task_list_valid() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    let ref_path = create_fasta_file(
        temp_path,
        "reference.fa",
        &[(
            "ref",
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
        )],
    );

    let reads_path = create_fastq_file(
        temp_path,
        "reads.fq",
        &[
            (
                "read1",
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
                "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
            ),
            (
                "read2",
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
                "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
            ),
            (
                "read3",
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
                "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
            ),
        ],
    );

    let output_path = temp_path.join("plan.phrayaplan");

    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "/home/crash/phraya/phraya-cli/Cargo.toml",
            "--",
            "plan",
            "--inputs",
            reads_path.to_str().unwrap(),
            "--reference",
            ref_path.to_str().unwrap(),
            "--output",
            output_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya plan");

    assert!(output.status.success(), "should succeed");

    let plan = phraya_io::plan::read_plan(&output_path).expect("plan file should be readable");

    // For Case 2 (reads + ref), each task should be (query_id, target_id)
    // where target_id is the reference (0) and query_id is each read (1, 2, 3, ...)
    assert_eq!(plan.task_list.len(), 3, "should have 3 tasks for 3 reads");

    for (query_id, target_id) in &plan.task_list {
        // Target should be reference (index 0)
        assert_eq!(*target_id, 0, "target should be reference");
        // Query IDs should be reads (index 1, 2, 3)
        assert!(*query_id >= 1, "query should be a read");
    }
}

/// Test: phraya plan returns non-zero exit code on success and writes file
#[test]
#[ignore = "test: implement phraya plan CLI"]
fn issue_68_plan_exit_code_on_success() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    let reads_path = create_fasta_file(
        temp_path,
        "reads.fa",
        &[(
            "read1",
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
        )],
    );

    let output_path = temp_path.join("plan.phrayaplan");

    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "/home/crash/phraya/phraya-cli/Cargo.toml",
            "--",
            "plan",
            "--inputs",
            reads_path.to_str().unwrap(),
            "--output",
            output_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya plan");

    // Success = exit code 0, not non-zero
    assert!(
        output.status.success(),
        "should exit with code 0 on success"
    );

    // Verify file was created
    assert!(output_path.exists(), "output file should exist");
}
