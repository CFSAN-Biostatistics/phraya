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

/// Helper to create a plan file programmatically using the library
fn create_plan_file(dir: &Path, filename: &str, task_list: Vec<(u32, u32)>) -> PathBuf {
    use phraya_io::plan::{PhrayaPlan, UseCase};
    use std::collections::HashMap;

    let plan = PhrayaPlan::new(
        UseCase::ReadsWithRef,
        vec![],
        "2026-05-31T12:00:00Z".to_string(),
        vec![],
        HashMap::new(),
        task_list,
    );

    let path = dir.join(filename);
    phraya_io::plan::write_plan(&path, &plan).unwrap();
    path
}

/// Test: phraya plan-tasks reads plan file and outputs TSV
#[test]
#[ignore = "test: implement phraya plan-tasks CLI"]
fn issue_69_plan_tasks_basic_output() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    // Create a plan file with 3 tasks
    let plan_path = create_plan_file(temp_path, "test.phrayaplan", vec![(1, 0), (2, 0), (3, 0)]);

    // Command: phraya plan-tasks test.phrayaplan
    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "phraya-cli/Cargo.toml",
            "--",
            "plan-tasks",
            plan_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya plan-tasks");

    // Verify successful execution
    assert!(
        output.status.success(),
        "phraya plan-tasks should succeed. stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Verify stdout contains TSV output
    let stdout = String::from_utf8_lossy(&output.stdout);

    // Should have header
    assert!(
        stdout.contains("query_id") && stdout.contains("target_id"),
        "output should contain header with query_id and target_id. got: {}",
        stdout
    );

    // Should have 3 data lines
    let lines: Vec<&str> = stdout.lines().collect();
    assert!(
        lines.len() >= 4,
        "output should have at least 4 lines (1 header + 3 data). got {} lines: {}",
        lines.len(),
        stdout
    );
}

/// Test: output TSV format is correct (tab-separated)
#[test]
#[ignore = "test: implement phraya plan-tasks CLI"]
fn issue_69_plan_tasks_tsv_format() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    // Create a plan with specific task IDs
    let plan_path = create_plan_file(temp_path, "test.phrayaplan", vec![(5, 10), (7, 10)]);

    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "phraya-cli/Cargo.toml",
            "--",
            "plan-tasks",
            plan_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya plan-tasks");

    assert!(output.status.success(), "should succeed");

    let stdout = String::from_utf8_lossy(&output.stdout);
    let lines: Vec<&str> = stdout.lines().collect();

    // First line is header
    assert_eq!(
        lines[0].trim(),
        "query_id\ttarget_id",
        "first line should be header"
    );

    // Data lines should be tab-separated
    let data_line_1 = lines[1];
    let parts: Vec<&str> = data_line_1.split('\t').collect();
    assert_eq!(
        parts.len(),
        2,
        "each line should have exactly 2 tab-separated fields"
    );

    // Verify values match input tasks
    assert!(
        (parts[0] == "5" && parts[1] == "10") || (parts[0] == "7" && parts[1] == "10"),
        "data line values should match task IDs"
    );
}

/// Test: row count matches task count
#[test]
#[ignore = "test: implement phraya plan-tasks CLI"]
fn issue_69_plan_tasks_row_count_matches_task_count() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    // Create plan with 10 tasks
    let task_count = 10;
    let tasks: Vec<(u32, u32)> = (0..task_count).map(|i| (i as u32, 100)).collect();
    let plan_path = create_plan_file(temp_path, "test.phrayaplan", tasks);

    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "phraya-cli/Cargo.toml",
            "--",
            "plan-tasks",
            plan_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya plan-tasks");

    assert!(output.status.success(), "should succeed");

    let stdout = String::from_utf8_lossy(&output.stdout);
    let lines: Vec<&str> = stdout.lines().collect();

    // Should have header + 10 data lines = 11 total
    // (lines may include trailing newline, so >= check)
    let data_lines = lines.len() - 1; // subtract header
    assert_eq!(
        data_lines, task_count,
        "should have {} data lines (1 per task), got {}",
        task_count, data_lines
    );
}

/// Test: can be piped to wc -l for task count
#[test]
#[ignore = "test: implement phraya plan-tasks CLI"]
fn issue_69_plan_tasks_pipeline_wc_l() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    // Create plan with 5 tasks
    let plan_path = create_plan_file(
        temp_path,
        "test.phrayaplan",
        vec![(1, 0), (2, 0), (3, 0), (4, 0), (5, 0)],
    );

    // Execute: phraya plan-tasks test.phrayaplan | wc -l
    let phraya_output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "phraya-cli/Cargo.toml",
            "--",
            "plan-tasks",
            plan_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya plan-tasks");

    assert!(
        phraya_output.status.success(),
        "phraya plan-tasks should succeed"
    );

    let stdout = String::from_utf8_lossy(&phraya_output.stdout);

    // Count lines in the output (wc -l equivalent)
    let line_count = stdout.lines().count();

    // Should have 6 lines: 1 header + 5 data
    assert_eq!(
        line_count, 6,
        "output should have 6 lines total (header + 5 data), got {}",
        line_count
    );
}

/// Test: missing plan file returns non-zero exit code
#[test]
#[ignore = "test: implement phraya plan-tasks CLI"]
fn issue_69_plan_tasks_missing_file_error() {
    let nonexistent_path = "/tmp/nonexistent_plan_69_12345.phrayaplan";

    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "phraya-cli/Cargo.toml",
            "--",
            "plan-tasks",
            nonexistent_path,
        ])
        .output()
        .expect("Failed to execute phraya plan-tasks");

    // Verify non-zero exit code
    assert!(
        !output.status.success(),
        "should return non-zero exit code for missing file"
    );

    // Verify stderr contains error message
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        !stderr.is_empty(),
        "should write error message to stderr: {}",
        stderr
    );
}

/// Test: corrupt plan file returns non-zero exit code
#[test]
#[ignore = "test: implement phraya plan-tasks CLI"]
fn issue_69_plan_tasks_corrupt_plan_file_error() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    // Create corrupted plan file (invalid binary data)
    let corrupt_path = temp_path.join("corrupt.phrayaplan");
    fs::write(&corrupt_path, b"this is not a valid phrayaplan file").unwrap();

    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "phraya-cli/Cargo.toml",
            "--",
            "plan-tasks",
            corrupt_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya plan-tasks");

    // Verify non-zero exit code
    assert!(
        !output.status.success(),
        "should return non-zero exit code for corrupt file"
    );

    // Verify stderr contains informative error message
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        !stderr.is_empty(),
        "should provide error message: {}",
        stderr
    );
}

/// Test: empty task list produces only header
#[test]
#[ignore = "test: implement phraya plan-tasks CLI"]
fn issue_69_plan_tasks_empty_task_list() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    // Create plan with empty task list
    let plan_path = create_plan_file(temp_path, "test.phrayaplan", vec![]);

    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "phraya-cli/Cargo.toml",
            "--",
            "plan-tasks",
            plan_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya plan-tasks");

    assert!(
        output.status.success(),
        "should succeed even with empty task list"
    );

    let stdout = String::from_utf8_lossy(&output.stdout);
    let lines: Vec<&str> = stdout.lines().collect();

    // Should have at least header line
    assert!(lines.len() >= 1, "output should have at least header line");

    // First line should be header
    assert_eq!(lines[0].trim(), "query_id\ttarget_id", "should have header");

    // Should have exactly 1 line (header only)
    assert_eq!(
        lines.len(),
        1,
        "should have only header line for empty task list"
    );
}

/// Test: large task list is handled correctly
#[test]
#[ignore = "test: implement phraya plan-tasks CLI"]
fn issue_69_plan_tasks_large_task_list() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    // Create plan with 1000 tasks
    let task_count = 1000;
    let tasks: Vec<(u32, u32)> = (0..task_count).map(|i| (i as u32, 999)).collect();
    let plan_path = create_plan_file(temp_path, "test.phrayaplan", tasks);

    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "phraya-cli/Cargo.toml",
            "--",
            "plan-tasks",
            plan_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya plan-tasks");

    assert!(output.status.success(), "should handle large task lists");

    let stdout = String::from_utf8_lossy(&output.stdout);
    let lines: Vec<&str> = stdout.lines().collect();

    // Should have header + 1000 data lines
    assert_eq!(
        lines.len(),
        task_count + 1,
        "should output all {} tasks plus header",
        task_count
    );
}

/// Test: integration with phraya plan (plan → plan-tasks → verify output)
#[test]
#[ignore = "test: implement phraya plan-tasks CLI"]
fn issue_69_plan_tasks_integration_with_plan_command() {
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
            (
                "read3",
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
                "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
            ),
        ],
    );

    let plan_path = temp_path.join("test.phrayaplan");

    // First: phraya plan to generate the plan file
    let plan_output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "phraya-cli/Cargo.toml",
            "--",
            "plan",
            "--inputs",
            reads_path.to_str().unwrap(),
            "--reference",
            ref_path.to_str().unwrap(),
            "--output",
            plan_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya plan");

    assert!(plan_output.status.success(), "phraya plan should succeed");

    // Verify plan file was created
    assert!(plan_path.exists(), "plan file should exist");

    // Second: phraya plan-tasks to extract task list
    let tasks_output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "phraya-cli/Cargo.toml",
            "--",
            "plan-tasks",
            plan_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya plan-tasks");

    assert!(
        tasks_output.status.success(),
        "phraya plan-tasks should succeed. stderr: {}",
        String::from_utf8_lossy(&tasks_output.stderr)
    );

    // Verify output is valid TSV
    let stdout = String::from_utf8_lossy(&tasks_output.stdout);
    let lines: Vec<&str> = stdout.lines().collect();

    // Should have header + tasks
    assert!(lines.len() > 1, "output should have header and task lines");

    // Verify header
    assert_eq!(lines[0].trim(), "query_id\ttarget_id");

    // Verify row count matches plan's task count
    let plan = phraya_io::plan::read_plan(&plan_path).expect("plan should be readable");
    let expected_data_lines = plan.task_list.len();
    let actual_data_lines = lines.len() - 1; // subtract header
    assert_eq!(
        actual_data_lines, expected_data_lines,
        "TSV row count should match plan's task count"
    );
}

/// Test: stdin cannot be used (file argument required)
#[test]
#[ignore = "test: implement phraya plan-tasks CLI"]
fn issue_69_plan_tasks_requires_file_argument() {
    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "phraya-cli/Cargo.toml",
            "--",
            "plan-tasks",
        ])
        .output()
        .expect("Failed to execute phraya plan-tasks");

    // Should fail without file argument
    assert!(
        !output.status.success(),
        "should require plan file argument"
    );

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(!stderr.is_empty(), "should provide error message");
}

/// Test: header format is exactly "query_id\ttarget_id" (case-sensitive)
#[test]
#[ignore = "test: implement phraya plan-tasks CLI"]
fn issue_69_plan_tasks_header_exact_format() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    let plan_path = create_plan_file(temp_path, "test.phrayaplan", vec![(1, 0)]);

    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "phraya-cli/Cargo.toml",
            "--",
            "plan-tasks",
            plan_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya plan-tasks");

    assert!(output.status.success(), "should succeed");

    let stdout = String::from_utf8_lossy(&output.stdout);
    let first_line = stdout
        .lines()
        .next()
        .expect("should have at least one line");

    // Header must be exactly this string
    assert_eq!(
        first_line.trim(),
        "query_id\ttarget_id",
        "header must be exactly 'query_id\\ttarget_id'"
    );
}

/// Test: numeric output is plain integers (no leading zeros or extra formatting)
#[test]
#[ignore = "test: implement phraya plan-tasks CLI"]
fn issue_69_plan_tasks_numeric_format() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    let plan_path = create_plan_file(temp_path, "test.phrayaplan", vec![(100, 50)]);

    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "phraya-cli/Cargo.toml",
            "--",
            "plan-tasks",
            plan_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya plan-tasks");

    assert!(output.status.success(), "should succeed");

    let stdout = String::from_utf8_lossy(&output.stdout);
    let lines: Vec<&str> = stdout.lines().collect();
    let data_line = lines[1].trim();

    // Should be exactly "100\t50"
    assert_eq!(
        data_line, "100\t50",
        "numeric output should be plain integers"
    );
}

/// Test: all tasks are output, none are skipped
#[test]
#[ignore = "test: implement phraya plan-tasks CLI"]
fn issue_69_plan_tasks_all_tasks_present() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    // Create plan with specific task IDs to verify order
    let tasks = vec![(1, 10), (2, 10), (3, 10), (4, 10), (5, 10)];
    let plan_path = create_plan_file(temp_path, "test.phrayaplan", tasks.clone());

    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            "phraya-cli/Cargo.toml",
            "--",
            "plan-tasks",
            plan_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya plan-tasks");

    assert!(output.status.success(), "should succeed");

    let stdout = String::from_utf8_lossy(&output.stdout);
    let lines: Vec<&str> = stdout.lines().collect();

    // Verify all tasks are present
    for (i, task) in tasks.iter().enumerate() {
        let line = lines[i + 1].trim(); // skip header
        let parts: Vec<&str> = line.split('\t').collect();
        let query_id: u32 = parts[0].parse().expect("should parse query_id");
        let target_id: u32 = parts[1].parse().expect("should parse target_id");

        assert_eq!(
            (query_id, target_id),
            *task,
            "task {} should match input",
            i
        );
    }
}
