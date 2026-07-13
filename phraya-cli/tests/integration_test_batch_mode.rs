use std::fs;
use std::process::Command;
use tempfile::TempDir;

#[test]
fn test_batch_mode_end_to_end() {
    // Create temp directory for test
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    // Generate synthetic 100-read FASTQ
    let reads_file = temp_path.join("reads.fq");
    let mut reads_content = String::new();
    for i in 0..100 {
        reads_content.push_str(&format!("@read_{}\n", i));
        reads_content.push_str("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n");
        reads_content.push_str("+\n");
        reads_content.push_str("IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n");
    }
    fs::write(&reads_file, reads_content).unwrap();

    // Generate reference
    let ref_file = temp_path.join("ref.fa");
    fs::write(
        &ref_file,
        ">reference\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n",
    )
    .unwrap();

    // Run phraya plan with batch-to 4
    let plan_file = temp_path.join("plan.phrayaplan");
    let output_pattern = temp_path
        .join("out_{worker}.phraya")
        .to_str()
        .unwrap()
        .to_string();

    let plan_status = Command::new("cargo")
        .args([
            "run",
            "--bin",
            "phraya",
            "--",
            "plan",
            "--inputs",
            reads_file.to_str().unwrap(),
            "--reference",
            ref_file.to_str().unwrap(),
            "--output",
            plan_file.to_str().unwrap(),
            "--batch-to",
            "4",
            "--batch-output-pattern",
            &output_pattern,
        ])
        .status()
        .unwrap();
    assert!(plan_status.success(), "phraya plan failed");

    // Verify plan file exists
    assert!(plan_file.exists(), "Plan file not created");

    // Run align for workers 0-3
    for worker_id in 0..4 {
        let align_status = Command::new("cargo")
            .args([
                "run",
                "--bin",
                "phraya",
                "--",
                "align",
                plan_file.to_str().unwrap(),
                "--worker",
                &worker_id.to_string(),
            ])
            .status()
            .unwrap();
        assert!(
            align_status.success(),
            "phraya align --worker {} failed",
            worker_id
        );
    }

    // Verify all worker outputs exist
    for worker_id in 0..4 {
        let output_file = temp_path.join(format!("out_{}.phraya", worker_id));
        assert!(
            output_file.exists(),
            "Worker {} output not created",
            worker_id
        );
    }

    // Merge via plan
    let merged_file = temp_path.join("merged.phraya");
    let merge_status = Command::new("cargo")
        .args([
            "run",
            "--bin",
            "phraya",
            "--",
            "merge",
            plan_file.to_str().unwrap(),
            "--output",
            merged_file.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(merge_status.success(), "phraya merge failed");

    // Verify merged file exists
    assert!(merged_file.exists(), "Merged file not created");

    println!("Batch mode end-to-end test passed");
}

#[test]
fn test_batch_mode_ensure() {
    // Create temp directory for test
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    // Generate synthetic 50-read FASTQ
    let reads_file = temp_path.join("reads.fq");
    let mut reads_content = String::new();
    for i in 0..50 {
        reads_content.push_str(&format!("@read_{}\n", i));
        reads_content.push_str("ACGTACGTACGTACGTACGTACGTACGTACGT\n");
        reads_content.push_str("+\n");
        reads_content.push_str("IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n");
    }
    fs::write(&reads_file, reads_content).unwrap();

    // Generate reference
    let ref_file = temp_path.join("ref.fa");
    fs::write(
        &ref_file,
        ">reference\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n",
    )
    .unwrap();

    // Run phraya plan with batch-by 20
    let plan_file = temp_path.join("plan.phrayaplan");
    let output_pattern = temp_path
        .join("chunk_{worker}.phraya")
        .to_str()
        .unwrap()
        .to_string();

    let plan_status = Command::new("cargo")
        .args([
            "run",
            "--bin",
            "phraya",
            "--",
            "plan",
            "--inputs",
            reads_file.to_str().unwrap(),
            "--reference",
            ref_file.to_str().unwrap(),
            "--output",
            plan_file.to_str().unwrap(),
            "--batch-by",
            "20",
            "--batch-output-pattern",
            &output_pattern,
        ])
        .status()
        .unwrap();
    assert!(plan_status.success(), "phraya plan failed");

    // Run only worker 0 and 2
    for worker_id in [0, 2] {
        let align_status = Command::new("cargo")
            .args([
                "run",
                "--bin",
                "phraya",
                "--",
                "align",
                plan_file.to_str().unwrap(),
                "--worker",
                &worker_id.to_string(),
            ])
            .status()
            .unwrap();
        assert!(
            align_status.success(),
            "phraya align --worker {} failed",
            worker_id
        );
    }

    // Run ensure mode to process missing chunks
    let ensure_status = Command::new("cargo")
        .args([
            "run",
            "--bin",
            "phraya",
            "--",
            "align",
            plan_file.to_str().unwrap(),
            "--ensure",
        ])
        .status()
        .unwrap();
    assert!(ensure_status.success(), "phraya align --ensure failed");

    // Verify all outputs exist (3 chunks: 0-19, 20-39, 40-49)
    for worker_id in 0..3 {
        let output_file = temp_path.join(format!("chunk_{}.phraya", worker_id));
        assert!(
            output_file.exists(),
            "Worker {} output not created",
            worker_id
        );
    }

    println!("Batch mode ensure test passed");
}
