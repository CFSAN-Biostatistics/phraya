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

/// Build a minimal plan file programmatically referencing the given FASTA.
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

// ── Tracer bullet ────────────────────────────────────────────────────────────

/// phraya align creates a .phraya output file for a known query/target pair
#[test]
fn issue_78_align_creates_phraya_file() {
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

    let status = std::process::Command::new("cargo")
        .args([
            "run",
            "--manifest-path",
            get_manifest_path().to_str().unwrap(),
            "--",
            "align",
            plan_path.to_str().unwrap(),
            "read1",
            "ref",
            "--output",
            output_path.to_str().unwrap(),
        ])
        .output()
        .expect("cargo run failed");

    assert!(
        status.status.success(),
        "phraya align should succeed.\nstderr: {}",
        String::from_utf8_lossy(&status.stderr)
    );
    assert!(
        output_path.exists(),
        ".phraya output file should be created"
    );
}

// ── Slice 2 ──────────────────────────────────────────────────────────────────

/// phraya align also writes the .phraya.queries sidecar
#[test]
fn issue_78_align_writes_queries_sidecar() {
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
    let queries_path = p.join("out.phraya.queries");

    let status = std::process::Command::new("cargo")
        .args([
            "run",
            "--manifest-path",
            get_manifest_path().to_str().unwrap(),
            "--",
            "align",
            plan_path.to_str().unwrap(),
            "read1",
            "ref",
            "--output",
            output_path.to_str().unwrap(),
        ])
        .output()
        .expect("cargo run failed");

    assert!(status.status.success());
    assert!(
        queries_path.exists(),
        ".phraya.queries sidecar should be written at {queries_path:?}"
    );
}

// ── Slice 3 ──────────────────────────────────────────────────────────────────

/// phraya align exits non-zero for an unknown query_id
#[test]
fn issue_78_align_unknown_query_id_fails() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();

    let fasta = create_fasta(
        p,
        "seqs.fa",
        &[("ref", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT")],
    );

    let plan_path = p.join("test.phrayaplan");
    write_test_plan(&plan_path, &fasta);

    let output_path = p.join("out.phraya");

    let status = std::process::Command::new("cargo")
        .args([
            "run",
            "--manifest-path",
            get_manifest_path().to_str().unwrap(),
            "--",
            "align",
            plan_path.to_str().unwrap(),
            "nonexistent_read",
            "ref",
            "--output",
            output_path.to_str().unwrap(),
        ])
        .output()
        .expect("cargo run failed");

    assert!(
        !status.status.success(),
        "phraya align should fail for unknown query_id"
    );

    let stderr = String::from_utf8_lossy(&status.stderr);
    assert!(
        stderr.contains("nonexistent_read") || stderr.contains("not found") || stderr.contains("unknown"),
        "error message should mention the bad ID: {stderr}"
    );
}
