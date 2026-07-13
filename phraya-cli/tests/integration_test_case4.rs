use std::path::{Path, PathBuf};
use tempfile::TempDir;

fn manifest() -> PathBuf {
    let d = std::env::var("CARGO_MANIFEST_DIR").unwrap();
    Path::new(&d).join("Cargo.toml")
}

/// Write a FASTA file where each sequence is `seq_len` bp (all 'A').
/// seq_len ≥ 5000 → classify_input treats them as contigs.
fn write_fasta(dir: &Path, name: &str, ids: &[&str], seq_len: usize) -> PathBuf {
    let path = dir.join(name);
    let seq = "ACGTACGTACGTACGT".repeat(seq_len / 16 + 1);
    let seq = &seq[..seq_len];
    let mut content = String::new();
    for id in ids {
        content.push_str(&format!(">{id}\n{seq}\n"));
    }
    std::fs::write(&path, content).unwrap();
    path
}

fn phraya(args: &[&str]) -> std::process::Output {
    std::process::Command::new("cargo")
        .arg("run")
        .arg("--manifest-path")
        .arg(manifest().to_str().unwrap())
        .arg("--")
        .args(args)
        .output()
        .expect("cargo run failed")
}

// ── Tracer bullet ────────────────────────────────────────────────────────────

/// phraya plan with contigs + reference stores UseCase::ContigsOnly, not ReadsWithRef
#[test]
fn issue_87_case4_with_ref_use_case_is_contigs_only() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();

    // 5 contig-length sequences (≥5kb each → classified as contigs)
    let contigs = write_fasta(
        p,
        "contigs.fa",
        &["ctg1", "ctg2", "ctg3", "ctg4", "ctg5"],
        5120,
    );
    let reference = write_fasta(p, "ref.fa", &["ref"], 5120);
    let plan_path = p.join("case4.phrayaplan");

    let out = phraya(&[
        "plan",
        "--inputs",
        contigs.to_str().unwrap(),
        "--reference",
        reference.to_str().unwrap(),
        "--output",
        plan_path.to_str().unwrap(),
    ]);
    assert!(
        out.status.success(),
        "phraya plan failed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );

    let plan = phraya_io::plan::read_plan(&plan_path).unwrap();
    assert_eq!(
        plan.use_case,
        phraya_io::plan::UseCase::ContigsOnly,
        "contigs + reference should produce ContigsOnly, got {:?}",
        plan.use_case
    );
}

// ── Slice 2 ──────────────────────────────────────────────────────────────────

/// phraya plan case4 with ref: task list has exactly M tasks, all targeting the reference
#[test]
fn issue_87_case4_with_ref_generates_m_tasks_targeting_ref() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();

    let contigs = write_fasta(
        p,
        "contigs.fa",
        &["ctg1", "ctg2", "ctg3", "ctg4", "ctg5"],
        5120,
    );
    let reference = write_fasta(p, "ref.fa", &["ref"], 5120);
    let plan_path = p.join("case4.phrayaplan");

    phraya(&[
        "plan",
        "--inputs",
        contigs.to_str().unwrap(),
        "--reference",
        reference.to_str().unwrap(),
        "--output",
        plan_path.to_str().unwrap(),
    ]);

    let plan = phraya_io::plan::read_plan(&plan_path).unwrap();

    assert_eq!(
        plan.task_list.len(),
        5,
        "5 contigs + reference → 5 tasks, got {}",
        plan.task_list.len()
    );

    // Reference is stored first (index 0); all tasks should target it
    for (query_id, target_id) in &plan.task_list {
        assert_eq!(
            *target_id, 0,
            "all tasks must target the reference (index 0), got target={target_id}"
        );
        assert_ne!(*query_id, 0, "query must not be the reference itself");
    }
}

// ── Slice 3 ──────────────────────────────────────────────────────────────────

/// phraya align runs successfully on a Case 4 (contig → reference) task
#[test]
fn issue_87_case4_align_contig_to_reference() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();

    // Use short sequences so the test runs fast; plan is written programmatically
    // so classify_input doesn't gate us here.
    let seqs = write_fasta(
        p,
        "seqs.fa",
        &["ref", "ctg1", "ctg2"],
        64, // short for speed; we're testing CLI wiring, not classification
    );
    let plan_path = p.join("case4.phrayaplan");

    // Write plan directly: ContigsOnly, tasks = (1,0) and (2,0)
    {
        use phraya_io::plan::{write_plan, PhrayaPlan, UseCase};
        use std::collections::HashMap;
        let plan = PhrayaPlan::new(
            UseCase::ContigsOnly,
            vec![seqs.to_string_lossy().to_string()],
            "2026-06-02T00:00:00Z".to_string(),
            HashMap::new(),
            HashMap::new(),
            vec![(1, 0), (2, 0)],
        );
        write_plan(&plan_path, &plan).unwrap();
    }

    let output_path = p.join("ctg1.phraya");

    let out = phraya(&[
        "align",
        plan_path.to_str().unwrap(),
        "ctg1",
        "ref",
        "--output",
        output_path.to_str().unwrap(),
    ]);

    assert!(
        out.status.success(),
        "phraya align (case 4 task) should succeed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(
        output_path.exists(),
        ".phraya output file should be created"
    );
}
