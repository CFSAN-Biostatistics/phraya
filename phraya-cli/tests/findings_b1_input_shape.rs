/// Regression tests for CSP2 spike finding B1: input-shape misclassification.
///
/// The input shape — reads vs contigs, presence or absence of a reference — is
/// derived at plan time from what is presented (file contents, record lengths)
/// and is not stored in the plan; the plan records only the task list. These
/// tests verify the correct task list is produced for each shape.
///
/// Two independent bugs were fixed:
/// 1. The with-reference path required *every* non-reference record to be
///    contig-length (`.all(|seq| seq.len() >= 5000)`). A single short contig (a
///    plasmid, an assembly tail) mixed into an otherwise-contig draft assembly
///    flipped the whole comparison to a reads-oriented task list instead of
///    the intended direct input×reference cross-product.
/// 2. The no-reference path treated "more than one input file" as sufficient
///    evidence for centroid selection — two draft assemblies passed as plain
///    positional inputs with neither flagged `--reference`, was always
///    misclassified as centroid selection instead of the direct pairwise
///    comparison it should get.
use std::path::{Path, PathBuf};
use tempfile::TempDir;

fn manifest() -> PathBuf {
    let d = std::env::var("CARGO_MANIFEST_DIR").unwrap();
    Path::new(&d).join("Cargo.toml")
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

/// Write a FASTA file where each record is `seq_len` bp long, named from `ids`.
fn write_fasta(dir: &Path, name: &str, ids: &[&str], seq_len: usize) -> PathBuf {
    let path = dir.join(name);
    let unit = "ACGTACGTACGTACGT".repeat(seq_len / 16 + 1);
    let seq = &unit[..seq_len];
    let mut content = String::new();
    for id in ids {
        content.push_str(&format!(">{id}\n{seq}\n"));
    }
    std::fs::write(&path, content).unwrap();
    path
}

/// Write a FASTA file mixing record lengths (one entry per (id, seq_len) pair).
fn write_fasta_mixed(dir: &Path, name: &str, records: &[(&str, usize)]) -> PathBuf {
    let path = dir.join(name);
    let mut content = String::new();
    for &(id, seq_len) in records {
        let unit = "ACGTACGTACGTACGT".repeat(seq_len / 16 + 1);
        let seq = &unit[..seq_len];
        content.push_str(&format!(">{id}\n{seq}\n"));
    }
    std::fs::write(&path, content).unwrap();
    path
}

/// Two draft assemblies (all contig-length records) passed as plain positional inputs,
/// with neither flagged `--reference`, must produce a direct all-pairs task list —
/// not centroid selection.
///
/// All contigs share the same 8192bp "ACGT" content, so all pairs pass the
/// default 0.1 homology gate: C(5,2) = 10 tasks.
#[test]
fn findings_b1_two_assemblies_no_ref_produces_all_pairs() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();

    let query_asm = write_fasta(p, "query.fa", &["qctg1", "qctg2", "qctg3"], 8192);
    let ref_asm = write_fasta(p, "ref_asm.fa", &["rctg1", "rctg2"], 8192);
    let plan_path = p.join("plan.phrayaplan");

    let out = phraya(&[
        "plan",
        "--inputs",
        query_asm.to_str().unwrap(),
        "--inputs",
        ref_asm.to_str().unwrap(),
        "--output",
        plan_path.to_str().unwrap(),
    ]);
    assert!(
        out.status.success(),
        "phraya plan failed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );

    let plan = phraya_io::plan::read_plan(&plan_path).unwrap();

    // 5 contig-length sequences, no reference → all-pairs task list.
    // Every pair shares the same 8192bp minimizers → Jaccard 1.0 ≥ 0.1 threshold.
    assert_eq!(
        plan.task_list.len(),
        10,
        "two 8192bp contig files should produce 10 all-pairs tasks (C(5,2)), got {:?}",
        plan.task_list
    );

    // No task should reference a reference index — there is no reference.
    let max_idx = plan.sequence_ids.len() as u32;
    assert!(
        plan.task_list.iter().all(|(q, _)| *q < max_idx),
        "tasks should only reference input sequences, not reference indices"
    );
}

/// A single short contig (e.g. a small plasmid, below the 5kb contig-length threshold)
/// mixed into an otherwise-contig `--reference`'d input must not prevent the full
/// input×reference cross-product. Majority-by-bases, not per-record unanimity.
#[test]
fn findings_b1_one_short_contig_does_not_block_reference_tasks() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();

    // 4 contig-length records (8192bp) + 1 short plasmid-like record (1500bp) — the
    // short record is a small fraction of total bases, so the majority-by-bases rule
    // should still treat this as contigs + reference.
    let contigs = write_fasta_mixed(
        p,
        "contigs.fa",
        &[
            ("ctg1", 8192),
            ("ctg2", 8192),
            ("ctg3", 8192),
            ("ctg4", 8192),
            ("plasmid1", 1500),
        ],
    );
    let reference = write_fasta(p, "ref.fa", &["ref"], 8192);
    let plan_path = p.join("plan.phrayaplan");

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

    // 5 input records × 1 reference record = 5 tasks. The short plasmid must
    // not be excluded — every input is tasked against every reference.
    assert_eq!(
        plan.task_list.len(),
        5,
        "5 inputs x 1 reference should produce 5 tasks, got {:?}",
        plan.task_list
    );

    // All tasks target the reference at index 0.
    assert!(
        plan.task_list.iter().all(|(_, t)| *t == 0),
        "all tasks should target the reference at index 0, got {:?}",
        plan.task_list
    );
}