//! Regression tests for issue #233: `--reference` silently truncated a multi-record
//! FASTA to its first sequence in `phraya plan`, `phraya align --reference` (ADR-0011
//! palette mode), and batch-mode workers. A multi-contig/multi-chromosome reference
//! (draft assembly, chromosome + plasmids — the normal case for CSP2) is N reference
//! spaces, not one (AGENTS.md: reference spaces are content-hashed and mechanical).

use std::path::{Path, PathBuf};
use tempfile::TempDir;

fn manifest() -> PathBuf {
    Path::new(&std::env::var("CARGO_MANIFEST_DIR").unwrap()).join("Cargo.toml")
}

fn phraya(args: &[&str]) -> std::process::Output {
    let m = manifest();
    let mut full = vec!["run", "--quiet", "--manifest-path", m.to_str().unwrap(), "--"];
    full.extend_from_slice(args);
    std::process::Command::new("cargo").args(&full).output().expect("cargo run failed")
}

fn write_fasta(dir: &Path, name: &str, records: &[(&str, &str)]) -> PathBuf {
    let path = dir.join(name);
    let mut s = String::new();
    for (id, seq) in records {
        s.push_str(&format!(">{id}\n{seq}\n"));
    }
    std::fs::write(&path, s).unwrap();
    path
}

// Two distinct "contigs" for a fragmented reference (e.g. chromosome + plasmid), long
// enough to seed against, plus a read that's an exact substring of the second contig.
const CHR: &str = "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCC";
const PLASMID: &str = "ACGTACGTGGCCTTAAGGCCTTAAGGCCACGTACGTTTGGCCAATTGGCCAATTACGTACGTGGCCTTAAGGCC";
const PLASMID_READ: &str = "GGCCTTAAGGCCTTAAGGCCACGTACGTTTGGCCAATT";

/// #233: `phraya plan --reference` with a 2-record FASTA must carry **both** records
/// into the plan's palette (`reference_space`), not just the first.
#[test]
fn plan_reference_space_has_one_entry_per_fasta_record() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();
    let reads = write_fasta(p, "reads.fa", &[("read1", PLASMID_READ)]);
    let reference = write_fasta(p, "ref.fa", &[("chr1", CHR), ("plasmid1", PLASMID)]);
    let plan_path = p.join("plan.phrayaplan");

    let out = phraya(&[
        "plan",
        "--inputs", reads.to_str().unwrap(),
        "--reference", reference.to_str().unwrap(),
        "--output", plan_path.to_str().unwrap(),
    ]);
    assert!(out.status.success(), "plan failed: {}", String::from_utf8_lossy(&out.stderr));

    let plan = phraya_io::plan::read_plan(&plan_path).unwrap();
    assert_eq!(
        plan.reference_space.len(),
        2,
        "a 2-record --reference FASTA must produce 2 reference spaces, got {}",
        plan.reference_space.len()
    );
    // Distinct content -> distinct content hashes; both records preserved, not deduped/collapsed.
    assert_ne!(plan.reference_space[0].content_hash, plan.reference_space[1].content_hash);
}

/// #233: Case 2 (reads + reference) task generation must target every reference record,
/// not just the first. With 2 reference contigs and 1 read, expect 2 tasks (one per
/// reference record), both with the same query (the read) and distinct targets.
#[test]
fn plan_case2_multi_contig_reference_tasks_cover_every_contig() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();
    let reads = write_fasta(p, "reads.fa", &[("read1", PLASMID_READ)]);
    let reference = write_fasta(p, "ref.fa", &[("chr1", CHR), ("plasmid1", PLASMID)]);
    let plan_path = p.join("plan.phrayaplan");

    let out = phraya(&[
        "plan",
        "--inputs", reads.to_str().unwrap(),
        "--reference", reference.to_str().unwrap(),
        "--output", plan_path.to_str().unwrap(),
    ]);
    assert!(out.status.success(), "plan failed: {}", String::from_utf8_lossy(&out.stderr));

    let plan = phraya_io::plan::read_plan(&plan_path).unwrap();
    assert_eq!(
        plan.task_list.len(),
        2,
        "1 read x 2 reference contigs must produce 2 tasks, got {:?}",
        plan.task_list
    );
    let mut targets: Vec<u32> = plan.task_list.iter().map(|(_, t)| *t).collect();
    targets.sort();
    assert_eq!(targets, vec![0, 1], "tasks must target both reference contigs (0 and 1)");
    for (query, _) in &plan.task_list {
        assert_eq!(*query, 2, "the read is sequence index 2 (after the 2 reference records)");
    }
}

/// #233: `phraya align --reference` (ADR-0011 palette mode) with a 2-record FASTA must
/// align against **both** records, writing one `.phraya` per record — not just the first.
#[test]
fn align_reference_mode_writes_one_phraya_per_fasta_record() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();
    let reads = write_fasta(p, "reads.fa", &[("read1", PLASMID_READ)]);
    let reference = write_fasta(p, "ref.fa", &[("chr1", CHR), ("plasmid1", PLASMID)]);
    let plan_path = p.join("plan.phrayaplan");

    let out = phraya(&[
        "plan",
        "--inputs", reads.to_str().unwrap(),
        "--output", plan_path.to_str().unwrap(),
    ]);
    assert!(out.status.success(), "plan failed: {}", String::from_utf8_lossy(&out.stderr));

    let out_dir = p.join("out");
    let out = phraya(&[
        "align",
        plan_path.to_str().unwrap(),
        "--reference", reference.to_str().unwrap(),
        "--output", out_dir.to_str().unwrap(),
    ]);
    assert!(out.status.success(), "align --reference failed: {}", String::from_utf8_lossy(&out.stderr));

    let phraya_files: Vec<_> = std::fs::read_dir(&out_dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| e.path().extension().map(|x| x == "phraya").unwrap_or(false))
        .collect();
    assert_eq!(
        phraya_files.len(),
        2,
        "a 2-record --reference FASTA must write 2 .phraya files (one per space), got {}",
        phraya_files.len()
    );
}

/// #233: batch mode (single-target-per-worker output) must hard-error on a multi-record
/// reference/centroid file instead of silently aligning against only the first record.
#[test]
fn batch_mode_errors_loudly_on_multi_contig_reference() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();
    let reads = write_fasta(p, "reads.fa", &[("read1", PLASMID_READ), ("read2", CHR)]);
    let reference = write_fasta(p, "ref.fa", &[("chr1", CHR), ("plasmid1", PLASMID)]);
    let plan_path = p.join("plan.phrayaplan");
    let pattern = p.join("out_{worker}.phraya").to_str().unwrap().to_string();

    let out = phraya(&[
        "plan",
        "--inputs", reads.to_str().unwrap(),
        "--reference", reference.to_str().unwrap(),
        "--output", plan_path.to_str().unwrap(),
        "--batch-to", "1",
        "--batch-output-pattern", &pattern,
    ]);
    assert!(out.status.success(), "plan failed: {}", String::from_utf8_lossy(&out.stderr));

    let out = phraya(&["align", plan_path.to_str().unwrap(), "--worker", "0"]);
    assert!(
        !out.status.success(),
        "batch-mode align against a multi-contig reference must fail loudly, not silently truncate"
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("multi-record") || stderr.contains("--reference"),
        "error should explain the multi-record limitation and point at --reference mode: {stderr}"
    );
}
