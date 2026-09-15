/// Regression tests for CSP2 spike finding B1: use-case misclassification with no
/// override. Two independent bugs, both in `detect_use_case` (`phraya-cli/src/main.rs`):
///
/// 1. The with-reference branch required *every* non-reference record to be
///    contig-length (`.all(|seq| seq.len() >= 5000)`). A single short contig (a
///    plasmid, an assembly tail) mixed into an otherwise-contig draft assembly flipped
///    the whole comparison to Case 2 (reads-with-ref) instead of Case 4.
/// 2. The no-reference branch treated "more than one input file" as sufficient evidence
///    for Case 3 (contigs + reads, centroid selection) — CSP2's actual invocation shape,
///    two draft assemblies passed as plain positional inputs with neither flagged
///    `--reference`, was always misclassified as Case 3 instead of the direct pairwise
///    Case 4 comparison it should be.
///
/// Also covers the `--use-case` override flag added as the escape hatch for whatever
/// heuristic edge case comes next.
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
/// with neither flagged `--reference`, must classify as Case 4 (ContigsOnly) — direct
/// pairwise comparison — not Case 3 (centroid selection).
#[test]
fn findings_b1_two_assemblies_no_ref_flag_classify_as_contigs_only() {
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
    assert_eq!(
        plan.use_case,
        phraya_io::plan::UseCase::ContigsOnly,
        "two contig-only assemblies with no --reference flag should classify as \
         ContigsOnly (Case 4), got {:?}",
        plan.use_case
    );
}

/// A single short contig (e.g. a small plasmid, below the 5kb contig-length threshold)
/// mixed into an otherwise-contig `--reference`'d input must not flip classification to
/// Case 2 (reads-with-ref). Majority-by-bases, not per-record unanimity.
#[test]
fn findings_b1_one_short_contig_does_not_flip_classification_to_reads() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();

    // 4 contig-length records (8192bp) + 1 short plasmid-like record (1500bp) — the
    // short record is a small fraction of total bases, so the majority-by-bases rule
    // should still classify this as contigs.
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
    assert_eq!(
        plan.use_case,
        phraya_io::plan::UseCase::ContigsOnly,
        "one short contig among four contig-length records should not flip \
         classification away from ContigsOnly, got {:?}",
        plan.use_case
    );
}

/// `--use-case contigs-only` forces Case 4 even on inputs the auto-detection heuristic
/// would otherwise route to Case 2 (short reads + reference).
#[test]
fn findings_b1_use_case_override_forces_contigs_only() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();

    // Short (150bp) records: auto-detection would classify this as ReadsWithRef.
    let reads = write_fasta(p, "reads.fa", &["read1", "read2"], 150);
    let reference = write_fasta(p, "ref.fa", &["ref"], 5120);
    let plan_path = p.join("plan.phrayaplan");

    let out = phraya(&[
        "plan",
        "--inputs",
        reads.to_str().unwrap(),
        "--reference",
        reference.to_str().unwrap(),
        "--use-case",
        "contigs-only",
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
        "--use-case contigs-only should override auto-detection"
    );
}

/// An invalid `--use-case` value is a clear CLI error, not a silent fallback.
#[test]
fn findings_b1_use_case_rejects_invalid_value() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();

    let reads = write_fasta(p, "reads.fa", &["read1"], 150);
    let reference = write_fasta(p, "ref.fa", &["ref"], 5120);
    let plan_path = p.join("plan.phrayaplan");

    let out = phraya(&[
        "plan",
        "--inputs",
        reads.to_str().unwrap(),
        "--reference",
        reference.to_str().unwrap(),
        "--use-case",
        "bogus-value",
        "--output",
        plan_path.to_str().unwrap(),
    ]);
    assert!(
        !out.status.success(),
        "an invalid --use-case value should be rejected, not silently accepted"
    );
}
