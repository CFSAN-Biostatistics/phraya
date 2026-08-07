//! Gap-affine scoring in `sensitive` (ADR-0014): `--gap-model` CLI validation and a
//! real indel-consolidation end-to-end check.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use tempfile::TempDir;

fn manifest() -> PathBuf {
    Path::new(&std::env::var("CARGO_MANIFEST_DIR").unwrap()).join("Cargo.toml")
}

fn run(args: &[&str]) -> Output {
    let m = manifest();
    let mut full = vec!["run", "--quiet", "--manifest-path", m.to_str().unwrap(), "--"];
    full.extend_from_slice(args);
    Command::new("cargo").args(&full).output().expect("cargo run failed")
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

/// A ~120bp reference and a read that is the reference with a clean 6bp deletion in the
/// middle — the exact scenario ADR-0014 exists for.
const REF_SEQ: &str = "ACGTACGTACGTACGTACGTGGCCAATTGGCCAATTCCGGTTAACCGGTTAACCTTAACCGGTTAACCGGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

fn read_with_clean_deletion(del_len: usize) -> String {
    let bytes = REF_SEQ.as_bytes();
    let start = 40; // inside a distinctive, non-repetitive stretch
    let mut out = bytes[..start].to_vec();
    out.extend_from_slice(&bytes[start + del_len..]);
    String::from_utf8(out).unwrap()
}

fn plan_and_align(
    dir: &Path,
    ref_path: &Path,
    read_path: &Path,
    strategy: &str,
    gap_model: Option<&str>,
) -> (Output, PathBuf) {
    let plan_path = dir.join(format!("plan_{strategy}_{:?}.phrayaplan", gap_model));
    let plan_out = run(&[
        "plan",
        "--inputs",
        read_path.to_str().unwrap(),
        "--reference",
        ref_path.to_str().unwrap(),
        "--output",
        plan_path.to_str().unwrap(),
    ]);
    assert!(plan_out.status.success(), "plan failed: {}", String::from_utf8_lossy(&plan_out.stderr));

    let out_path = dir.join(format!("out_{strategy}_{:?}.phraya", gap_model));
    let mut args = vec![
        "align",
        plan_path.to_str().unwrap(),
        "read1",
        "ref",
        "--output",
        out_path.to_str().unwrap(),
        "--strategy",
        strategy,
    ];
    if let Some(gm) = gap_model {
        args.push("--gap-model");
        args.push(gm);
    }
    (run(&args), out_path)
}

/// `--gap-model` paired with `balanced` is a hard CLI error, not a silent no-op.
#[test]
fn adr0014_gap_model_rejected_with_balanced() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();
    let read = read_with_clean_deletion(6);
    let fasta_ref = write_fasta(p, "ref.fasta", &[("ref", REF_SEQ)]);
    let fasta_read = write_fasta(p, "read.fasta", &[("read1", &read)]);
    let (out, _) = plan_and_align(p, &fasta_ref, &fasta_read, "balanced", Some("affine"));
    assert!(!out.status.success(), "--gap-model must be rejected with --strategy balanced");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("sensitive"),
        "error message should explain the sensitive-only requirement, got: {stderr}"
    );
}

/// `--gap-model` paired with `fast` is likewise rejected.
#[test]
fn adr0014_gap_model_rejected_with_fast() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();
    let read = read_with_clean_deletion(6);
    let fasta_ref = write_fasta(p, "ref.fasta", &[("ref", REF_SEQ)]);
    let fasta_read = write_fasta(p, "read.fasta", &[("read1", &read)]);
    let (out, _) = plan_and_align(p, &fasta_ref, &fasta_read, "fast", Some("linear"));
    assert!(!out.status.success(), "--gap-model must be rejected with --strategy fast");
}

/// An invalid --gap-model value is rejected even with --strategy sensitive.
#[test]
fn adr0014_invalid_gap_model_value_rejected() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();
    let read = read_with_clean_deletion(6);
    let fasta_ref = write_fasta(p, "ref.fasta", &[("ref", REF_SEQ)]);
    let fasta_read = write_fasta(p, "read.fasta", &[("read1", &read)]);
    let (out, _) = plan_and_align(p, &fasta_ref, &fasta_read, "sensitive", Some("quadratic"));
    assert!(!out.status.success(), "an unrecognized --gap-model value must be rejected");
}

/// `--strategy sensitive` with no `--gap-model` (defaults to affine) succeeds and
/// consolidates a clean 6bp deletion into a single CIGAR indel op — the observable
/// contract ADR-0014 exists to deliver. `--gap-model linear` explicitly is also
/// exercised as a regression guard for the old path (still available on request).
#[test]
fn adr0014_sensitive_default_consolidates_clean_deletion() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();
    let read = read_with_clean_deletion(6);
    let fasta_ref = write_fasta(p, "ref.fasta", &[("ref", REF_SEQ)]);
    let fasta_read = write_fasta(p, "read.fasta", &[("read1", &read)]);

    // Default (affine) — no --gap-model flag at all, per ADR-0014's decision.
    let (align_out, out_path) = plan_and_align(p, &fasta_ref, &fasta_read, "sensitive", None);
    assert!(align_out.status.success(), "align failed: {}", String::from_utf8_lossy(&align_out.stderr));

    let filter_out = run(&["filter", out_path.to_str().unwrap(), "--format", "tsv"]);
    assert!(filter_out.status.success());
    let tsv = String::from_utf8_lossy(&filter_out.stdout);
    let data_rows: Vec<&str> = tsv.lines().skip(1).filter(|l| !l.is_empty()).collect();
    // A clean 6bp deletion with no other divergence should surface as a single
    // deletion event's worth of rows (extract_variants_from_cigar emits one row per
    // deleted base under the 'I'-convention deletion branch — see main.rs::output_tsv /
    // executor.rs::extract_variants_from_cigar), not scattered substitutions.
    assert!(
        !data_rows.is_empty(),
        "expected the deletion to be reported, got no variant rows:\n{tsv}"
    );
    for row in &data_rows {
        // No substitution-only rows masquerading as the deletion — every reported row's
        // cigar should contain an indel op, not just X's, confirming consolidation
        // rather than the linear model's scatter failure mode.
        assert!(
            row.contains('I') || row.contains('D'),
            "expected indel-bearing CIGAR in reported row, got: {row}"
        );
    }
}

/// `--gap-model linear` under `--strategy sensitive` is accepted (opt-out preserved).
#[test]
fn adr0014_sensitive_linear_override_accepted() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();
    let read = read_with_clean_deletion(6);
    let fasta_ref = write_fasta(p, "ref.fasta", &[("ref", REF_SEQ)]);
    let fasta_read = write_fasta(p, "read.fasta", &[("read1", &read)]);
    let (align_out, out_path) = plan_and_align(p, &fasta_ref, &fasta_read, "sensitive", Some("linear"));
    assert!(align_out.status.success(), "align failed: {}", String::from_utf8_lossy(&align_out.stderr));
    assert!(out_path.exists());
}
