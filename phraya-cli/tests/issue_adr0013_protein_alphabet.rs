//! Protein-space alignment (ADR-0013): auto-detection, --alphabet override, and a full
//! plan+align+filter round trip over amino-acid sequences.

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

// 100 aa, deliberately containing E/F/I/L/P/Q (protein-exclusive letters — no IUPAC
// nucleotide meaning) throughout, so content-based alphabet auto-detection reliably
// classifies it as Protein.
const PROTEOME_SEQ: &str = "MEEPILQFGHKLNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTV";

/// Same sequence with 2 substitutions at known offsets (10 and 40), for a query homolog.
fn mutated_query() -> String {
    let mut bytes = PROTEOME_SEQ.as_bytes().to_vec();
    assert_eq!(bytes[10], b'K');
    bytes[10] = b'W';
    assert_eq!(bytes[40], b'Y');
    bytes[40] = b'A';
    String::from_utf8(bytes).unwrap()
}

/// ADR-0013: `phraya plan` auto-detects Protein alphabet from content (no --alphabet
/// flag needed), and the resulting plan aligns correctly — a full plan+align+filter
/// round trip over amino-acid sequences, with the expected 2 substitutions reported.
#[test]
fn adr0013_protein_alphabet_auto_detected_and_aligns() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();

    let query_seq = mutated_query();
    let proteome = write_fasta(p, "proteome.fasta", &[("target1", PROTEOME_SEQ)]);
    let query = write_fasta(p, "query.fasta", &[("query1", &query_seq)]);
    let plan_path = p.join("plan.phrayaplan");

    let plan_out = run(&[
        "plan",
        "--inputs",
        query.to_str().unwrap(),
        "--reference",
        proteome.to_str().unwrap(),
        "--output",
        plan_path.to_str().unwrap(),
    ]);
    assert!(
        plan_out.status.success(),
        "plan failed: {}",
        String::from_utf8_lossy(&plan_out.stderr)
    );
    let plan_stderr = String::from_utf8_lossy(&plan_out.stderr);
    assert!(
        plan_stderr.contains("Detected alphabet: Protein"),
        "expected auto-detected Protein alphabet in stderr, got: {plan_stderr}"
    );

    let out_path = p.join("out.phraya");
    let align_out = run(&[
        "align",
        plan_path.to_str().unwrap(),
        "query1",
        "target1",
        "--output",
        out_path.to_str().unwrap(),
        "--strategy",
        "sensitive",
    ]);
    assert!(
        align_out.status.success(),
        "align failed: {}",
        String::from_utf8_lossy(&align_out.stderr)
    );
    assert!(out_path.exists(), ".phraya output should be created");

    let filter_out = run(&["filter", out_path.to_str().unwrap(), "--format", "tsv"]);
    assert!(
        filter_out.status.success(),
        "filter failed: {}",
        String::from_utf8_lossy(&filter_out.stderr)
    );
    let tsv = String::from_utf8_lossy(&filter_out.stdout);
    let data_rows: Vec<&str> = tsv.lines().skip(1).filter(|l| !l.is_empty()).collect();
    assert_eq!(
        data_rows.len(),
        2,
        "expected exactly 2 substitution rows (the 2 introduced mutations), got:\n{tsv}"
    );
}

/// ADR-0013: `--alphabet protein` explicit override works (redundant with auto-detection
/// here, but exercises the override path itself, not just the default).
#[test]
fn adr0013_explicit_alphabet_protein_override() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();

    let query_seq = mutated_query();
    let proteome = write_fasta(p, "proteome.fasta", &[("target1", PROTEOME_SEQ)]);
    let query = write_fasta(p, "query.fasta", &[("query1", &query_seq)]);
    let plan_path = p.join("plan.phrayaplan");

    let plan_out = run(&[
        "plan",
        "--inputs",
        query.to_str().unwrap(),
        "--reference",
        proteome.to_str().unwrap(),
        "--output",
        plan_path.to_str().unwrap(),
        "--alphabet",
        "protein",
    ]);
    assert!(
        plan_out.status.success(),
        "plan failed: {}",
        String::from_utf8_lossy(&plan_out.stderr)
    );
    assert!(String::from_utf8_lossy(&plan_out.stderr).contains("Detected alphabet: Protein"));
}

/// ADR-0013: `--alphabet dna` forces DNA interpretation even on protein-looking content —
/// the override always wins over auto-detection.
#[test]
fn adr0013_explicit_alphabet_dna_override_wins_over_detection() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();

    let query_seq = mutated_query();
    let proteome = write_fasta(p, "proteome.fasta", &[("target1", PROTEOME_SEQ)]);
    let query = write_fasta(p, "query.fasta", &[("query1", &query_seq)]);
    let plan_path = p.join("plan.phrayaplan");

    let plan_out = run(&[
        "plan",
        "--inputs",
        query.to_str().unwrap(),
        "--reference",
        proteome.to_str().unwrap(),
        "--output",
        plan_path.to_str().unwrap(),
        "--alphabet",
        "dna",
    ]);
    assert!(
        plan_out.status.success(),
        "plan failed: {}",
        String::from_utf8_lossy(&plan_out.stderr)
    );
    assert!(
        String::from_utf8_lossy(&plan_out.stderr).contains("Detected alphabet: Dna"),
        "explicit --alphabet dna must override content auto-detection"
    );
}

/// `--alphabet` rejects any value outside {auto, dna, protein}.
#[test]
fn adr0013_invalid_alphabet_value_rejected() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();

    let proteome = write_fasta(p, "proteome.fasta", &[("target1", PROTEOME_SEQ)]);
    let query = write_fasta(p, "query.fasta", &[("query1", PROTEOME_SEQ)]);
    let plan_path = p.join("plan.phrayaplan");

    let plan_out = run(&[
        "plan",
        "--inputs",
        query.to_str().unwrap(),
        "--reference",
        proteome.to_str().unwrap(),
        "--output",
        plan_path.to_str().unwrap(),
        "--alphabet",
        "rna",
    ]);
    assert!(
        !plan_out.status.success(),
        "plan should reject an unrecognized --alphabet value"
    );
}
