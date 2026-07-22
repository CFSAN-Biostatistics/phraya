//! Reference-palette alignment mode (ADR-0011): issues #226 (--reference wiring +
//! hit/miss resolution + composability), #199 (--sealed), and #198 (cross-space sidecar).

use std::path::{Path, PathBuf};
use std::process::Command;
use tempfile::TempDir;

fn manifest() -> PathBuf {
    Path::new(&std::env::var("CARGO_MANIFEST_DIR").unwrap()).join("Cargo.toml")
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

fn run(args: &[&str]) -> std::process::Output {
    let m = manifest();
    let mut full = vec!["run", "--quiet", "--manifest-path", m.to_str().unwrap(), "--"];
    full.extend_from_slice(args);
    Command::new("cargo").args(&full).output().expect("cargo run failed")
}

/// A reference long enough to seed against, plus a read that is an exact substring.
const REF_A: &str = "ACGTACGTGGCCTTAAGGCCTTAAGGCCACGTACGTTTGGCCAATTGGCCAATTACGTACGTGGCCTTAAGGCC";
const READ: &str = "GGCCTTAAGGCCTTAAGGCCACGTACGTTTGGCCAATT";
// A second, distinct reference space (different content -> different hash).
const REF_B: &str = "TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAAACGTACGTTTGGCCAATTGGCCAATTTTTTGGGGCCCCAAAA";

fn plan_with_ref(dir: &Path, plan_name: &str, reads: &Path, reference: &Path) {
    let out = run(&[
        "plan",
        "--inputs",
        reads.to_str().unwrap(),
        "--reference",
        reference.to_str().unwrap(),
        "--output",
        dir.join(plan_name).to_str().unwrap(),
    ]);
    assert!(out.status.success(), "plan failed: {}", String::from_utf8_lossy(&out.stderr));
}

/// #226: a planned reference is a palette **hit** — align reuses it and emits one .phraya.
#[test]
fn reference_hit_emits_per_space_phraya() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();
    let reads = write_fasta(p, "reads.fa", &[("read1", READ)]);
    let refa = write_fasta(p, "refA.fa", &[("refA", REF_A)]);
    plan_with_ref(p, "plan.phrayaplan", &reads, &refa);

    let out_dir = p.join("out");
    let out = run(&[
        "align",
        p.join("plan.phrayaplan").to_str().unwrap(),
        "--reference",
        refa.to_str().unwrap(),
        "--output",
        out_dir.to_str().unwrap(),
    ]);
    assert!(out.status.success(), "align failed: {}", String::from_utf8_lossy(&out.stderr));
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("palette hit"), "expected palette hit, got: {stderr}");

    // One .phraya per space, plus the cross-space sidecar.
    let phrayas: Vec<_> = std::fs::read_dir(&out_dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| e.path().extension().map(|x| x == "phraya").unwrap_or(false))
        .collect();
    assert_eq!(phrayas.len(), 1, "exactly one .phraya per reference space");
    assert!(out_dir.join("cross_space.phraya.queries").exists(), "cross-space sidecar written");
}

/// #199: a reference not in the palette is a **miss**.
/// Tolerant default -> warn + proceed; --sealed -> hard error naming the reference.
#[test]
fn reference_miss_tolerant_vs_sealed() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();
    let reads = write_fasta(p, "reads.fa", &[("read1", READ)]);
    let refa = write_fasta(p, "refA.fa", &[("refA", REF_A)]);
    let refb = write_fasta(p, "refB.fa", &[("refB", REF_B)]);
    // Plan knows only refA; refB will miss.
    plan_with_ref(p, "plan.phrayaplan", &reads, &refa);

    // Tolerant: warns and proceeds.
    let out_dir = p.join("out_tolerant");
    let tolerant = run(&[
        "align",
        p.join("plan.phrayaplan").to_str().unwrap(),
        "--reference",
        refb.to_str().unwrap(),
        "--output",
        out_dir.to_str().unwrap(),
    ]);
    assert!(tolerant.status.success(), "tolerant miss should succeed: {}", String::from_utf8_lossy(&tolerant.stderr));
    let terr = String::from_utf8_lossy(&tolerant.stderr);
    assert!(terr.contains("not in the plan's palette"), "tolerant miss should warn: {terr}");
    assert!(terr.to_lowercase().contains("sketching on the fly"), "tolerant miss should sketch on the fly: {terr}");

    // Sealed: hard error, names the offending reference, and writes nothing.
    let sealed_dir = p.join("out_sealed");
    let sealed = run(&[
        "align",
        p.join("plan.phrayaplan").to_str().unwrap(),
        "--reference",
        refb.to_str().unwrap(),
        "--output",
        sealed_dir.to_str().unwrap(),
        "--sealed",
    ]);
    assert!(!sealed.status.success(), "sealed miss must fail");
    let serr = String::from_utf8_lossy(&sealed.stderr);
    assert!(serr.contains("sealed"), "error should mention sealed: {serr}");
    assert!(serr.contains("refB.fa") || serr.contains("not in the plan's palette"), "error names the offending reference: {serr}");
}

/// #226 composability: align({A,B}) yields the same per-space .phraya files as
/// align({A}) ∪ align({B}). We check byte-identity of each space's .phraya, with the
/// output timestamp pinned so the only run-to-run variation is content.
#[test]
fn multi_reference_is_composable() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();
    let reads = write_fasta(p, "reads.fa", &[("read1", READ)]);
    let refa = write_fasta(p, "refA.fa", &[("refA", REF_A)]);
    let refb = write_fasta(p, "refB.fa", &[("refB", REF_B)]);
    // Plan carrying both spaces in its palette (repeatable --reference on plan).
    let out = run(&[
        "plan",
        "--inputs",
        reads.to_str().unwrap(),
        "--reference",
        refa.to_str().unwrap(),
        "--reference",
        refb.to_str().unwrap(),
        "--output",
        p.join("plan.phrayaplan").to_str().unwrap(),
    ]);
    assert!(out.status.success(), "plan failed: {}", String::from_utf8_lossy(&out.stderr));

    let env_ts = ("PHRAYA_SOURCE_DATE", "2026-01-01T00:00:00Z");

    let run_ts = |args: &[&str]| -> std::process::Output {
        let m = manifest();
        let mut full = vec!["run", "--quiet", "--manifest-path", m.to_str().unwrap(), "--"];
        full.extend_from_slice(args);
        Command::new("cargo").args(&full).env(env_ts.0, env_ts.1).output().unwrap()
    };

    // Combined {A,B}
    let both = p.join("both");
    let o = run_ts(&[
        "align", p.join("plan.phrayaplan").to_str().unwrap(),
        "--reference", refa.to_str().unwrap(),
        "--reference", refb.to_str().unwrap(),
        "--output", both.to_str().unwrap(),
    ]);
    assert!(o.status.success(), "combined align failed: {}", String::from_utf8_lossy(&o.stderr));

    // Separate {A} and {B}
    let only_a = p.join("only_a");
    let o = run_ts(&["align", p.join("plan.phrayaplan").to_str().unwrap(), "--reference", refa.to_str().unwrap(), "--output", only_a.to_str().unwrap()]);
    assert!(o.status.success(), "align A failed: {}", String::from_utf8_lossy(&o.stderr));
    let only_b = p.join("only_b");
    let o = run_ts(&["align", p.join("plan.phrayaplan").to_str().unwrap(), "--reference", refb.to_str().unwrap(), "--output", only_b.to_str().unwrap()]);
    assert!(o.status.success(), "align B failed: {}", String::from_utf8_lossy(&o.stderr));

    // Each space's .phraya must be byte-identical whether produced alone or alongside the other.
    for space_dir in [&only_a, &only_b] {
        for entry in std::fs::read_dir(space_dir).unwrap().filter_map(|e| e.ok()) {
            let path = entry.path();
            if path.extension().map(|x| x == "phraya").unwrap_or(false) {
                let fname = path.file_name().unwrap();
                let combined = both.join(fname);
                assert!(combined.exists(), "combined run missing {:?}", fname);
                assert_eq!(
                    std::fs::read(&path).unwrap(),
                    std::fs::read(&combined).unwrap(),
                    "per-space .phraya {:?} must be identical alone vs combined (composability)",
                    fname
                );
            }
        }
    }
}
