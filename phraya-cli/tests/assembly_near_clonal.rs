//! End-to-end smoke for `--strategy assembly` on near-clonal synthetic pairs,
//! exercising §5 (identical contigs retained) + §6 (no false variants).
//!
//! Builds on the `palette_fallback_noise_suppression.rs` harness shape: a small
//! multi-contig reference plus near-clonal query contigs, run through the full
//! `plan -> align -> filter` CLI and asserted on the `.phraya` observations and
//! the snpdiffs header.

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

fn dna(len: usize, seed: u64) -> String {
    let mut x = seed;
    (0..len)
        .map(|_| {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            b"ACGT"[(x >> 33 & 3) as usize] as char
        })
        .collect()
}

fn flip(c: char) -> char {
    match c { 'A' => 'C', 'C' => 'G', 'G' => 'T', _ => 'A' }
}

/// Near-clonal pair with mixed contig fates:
///  - `ref1` (400 bp): query copy has one planted SNP -> exercises §6 (no false variants at 400 bp).
///  - `ref2` (1200 bp): query copy is byte-identical, header renamed -> exercises §5 (must be
///    retained and aligned, not dropped by content-hash exclusion).
#[test]
fn assembly_near_clonal_mixed_contigs() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();
    let ref1 = dna(400, 11);
    let ref2 = dna(1200, 99);
    let ref_fasta = write_fasta(p, "ref.fasta", &[("ref1", ref1.as_str()), ("ref2", ref2.as_str())]);

    let q1 = {
        let mut b: Vec<char> = ref1.chars().collect();
        b[200] = flip(b[200]);
        b.into_iter().collect::<String>()
    };
    let q2 = ref2.clone(); // identical to ref2, renamed header
    let query_fasta = write_fasta(p, "query.fasta", &[("q1", q1.as_str()), ("q2", q2.as_str())]);

    let m = manifest();
    let mut full = |args: &[&str]| -> std::process::Output {
        let mut a = vec!["run", "--quiet", "--bin", "phraya", "--manifest-path", m.to_str().unwrap(), "--"];
        a.extend_from_slice(args);
        Command::new("cargo").args(&a).output().expect("cargo run failed")
    };

    let plan_out = full(&["plan", "--inputs", query_fasta.to_str().unwrap(),
        "--reference", ref_fasta.to_str().unwrap(),
        "--output", p.join("plan.phrayaplan").to_str().unwrap()]);
    assert!(plan_out.status.success(), "plan failed: {}", String::from_utf8_lossy(&plan_out.stderr));

    let out_dir = p.join("out");
    let align_out = full(&["align", "--reference", ref_fasta.to_str().unwrap(),
        p.join("plan.phrayaplan").to_str().unwrap(),
        "--output", out_dir.to_str().unwrap(), "--strategy", "assembly"]);
    assert!(align_out.status.success(), "align failed: {}", String::from_utf8_lossy(&align_out.stderr));

    // §5: q2 (byte-identical to ref2, renamed) MUST be placed and yield zero variants.
    let f2 = phraya_io::phraya::read_phraya(&out_dir.join("ref2.phraya")).unwrap();
    assert!(
        !f2.alignments.is_empty(),
        "byte-identical renamed contig q2 was dropped from the read pool (§5 not fixed)"
    );
    assert_eq!(f2.observations.len(), 0, "identical contig q2 must yield 0 variants, got {}", f2.observations.len());

    // §6: q1 vs ref1 carries exactly one planted SNP; no false variants.
    let f1 = phraya_io::phraya::read_phraya(&out_dir.join("ref1.phraya")).unwrap();
    assert_eq!(
        f1.observations.len(),
        1,
        "q1 vs ref1 should have exactly the one planted SNP, got {} observations",
        f1.observations.len()
    );
    let obs = &f1.observations[0];
    assert_eq!(obs.position(), 200, "planted SNP should be at position 200");

    // End-to-end snpdiffs across both contigs: exactly 1 SNP, 0 indels.
    let snpdiffs_out = p.join("merged.snpdiffs");
    let filter_out = full(&[
        "filter", out_dir.join("ref1.phraya").to_str().unwrap(), out_dir.join("ref2.phraya").to_str().unwrap(),
        "--format", "snpdiffs",
        "--reference-fasta", ref_fasta.to_str().unwrap(),
        "--query-fasta", query_fasta.to_str().unwrap(),
        "--reference-id", "REF", "--query-id", "QRY",
        "--output", snpdiffs_out.to_str().unwrap(),
    ]);
    assert!(filter_out.status.success(), "filter failed: {}", String::from_utf8_lossy(&filter_out.stderr));
    let out = std::fs::read_to_string(&snpdiffs_out).unwrap();
    let header = out.lines().find(|l| l.starts_with("#\t")).unwrap();
    let snps = header.split('\t').find(|k| k.starts_with("SNPs:")).unwrap();
    assert_eq!(snps, "SNPs:1", "end-to-end snpdiffs should report exactly 1 SNP; header: {header}");
}
