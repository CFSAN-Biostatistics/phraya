//! Integration test for the CSP2 snpdiffs escalation.
//!
//! Reproduces the scenario from `phraya-csp2-escalation.md`: a multi-contig
//! bacterial assembly (two ~400 bp random contigs) aligned against another
//! assembly whose contigs are non-homologous to each other.  Without the
//! `(0,0)` fallback noise gate each query contig would deposit ~200 fabricated
//! `VariantObservation`s into every non-homologous reference space; with the
//! gate only the two real SNPs (planted at fixed offsets) survive.
//!
//! Also exercises `phraya filter --format snpdiffs` with multiple input
//! `.phraya` files (multi-input merge).

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


/// Pseudo-random DNA via a 64-bit LCG (same algorithm as the aligner internals).
fn dna(len: usize, seed: u64) -> String {
    let mut x = seed;
    (0..len)
        .map(|_| {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"ACGT"[(x >> 33 & 3) as usize] as char
        })
        .collect()
}

/// Flip a single base to a different one (A↔C↔G↔T cycle).
fn flip(c: char) -> char {
    match c {
        'A' => 'C',
        'C' => 'G',
        'G' => 'T',
        _ => 'A',
    }
}

/// Embed `n_snps` substitutions into `src` at positions `positions`.
fn with_snps(src: &str, positions: &[usize]) -> String {
    let bytes: Vec<char> = src.chars().collect();
    let mut out = bytes.clone();
    for &pos in positions {
        out[pos] = flip(out[pos]);
    }
    out.into_iter().collect()
}

/// Contig length chosen so that `fallback_anchor_applies(len, len)` is true
/// (target ≤ 10× query) while staying within Myers range (≤ 500 bp) to avoid
/// WFA cap issues on CI machines.
const CONTIG_LEN: usize = 400;

/// Two planted SNPs per homologous contig pair.
const SNP_POSITIONS: &[usize] = &[50, 150];

#[test]
fn palette_comparable_contigs_no_fallback_noise() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();

    // ── Reference: two random contigs that share no minimizer ──
    let ref_a = dna(CONTIG_LEN, 11);
    let ref_b = dna(CONTIG_LEN, 99);
    let ref_fasta = write_fasta(p, "ref.fasta", &[("ref1", ref_a.as_str()), ("ref2", ref_b.as_str())]);

    // ── Query: each contig is a copy of one reference + 2 SNPs ──
    let query1 = with_snps(&ref_a, SNP_POSITIONS); // homologous to ref1
    let query2 = with_snps(&ref_b, SNP_POSITIONS); // homologous to ref2
    let query_fasta = write_fasta(p, "query.fasta", &[("query1", query1.as_str()), ("query2", query2.as_str())]);

    // ── Plan ──
    let plan_out = run(&[
        "plan",
        "--inputs",
        query_fasta.to_str().unwrap(),
        "--reference",
        ref_fasta.to_str().unwrap(),
        "--output",
        p.join("plan.phrayaplan").to_str().unwrap(),
    ]);
    assert!(plan_out.status.success(), "plan failed: {}", String::from_utf8_lossy(&plan_out.stderr));

    // ── Align (palette mode) ──
    let out_dir = p.join("out");
    let align_out = run(&[
        "align",
        "--reference",
        ref_fasta.to_str().unwrap(),
        p.join("plan.phrayaplan").to_str().unwrap(),
        "--output",
        out_dir.to_str().unwrap(),
        "--strategy",
        "balanced",
    ]);
    assert!(align_out.status.success(), "align failed: {}", String::from_utf8_lossy(&align_out.stderr));

    let file1 = phraya_io::phraya::read_phraya(&out_dir.join("ref1.phraya")).unwrap();
    let file2 = phraya_io::phraya::read_phraya(&out_dir.join("ref2.phraya")).unwrap();

    // ── Assert: noise suppressed, real SNPs retained ──
    // ref1.phraya: Q1 aligns here (homologous, 2 real SNPs).
    //              Q2 aligns here (non-homologous, (0,0) fallback at ~25% identity → suppressed).
    assert_eq!(
        file1.observations.len(),
        2,
        "ref1.phraya should have exactly 2 variant observations (real SNPs from query1), \
         got {} — fallback noise not suppressed",
        file1.observations.len()
    );
    assert_eq!(
        file2.observations.len(),
        2,
        "ref2.phraya should have exactly 2 variant observations (real SNPs from query2), \
         got {} — fallback noise not suppressed",
        file2.observations.len()
    );

    // Verify the SNPs are at the expected positions.
    let positions1: Vec<u32> = file1.observations.iter().map(|o| o.position()).collect();
    let positions2: Vec<u32> = file2.observations.iter().map(|o| o.position()).collect();
    assert!(positions1.contains(&(SNP_POSITIONS[0] as u32)), "ref1 should have SNP at position {}", SNP_POSITIONS[0]);
    assert!(positions1.contains(&(SNP_POSITIONS[1] as u32)), "ref1 should have SNP at position {}", SNP_POSITIONS[1]);
    assert!(positions2.contains(&(SNP_POSITIONS[0] as u32)), "ref2 should have SNP at position {}", SNP_POSITIONS[0]);
    assert!(positions2.contains(&(SNP_POSITIONS[1] as u32)), "ref2 should have SNP at position {}", SNP_POSITIONS[1]);
}

#[test]
fn multi_input_snpdiffs_merges_contig_files() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();

    // ── Same setup as the noise suppression test ──
    let ref_a = dna(CONTIG_LEN, 11);
    let ref_b = dna(CONTIG_LEN, 99);
    let ref_fasta = write_fasta(p, "ref.fasta", &[("ref1", ref_a.as_str()), ("ref2", ref_b.as_str())]);

    let query1 = with_snps(&ref_a, SNP_POSITIONS);
    let query2 = with_snps(&ref_b, SNP_POSITIONS);
    let query_fasta = write_fasta(p, "query.fasta", &[("query1", query1.as_str()), ("query2", query2.as_str())]);

    // Plan + align
    let plan_out = run(&[
        "plan",
        "--inputs", query_fasta.to_str().unwrap(),
        "--reference", ref_fasta.to_str().unwrap(),
        "--output", p.join("plan.phrayaplan").to_str().unwrap(),
    ]);
    assert!(plan_out.status.success(), "plan failed: {}", String::from_utf8_lossy(&plan_out.stderr));

    let out_dir = p.join("out");
    let align_out = run(&[
        "align",
        "--reference", ref_fasta.to_str().unwrap(),
        p.join("plan.phrayaplan").to_str().unwrap(),
        "--output", out_dir.to_str().unwrap(),
    ]);
    assert!(align_out.status.success(), "align failed: {}", String::from_utf8_lossy(&align_out.stderr));

    // ── Multi-input snpdiffs filter ──
    let snpdiffs_out = p.join("merged.snpdiffs");
    let filter_out = run(&[
        "filter",
        out_dir.join("ref1.phraya").to_str().unwrap(),
        out_dir.join("ref2.phraya").to_str().unwrap(),
        "--format", "snpdiffs",
        "--reference-fasta", ref_fasta.to_str().unwrap(),
        "--query-fasta", query_fasta.to_str().unwrap(),
        "--reference-id", "assembly_ref",
        "--query-id", "assembly_query",
        "--output", snpdiffs_out.to_str().unwrap(),
    ]);
    assert!(filter_out.status.success(), "filter failed: {}", String::from_utf8_lossy(&filter_out.stderr));

    let output = std::fs::read_to_string(&snpdiffs_out).unwrap();
    let lines: Vec<&str> = output.lines().collect();

    // ── Header ──
    let header: Vec<&str> = lines[0].split('\t').collect();
    assert_eq!(header[0], "#");
    assert!(header.iter().any(|f| f == &"Reference_ID:assembly_ref"));
    assert!(header.iter().any(|f| f == &"Query_ID:assembly_query"));
    assert!(header.iter().any(|f| f == &"SNPs:4"), "header SNPs should be 4 (2 per contig), got: {:?}",
        header.iter().find(|f| f.starts_with("SNPs")));

    // ── Two BED rows (one per reference contig) ──
    let bed1: Vec<&str> = lines[1].split('\t').collect();
    assert_eq!(bed1[0], "##");
    assert_eq!(bed1[1], "ref1");

    let bed2: Vec<&str> = lines[2].split('\t').collect();
    assert_eq!(bed2[0], "##");
    assert_eq!(bed2[1], "ref2");

    // ── Variant rows: 2 per contig, with correct Ref_Contig ──
    let var_lines: Vec<&str> = lines[3..]
        .iter()
        .filter(|l| !l.starts_with("#\t") && !l.starts_with("##"))
        .copied()
        .collect();
    assert_eq!(var_lines.len(), 4, "should have 4 variant rows (2 per contig)");

    for line in &var_lines {
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), 21, "variant row must have 21 columns");
        assert_eq!(fields[20], "SNP", "all variants should be SNPs");
    }

    // Issue 4 regression: Ref_Contig must be the FASTA header ID (string), not a
    // numeric index (e.g. "0", "1"). If numeric IDs regressed, these would fail.
    assert_eq!(var_lines[0].split('\t').next().unwrap(), "ref1");
    assert_eq!(var_lines[1].split('\t').next().unwrap(), "ref1");
    assert_eq!(var_lines[2].split('\t').next().unwrap(), "ref2");
    assert_eq!(var_lines[3].split('\t').next().unwrap(), "ref2");
}

