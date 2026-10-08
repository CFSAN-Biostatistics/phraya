//! Assembly-vs-assembly regression tests for `--strategy assembly` (ADR-0015; CSP2
//! handoff `phraya-csp2-snp-discrepancy.md`, failures A and B).
//!
//! Truth is exact: a random genome, and a query that is the genome with one SNP every
//! `SNP_PERIOD` bp. Both are cut into contigs. SNPs within `EDGE` bp of any contig end in
//! either assembly are excluded from the recall denominator ("callable" truth): a local
//! aligner cannot place a mismatch at an alignment terminus, and dnadiff does not either.
//!
//! Fixture A keeps identical breakpoints but includes a 70 kb contig, which a whole-contig
//! WFA abandons at 0.5% divergence (failure A1). Fixture B re-cuts the query at different
//! offsets so query contigs span several reference contigs, reverse-complements one of
//! them, and adds an unrelated decoy contig to each assembly (failures A2, A3 and B).

use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::LazyLock;
use tempfile::TempDir;

const GENOME_LEN: usize = 150_000;
const SNP_PERIOD: usize = 200;
const SNP_OFFSET: usize = 101;
const EDGE: usize = 50;
const DECOY_LEN: usize = 20_000;

fn lcg_dna(len: usize, seed: u64) -> Vec<u8> {
    let mut x = seed;
    (0..len)
        .map(|_| {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            b"ACGT"[((x >> 33) & 3) as usize]
        })
        .collect()
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        other => other,
    }
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&b| complement(b)).collect()
}

fn alt_base(b: u8) -> u8 {
    match b {
        b'A' => b'C',
        b'C' => b'G',
        b'G' => b'T',
        _ => b'A',
    }
}

/// A contig cut from the genome interval `[start, end)`, optionally reverse-complemented.
struct Cut {
    id: String,
    start: usize,
    end: usize,
    reverse: bool,
}

fn cuts(prefix: &str, points: &[usize], reversed: &[usize]) -> Vec<Cut> {
    points
        .windows(2)
        .enumerate()
        .map(|(i, w)| Cut {
            id: format!("{prefix}{i}"),
            start: w[0],
            end: w[1],
            reverse: reversed.contains(&i),
        })
        .collect()
}

struct Truth {
    alt: u8,
    callable: bool,
}

struct Fixture {
    /// (id, sequence) records exactly as written to the FASTA files.
    ref_records: Vec<(String, Vec<u8>)>,
    qry_records: Vec<(String, Vec<u8>)>,
    /// Keyed by (reference contig, 1-based reference position).
    truth: HashMap<(String, u32), Truth>,
    decoy_ref: Option<String>,
    decoy_qry: Option<String>,
}

fn build_fixture(seed: u64, ref_points: &[usize], qry_points: &[usize], reversed: &[usize], decoys: bool) -> Fixture {
    let genome = lcg_dna(GENOME_LEN, seed);
    let mut mutated = genome.clone();
    let snps: Vec<usize> = (SNP_OFFSET..GENOME_LEN).step_by(SNP_PERIOD).collect();
    for &g in &snps {
        mutated[g] = alt_base(genome[g]);
    }

    let ref_cuts = cuts("r", ref_points, &[]);
    let qry_cuts = cuts("q", qry_points, reversed);
    let containing = |set: &[Cut], g: usize| -> usize {
        set.iter().position(|c| c.start <= g && g < c.end).expect("genome position outside every contig")
    };
    let off_edge = |c: &Cut, g: usize| g - c.start >= EDGE && c.end - 1 - g >= EDGE;

    let mut truth = HashMap::new();
    for &g in &snps {
        let r = &ref_cuts[containing(&ref_cuts, g)];
        let q = &qry_cuts[containing(&qry_cuts, g)];
        truth.insert(
            (r.id.clone(), (g - r.start + 1) as u32),
            Truth { alt: mutated[g], callable: off_edge(r, g) && off_edge(q, g) },
        );
    }

    let mut ref_records: Vec<(String, Vec<u8>)> =
        ref_cuts.iter().map(|c| (c.id.clone(), genome[c.start..c.end].to_vec())).collect();
    let mut qry_records: Vec<(String, Vec<u8>)> = qry_cuts
        .iter()
        .map(|c| {
            let seq = &mutated[c.start..c.end];
            (c.id.clone(), if c.reverse { revcomp(seq) } else { seq.to_vec() })
        })
        .collect();

    let (mut decoy_ref, mut decoy_qry) = (None, None);
    if decoys {
        ref_records.push(("decoy_ref".to_string(), lcg_dna(DECOY_LEN, seed ^ 0xDEC0_1)));
        qry_records.push(("decoy_qry".to_string(), lcg_dna(DECOY_LEN, seed ^ 0xDEC0_2)));
        decoy_ref = Some("decoy_ref".to_string());
        decoy_qry = Some("decoy_qry".to_string());
    }
    Fixture { ref_records, qry_records, truth, decoy_ref, decoy_qry }
}

fn write_fasta(path: &Path, records: &[(String, Vec<u8>)]) {
    let mut s = String::new();
    for (id, seq) in records {
        s.push('>');
        s.push_str(id);
        s.push('\n');
        s.push_str(std::str::from_utf8(seq).unwrap());
        s.push('\n');
    }
    std::fs::write(path, s).unwrap();
}

fn phraya(args: &[&str]) {
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("Cargo.toml");
    let out = Command::new("cargo")
        .args(["run", "--quiet", "--manifest-path", manifest.to_str().unwrap(), "--"])
        .args(args)
        .output()
        .expect("failed to spawn cargo run");
    assert!(
        out.status.success(),
        "phraya {} failed:\nstdout:\n{}\nstderr:\n{}",
        args.join(" "),
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );
}

struct Snpdiffs {
    header: HashMap<String, String>,
    /// BED rows without the `##` prefix: 11 columns.
    bed: Vec<Vec<String>>,
    /// Variant rows: 21 columns.
    variants: Vec<Vec<String>>,
}

impl Snpdiffs {
    fn header_f64(&self, key: &str) -> f64 {
        let raw = self.header.get(key).unwrap_or_else(|| panic!("header is missing {key}"));
        raw.parse().unwrap_or_else(|_| panic!("header {key} is not numeric: {raw}"))
    }

    fn rows_of(&self, cat: &str) -> Vec<&Vec<String>> {
        self.variants.iter().filter(|r| r[20] == cat).collect()
    }

    /// BED rows describing an alignment block (not an uncovered `.` region).
    fn covered_bed(&self) -> Vec<&Vec<String>> {
        self.bed.iter().filter(|r| r[1] != "." && r[6] != ".").collect()
    }
}

fn parse_snpdiffs(text: &str) -> Snpdiffs {
    let mut header = HashMap::new();
    let mut bed = Vec::new();
    let mut variants = Vec::new();
    for line in text.lines().filter(|l| !l.is_empty()) {
        if let Some(rest) = line.strip_prefix("##\t") {
            let row: Vec<String> = rest.split('\t').map(str::to_string).collect();
            assert_eq!(row.len(), 11, "BED row must have 11 columns: {line}");
            bed.push(row);
        } else if let Some(rest) = line.strip_prefix("#\t") {
            for field in rest.split('\t') {
                let (k, v) = field.split_once(':').unwrap_or_else(|| panic!("header field without ':' : {field}"));
                header.insert(k.to_string(), v.to_string());
            }
        } else {
            let row: Vec<String> = line.split('\t').map(str::to_string).collect();
            assert_eq!(row.len(), 21, "variant row must have 21 columns: {line}");
            variants.push(row);
        }
    }
    Snpdiffs { header, bed, variants }
}

fn run_pipeline(fx: &Fixture) -> Snpdiffs {
    let dir = TempDir::new().unwrap();
    let p = dir.path();
    let (ref_fa, qry_fa) = (p.join("ref.fasta"), p.join("qry.fasta"));
    write_fasta(&ref_fa, &fx.ref_records);
    write_fasta(&qry_fa, &fx.qry_records);
    let (plan, out_dir, snpdiffs) = (p.join("plan.phrayaplan"), p.join("aln"), p.join("out.snpdiffs"));
    let s = |path: &Path| path.to_str().unwrap().to_string();

    phraya(&["plan", "--reference", &s(&ref_fa), "--inputs", &s(&qry_fa), "--output", &s(&plan)]);
    phraya(&["align", "--reference", &s(&ref_fa), &s(&plan), "--output", &s(&out_dir), "--strategy", "assembly"]);

    let mut files: Vec<String> = std::fs::read_dir(&out_dir)
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|path| path.extension().is_some_and(|ext| ext == "phraya"))
        .map(|path| s(&path))
        .collect();
    files.sort();
    assert!(!files.is_empty(), "align wrote no .phraya files");
    let mut args = vec!["filter".to_string()];
    args.extend(files);
    for flag in [
        ["--format", "snpdiffs"],
        ["--reference-fasta", &s(&ref_fa)],
        ["--query-fasta", &s(&qry_fa)],
        ["--reference-id", "REF"],
        ["--query-id", "QRY"],
        ["--output", &s(&snpdiffs)],
    ] {
        args.extend(flag.iter().map(|a| a.to_string()));
    }
    phraya(&args.iter().map(String::as_str).collect::<Vec<_>>());
    parse_snpdiffs(&std::fs::read_to_string(&snpdiffs).unwrap())
}

/// Identical breakpoints, including a 70 kb contig (> the ~60 kb whole-contig WFA limit at
/// 0.5% divergence).
static FIXTURE_A: LazyLock<(Fixture, Snpdiffs)> = LazyLock::new(|| {
    let points = [0, 70_013, 110_029, GENOME_LEN];
    let fx = build_fixture(0x5EED_A, &points, &points, &[], false);
    let out = run_pipeline(&fx);
    (fx, out)
});

/// Shifted breakpoints (query contigs span 2-3 reference contigs), query contig `q1`
/// reverse-complemented, plus one unrelated decoy contig in each assembly.
static FIXTURE_B: LazyLock<(Fixture, Snpdiffs)> = LazyLock::new(|| {
    let ref_points = [0, 30_011, 61_007, 95_003, GENOME_LEN];
    let qry_points = [0, 18_017, 52_021, 99_013, GENOME_LEN];
    let fx = build_fixture(0x5EED_B, &ref_points, &qry_points, &[1], true);
    let out = run_pipeline(&fx);
    (fx, out)
});

/// Every planted SNP that is callable is reported with the planted allele, and nothing
/// else is reported.
fn assert_variant_truth(fx: &Fixture, out: &Snpdiffs) {
    let indels = out.rows_of("Indel");
    assert!(indels.is_empty(), "no indels were planted, got {} Indel rows, e.g. {:?}", indels.len(), indels.first());

    let snps = out.rows_of("SNP");
    let found: HashMap<(String, u32), &Vec<String>> =
        snps.iter().map(|r| ((r[0].clone(), r[2].parse().unwrap()), *r)).collect();

    let false_positives: Vec<_> = found.keys().filter(|k| !fx.truth.contains_key(*k)).take(10).collect();
    assert!(false_positives.is_empty(), "SNP rows at unplanted positions: {false_positives:?}");

    let wrong_allele: Vec<_> = found
        .iter()
        .filter(|(k, r)| r[13].as_bytes() != [fx.truth[*k].alt])
        .map(|(k, r)| (k.clone(), r[13].clone()))
        .take(10)
        .collect();
    assert!(wrong_allele.is_empty(), "SNP rows with the wrong Query_Base: {wrong_allele:?}");

    let callable: Vec<_> = fx.truth.iter().filter(|(_, t)| t.callable).map(|(k, _)| k).collect();
    let recalled = callable.iter().filter(|k| found.contains_key(**k)).count();
    let recall = recalled as f64 / callable.len() as f64;
    let mut missed: Vec<_> = callable.iter().filter(|k| !found.contains_key(**k)).take(10).collect();
    missed.sort();
    assert!(
        recall >= 0.99,
        "callable SNP recall {:.3} ({recalled}/{}) below 0.99; first missed: {missed:?}",
        recall,
        callable.len()
    );
}

#[test]
fn assembly_pair_variants_identical_breakpoints() {
    let (fx, out) = &*FIXTURE_A;
    assert_variant_truth(fx, out);
}

#[test]
fn assembly_pair_variants_shifted_and_reverse_complement() {
    let (fx, out) = &*FIXTURE_B;
    assert_variant_truth(fx, out);
    let decoy = fx.decoy_ref.as_deref().unwrap();
    assert!(
        out.variants.iter().all(|r| r[0] != decoy),
        "variant rows reported on the unrelated decoy contig"
    );
    assert!(
        out.rows_of("SNP").iter().any(|r| r[3] == "q1" && r[18] == "-1"),
        "SNPs from the reverse-complemented query contig q1 must carry Query_Direction -1"
    );
}

/// The `.snpdiffs` contract CSP2 consumes: header invariants, one BED row per kept
/// alignment block with real query coordinates, coverage computed from those blocks,
/// forward query coordinates on reverse-strand rows, and no duplicated positions.
#[test]
fn assembly_pair_snpdiffs_contract() {
    let (fx, out) = &*FIXTURE_B;
    let ref_seqs: HashMap<&str, &[u8]> = fx.ref_records.iter().map(|(id, s)| (id.as_str(), s.as_slice())).collect();
    let qry_seqs: HashMap<&str, &[u8]> = fx.qry_records.iter().map(|(id, s)| (id.as_str(), s.as_slice())).collect();
    let ref_total: usize = ref_seqs.values().map(|s| s.len()).sum();
    let qry_total: usize = qry_seqs.values().map(|s| s.len()).sum();
    let (decoy_ref, decoy_qry) = (fx.decoy_ref.as_deref().unwrap(), fx.decoy_qry.as_deref().unwrap());

    // Header counts match the rows.
    assert_eq!(out.header_f64("SNPs") as usize, out.rows_of("SNP").len(), "header SNPs vs SNP rows");
    assert_eq!(out.header_f64("Indels") as usize, out.rows_of("Indel").len(), "header Indels vs Indel rows");
    for key in ["Reference_Percent_Aligned", "Query_Percent_Aligned", "Median_Percent_Identity"] {
        let v = out.header_f64(key);
        assert!((0.0..=100.0).contains(&v), "{key} out of [0, 100]: {v}");
    }

    // BED rows are alignment blocks with real, self-consistent coordinates.
    let blocks = out.covered_bed();
    assert!(!blocks.is_empty(), "no covered BED rows");
    let num = |r: &Vec<String>, i: usize| -> usize {
        r[i].parse().unwrap_or_else(|_| panic!("BED column {i} not numeric: {r:?}"))
    };
    let mut ref_cov: HashMap<&str, HashSet<usize>> = HashMap::new();
    let mut qry_cov: HashMap<&str, HashSet<usize>> = HashMap::new();
    for r in &blocks {
        let (rc, qc) = (r[0].as_str(), r[5].as_str());
        assert!(rc != decoy_ref && qc != decoy_qry, "decoy contig in a covered BED row: {r:?}");
        let rseq = ref_seqs.get(rc).unwrap_or_else(|| panic!("BED Ref_Contig not in reference FASTA: {r:?}"));
        let qseq = qry_seqs.get(qc).unwrap_or_else(|| panic!("BED Query_Contig not in query FASTA: {r:?}"));
        let (rs, re, rlen, ralen) = (num(r, 1), num(r, 2), num(r, 3), num(r, 4));
        let (qs, qe, qlen, qalen) = (num(r, 6), num(r, 7), num(r, 8), num(r, 9));
        assert!(rs < re && re <= rlen && rlen == rseq.len(), "bad reference interval: {r:?}");
        assert!(qs < qe && qe <= qlen && qlen == qseq.len(), "bad query interval: {r:?}");
        assert_eq!(ralen, re - rs, "Ref_Aligned must equal Ref_End - Ref_Start: {r:?}");
        assert_eq!(qalen, qe - qs, "Query_Aligned must equal Query_End - Query_Start: {r:?}");
        let iden: f64 = r[10].parse().unwrap();
        assert!((80.0..=100.0).contains(&iden), "block Perc_Iden out of range: {r:?}");
        let rset = ref_cov.entry(rc).or_default();
        for pos in rs..re {
            assert!(rset.insert(pos), "reference base {rc}:{pos} covered by two kept blocks");
        }
        qry_cov.entry(qc).or_default().extend(qs..qe);
    }

    // Header coverage equals the union of the blocks, and the blocks cover nearly all of the
    // homologous sequence (everything except the decoys).
    let ref_pct = 100.0 * ref_cov.values().map(HashSet::len).sum::<usize>() as f64 / ref_total as f64;
    let qry_pct = 100.0 * qry_cov.values().map(HashSet::len).sum::<usize>() as f64 / qry_total as f64;
    assert!((ref_pct - out.header_f64("Reference_Percent_Aligned")).abs() <= 0.01, "Reference_Percent_Aligned {} vs BED union {ref_pct:.2}", out.header_f64("Reference_Percent_Aligned"));
    assert!((qry_pct - out.header_f64("Query_Percent_Aligned")).abs() <= 0.01, "Query_Percent_Aligned {} vs BED union {qry_pct:.2}", out.header_f64("Query_Percent_Aligned"));
    let homologous_ref = 100.0 * GENOME_LEN as f64 / ref_total as f64;
    let homologous_qry = 100.0 * GENOME_LEN as f64 / qry_total as f64;
    assert!(ref_pct >= 0.98 * homologous_ref, "reference breadth {ref_pct:.2}% below 98% of {homologous_ref:.2}%");
    assert!(qry_pct >= 0.98 * homologous_qry, "query breadth {qry_pct:.2}% below 98% of {homologous_qry:.2}%");

    // Uncovered regions are reported as `.` rows, like compileMUMmer.py.
    assert!(
        out.bed.iter().any(|r| r[0] == decoy_ref && r[1] == "0" && r[2] == DECOY_LEN.to_string() && r[5] == "."),
        "missing the uncovered-reference `.` row for the decoy contig"
    );

    // Variant rows: one per reference position, bases and coordinates consistent with the
    // FASTAs, and alignment fields taken from the containing BED block.
    let mut seen = HashSet::new();
    for r in out.rows_of("SNP") {
        assert!(seen.insert((r[0].clone(), r[2].clone())), "duplicate reference position {}:{}", r[0], r[2]);
        let ref_pos: usize = r[2].parse().unwrap();
        let qry_pos: usize = r[5].parse().unwrap();
        assert_eq!(r[1], (ref_pos - 1).to_string(), "Start_Ref must be Ref_Pos - 1: {r:?}");
        assert_eq!(r[4], (qry_pos - 1).to_string(), "Start_Query must be Query_Pos - 1: {r:?}");
        let rseq = ref_seqs[r[0].as_str()];
        assert_eq!(r[12].as_bytes(), [rseq[ref_pos - 1]], "Ref_Base does not match the reference FASTA: {r:?}");
        let qbase = qry_seqs[r[3].as_str()][qry_pos - 1];
        let expected = match r[18].as_str() {
            "1" => qbase,
            "-1" => complement(qbase),
            other => panic!("Query_Direction must be 1 or -1, got {other}: {r:?}"),
        };
        assert_eq!(r[13].as_bytes(), [expected], "Query_Pos is not a forward query coordinate: {r:?}");
        assert!(
            blocks.iter().any(|b| b[0] == r[0] && b[5] == r[3] && b[1] == r[8] && b[2] == r[9] && b[6] == r[10] && b[7] == r[11]),
            "variant row's block fields do not match any covered BED row: {r:?}"
        );
    }
}
