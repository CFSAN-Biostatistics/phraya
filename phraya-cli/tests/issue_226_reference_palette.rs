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
    let mut full = vec!["run", "--quiet", "--bin", "phraya", "--manifest-path", m.to_str().unwrap(), "--"];
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

fn plan_with_refs(dir: &Path, plan_name: &str, reads: &Path, references: &[&Path]) {
    let mut args: Vec<String> = vec![
        "plan".to_string(),
        "--inputs".to_string(),
        reads.to_str().unwrap().to_string(),
    ];
    for r in references {
        args.push("--reference".to_string());
        args.push(r.to_str().unwrap().to_string());
    }
    args.push("--output".to_string());
    args.push(dir.join(plan_name).to_str().unwrap().to_string());
    let out = run(&args.iter().map(|s| s.as_str()).collect::<Vec<_>>());
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

/// palette-perf-and-pa-scoring plan, step 1: a reference space takes its name from the
/// FASTA record ID at plan time, so palette output filenames are the (sanitized) record ID
/// instead of a content-hash prefix. This is the observable contract the benchmark's PA
/// scorer depends on to match a cross-space placement's space label to a read's true
/// chromosome.
#[test]
fn palette_output_is_named_by_fasta_record_id() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();
    let reads = write_fasta(p, "reads.fa", &[("read1", READ)]);
    let two_chrom = write_fasta(p, "two_chrom.fa", &[("chr1", REF_A), ("chr2", REF_B)]);
    plan_with_ref(p, "plan.phrayaplan", &reads, &two_chrom);

    let out_dir = p.join("out");
    let out = run(&[
        "align",
        p.join("plan.phrayaplan").to_str().unwrap(),
        "--reference",
        two_chrom.to_str().unwrap(),
        "--output",
        out_dir.to_str().unwrap(),
    ]);
    assert!(out.status.success(), "align failed: {}", String::from_utf8_lossy(&out.stderr));

    let names: Vec<_> = std::fs::read_dir(&out_dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .map(|e| e.file_name())
        .collect();
    assert!(out_dir.join("chr1.phraya").exists(), "expected chr1.phraya, found: {names:?}");
    assert!(out_dir.join("chr2.phraya").exists(), "expected chr2.phraya, found: {names:?}");
    assert!(out_dir.join("cross_space.phraya.queries").exists(), "cross-space sidecar written");
}

fn dna(len: usize, seed: u64) -> String {
    let mut x = seed;
    (0..len)
        .map(|_| {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"ACGT"[((x >> 33) & 3) as usize] as char
        })
        .collect()
}

/// A space no read comes from must be *empty*, not full of fabricated variants.
///
/// Palette mode aligns every read against every space, so the `(0,0)` fallback anchor used
/// to force a ~50-60% identity alignment of every read against every non-homologous space,
/// each depositing ~`read_len/2` `VariantObservation`s. The cross-space sidecar filters
/// those out at its 0.95 identity threshold, but the per-space `.phraya` does not — measured
/// at 1.87 MB of fabricated calls per non-homologous space against 42 KB of real signal in
/// the one space the reads actually came from.
#[test]
fn non_homologous_space_yields_no_variants() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();

    // Target/read ratio must exceed the fallback scope (10x) for the guard to engage; the
    // shared REF_A/READ consts are deliberately tiny and comparable-length.
    let chr_x = dna(3_000, 11);
    let chr_y = dna(3_000, 99);
    let flip = |c: char| match c {
        'A' => 'C',
        'C' => 'G',
        'G' => 'T',
        _ => 'A',
    };
    let read_records: Vec<(String, String)> = (0..5)
        .map(|i| {
            let off = 100 + i * 400;
            let mut b: Vec<char> = chr_x[off..off + 150].chars().collect();
            b[40] = flip(b[40]);
            b[90] = flip(b[90]);
            (format!("read_{i}"), b.into_iter().collect())
        })
        .collect();

    let reads = write_fasta(
        p,
        "reads.fa",
        &read_records
            .iter()
            .map(|(id, s)| (id.as_str(), s.as_str()))
            .collect::<Vec<_>>(),
    );
    let two_chrom = write_fasta(
        p,
        "two_chrom.fa",
        &[("chrX", chr_x.as_str()), ("chrY", chr_y.as_str())],
    );
    plan_with_ref(p, "plan.phrayaplan", &reads, &two_chrom);

    let out_dir = p.join("out");
    let out = run(&[
        "align",
        p.join("plan.phrayaplan").to_str().unwrap(),
        "--reference",
        two_chrom.to_str().unwrap(),
        "--output",
        out_dir.to_str().unwrap(),
    ]);
    assert!(out.status.success(), "align failed: {}", String::from_utf8_lossy(&out.stderr));

    let x = phraya_io::phraya::read_phraya(&out_dir.join("chrX.phraya")).unwrap();
    let y = phraya_io::phraya::read_phraya(&out_dir.join("chrY.phraya")).unwrap();

    // chrX is the control: the reads do come from here, carrying two substitutions each, so
    // real coverage and real calls must survive. Without this the test would pass trivially
    // if alignment stopped emitting anything at all.
    assert!(
        x.header.covered_positions.unwrap_or(0) > 0,
        "reads must still place in the space they came from"
    );
    assert!(
        !x.observations.is_empty(),
        "the planted substitutions must still be called in chrX"
    );

    assert!(
        y.observations.is_empty(),
        "chrY shares no read; it must carry zero variant observations, got {}",
        y.observations.len()
    );
    assert_eq!(
        y.header.covered_positions,
        Some(0),
        "chrY must have zero covered positions"
    );
    assert_eq!(
        y.header.reference_length, 3_000,
        "a skipped space still reports its true reference length"
    );
    }

    /// An input contig whose bytes are identical to a reference contig (and whose header is
    /// renamed) must NOT be dropped from the read pool. This is the assembly-mode form of
    /// issue-17's "identical contig" failure. With the old content-hash exclusion, such a contig
    /// was dropped (log: `0 reads`); with source-based exclusion (`phraya plan` folds the first
    /// reference to `input_files[0]`), only that folded file is skipped.
    #[test]
    fn assembly_identical_contig_is_aligned_not_dropped() {
        let dir = TempDir::new().unwrap();
        let p = dir.path();
        // Reference: a ~4 kb seeded contig.
        let ref_seq = dna(4_000, 7);
        let refa = write_fasta(p, "refA.fa", &[("refA", ref_seq.as_str())]);
        // Query is a byte-identical copy of a reference contig, header renamed.
        let reads = write_fasta(p, "reads.fa", &[("identical_contig", ref_seq.as_str())]);
        plan_with_ref(p, "plan.phrayaplan", &reads, &refa);

        let out_dir = p.join("out");
        let out = run(&[
            "align",
            "--reference",
            refa.to_str().unwrap(),
            p.join("plan.phrayaplan").to_str().unwrap(),
            "--strategy",
            "assembly",
            "--output",
            out_dir.to_str().unwrap(),
        ]);
        assert!(out.status.success(), "align failed: {}", String::from_utf8_lossy(&out.stderr));
        let stderr = String::from_utf8_lossy(&out.stderr);
        // Source exclusion keeps the renamed contig; the log must count it as a read.
        assert!(
            stderr.contains("Reference-palette mode: 1 reference space(s), 1 reads"),
            "byte-identical contig must enter the read pool; log was:\n{stderr}"
        );
    }

    /// A query contig byte-identical to a contig in a *second* `--reference` file (which is
    /// never folded into `input_files[0]`) must also be aligned, not dropped. Under the old
    /// content-hash exclusion it matched the palette and was dropped; source exclusion skips
    /// only `input_files[0]`, so it is a read and is placed in that space.
    #[test]
    fn identical_contig_matching_second_reference_is_aligned() {
        let dir = TempDir::new().unwrap();
        let p = dir.path();
        let seq_b = dna(4_000, 5);
        let refa = write_fasta(p, "refA.fa", &[("refA", dna(6_000, 7).as_str())]);
        let refb = write_fasta(p, "refB.fa", &[("refB", seq_b.as_str())]);
        // Query == refB exactly (renamed). refB is a space but not in input_files[0].
        let reads = write_fasta(p, "reads.fa", &[("identical_to_B", seq_b.as_str())]);
        plan_with_refs(p, "plan.phrayaplan", &reads, &[&refa, &refb]);

        let out_dir = p.join("out");
        let out = run(&[
            "align",
            p.join("plan.phrayaplan").to_str().unwrap(),
            "--reference",
            refa.to_str().unwrap(),
            "--reference",
            refb.to_str().unwrap(),
            "--strategy",
            "assembly",
            "--output",
            out_dir.to_str().unwrap(),
        ]);
        let stderr = String::from_utf8_lossy(&out.stderr);
        // The query is NOT a present reference this invocation, so the read-path hash guard
        // does not drop it; source exclusion skips only refA (input_files[0]), keeping the query.
        assert!(
            stderr.contains("Reference-palette mode: 2 reference space(s), 1 reads"),
            "query identical to a 2nd reference must be a read; log was:\n{stderr}"
        );
        // And it must place in the refB space with zero variants (100% identity).
        let ph = phraya_io::phraya::read_phraya(&out_dir.join("refB.phraya")).unwrap();
        assert!(
            !ph.alignments.is_empty(),
            "identical-contig-to-B must place in the refB space"
        );
        assert_eq!(ph.observations.len(), 0, "100%-identical contig yields zero variants, got {}", ph.observations.len());
    }
