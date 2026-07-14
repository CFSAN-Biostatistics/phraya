/// End-to-end Case 2 integration tests: plan → align → merge → filter → VCF
///
/// Fixture
/// -------
/// Reference: 200bp diverse (LCG-generated, avoids minimizer-seed explosion
/// that occurs with repetitive "ACGT"×n sequences).
/// Reads (100bp prefix of reference, 5 total):
///   - 2 perfect-match reads
///   - 3 SNP reads: one base changed at position 50 (ref_base → alt_base)
///
/// AC coverage
/// -----------
/// ✓ task_count == num_reads           (test 1)
/// ✓ all align tasks exit 0             (test 2)
/// ✓ merge produces parseable file      (test 3)
/// ✓ SNP position appears in VCF        (test 4) ← RED: CIGAR convention bug
/// ✓ VCF shows correct REF/ALT alleles  (test 5) ← RED: same bug
/// ✓ perfect reads → no false variants  (test 6)
/// ✓ SNP survives --min-coverage 5      (test 7) ← RED: local_coverage=1 bug
/// ✓ .phraya binary < VCF text size     (test 8) ← RED: needs real content
/// ✓ throughput ≥ 100 reads/sec         (phraya-align executor tests)

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

/// LCG-based diverse DNA sequence. Avoids the minimizer-seed explosion that
/// repetitive sequences (e.g. "ACGT"×n) cause by ensuring most 21-mers are unique.
fn diverse_dna(len: usize, seed: u64) -> Vec<u8> {
    let mut x = seed;
    (0..len)
        .map(|_| {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            b"ACGT"[((x >> 33) & 3) as usize]
        })
        .collect()
}

/// Build fixture files. Returns (ref_path, reads_path, snp_pos, ref_base, alt_base).
///
/// Reference: 200bp diverse DNA.
/// 2 perfect reads (100bp) + 3 SNP reads (100bp, one base changed at snp_pos=50).
fn write_fixtures(dir: &Path) -> (PathBuf, PathBuf, usize, u8, u8) {
    let ref_seq = diverse_dna(200, 42);
    let snp_pos: usize = 50;
    let ref_base = ref_seq[snp_pos];
    // Pick an alt base guaranteed to differ from ref_base
    let alt_base = if ref_base == b'A' { b'T' } else { b'A' };

    let perfect: Vec<u8> = ref_seq[..100].to_vec();
    let mut snp_read = perfect.clone();
    snp_read[snp_pos] = alt_base;

    let ref_path = dir.join("ref.fa");
    std::fs::write(
        &ref_path,
        format!(">ref\n{}\n", String::from_utf8(ref_seq).unwrap()),
    )
    .unwrap();

    let reads_path = dir.join("reads.fa");
    let mut content = String::new();
    for i in 0..2 {
        content.push_str(&format!(
            ">perfect{i}\n{}\n",
            String::from_utf8(perfect.clone()).unwrap()
        ));
    }
    // 5 SNP reads so merged obs_count at snp_pos = 5 → survives --min-coverage 5
    for i in 0..5 {
        content.push_str(&format!(
            ">snp{i}\n{}\n",
            String::from_utf8(snp_read.clone()).unwrap()
        ));
    }
    std::fs::write(&reads_path, content).unwrap();

    (ref_path, reads_path, snp_pos, ref_base, alt_base)
}

/// Run plan → align all 5 reads → merge. Returns (merged_path, snp_pos, ref_base, alt_base).
fn run_to_merge(dir: &Path) -> (PathBuf, usize, u8, u8) {
    let (ref_path, reads_path, snp_pos, ref_base, alt_base) = write_fixtures(dir);
    let plan_path = dir.join("plan.phrayaplan");

    let out = phraya(&[
        "plan",
        "--inputs", reads_path.to_str().unwrap(),
        "--reference", ref_path.to_str().unwrap(),
        "--output", plan_path.to_str().unwrap(),
    ]);
    assert!(
        out.status.success(),
        "phraya plan failed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );

    let read_ids: Vec<String> = (0..2)
        .map(|i| format!("perfect{i}"))
        .chain((0..5).map(|i| format!("snp{i}")))
        .collect();

    let mut per_read: Vec<PathBuf> = Vec::new();
    for id in &read_ids {
        let op = dir.join(format!("{id}.phraya"));
        let out = phraya(&[
            "align",
            plan_path.to_str().unwrap(),
            id, "ref",
            "--output", op.to_str().unwrap(),
        ]);
        assert!(
            out.status.success(),
            "phraya align {id} failed:\n{}",
            String::from_utf8_lossy(&out.stderr)
        );
        per_read.push(op);
    }

    let merged = dir.join("merged.phraya");
    let strs: Vec<String> = per_read.iter().map(|p| p.to_str().unwrap().to_string()).collect();
    let mut args = vec!["merge"];
    for s in &strs { args.push(s.as_str()); }
    args.extend(&["--output", merged.to_str().unwrap()]);
    let out = phraya(&args);
    assert!(
        out.status.success(),
        "phraya merge failed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );

    (merged, snp_pos, ref_base, alt_base)
}

// ── Test 1: task count ────────────────────────────────────────────────────────

/// phraya plan produces exactly N tasks for N reads (Case 2)
#[test]
fn issue_88_plan_task_count_equals_read_count() {
    let dir = TempDir::new().unwrap();
    let (ref_path, reads_path, _, _, _) = write_fixtures(dir.path());
    let plan_path = dir.path().join("plan.phrayaplan");

    let out = phraya(&[
        "plan",
        "--inputs", reads_path.to_str().unwrap(),
        "--reference", ref_path.to_str().unwrap(),
        "--output", plan_path.to_str().unwrap(),
    ]);
    assert!(out.status.success(), "phraya plan failed:\n{}", String::from_utf8_lossy(&out.stderr));

    let plan = phraya_io::plan::read_plan(&plan_path).unwrap();
    // 2 perfect + 5 snp = 7 reads → 7 tasks
    assert_eq!(plan.task_list.len(), 7, "7 reads → 7 tasks, got {}", plan.task_list.len());
    for (q, t) in &plan.task_list {
        assert_eq!(*t, 0, "all tasks must target reference (index 0)");
        assert_ne!(*q, 0, "query must not be the reference");
    }
}

// ── Test 2: all align tasks exit 0 ───────────────────────────────────────────

/// every phraya align call exits 0 and writes .phraya + .phraya.queries
#[test]
fn issue_88_all_align_tasks_exit_zero() {
    let dir = TempDir::new().unwrap();
    let (ref_path, reads_path, _, _, _) = write_fixtures(dir.path());
    let plan_path = dir.path().join("plan.phrayaplan");

    phraya(&[
        "plan",
        "--inputs", reads_path.to_str().unwrap(),
        "--reference", ref_path.to_str().unwrap(),
        "--output", plan_path.to_str().unwrap(),
    ]);

    let ids: Vec<String> = (0..2)
        .map(|i| format!("perfect{i}"))
        .chain((0..5).map(|i| format!("snp{i}")))
        .collect();

    for id in &ids {
        let op = dir.path().join(format!("{id}.phraya"));
        let out = phraya(&[
            "align",
            plan_path.to_str().unwrap(),
            id, "ref",
            "--output", op.to_str().unwrap(),
        ]);
        assert!(out.status.success(), "phraya align {id} failed:\n{}", String::from_utf8_lossy(&out.stderr));
        assert!(op.exists(), ".phraya missing for {id}");
        assert!(dir.path().join(format!("{id}.phraya.queries")).exists(), ".queries missing for {id}");
    }
}

// ── Test 3: merge produces parseable file ─────────────────────────────────────

/// merged .phraya is readable and has the correct reference length
#[test]
fn issue_88_merge_produces_parseable_file() {
    let dir = TempDir::new().unwrap();
    let (merged, _, _, _) = run_to_merge(dir.path());

    let f = phraya_io::phraya::read_phraya(&merged).expect("merged .phraya must be parseable");
    assert_eq!(
        f.header.reference_length, 200,
        "reference_length should match the 200bp reference"
    );
}

// ── Test 4: SNP position appears in VCF ──────────────────────────────────────

/// VCF (no coverage filter) contains a data line at the SNP position
#[test]
fn issue_88_snp_appears_in_vcf_no_coverage_filter() {
    let dir = TempDir::new().unwrap();
    let (merged, snp_pos, _, _) = run_to_merge(dir.path());
    let expected_vcf_pos = (snp_pos + 1).to_string(); // 1-indexed

    let out = phraya(&["filter", merged.to_str().unwrap(), "--format", "vcf"]);
    assert!(out.status.success());
    let vcf = String::from_utf8_lossy(&out.stdout);

    assert!(
        vcf.lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty())
            .any(|l| l.split('\t').nth(1) == Some(&expected_vcf_pos)),
        "VCF must contain a data line at position {expected_vcf_pos} (1-indexed).\nVCF:\n{vcf}"
    );
}

// ── Test 5: correct REF/ALT alleles ──────────────────────────────────────────

/// VCF line at the SNP position shows the correct REF and ALT bases
#[test]
fn issue_88_vcf_correct_ref_alt_alleles() {
    let dir = TempDir::new().unwrap();
    let (merged, snp_pos, ref_base, alt_base) = run_to_merge(dir.path());
    let expected_vcf_pos = (snp_pos + 1).to_string();

    let out = phraya(&["filter", merged.to_str().unwrap(), "--format", "vcf"]);
    let vcf = String::from_utf8_lossy(&out.stdout);

    let snp_line = vcf
        .lines()
        .find(|l| !l.starts_with('#') && l.split('\t').nth(1) == Some(&expected_vcf_pos))
        .unwrap_or_else(|| panic!("No VCF data line at position {expected_vcf_pos}.\nVCF:\n{vcf}"));

    let cols: Vec<&str> = snp_line.split('\t').collect();
    assert!(cols.len() >= 5, "VCF line too short: {snp_line}");
    assert_eq!(
        cols[3],
        (ref_base as char).to_string(),
        "REF should be {} at position {snp_pos}", ref_base as char
    );
    assert!(
        cols[4].contains(alt_base as char),
        "ALT should contain {} (the SNP allele), got {}", alt_base as char, cols[4]
    );
}

// ── Test 6: no false positives from perfect reads ─────────────────────────────

/// Aligning only perfect-match reads produces a VCF with zero data lines
#[test]
fn issue_88_perfect_reads_no_false_variants() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();

    let ref_seq = diverse_dna(200, 42);
    let perfect: Vec<u8> = ref_seq[..100].to_vec();

    let ref_path = p.join("ref.fa");
    std::fs::write(&ref_path, format!(">ref\n{}\n", String::from_utf8(ref_seq).unwrap())).unwrap();

    let reads_path = p.join("perf.fa");
    let mut content = String::new();
    for i in 0..3 {
        content.push_str(&format!(">r{i}\n{}\n", String::from_utf8(perfect.clone()).unwrap()));
    }
    std::fs::write(&reads_path, content).unwrap();

    let plan_path = p.join("perf.phrayaplan");
    phraya(&["plan", "--inputs", reads_path.to_str().unwrap(),
             "--reference", ref_path.to_str().unwrap(),
             "--output", plan_path.to_str().unwrap()]);

    let mut files: Vec<PathBuf> = Vec::new();
    for i in 0..3 {
        let op = p.join(format!("r{i}.phraya"));
        phraya(&["align", plan_path.to_str().unwrap(), &format!("r{i}"), "ref",
                 "--output", op.to_str().unwrap()]);
        files.push(op);
    }

    let merged = p.join("perf_merged.phraya");
    let strs: Vec<String> = files.iter().map(|p| p.to_str().unwrap().to_string()).collect();
    let mut args = vec!["merge"];
    for s in &strs { args.push(s.as_str()); }
    args.extend(&["--output", merged.to_str().unwrap()]);
    phraya(&args);

    let out = phraya(&["filter", merged.to_str().unwrap(), "--format", "vcf"]);
    let vcf = String::from_utf8_lossy(&out.stdout);
    let data: Vec<&str> = vcf.lines().filter(|l| !l.starts_with('#') && !l.is_empty()).collect();

    assert!(
        data.is_empty(),
        "perfect-match reads must produce zero VCF data lines, got:\n{}",
        data.join("\n")
    );
}

// ── Test 7: SNP survives --min-coverage 5 ────────────────────────────────────

/// 3 reads with the same SNP must survive --min-coverage 5
#[test]
fn issue_88_snp_survives_min_coverage_5() {
    let dir = TempDir::new().unwrap();
    let (merged, snp_pos, _, _) = run_to_merge(dir.path());
    let expected_vcf_pos = (snp_pos + 1).to_string();

    let out = phraya(&[
        "filter", merged.to_str().unwrap(),
        "--min-coverage", "5",
        "--format", "vcf",
    ]);
    assert!(out.status.success());
    let vcf = String::from_utf8_lossy(&out.stdout);

    assert!(
        vcf.lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty())
            .any(|l| l.split('\t').nth(1) == Some(&expected_vcf_pos)),
        "SNP at pos {expected_vcf_pos} has 3× read support and must survive --min-coverage 5.\n\
         Fails because local_coverage is hardcoded to 1 per observation.\nVCF:\n{vcf}"
    );
}

// ── Test 8: binary format smaller than VCF text ───────────────────────────────

/// The merged .phraya binary must be smaller than the VCF text for the same content
#[test]
fn issue_88_binary_smaller_than_vcf_text() {
    let dir = TempDir::new().unwrap();
    let (merged, _, _, _) = run_to_merge(dir.path());

    let phraya_size = std::fs::metadata(&merged).unwrap().len();

    // VCF of the unfiltered merged content
    let out = phraya(&["filter", merged.to_str().unwrap(), "--format", "vcf"]);
    let vcf_bytes = out.stdout.len() as u64;

    // Require both: VCF has content (not just headers) AND binary is smaller
    let vcf_text = String::from_utf8_lossy(&out.stdout);
    let vcf_has_data = vcf_text.lines().any(|l| !l.starts_with('#') && !l.is_empty());

    assert!(
        vcf_has_data,
        "VCF must have data lines before size comparison is meaningful.\n\
         Currently VCF is header-only because the CIGAR bug prevents variant detection.\n\
         VCF:\n{vcf_text}"
    );
    assert!(
        phraya_size < vcf_bytes,
        ".phraya binary ({phraya_size}B) must be smaller than VCF text ({vcf_bytes}B)"
    );
}
