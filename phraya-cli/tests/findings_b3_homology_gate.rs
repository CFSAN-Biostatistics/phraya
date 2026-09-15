/// Regression test for CSP2 spike finding B3: Case 4 (contigs-only, no reference)
/// generated a dense `i<j` all-pairs task list with no homology pre-filter — two draft
/// assemblies with a few hundred contigs each produce tens of thousands of tasks, almost
/// all between contigs that share no real homology.
///
/// Fixture: 5 contigs, no reference, all in one input file (single-file Case 4 MSA).
/// Two genuinely related pairs (mutated copies of each other, ~95% similar — the
/// same fixture convention `issue_210_fixture_tests.rs` documents) plus one unrelated
/// singleton. Independent-seed LCG sequences are ~75% divergent from each other per
/// that same documentation, well below any reasonable homology threshold.
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

/// LCG-based diverse DNA sequence. Independent seeds are ~75% divergent from each
/// other (4-letter alphabet, uncorrelated streams).
fn diverse_dna(len: usize, seed: u64) -> Vec<u8> {
    let mut state = seed;
    let bases = [b'A', b'C', b'G', b'T'];
    (0..len)
        .map(|_| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            bases[((state >> 32) as usize) % 4]
        })
        .collect()
}

/// Mutated copy of `base` with `mutation_rate` fraction of positions substituted.
fn mutate_sequence(base: &[u8], mutation_rate: f64, seed: u64) -> Vec<u8> {
    let mut x = seed;
    let bases = [b'A', b'C', b'G', b'T'];
    base.iter()
        .map(|&b| {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            let rand_01 = ((x >> 32) as u32) as f64 / (u32::MAX as f64 + 1.0);
            if rand_01 < mutation_rate {
                // Filter out the current base *before* indexing, matching
                // `issue_210_fixture_tests.rs`'s convention — a "keep retrying" loop
                // reusing an already-consumed random draw can pick the same excluded
                // value forever (found via B4's regression test hanging at 40kb+).
                let candidates: Vec<u8> = bases.iter().copied().filter(|&c| c != b).collect();
                let idx = ((x >> 16) as usize) % candidates.len();
                candidates[idx]
            } else {
                b
            }
        })
        .collect()
}

/// Write the 5-contig, no-reference, single-file fixture. Returns the FASTA path and
/// the ordered list of contig IDs (matching plan-time index order).
fn write_fixture(dir: &Path) -> (PathBuf, Vec<&'static str>) {
    let contig_a = diverse_dna(8000, 1);
    let contig_a_mut = mutate_sequence(&contig_a, 0.05, 100);
    let contig_b = diverse_dna(8000, 2);
    let contig_b_mut = mutate_sequence(&contig_b, 0.05, 200);
    let contig_c = diverse_dna(8000, 3);

    let ids = vec!["contigA", "contigA_mut", "contigB", "contigB_mut", "contigC"];
    let seqs = [contig_a, contig_a_mut, contig_b, contig_b_mut, contig_c];

    let path = dir.join("assembly.fa");
    let mut content = String::new();
    for (id, seq) in ids.iter().zip(seqs.iter()) {
        content.push_str(&format!(">{id}\n{}\n", String::from_utf8(seq.clone()).unwrap()));
    }
    std::fs::write(&path, content).unwrap();
    (path, ids)
}

/// Default homology gate keeps the two genuinely related pairs and drops (most/all of)
/// the unrelated ones — task count well below the full 10-pair dense all-pairs list.
#[test]
fn findings_b3_default_gate_keeps_related_pairs_drops_unrelated() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();
    let (fasta_path, ids) = write_fixture(p);
    let plan_path = p.join("plan.phrayaplan");

    let out = phraya(&[
        "plan",
        "--inputs",
        fasta_path.to_str().unwrap(),
        "--output",
        plan_path.to_str().unwrap(),
    ]);
    assert!(
        out.status.success(),
        "phraya plan failed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );

    let plan = phraya_io::plan::read_plan(&plan_path).unwrap();
    assert_eq!(plan.use_case, phraya_io::plan::UseCase::ContigsOnly);

    let a = ids.iter().position(|&x| x == "contigA").unwrap() as u32;
    let a_mut = ids.iter().position(|&x| x == "contigA_mut").unwrap() as u32;
    let b = ids.iter().position(|&x| x == "contigB").unwrap() as u32;
    let b_mut = ids.iter().position(|&x| x == "contigB_mut").unwrap() as u32;

    let pair = |x: u32, y: u32| if x < y { (x, y) } else { (y, x) };

    assert!(
        plan.task_list.contains(&pair(a, a_mut)),
        "gated task list should keep the genuinely related (contigA, contigA_mut) pair, \
         got: {:?}",
        plan.task_list
    );
    assert!(
        plan.task_list.contains(&pair(b, b_mut)),
        "gated task list should keep the genuinely related (contigB, contigB_mut) pair, \
         got: {:?}",
        plan.task_list
    );

    // Full dense all-pairs over 5 contigs is C(5,2) = 10. The gate must strictly
    // reduce this — independent-seed contigs share no real homology.
    assert!(
        plan.task_list.len() < 10,
        "homology gate should reduce task count below the full dense all-pairs list \
         (10), got {} tasks: {:?}",
        plan.task_list.len(),
        plan.task_list
    );
}

/// `--no-homology-gate` restores the full dense all-pairs list (all 10 pairs for 5
/// contigs), overriding the default gate.
#[test]
fn findings_b3_no_homology_gate_restores_full_all_pairs() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();
    let (fasta_path, _ids) = write_fixture(p);
    let plan_path = p.join("plan.phrayaplan");

    let out = phraya(&[
        "plan",
        "--inputs",
        fasta_path.to_str().unwrap(),
        "--no-homology-gate",
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
        plan.task_list.len(),
        10,
        "--no-homology-gate should restore the full C(5,2)=10 dense all-pairs list, \
         got {} tasks: {:?}",
        plan.task_list.len(),
        plan.task_list
    );
}

/// A very high `--min-homology` threshold drops even the genuinely related pairs —
/// confirms the flag actually changes gating behavior, not just accepted-and-ignored.
#[test]
fn findings_b3_min_homology_threshold_is_respected() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();
    let (fasta_path, _ids) = write_fixture(p);
    let plan_path = p.join("plan.phrayaplan");

    let out = phraya(&[
        "plan",
        "--inputs",
        fasta_path.to_str().unwrap(),
        "--min-homology",
        "0.999",
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
        plan.task_list.len(),
        0,
        "an unreachably high --min-homology threshold should drop every pair, \
         got {} tasks: {:?}",
        plan.task_list.len(),
        plan.task_list
    );
}
