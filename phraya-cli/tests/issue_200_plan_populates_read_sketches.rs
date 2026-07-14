/// RED acceptance test for issue #200: `phraya plan` populates `read_sketches`
/// (content-hash-keyed sketch cache for reads, added in #200's first slice, PR #225)
/// during its existing read pass, rather than leaving the field permanently empty.
///
/// `run_plan` (phraya-cli/src/main.rs) already computes a `MinimizerSketch` per input
/// sequence and stores it in `kmer_index` (keyed by sequence ID, covers both reference
/// and reads). This test checks that reads' sketches are *also* stored in
/// `read_sketches`, keyed by `read_content_hash(bases)` — the amortization #200 asks
/// for, so a later pipeline stage can look a read's sketch up by content even if its ID
/// changed or it was re-batched.
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

fn diverse_dna(len: usize, seed: u64) -> Vec<u8> {
    let mut x = seed;
    (0..len)
        .map(|_| {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"ACGT"[((x >> 33) & 3) as usize]
        })
        .collect()
}

fn write_fixtures(dir: &Path) -> (PathBuf, PathBuf, Vec<u8>) {
    let ref_seq = diverse_dna(200, 42);
    let read_bases: Vec<u8> = ref_seq[..100].to_vec();

    let ref_path = dir.join("ref.fa");
    std::fs::write(
        &ref_path,
        format!(">ref\n{}\n", String::from_utf8(ref_seq).unwrap()),
    )
    .unwrap();

    let reads_path = dir.join("reads.fa");
    std::fs::write(
        &reads_path,
        format!(
            ">read0\n{}\n",
            String::from_utf8(read_bases.clone()).unwrap()
        ),
    )
    .unwrap();

    (ref_path, reads_path, read_bases)
}

#[test]
fn plan_stores_read_sketch_keyed_by_content_hash() {
    let dir = TempDir::new().unwrap();
    let (ref_path, reads_path, read_bases) = write_fixtures(dir.path());
    let plan_path = dir.path().join("plan.phrayaplan");

    let out = phraya(&[
        "plan",
        "--inputs",
        reads_path.to_str().unwrap(),
        "--reference",
        ref_path.to_str().unwrap(),
        "--output",
        plan_path.to_str().unwrap(),
    ]);
    assert!(
        out.status.success(),
        "phraya plan failed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );

    let plan = phraya_io::plan::read_plan(&plan_path).unwrap();
    let hash = phraya_io::plan::read_content_hash(&read_bases);

    assert!(
        !plan.read_sketches.is_empty(),
        "phraya plan should populate read_sketches during its existing read pass, \
         but it is empty — reads' sketches are only being stored in kmer_index (by ID), \
         not read_sketches (by content hash)"
    );
    assert!(
        plan.get_read_sketch(hash).is_some(),
        "the read's sketch should be retrievable by its own content hash"
    );
}

/// The reference's sketch must NOT appear in `read_sketches` — that map is for reads
/// only; the reference already has its own storage (`kmer_index`/`reference_space`).
#[test]
fn plan_does_not_store_reference_sketch_in_read_sketches() {
    let dir = TempDir::new().unwrap();
    let (ref_path, reads_path, _read_bases) = write_fixtures(dir.path());
    let plan_path = dir.path().join("plan.phrayaplan");

    let out = phraya(&[
        "plan",
        "--inputs",
        reads_path.to_str().unwrap(),
        "--reference",
        ref_path.to_str().unwrap(),
        "--output",
        plan_path.to_str().unwrap(),
    ]);
    assert!(out.status.success());

    let plan = phraya_io::plan::read_plan(&plan_path).unwrap();
    let ref_seq = diverse_dna(200, 42);
    let ref_hash = phraya_io::plan::read_content_hash(&ref_seq);

    assert!(
        plan.get_read_sketch(ref_hash).is_none(),
        "the reference's own content hash should not resolve to an entry in \
         read_sketches — that map is scoped to reads, not the reference"
    );
}
