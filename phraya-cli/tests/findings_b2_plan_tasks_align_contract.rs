/// Regression test for CSP2 spike finding B2: `phraya plan-tasks` must emit
/// sequence IDs that `phraya align`'s traditional mode (`align PLAN QUERY_ID TARGET_ID`)
/// can look up directly. Before the fix, `plan-tasks` printed raw `u32` plan-internal
/// array indices while `align` resolved its `QUERY_ID`/`TARGET_ID` arguments by FASTA/
/// FASTQ header string — the documented `plan-tasks | parallel align` pipeline in
/// README.md could never work for any use case, not just contigs.
///
/// This test exercises the full CLI pipeline end-to-end (not the library API directly,
/// and not a hand-decoded index→ID mapping) — `plan-tasks`'s own stdout, unmodified,
/// is fed straight into `align`.
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

fn write_fasta(dir: &Path, name: &str, records: &[(&str, &str)]) -> PathBuf {
    let path = dir.join(name);
    let mut content = String::new();
    for (id, seq) in records {
        content.push_str(&format!(">{id}\n{seq}\n"));
    }
    std::fs::write(&path, content).unwrap();
    path
}

/// `phraya plan-tasks`'s raw stdout, piped directly into `phraya align`, succeeds for
/// every task — no manual index→ID translation.
#[test]
fn findings_b2_plan_tasks_output_is_directly_alignable() {
    let temp_dir = TempDir::new().unwrap();
    let dir = temp_dir.path();

    let reference = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    let ref_path = write_fasta(dir, "ref.fa", &[("chrom1", reference)]);
    let reads_path = write_fasta(
        dir,
        "reads.fa",
        &[
            ("readA", reference),
            ("readB", reference),
            ("readC", reference),
        ],
    );

    let plan_path = dir.join("plan.phrayaplan");
    let plan_out = phraya(&[
        "plan",
        "--inputs",
        reads_path.to_str().unwrap(),
        "--reference",
        ref_path.to_str().unwrap(),
        "--output",
        plan_path.to_str().unwrap(),
    ]);
    assert!(
        plan_out.status.success(),
        "phraya plan failed:\n{}",
        String::from_utf8_lossy(&plan_out.stderr)
    );

    let tasks_out = phraya(&["plan-tasks", plan_path.to_str().unwrap()]);
    assert!(
        tasks_out.status.success(),
        "phraya plan-tasks failed:\n{}",
        String::from_utf8_lossy(&tasks_out.stderr)
    );

    let stdout = String::from_utf8_lossy(&tasks_out.stdout).into_owned();
    let mut lines = stdout.lines();
    assert_eq!(
        lines.next().map(str::trim),
        Some("query_id\ttarget_id"),
        "plan-tasks header"
    );

    let mut task_count = 0;
    for line in lines {
        let mut parts = line.split('\t');
        let query_id = parts.next().expect("query_id column");
        let target_id = parts.next().expect("target_id column");

        // The B2 bug: these would be plan-internal indices like "0"/"1", which
        // `align` cannot resolve against the FASTA/FASTQ header IDs it reads from
        // the plan's input files.
        assert!(
            query_id == "readA" || query_id == "readB" || query_id == "readC",
            "query_id should be a real sequence ID, got: {query_id}"
        );
        assert_eq!(target_id, "chrom1", "target_id should be the reference ID");

        let out_path = dir.join(format!("{query_id}_{target_id}.phraya"));
        let align_out = phraya(&[
            "align",
            plan_path.to_str().unwrap(),
            query_id,
            target_id,
            "--output",
            out_path.to_str().unwrap(),
        ]);
        assert!(
            align_out.status.success(),
            "phraya align {query_id} {target_id} failed:\n{}",
            String::from_utf8_lossy(&align_out.stderr)
        );
        assert!(out_path.exists(), "align should write {out_path:?}");
        task_count += 1;
    }

    assert_eq!(task_count, 3, "expected one task per read against the reference");
}
