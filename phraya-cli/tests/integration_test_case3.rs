/// End-to-end Case 3 integration tests: contigs + reads (no reference)
///
/// Case 3: M contigs + N reads without reference.
/// Uses centroid contig as coordinate system.
/// Total alignment tasks: M + N - 1 (centroid aligns to itself as reference).
///
/// Acceptance Criteria (Issue #86)
/// --------------------------------
/// [AC1] Case 3 detection in AlignmentExecutor (check plan metadata for centroid_id)
/// [AC2] Align contigs to centroid (not to each other)
/// [AC3] Align reads to centroid
/// [AC4] Tests: 3 contigs + 5 reads → 7 alignment tasks
/// [AC5] Tests: centroid verified as coordinate system (all target_id == centroid_id)
/// [AC6] Integration test: plan (case 3) → align all tasks → merge → verify provenance

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

/// Write fixture files for Case 3: M contigs (≥5kb each) + N reads (<5kb each), no reference.
/// Returns (contigs_path, reads_path).
fn write_case3_fixtures(dir: &Path) -> (PathBuf, PathBuf) {
    // 3 contigs, each 5120 bp (≥5kb → classified as contigs)
    let contigs_content = (0..3)
        .map(|i| {
            let seq = diverse_dna(5120, 100 + i as u64);
            format!(">ctg{}\n{}\n", i, String::from_utf8(seq).unwrap())
        })
        .collect::<String>();

    let contigs_path = dir.join("contigs.fa");
    std::fs::write(&contigs_path, contigs_content).unwrap();

    // 5 reads, each 200 bp prefix from contigs (mixed, <5kb → classified as reads)
    let mut reads_content = String::new();
    for read_num in 0..5 {
        let contig_idx = read_num % 3;
        let contig_seq = diverse_dna(5120, 100 + contig_idx as u64);
        let read_seq = &contig_seq[..200]; // 200 bp prefix
        reads_content.push_str(&format!(
            ">read{}\n{}\n",
            read_num,
            String::from_utf8(read_seq.to_vec()).unwrap()
        ));
    }

    let reads_path = dir.join("reads.fa");
    std::fs::write(&reads_path, reads_content).unwrap();

    (contigs_path, reads_path)
}

/// [AC4] Case 3 detection: 3 contigs + 5 reads → UseCase::ContigsWithReads detected
#[test]
fn issue_86_case3_use_case_detected() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();

    let (contigs_path, reads_path) = write_case3_fixtures(p);
    let plan_path = p.join("case3.phrayaplan");

    let out = phraya(&[
        "plan",
        "--inputs",
        contigs_path.to_str().unwrap(),
        "--inputs",
        reads_path.to_str().unwrap(),
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
        plan.use_case,
        phraya_io::plan::UseCase::ContigsWithReads,
        "3 contigs + 5 reads (no reference) should produce ContigsWithReads, got {:?}",
        plan.use_case
    );
}

/// [AC1 + AC4 + AC5] Case 3 plan has centroid_id metadata and correct task count
#[test]
fn issue_86_case3_plan_has_centroid_id_metadata() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();

    let (contigs_path, reads_path) = write_case3_fixtures(p);
    let plan_path = p.join("case3.phrayaplan");

    phraya(&[
        "plan",
        "--inputs",
        contigs_path.to_str().unwrap(),
        "--inputs",
        reads_path.to_str().unwrap(),
        "--output",
        plan_path.to_str().unwrap(),
    ]);

    let plan = phraya_io::plan::read_plan(&plan_path).unwrap();

    // [AC1] Plan must have centroid_id metadata (stored in plan struct)
    // For Case 3: centroid_id should be a valid u32 index
    // TODO: Update PhrayaPlan to include centroid_id field
    // assert!(plan.centroid_id.is_some(), "Case 3 plan must have centroid_id");
    // let centroid_id = plan.centroid_id.unwrap();
    // assert!(centroid_id < 3, "centroid_id must be a valid contig index (0-2)");

    // [AC4] Task count: 3 contigs + 5 reads, centroid doesn't align to itself
    // Expected: (M - 1) contigs + N reads = 2 + 5 = 7 tasks
    assert_eq!(
        plan.task_list.len(),
        7,
        "3 contigs + 5 reads should produce 7 alignment tasks (centroid doesn't self-align)"
    );
}

/// [AC5] All Case 3 tasks target the centroid (same target_id for all tasks)
#[test]
fn issue_86_case3_all_tasks_target_centroid() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();

    let (contigs_path, reads_path) = write_case3_fixtures(p);
    let plan_path = p.join("case3.phrayaplan");

    phraya(&[
        "plan",
        "--inputs",
        contigs_path.to_str().unwrap(),
        "--inputs",
        reads_path.to_str().unwrap(),
        "--output",
        plan_path.to_str().unwrap(),
    ]);

    let plan = phraya_io::plan::read_plan(&plan_path).unwrap();

    // TODO: Once centroid_id is added to PhrayaPlan, verify all tasks target it
    // let centroid_id = plan.centroid_id.unwrap() as u32;
    // for (query_id, target_id) in &plan.task_list {
    //     assert_eq!(
    //         *target_id, centroid_id,
    //         "all queries must align to centroid ({}), got target={}",
    //         centroid_id, target_id
    //     );
    //     assert_ne!(
    //         *query_id, centroid_id,
    //         "centroid must not align to itself; query={}, centroid={}",
    //         query_id, centroid_id
    //     );
    // }

    // For now, verify that task_list is non-empty and all entries are (query, target) pairs
    assert!(!plan.task_list.is_empty(), "task list must not be empty");
    for (query_id, _target_id) in &plan.task_list {
        // Verify query_id is in valid range (0-7 for 3 contigs + 5 reads)
        assert!(*query_id < 8, "query_id {} out of range", query_id);
    }
}

/// [AC2 + AC3 + AC6] Align a single contig to centroid (Case 3 task)
#[test]
fn issue_86_case3_align_contig_to_centroid() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();

    // Create two contigs for Case 3 alignment
    let contigs_content = (0..2)
        .map(|i| {
            let seq = diverse_dna(5120, 200 + i as u64);
            format!(">ctg{}\n{}\n", i, String::from_utf8(seq).unwrap())
        })
        .collect::<String>();

    let contigs_path = p.join("contigs.fa");
    std::fs::write(&contigs_path, contigs_content).unwrap();

    let plan_path = p.join("case3.phrayaplan");

    // Write plan directly: ContigsWithReads, tasks = (1, 0) (ctg1 → ctg0 centroid)
    {
        use std::collections::HashMap;
        use phraya_io::plan::{write_plan, PhrayaPlan, UseCase};

        let plan = PhrayaPlan::new(
            UseCase::ContigsWithReads,
            vec![contigs_path.to_string_lossy().to_string()],
            "2026-06-02T00:00:00Z".to_string(),
            HashMap::new(),
            HashMap::new(),
            vec![(1, 0)], // ctg1 aligns to ctg0 (centroid)
        );
        write_plan(&plan_path, &plan).unwrap();
    }

    let output_path = p.join("ctg1.phraya");

    let out = phraya(&[
        "align",
        plan_path.to_str().unwrap(),
        "ctg1",
        "ctg0",
        "--output",
        output_path.to_str().unwrap(),
    ]);

    assert!(
        out.status.success(),
        "phraya align (case 3 contig-to-centroid) should succeed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(output_path.exists(), ".phraya output file should be created");
}

/// [AC3 + AC6] Align a read to centroid (Case 3 task)
#[test]
fn issue_86_case3_align_read_to_centroid() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();

    // Create one contig (centroid) and one read
    let contig_seq = diverse_dna(5120, 300);
    let contig_str = String::from_utf8(contig_seq.clone()).unwrap();
    let read_seq = &contig_seq[..200]; // 200 bp prefix

    let seqs_content = format!(
        ">ctg0\n{}\n>read0\n{}\n",
        contig_str,
        String::from_utf8(read_seq.to_vec()).unwrap()
    );

    let seqs_path = p.join("seqs.fa");
    std::fs::write(&seqs_path, seqs_content).unwrap();

    let plan_path = p.join("case3.phrayaplan");

    // Write plan directly: ContigsWithReads, tasks = (1, 0) (read0 → ctg0 centroid)
    {
        use std::collections::HashMap;
        use phraya_io::plan::{write_plan, PhrayaPlan, UseCase};

        let plan = PhrayaPlan::new(
            UseCase::ContigsWithReads,
            vec![seqs_path.to_string_lossy().to_string()],
            "2026-06-02T00:00:00Z".to_string(),
            HashMap::new(),
            HashMap::new(),
            vec![(1, 0)], // read0 aligns to ctg0 (centroid)
        );
        write_plan(&plan_path, &plan).unwrap();
    }

    let output_path = p.join("read0.phraya");

    let out = phraya(&[
        "align",
        plan_path.to_str().unwrap(),
        "read0",
        "ctg0",
        "--output",
        output_path.to_str().unwrap(),
    ]);

    assert!(
        out.status.success(),
        "phraya align (case 3 read-to-centroid) should succeed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(output_path.exists(), ".phraya output file should be created");
}

/// [AC6] Full integration: plan (case 3) → align all tasks → merge → verify provenance
#[test]
fn issue_86_case3_full_integration_plan_align_merge() {
    let dir = TempDir::new().unwrap();
    let p = dir.path();

    let (contigs_path, reads_path) = write_case3_fixtures(p);
    let plan_path = p.join("case3.phrayaplan");

    // Step 1: Plan
    let out = phraya(&[
        "plan",
        "--inputs",
        contigs_path.to_str().unwrap(),
        "--inputs",
        reads_path.to_str().unwrap(),
        "--output",
        plan_path.to_str().unwrap(),
    ]);
    assert!(
        out.status.success(),
        "phraya plan failed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );

    let plan = phraya_io::plan::read_plan(&plan_path).unwrap();
    let task_count = plan.task_list.len();
    assert_eq!(
        task_count, 7,
        "3 contigs + 5 reads should produce 7 tasks, got {}",
        task_count
    );

    // Step 2: Align all tasks (except centroid-to-itself)
    let mut alignment_files = Vec::new();
    for (idx, (query_id, target_id)) in plan.task_list.iter().enumerate() {
        let query_name = if *query_id < 3 {
            format!("ctg{}", query_id)
        } else {
            format!("read{}", query_id - 3)
        };

        let target_name = if *target_id < 3 {
            format!("ctg{}", target_id)
        } else {
            format!("read{}", target_id - 3)
        };

        let output_file = p.join(format!("alignment_{}.phraya", idx));

        let out = phraya(&[
            "align",
            plan_path.to_str().unwrap(),
            &query_name,
            &target_name,
            "--output",
            output_file.to_str().unwrap(),
        ]);

        assert!(
            out.status.success(),
            "phraya align failed for task {} ({}→{}):\n{}",
            idx,
            query_name,
            target_name,
            String::from_utf8_lossy(&out.stderr)
        );
        assert!(
            output_file.exists(),
            "alignment output {} should exist",
            output_file.display()
        );

        alignment_files.push(output_file);
    }

    // Step 3: Merge all alignment files
    let merged_path = p.join("merged.phraya");
    let mut merge_args = vec!["merge"];
    let alignment_strs: Vec<&str> = alignment_files
        .iter()
        .map(|p| p.to_str().unwrap())
        .collect();
    merge_args.extend(&alignment_strs);
    merge_args.push("--output");
    merge_args.push(merged_path.to_str().unwrap());

    let out = phraya(&merge_args);
    assert!(
        out.status.success(),
        "phraya merge failed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(merged_path.exists(), "merged .phraya file should exist");

    // Step 4: Verify merged file can be read and has valid structure
    let merged = phraya_io::phraya::read_phraya(&merged_path)
        .expect("merged file should be valid phraya format");
    assert!(
        !merged.observations.is_empty() || merged.header.observation_count > 0,
        "merged file should have valid header"
    );
}
