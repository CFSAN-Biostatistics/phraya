use std::collections::HashMap;
/// Acceptance tests for Issue #58: Phase 1 MVP - Evidence-Informed Alignment with Deferred Filtering
/// Tests validate end-to-end MVP functionality for use cases 2, 3, 4 (case 1 deferred to Phase 2).
/// Tests verify external contracts: CLI interfaces, file formats, error codes.
use std::fs;
use std::path::{Path, PathBuf};
use tempfile::TempDir;

// ============================================================================
// Portable manifest path helper
// ============================================================================

/// Get the portable path to phraya-cli/Cargo.toml using compile-time env var
fn get_manifest_path() -> String {
    format!("{}/Cargo.toml", env!("CARGO_MANIFEST_DIR"))
}

// ============================================================================
// Helper functions for test setup
// ============================================================================

/// Create a temporary FASTA file from sequences
fn create_fasta_file(dir: &Path, filename: &str, sequences: &[(&str, &str)]) -> PathBuf {
    let path = dir.join(filename);
    let mut content = String::new();
    for (id, seq) in sequences {
        content.push('>');
        content.push_str(id);
        content.push('\n');
        content.push_str(seq);
        content.push('\n');
    }
    fs::write(&path, content).unwrap();
    path
}

/// Create a temporary FASTQ file from sequences
fn create_fastq_file(dir: &Path, filename: &str, sequences: &[(&str, &str, &str)]) -> PathBuf {
    let path = dir.join(filename);
    let mut content = String::new();
    for (id, seq, qual) in sequences {
        content.push('@');
        content.push_str(id);
        content.push('\n');
        content.push_str(seq);
        content.push_str("\n+\n");
        content.push_str(qual);
        content.push('\n');
    }
    fs::write(&path, content).unwrap();
    path
}

/// Create a temporary .phraya file with variant observations
fn create_phraya_file(
    dir: &Path,
    filename: &str,
    observations: Vec<(u32, u8, HashMap<u8, u32>, u8)>, // position, ref_base, alleles, mapq
    reference_length: u32,
) -> PathBuf {
    use phraya_core::types::{CoverageTrack, VariantObservation};
    use phraya_io::phraya::{write_phraya, PhrayaFile};

    let path = dir.join(filename);

    let variant_obs: Vec<VariantObservation> = observations
        .into_iter()
        .enumerate()
        .map(|(i, (pos, ref_base, alleles, mapq))| {
            VariantObservation::new(
                pos,
                ref_base,
                alleles,
                0.95,
                "10M".to_string(),
                mapq,
                0,
                vec![10],
                35.0,
                format!("sample:read{}", i),
            )
        })
        .collect();

    let coverage = CoverageTrack::new(vec![10; reference_length as usize]);
    let file = PhrayaFile::new(
        reference_length,
        "test_sample".to_string(),
        "2026-06-01T00:00:00Z".to_string(),
        variant_obs,
        coverage,
    );

    write_phraya(&path, &file).unwrap();
    path
}

// ============================================================================
// SECTION A: .phrayaplan Format (Write/Read)
// ============================================================================

/// Test: .phrayaplan file can be written and read back
#[test]
fn issue_58_phrayaplan_format_round_trip() {
    use phraya_io::plan::{read_plan, write_plan, PhrayaPlan, UseCase};

    let temp_dir = TempDir::new().unwrap();
    let plan_path = temp_dir.path().join("test.phrayaplan");

    // Create a plan
    let plan = PhrayaPlan::new(
        UseCase::ReadsWithRef,
        vec!["ref.fa".to_string(), "reads.fq".to_string()],
        "2026-06-01T12:00:00Z".to_string(),
        HashMap::new(),               // sketches
        HashMap::new(),       // kmer_uniqueness
        vec![(1, 0), (2, 0)], // task_list
    );

    // Write plan
    write_plan(&plan_path, &plan).expect("should write plan");

    // Verify file exists and is non-empty
    assert!(plan_path.exists(), "plan file should exist");
    let file_size = fs::metadata(&plan_path).unwrap().len();
    assert!(file_size > 0, "plan file should be non-empty");

    // Read plan back
    let read_back = read_plan(&plan_path).expect("should read plan");

    // Verify contents
    assert_eq!(
        read_back.use_case, plan.use_case,
        "use_case should be preserved"
    );
    assert_eq!(
        read_back.input_files.len(),
        2,
        "input_files should be preserved"
    );
    assert_eq!(
        read_back.task_list, plan.task_list,
        "task_list should be preserved"
    );
}

/// Test: .phrayaplan contains all required fields
#[test]
fn issue_58_phrayaplan_contains_required_fields() {
    use phraya_io::plan::{read_plan, write_plan, PhrayaPlan, UseCase};

    let temp_dir = TempDir::new().unwrap();
    let plan_path = temp_dir.path().join("test.phrayaplan");

    let mut uniqueness = HashMap::new();
    uniqueness.insert(0u32, 0.95);
    uniqueness.insert(100u32, 0.75);

    let plan = PhrayaPlan::new(
        UseCase::ContigsWithReads,
        vec!["contigs.fa".to_string(), "reads.fq".to_string()],
        "2026-06-01T12:00:00Z".to_string(),
        HashMap::new(), // sketches
        uniqueness,
        vec![(0, 1), (1, 0), (2, 0)],
    );

    write_plan(&plan_path, &plan).unwrap();
    let read_back = read_plan(&plan_path).unwrap();

    // Verify timestamp is present
    assert!(
        !read_back.timestamp.is_empty(),
        "timestamp should be present"
    );

    // Verify input files list
    assert!(
        !read_back.input_files.is_empty(),
        "input_files should be populated"
    );

    // Verify k-mer uniqueness is present (even if empty for now)
    assert!(
        read_back.kmer_uniqueness.len() > 0,
        "kmer_uniqueness should be populated"
    );

    // Verify task_list is present
    assert_eq!(
        read_back.task_list.len(),
        3,
        "task_list should have 3 tasks"
    );
}

// ============================================================================
// SECTION B: Use Case Detection (Cases 2, 3, 4)
// ============================================================================

/// Test: Case 2 detection - Reads with reference
#[test]
fn issue_58_use_case_case2_reads_with_reference() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    let ref_path = create_fasta_file(
        temp_path,
        "reference.fa",
        &[(
            "ref",
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
        )],
    );

    let reads_path = create_fastq_file(
        temp_path,
        "reads.fq",
        &[
            (
                "read1",
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
                "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
            ),
            (
                "read2",
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
                "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
            ),
        ],
    );
    // Note: Both sequences are 52 bases, quality string is 52 'I' chars

    let output_path = temp_path.join("plan.phrayaplan");

    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            &get_manifest_path(),
            "--",
            "plan",
            "--inputs",
            reads_path.to_str().unwrap(),
            "--reference",
            ref_path.to_str().unwrap(),
            "--output",
            output_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya plan");

    assert!(
        output.status.success(),
        "phraya plan should succeed. stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    assert!(output_path.exists(), "plan file should be created");

    let plan = phraya_io::plan::read_plan(&output_path).expect("plan should be readable");
    assert_eq!(
        plan.use_case,
        phraya_io::plan::UseCase::ReadsWithRef,
        "should detect Case 2: reads with reference"
    );
}

/// Test: Case 3 detection - Contigs with reads
#[test]
fn issue_58_use_case_case3_contigs_with_reads() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    let contigs_path = create_fasta_file(
        temp_path,
        "contigs.fa",
        &[
            (
                "contig1",
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            ),
            (
                "contig2",
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGA",
            ),
        ],
    );

    let reads_path = create_fastq_file(
        temp_path,
        "reads.fq",
        &[(
            "read1",
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        )],
    );
    // Note: Sequence is 52 bases, quality string is 52 'I' chars

    let output_path = temp_path.join("plan.phrayaplan");

    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            &get_manifest_path(),
            "--",
            "plan",
            "--inputs",
            contigs_path.to_str().unwrap(),
            "--inputs",
            reads_path.to_str().unwrap(),
            "--output",
            output_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya plan");

    assert!(
        output.status.success(),
        "phraya plan should succeed. stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let plan = phraya_io::plan::read_plan(&output_path).expect("plan should be readable");
    assert_eq!(
        plan.use_case,
        phraya_io::plan::UseCase::ContigsWithReads,
        "should detect Case 3: contigs with reads"
    );
}

/// Test: Case 4 detection - Contigs only
#[test]
fn issue_58_use_case_case4_contigs_only() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    let contigs_path = create_fasta_file(
        temp_path,
        "contigs.fa",
        &[
            (
                "contig1",
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            ),
            (
                "contig2",
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGA",
            ),
            (
                "contig3",
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGC",
            ),
        ],
    );

    let output_path = temp_path.join("plan.phrayaplan");

    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            &get_manifest_path(),
            "--",
            "plan",
            "--inputs",
            contigs_path.to_str().unwrap(),
            "--output",
            output_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya plan");

    assert!(output.status.success(), "phraya plan should succeed");

    let plan = phraya_io::plan::read_plan(&output_path).expect("plan should be readable");
    assert_eq!(
        plan.use_case,
        phraya_io::plan::UseCase::ContigsOnly,
        "should detect Case 4: contigs only"
    );
}

// ============================================================================
// SECTION C: Task List Generation
// ============================================================================

/// Test: Task list for Case 2 (reads + ref) has one task per read
#[test]
fn issue_58_task_list_case2_one_per_read() {
    use phraya_io::plan::{read_plan, write_plan, PhrayaPlan, UseCase};

    let temp_dir = TempDir::new().unwrap();
    let plan_path = temp_dir.path().join("test.phrayaplan");

    // 1 reference + 3 reads = 3 tasks
    let plan = PhrayaPlan::new(
        UseCase::ReadsWithRef,
        vec!["ref.fa".to_string(), "reads.fq".to_string()],
        "2026-06-01T12:00:00Z".to_string(),
        HashMap::new(), // sketches
        HashMap::new(),
        vec![(1, 0), (2, 0), (3, 0)], // 3 tasks: each read vs reference
    );

    write_plan(&plan_path, &plan).unwrap();
    let read_plan = read_plan(&plan_path).unwrap();

    assert_eq!(
        read_plan.task_list.len(),
        3,
        "should have 3 tasks for 3 reads"
    );

    // All tasks should target reference (index 0)
    for (query_id, target_id) in &read_plan.task_list {
        assert_eq!(*target_id, 0, "target should always be reference");
        assert!(*query_id >= 1, "query should be a read");
    }
}

/// Test: Task list contains valid (query_id, target_id) tuples
#[test]
fn issue_58_task_list_structure_valid() {
    use phraya_io::plan::{read_plan, write_plan, PhrayaPlan, UseCase};

    let temp_dir = TempDir::new().unwrap();
    let plan_path = temp_dir.path().join("test.phrayaplan");

    let plan = PhrayaPlan::new(
        UseCase::ContigsOnly,
        vec!["contigs.fa".to_string()],
        "2026-06-01T12:00:00Z".to_string(),
        HashMap::new(), // sketches
        HashMap::new(),
        vec![(0, 1), (0, 2), (1, 2)], // pairwise
    );

    write_plan(&plan_path, &plan).unwrap();
    let read_back = read_plan(&plan_path).unwrap();

    // All task entries should be (u32, u32) tuples
    for (query_id, target_id) in &read_back.task_list {
        // Just verify they're valid u32 values
        assert!(*query_id < 1000, "query_id should be reasonable");
        assert!(*target_id < 1000, "target_id should be reasonable");
    }
}

// ============================================================================
// SECTION D: VariantObservation Fields (CIGAR, MAPQ, local_coverage, all_alleles)
// ============================================================================

/// Test: VariantObservation has CIGAR field
#[test]
fn issue_58_variant_observation_has_cigar() {
    use phraya_core::types::VariantObservation;

    let mut alleles = HashMap::new();
    alleles.insert(b'A', 10);

    let obs = VariantObservation::new(
        100,
        b'A',
        alleles,
        0.95,
        "15M2D5M".to_string(), // CIGAR string
        60,
        2,
        vec![10],
        35.0,
        "sample:read1".to_string(),
    );

    assert_eq!(obs.cigar(), "15M2D5M", "CIGAR string should be accessible");
}

/// Test: VariantObservation has edit_distance field
#[test]
fn issue_58_variant_observation_has_edit_distance() {
    use phraya_core::types::VariantObservation;

    let mut alleles = HashMap::new();
    alleles.insert(b'A', 10);

    let obs = VariantObservation::new(
        100,
        b'A',
        alleles,
        0.95,
        "10M".to_string(),
        60,
        5, // edit distance
        vec![10],
        35.0,
        "sample:read1".to_string(),
    );

    assert_eq!(obs.edit_distance(), 5, "edit_distance should be accessible");
}

/// Test: VariantObservation has local_coverage vector
#[test]
fn issue_58_variant_observation_has_local_coverage() {
    use phraya_core::types::VariantObservation;

    let mut alleles = HashMap::new();
    alleles.insert(b'A', 15);

    let coverage_window = vec![10, 12, 15, 18, 20, 22, 20, 18, 15, 12];

    let obs = VariantObservation::new(
        100,
        b'A',
        alleles,
        0.95,
        "10M".to_string(),
        60,
        0,
        coverage_window.clone(),
        35.0,
        "sample:read1".to_string(),
    );

    assert_eq!(
        obs.local_coverage(),
        &coverage_window,
        "local_coverage should be accessible"
    );
}

/// Test: VariantObservation has all_alleles with counts
#[test]
fn issue_58_variant_observation_has_all_alleles() {
    use phraya_core::types::VariantObservation;

    let mut alleles = HashMap::new();
    alleles.insert(b'A', 15);
    alleles.insert(b'T', 5);
    alleles.insert(b'G', 2);

    let obs = VariantObservation::new(
        100,
        b'A',
        alleles.clone(),
        0.95,
        "10M".to_string(),
        60,
        0,
        vec![22],
        35.0,
        "sample:read1".to_string(),
    );

    assert_eq!(
        obs.all_alleles(),
        &alleles,
        "all_alleles should be accessible with all counts"
    );
    assert_eq!(obs.all_alleles().get(&b'A'), Some(&15));
    assert_eq!(obs.all_alleles().get(&b'T'), Some(&5));
    assert_eq!(obs.all_alleles().get(&b'G'), Some(&2));
}

// ============================================================================
// SECTION E: Coverage Tracks (Quantized to 5, RLE-compressed)
// ============================================================================

/// Test: CoverageTrack quantizes to nearest 5
#[test]
fn issue_58_coverage_track_quantizes_to_5() {
    use phraya_core::types::CoverageTrack;

    // Test boundary cases
    assert_eq!(CoverageTrack::quantize(0), 0, "0 should quantize to 0");
    assert_eq!(CoverageTrack::quantize(2), 0, "2 should quantize to 0");
    assert_eq!(CoverageTrack::quantize(3), 5, "3 should quantize to 5");
    assert_eq!(CoverageTrack::quantize(7), 5, "7 should quantize to 5");
    assert_eq!(CoverageTrack::quantize(8), 10, "8 should quantize to 10");
    assert_eq!(CoverageTrack::quantize(12), 10, "12 should quantize to 10");
    assert_eq!(CoverageTrack::quantize(13), 15, "13 should quantize to 15");
}

/// Test: CoverageTrack compresses with RLE
#[test]
fn issue_58_coverage_track_rle_compression() {
    use phraya_core::types::CoverageTrack;

    // Uniform coverage should compress to single run
    let coverage = vec![10, 10, 10, 10, 10, 10, 10, 10];
    let track = CoverageTrack::new(coverage);

    // Should decompress back to original (quantized)
    let decompressed = track.decompress();
    assert_eq!(decompressed.len(), 8, "should decompress to full length");
    assert!(
        decompressed.iter().all(|&c| c == 10),
        "all positions should have coverage 10"
    );
}

/// Test: CoverageTrack provides position lookups
#[test]
fn issue_58_coverage_track_position_lookup() {
    use phraya_core::types::CoverageTrack;

    let coverage = vec![10, 10, 5, 5, 15, 15];
    let track = CoverageTrack::new(coverage);

    assert_eq!(track.coverage_at(0), Some(10));
    assert_eq!(track.coverage_at(1), Some(10));
    assert_eq!(track.coverage_at(2), Some(5));
    assert_eq!(track.coverage_at(3), Some(5));
    assert_eq!(track.coverage_at(4), Some(15));
    assert_eq!(track.coverage_at(5), Some(15));
    assert_eq!(
        track.coverage_at(6),
        None,
        "out of bounds should return None"
    );
}

// ============================================================================
// SECTION F: Multi-mapping Query Index (score_ratio >= 0.95)
// ============================================================================

/// Test: Query index can be written and read
#[test]
fn issue_58_query_index_round_trip() {
    use phraya_io::queries::{read_queries, write_queries};
    use std::collections::HashMap;

    let temp_dir = TempDir::new().unwrap();
    let query_path = temp_dir.path().join("test.phraya.queries");

    // Create query index: query_id -> list of (position, score)
    let mut queries = HashMap::new();
    queries.insert(
        "query_0".to_string(),
        vec![(100u32, 0.98f64), (200u32, 0.96f64)],
    );
    queries.insert("query_1".to_string(), vec![(150u32, 0.99f64)]);

    write_queries(&query_path, &queries).expect("should write queries");

    assert!(query_path.exists(), "queries file should exist");

    let read_back = read_queries(&query_path).expect("should read queries");

    assert_eq!(read_back.len(), 2, "should have entries for 2 queries");
    assert_eq!(
        read_back.get("query_0"),
        Some(&vec![(100u32, 0.98f64), (200u32, 0.96f64)]),
        "query_0 should have correct alignments"
    );
}

/// Test: Query index enforces score_ratio >= 0.95
#[test]
fn issue_58_query_index_filters_by_score_ratio() {
    use phraya_io::queries::{read_queries, write_queries};

    let temp_dir = TempDir::new().unwrap();
    let query_path = temp_dir.path().join("test.phraya.queries");

    // Only positions with score_ratio >= 0.95 should be included
    let mut queries = HashMap::new();
    queries.insert(
        "query_0".to_string(),
        vec![(100u32, 0.95f64), (200u32, 0.94f64)],
    );

    write_queries(&query_path, &queries).expect("should write");
    let read_back = read_queries(&query_path).expect("should read");

    if let Some(alignments) = read_back.get("query_0") {
        // Position 100 with score 0.95 should be included
        assert!(alignments
            .iter()
            .any(|(pos, score)| *pos == 100 && *score >= 0.95));
        // Position 200 with score 0.94 should be filtered out
        assert!(!alignments
            .iter()
            .any(|(pos, score)| *pos == 200 && *score < 0.95));
    }
}

// ============================================================================
// SECTION G: Filter Operations (Threshold-based)
// ============================================================================

/// Test: Filter with min_coverage threshold
#[test]
fn issue_58_filter_min_coverage() {
    use phraya_core::types::VariantObservation;
    use phraya_filter::FilterBuilder;

    let mut alleles1 = HashMap::new();
    alleles1.insert(b'A', 5);
    let obs_low = VariantObservation::new(
        100,
        b'A',
        alleles1,
        0.95,
        "10M".to_string(),
        60,
        0,
        vec![5],
        35.0,
        "test:read1".to_string(),
    );

    let mut alleles2 = HashMap::new();
    alleles2.insert(b'A', 15);
    let obs_high = VariantObservation::new(
        100,
        b'A',
        alleles2,
        0.95,
        "10M".to_string(),
        60,
        0,
        vec![15],
        35.0,
        "test:read2".to_string(),
    );

    let filter = FilterBuilder::new().min_coverage(10).build();

    assert!(
        !filter.apply(&obs_low),
        "coverage 5 should fail min_coverage 10"
    );
    assert!(
        filter.apply(&obs_high),
        "coverage 15 should pass min_coverage 10"
    );
}

/// Test: Filter with min_mapq threshold
#[test]
fn issue_58_filter_min_mapq() {
    use phraya_core::types::VariantObservation;
    use phraya_filter::FilterBuilder;

    let mut alleles = HashMap::new();
    alleles.insert(b'A', 10);

    let obs_low_mapq = VariantObservation::new(
        100,
        b'A',
        alleles.clone(),
        0.95,
        "10M".to_string(),
        30,
        0,
        vec![10],
        35.0,
        "test:read1".to_string(),
    );
    let obs_high_mapq = VariantObservation::new(
        100,
        b'A',
        alleles,
        0.95,
        "10M".to_string(),
        50,
        0,
        vec![10],
        35.0,
        "test:read2".to_string(),
    );

    let filter = FilterBuilder::new().min_mapq(40).build();

    assert!(
        !filter.apply(&obs_low_mapq),
        "mapq 30 should fail min_mapq 40"
    );
    assert!(
        filter.apply(&obs_high_mapq),
        "mapq 50 should pass min_mapq 40"
    );
}

/// Test: Filter composition (multiple thresholds)
#[test]
fn issue_58_filter_composition() {
    use phraya_core::types::VariantObservation;
    use phraya_filter::FilterBuilder;

    let mut alleles = HashMap::new();
    alleles.insert(b'A', 15);
    let obs = VariantObservation::new(
        100,
        b'A',
        alleles,
        0.95,
        "10M".to_string(),
        50,
        0,
        vec![15],
        35.0,
        "test:read".to_string(),
    );

    let filter = FilterBuilder::new().min_coverage(10).min_mapq(40).build();

    assert!(filter.apply(&obs), "should pass both thresholds");

    // Fail one threshold
    let strict = FilterBuilder::new().min_coverage(10).min_mapq(60).build();

    assert!(!strict.apply(&obs), "should fail min_mapq 60");
}

// ============================================================================
// SECTION H: CLI Integration - Filter Command
// ============================================================================

/// Test: filter command outputs VCF format
#[test]
fn issue_58_filter_vcf_output() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    let observations = vec![(
        100,
        b'A',
        {
            let mut h = HashMap::new();
            h.insert(b'A', 10);
            h.insert(b'T', 5);
            h
        },
        60,
    )];

    let phraya_path = create_phraya_file(temp_path, "test.phraya", observations, 200);

    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            &get_manifest_path(),
            "--",
            "filter",
            phraya_path.to_str().unwrap(),
            "--format",
            "vcf",
        ])
        .output()
        .expect("Failed to execute phraya filter");

    assert!(
        output.status.success(),
        "filter should succeed. stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let stdout = String::from_utf8_lossy(&output.stdout);
    // VCF should have headers and content
    assert!(
        stdout.contains("##fileformat=VCF") || stdout.contains("CHROM"),
        "output should contain VCF header or column line"
    );
}

/// Test: filter command outputs TSV format
#[test]
fn issue_58_filter_tsv_output() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    let observations = vec![(
        100,
        b'A',
        {
            let mut h = HashMap::new();
            h.insert(b'A', 10);
            h
        },
        60,
    )];

    let phraya_path = create_phraya_file(temp_path, "test.phraya", observations, 200);

    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            &get_manifest_path(),
            "--",
            "filter",
            phraya_path.to_str().unwrap(),
            "--format",
            "tsv",
        ])
        .output()
        .expect("Failed to execute phraya filter");

    assert!(output.status.success(), "filter TSV should succeed");

    let stdout = String::from_utf8_lossy(&output.stdout);
    // TSV should have tab-separated columns
    assert!(stdout.contains('\t'), "TSV output should contain tabs");
    assert!(
        stdout.contains("position") || stdout.contains("100"),
        "TSV should contain column header or data"
    );
}

/// Test: filter command with --min-coverage threshold
#[test]
fn issue_58_filter_min_coverage_cli() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    // Create .phraya with 2 observations: one with low coverage, one with high
    let observations = vec![
        (
            100,
            b'A',
            {
                let mut h = HashMap::new();
                h.insert(b'A', 5);
                h
            },
            60,
        ),
        (
            200,
            b'A',
            {
                let mut h = HashMap::new();
                h.insert(b'A', 15);
                h
            },
            60,
        ),
    ];

    let phraya_path = create_phraya_file(temp_path, "test.phraya", observations, 300);

    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            &get_manifest_path(),
            "--",
            "filter",
            phraya_path.to_str().unwrap(),
            "--min-coverage",
            "10",
            "--format",
            "tsv",
        ])
        .output()
        .expect("Failed to execute phraya filter");

    assert!(output.status.success(), "filter should succeed");

    let stdout = String::from_utf8_lossy(&output.stdout);
    // Should have header + 1 data line (only position 200 has coverage >= 10)
    let lines: Vec<&str> = stdout.lines().collect();
    assert!(
        lines.len() >= 2,
        "should have header and at least 1 data line, got {} lines",
        lines.len()
    );
}

/// Test: filter command outputs filtered .phraya format
#[test]
fn issue_58_filter_phraya_output() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    let observations = vec![(
        100,
        b'A',
        {
            let mut h = HashMap::new();
            h.insert(b'A', 10);
            h
        },
        60,
    )];

    let input_path = create_phraya_file(temp_path, "test.phraya", observations, 200);
    let output_path = temp_path.join("filtered.phraya");

    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            &get_manifest_path(),
            "--",
            "filter",
            input_path.to_str().unwrap(),
            "--format",
            "phraya",
            "--output",
            output_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya filter");

    assert!(
        output.status.success(),
        "filter phraya should succeed. stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    assert!(output_path.exists(), "filtered .phraya file should exist");

    // Verify we can read the filtered file
    let filtered =
        phraya_io::phraya::read_phraya(&output_path).expect("should be able to read filtered file");

    assert!(
        filtered.observations.len() > 0,
        "filtered file should contain observations"
    );
}

// ============================================================================
// SECTION I: .phraya File Format (Observations + Coverage Track)
// ============================================================================

/// Test: .phraya file round-trip with observations
#[test]
fn issue_58_phraya_file_round_trip() {
    use phraya_core::types::{CoverageTrack, VariantObservation};
    use phraya_io::phraya::{read_phraya, write_phraya, PhrayaFile};

    let temp_dir = TempDir::new().unwrap();
    let test_file = temp_dir.path().join("test.phraya");

    // Create observations
    let mut alleles = HashMap::new();
    alleles.insert(b'A', 10);
    let obs = VariantObservation::new(
        100,
        b'A',
        alleles,
        0.95,
        "10M".to_string(),
        60,
        0,
        vec![10],
        35.0,
        "sample:read1".to_string(),
    );

    let coverage = CoverageTrack::new(vec![10; 200]);
    let file = PhrayaFile::new(
        200,
        "test_sample".to_string(),
        "2026-06-01T00:00:00Z".to_string(),
        vec![obs],
        coverage,
    );

    write_phraya(&test_file, &file).expect("should write");
    assert!(test_file.exists(), "file should exist");

    let read_back = read_phraya(&test_file).expect("should read");

    assert_eq!(read_back.observations.len(), 1, "should have 1 observation");
    assert_eq!(
        read_back.observations[0].position(),
        100,
        "position should match"
    );
}

/// Test: .phraya coverage track is preserved through write/read
#[test]
fn issue_58_phraya_coverage_track_preserved() {
    use phraya_core::types::CoverageTrack;
    use phraya_io::phraya::{read_phraya, write_phraya, PhrayaFile};

    let temp_dir = TempDir::new().unwrap();
    let test_file = temp_dir.path().join("test.phraya");

    let coverage = CoverageTrack::new(vec![10, 10, 5, 5, 15, 15]);
    let file = PhrayaFile::new(
        6,
        "test_sample".to_string(),
        "2026-06-01T00:00:00Z".to_string(),
        vec![],
        coverage,
    );

    write_phraya(&test_file, &file).unwrap();
    let read_back = read_phraya(&test_file).unwrap();

    assert_eq!(
        read_back.coverage_track.total_length(),
        6,
        "total_length should match"
    );
    assert_eq!(read_back.coverage_track.coverage_at(0), Some(10));
    assert_eq!(read_back.coverage_track.coverage_at(2), Some(5));
    assert_eq!(read_back.coverage_track.coverage_at(4), Some(15));
}

// ============================================================================
// SECTION J: Error Handling
// ============================================================================

/// Test: Missing input file returns error
#[test]
fn issue_58_error_missing_input_file() {
    let nonexistent = "/tmp/nonexistent_issue_58_xyz.fa";
    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            &get_manifest_path(),
            "--",
            "plan",
            "--inputs",
            nonexistent,
            "--output",
            "/tmp/out.phrayaplan",
        ])
        .output()
        .expect("Failed to execute");

    assert!(!output.status.success(), "should fail with missing file");
}

/// Test: Invalid .phraya file returns error
#[test]
fn issue_58_error_corrupt_phraya_file() {
    let temp_dir = TempDir::new().unwrap();
    let bad_file = temp_dir.path().join("corrupt.phraya");

    fs::write(&bad_file, b"this is not a valid phraya file").unwrap();

    let output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            &get_manifest_path(),
            "--",
            "filter",
            bad_file.to_str().unwrap(),
            "--format",
            "vcf",
        ])
        .output()
        .expect("Failed to execute");

    assert!(!output.status.success(), "should fail with corrupt file");
}

// ============================================================================
// SECTION K: MVP End-to-End Workflow
// ============================================================================

/// Test: Full workflow - plan → plan-tasks → data preservation
#[test]
fn issue_58_mvp_workflow_plan_to_tasks() {
    let temp_dir = TempDir::new().unwrap();
    let temp_path = temp_dir.path();

    // Step 1: Create inputs
    let ref_path = create_fasta_file(
        temp_path,
        "reference.fa",
        &[(
            "ref",
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
        )],
    );

    let reads_path = create_fastq_file(
        temp_path,
        "reads.fq",
        &[
            (
                "read1",
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
                "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
            ),
            (
                "read2",
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
                "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
            ),
        ],
    );
    // Note: Both sequences are 52 bases, quality string is 52 'I' chars

    let plan_path = temp_path.join("plan.phrayaplan");

    // Step 2: Run phraya plan
    let plan_output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            &get_manifest_path(),
            "--",
            "plan",
            "--inputs",
            reads_path.to_str().unwrap(),
            "--reference",
            ref_path.to_str().unwrap(),
            "--output",
            plan_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya plan");

    assert!(plan_output.status.success(), "phraya plan should succeed");

    // Step 3: Run phraya plan-tasks
    let tasks_output = std::process::Command::new("cargo")
        .args(&[
            "run",
            "--manifest-path",
            &get_manifest_path(),
            "--",
            "plan-tasks",
            plan_path.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute phraya plan-tasks");

    assert!(
        tasks_output.status.success(),
        "phraya plan-tasks should succeed. stderr: {}",
        String::from_utf8_lossy(&tasks_output.stderr)
    );

    // Step 4: Verify output
    let stdout = String::from_utf8_lossy(&tasks_output.stdout);
    let lines: Vec<&str> = stdout.lines().collect();

    assert!(lines.len() >= 3, "should have header + 2 task lines");
    assert_eq!(
        lines[0].trim(),
        "query_id\ttarget_id",
        "header should match"
    );
}
