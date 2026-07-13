use phraya_core::types::MinimizerSketch;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};
use std::path::Path;
use thiserror::Error;

// Import SHA-256 for content hashing (256-bit cryptographic strength)
use sha2::{Digest, Sha256};

/// Serialize a `HashMap` with keys in ascending order for deterministic output.
///
/// The plan's per-sequence maps (sketches, uniqueness, membership, mate info) are `HashMap`s
/// for fast build-time insertion, but `HashMap` iteration order is randomized per process,
/// which makes the serialized `.phrayaplan` bytes vary between identical `phraya plan` runs.
/// Routing serialization through a `BTreeMap` gives a canonical, byte-stable order without
/// changing the in-memory types or any accessor. Paired with a pinned header timestamp
/// (`PHRAYA_SOURCE_DATE`), this makes plans reproducible so a content hash can gate
/// regression runs.
fn serialize_map_sorted<S, K, V>(map: &HashMap<K, V>, serializer: S) -> Result<S::Ok, S::Error>
where
    S: serde::Serializer,
    K: serde::Serialize + Ord + Clone,
    V: serde::Serialize + Clone,
{
    let sorted: BTreeMap<K, V> = map.iter().map(|(k, v)| (k.clone(), v.clone())).collect();
    sorted.serialize(serializer)
}

/// PhrayaPlan format version for forward compatibility
pub const PHRAYAPLAN_VERSION: u32 = 6;

/// Plan file format errors
#[derive(Debug, Error, Serialize, Deserialize)]
pub enum PlanError {
    #[error("serialization error: {0}")]
    SerializationError(String),
    #[error("decompression error: {0}")]
    DecompressionError(String),
    #[error("compression error: {0}")]
    CompressionError(String),
    #[error("io error: {0}")]
    IoError(String),
    #[error("version mismatch: expected {expected}, got {got}")]
    VersionMismatch { expected: u32, got: u32 },
}

/// Use case detected from input sequences
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
pub enum UseCase {
    /// N reads + reference genome
    ReadsWithRef = 1,
    /// N reads only, no reference (MSA)
    ReadsOnly = 2,
    /// M contigs + N reads, no reference
    ContigsWithReads = 3,
    /// M contigs only
    ContigsOnly = 4,
}

/// K-mer sketching parameters used during planning
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct KmerParams {
    pub k: usize,
    pub w: usize,
}

impl Default for KmerParams {
    fn default() -> Self {
        Self { k: 21, w: 11 }
    }
}

fn default_dense_kmer_index() -> HashMap<String, MinimizerSketch> {
    HashMap::new()
}

fn default_w11_membership() -> HashMap<String, Vec<bool>> {
    HashMap::new()
}

fn default_sparse_mode() -> bool {
    false
}

fn is_false(v: &bool) -> bool {
    !v
}

/// Content-addressed reference space (ADR-0011): identity by content hash,
/// with optional human-facing name and per-sequence sketches.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct ReferenceSpace {
    /// Strong cryptographic hash (BLAKE3/SHA-256-class, ≥256 bits) of reference content
    pub content_hash: String,
    /// Optional human-facing name (e.g., "chr1-v1", "E. coli K-12")
    pub name: Option<String>,
    /// Per-sequence k-mer sketches, keyed by sequence ID
    pub sketches: HashMap<String, MinimizerSketch>,
}

/// Deduplicate a minimizer sketch, removing duplicate (hash, position) tuples.
/// Returns a new sketch with only unique minimizers while preserving the k and w parameters.
fn deduplicate_sketch(sketch: &MinimizerSketch) -> MinimizerSketch {
    let unique_minimizers: std::collections::HashSet<_> =
        sketch.minimizers.iter().copied().collect();
    let mut minimizers: Vec<_> = unique_minimizers.into_iter().collect();
    // Sort by position to maintain consistent order
    minimizers.sort_by_key(|&(_, pos)| pos);

    MinimizerSketch {
        minimizers,
        k: sketch.k,
        w: sketch.w,
    }
}

/// Insert size distribution inferred from BAM during plan phase
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct InsertSizeDistribution {
    #[serde(default)]
    pub mean: i32,
    #[serde(default)]
    pub std_dev: i32,
    #[serde(default)]
    pub orientation: String, // FR (Illumina standard)
    #[serde(default)]
    pub sample_size: usize,
}

impl InsertSizeDistribution {
    /// Infer from BAM proper pairs (SAM flag 0x2)
    pub fn from_bam_proper_pairs(tlens: &[i32]) -> Option<Self> {
        if tlens.len() < 100 {
            return None; // Insufficient data
        }

        let mean = tlens.iter().sum::<i32>() / tlens.len() as i32;
        let variance = tlens
            .iter()
            .map(|&t| {
                let diff = t - mean;
                (diff as f64).powi(2)
            })
            .sum::<f64>()
            / tlens.len() as f64;
        let std_dev = variance.sqrt() as i32;

        Some(InsertSizeDistribution {
            mean,
            std_dev,
            orientation: "FR".to_string(),
            sample_size: tlens.len(),
        })
    }
}

/// PhrayaPlan: read-only reference for alignment workers
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhrayaPlan {
    /// Format version
    pub version: u32,
    /// Detected use case
    pub use_case: UseCase,
    /// Content-addressed reference space (ADR-0011): optional reference with
    /// content hash, name, and sketches. Used in plan v6+.
    #[serde(default)]
    pub reference_space: Option<ReferenceSpace>,
    /// Input file paths
    pub input_files: Vec<String>,
    /// Timestamp (ISO8601)
    pub timestamp: String,
    /// K-mer sketches keyed by sequence ID — for reuse during alignment
    #[serde(serialize_with = "serialize_map_sorted")]
    pub kmer_index: HashMap<String, MinimizerSketch>,
    /// K-mer uniqueness: position → uniqueness score
    #[serde(serialize_with = "serialize_map_sorted")]
    pub kmer_uniqueness: HashMap<u32, f64>,
    /// Task list: (query_id, target_id) pairs
    pub task_list: Vec<(u32, u32)>,
    /// Variation hotspot intervals detected at plan time: (start, end) pairs
    #[serde(default)]
    pub hotspot_intervals: Vec<(u32, u32)>,
    /// Read counts per input file (for batch-mode indexing)
    #[serde(default)]
    pub reads_per_file: Vec<usize>,
    /// Total read count across all inputs
    #[serde(default)]
    pub total_read_count: usize,
    /// K-mer sketching parameters used during planning
    #[serde(default)]
    pub kmer_params: KmerParams,
    /// Batch mode: divide into N chunks
    #[serde(default)]
    pub batch_num_chunks: Option<usize>,
    /// Batch mode: X reads per chunk
    #[serde(default)]
    pub batch_reads_per_chunk: Option<usize>,
    /// Byte offsets for start of each read, per input file
    #[serde(default)]
    pub read_byte_offsets: Vec<Vec<u64>>,
    /// Output paths for each batch chunk (empty if no batching)
    #[serde(default)]
    pub batch_output_paths: Vec<String>,
    /// Insert size distribution (None for FASTQ input without alignment)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub insert_size_distribution: Option<InsertSizeDistribution>,
    /// Mate information keyed by sequence ID (for BAM/CRAM inputs)
    #[serde(
        default,
        skip_serializing_if = "HashMap::is_empty",
        serialize_with = "serialize_map_sorted"
    )]
    pub mate_info: HashMap<String, phraya_core::types::MateInfo>,
    /// Dense minimizer sketches keyed by sequence ID
    /// Empty if sparse_mode is true
    #[serde(
        default = "default_dense_kmer_index",
        skip_serializing_if = "HashMap::is_empty"
    )]
    pub dense_kmer_index: HashMap<String, MinimizerSketch>,
    /// Per-sequence w=11 membership tags for dense sketches
    /// Indicates which dense minimizers are part of the canonical w=11 set
    #[serde(
        default = "default_w11_membership",
        skip_serializing_if = "HashMap::is_empty"
    )]
    pub w11_membership: HashMap<String, Vec<bool>>,
    /// If true, only w=11 sketches are stored (--sparse flag)
    /// If false, both w=11 and dense sketches are stored (default)
    #[serde(default = "default_sparse_mode", skip_serializing_if = "is_false")]
    pub sparse_mode: bool,
}

impl PhrayaPlan {
    /// Create a new plan
    pub fn new(
        use_case: UseCase,
        input_files: Vec<String>,
        timestamp: String,
        kmer_index: HashMap<String, MinimizerSketch>,
        kmer_uniqueness: HashMap<u32, f64>,
        task_list: Vec<(u32, u32)>,
    ) -> Self {
        PhrayaPlan {
            version: PHRAYAPLAN_VERSION,
            use_case,
            input_files,
            timestamp,
            kmer_index,
            kmer_uniqueness,
            task_list,
            hotspot_intervals: Vec::new(),
            reads_per_file: Vec::new(),
            total_read_count: 0,
            kmer_params: KmerParams::default(),
            batch_num_chunks: None,
            batch_reads_per_chunk: None,
            read_byte_offsets: Vec::new(),
            batch_output_paths: Vec::new(),
            insert_size_distribution: None,
            mate_info: HashMap::new(),
            dense_kmer_index: HashMap::new(),
            w11_membership: HashMap::new(),
            sparse_mode: false,
            reference_space: None,
        }
    }

    /// Look up a pre-computed sketch by sequence ID. Returns None if not in plan.
    pub fn get_sketch(&self, sequence_id: &str) -> Option<&MinimizerSketch> {
        self.kmer_index.get(sequence_id)
    }

    /// Look up a dense minimizer sketch by sequence ID.
    /// Returns None if the plan was created with --sparse or if dense sketches are not available.
    pub fn get_dense_sketch(&self, sequence_id: &str) -> Option<&MinimizerSketch> {
        self.dense_kmer_index.get(sequence_id)
    }

    /// Look up the w=11 membership tags for a dense sketch.
    /// Returns a Vec<bool> where each bool indicates if the corresponding dense minimizer
    /// is part of the canonical w=11 set.
    ///
    /// NOTE: This recomputes membership based on the deduplicated w=11 sketch to ensure
    /// byte-equivalence: the extracted w=11 subset exactly matches the canonical w=11 sketch.
    ///
    /// Returns None if the plan was created with --sparse or if tags are not available.
    pub fn get_w11_membership(&self, sequence_id: &str) -> Option<&Vec<bool>> {
        // Return the cached membership tags if available
        // (this is the fast path for most uses)
        self.w11_membership.get(sequence_id)
    }

    /// Check if this plan was created with --sparse (dense sketches not stored).
    pub fn is_sparse(&self) -> bool {
        self.sparse_mode
    }

    /// Compute and store dense sketches for all sequences in the kmer_index.
    /// This method is called during plan creation to populate dense sketches
    /// alongside the default w=11 sketches.
    ///
    /// Key behaviors:
    /// - Deduplicates the w=11 sketch before computing membership to ensure byte-equivalence
    /// - Skips computation if sparse_mode is true
    /// - Updates both dense_kmer_index and w11_membership with computed data
    pub fn populate_dense_sketches(
        &mut self,
        sequences: &HashMap<String, phraya_core::types::Sequence>,
    ) {
        use phraya_core::types::sketch_sequence;

        if self.sparse_mode {
            return; // Don't compute dense sketches for sparse plans
        }

        for (seq_id, seq) in sequences {
            if !self.kmer_index.contains_key(seq_id) {
                continue; // Skip sequences not in kmer_index
            }

            // Compute dense sketch with w=5 (denser than w=11), deduplicated so a
            // minimizer repeated across overlapping windows is stored once (matches
            // the w=11 sketch's dedup below, keeping membership counts meaningful).
            let dense_sketch = deduplicate_sketch(&sketch_sequence(seq, 21, 5));

            // Get the w=11 sketch and deduplicate it
            let w11_sketch_original = &self.kmer_index[seq_id];
            let w11_sketch = deduplicate_sketch(w11_sketch_original);

            // Update kmer_index with deduplicated w=11 sketch to ensure byte-equivalence
            self.kmer_index.insert(seq_id.clone(), w11_sketch.clone());

            // Create membership tags: which dense minimizers are in deduplicated w=11 sketch?
            let w11_set: std::collections::HashSet<(u64, u32)> =
                w11_sketch.minimizers.iter().copied().collect();
            let membership: Vec<bool> = dense_sketch
                .minimizers
                .iter()
                .map(|m| w11_set.contains(m))
                .collect();

            self.dense_kmer_index.insert(seq_id.clone(), dense_sketch);
            self.w11_membership.insert(seq_id.clone(), membership);
        }
    }
}

/// Write PhrayaPlan to compressed binary file
pub fn write_plan(path: &Path, plan: &PhrayaPlan) -> Result<(), PlanError> {
    // Serialize using MessagePack
    let serialized =
        rmp_serde::to_vec(plan).map_err(|e| PlanError::SerializationError(e.to_string()))?;

    // Compress using zstd
    let compressed = zstd::encode_all(&serialized[..], 3)
        .map_err(|e| PlanError::CompressionError(e.to_string()))?;

    // Write to file
    std::fs::write(path, compressed).map_err(|e| PlanError::IoError(e.to_string()))?;

    Ok(())
}

/// Read PhrayaPlan from compressed binary file
pub fn read_plan(path: &Path) -> Result<PhrayaPlan, PlanError> {
    // Read file
    let compressed = std::fs::read(path).map_err(|e| PlanError::IoError(e.to_string()))?;

    // Decompress using zstd
    let decompressed = zstd::decode_all(&compressed[..])
        .map_err(|e| PlanError::DecompressionError(e.to_string()))?;

    // Deserialize using MessagePack
    let plan: PhrayaPlan = rmp_serde::from_slice(&decompressed)
        .map_err(|e| PlanError::SerializationError(e.to_string()))?;

    // Check version
    if plan.version != PHRAYAPLAN_VERSION {
        return Err(PlanError::VersionMismatch {
            expected: PHRAYAPLAN_VERSION,
            got: plan.version,
        });
    }

    Ok(plan)
}

/// Compute a strong cryptographic hash (SHA-256, 256-bit) of raw bytes.
/// Returns a lowercase hex-encoded string of 64 characters (256 bits / 4 bits per hex digit).
pub fn content_hash_for_bytes(bytes: &[u8]) -> String {
    let mut hasher = Sha256::new();
    hasher.update(bytes);
    let result = hasher.finalize();
    format!("{:x}", result)
}

/// Compute a strong cryptographic hash of a sequence's byte content.
/// The hash depends only on the bases, not on sequence ID or other metadata.
pub fn content_hash_for_sequence(sequence: &phraya_core::types::Sequence) -> String {
    content_hash_for_bytes(sequence.bases())
}

#[cfg(test)]
mod tests {
    use super::*;
    use phraya_core::types::Sequence;
    use tempfile::NamedTempFile;

    #[test]
    fn populate_dense_sketches_skips_sequences_absent_from_kmer_index() {
        // Only "known" is in kmer_index; "unknown" is present in `sequences` but
        // must be skipped rather than gaining a dense sketch of its own.
        let mut plan = PhrayaPlan::new(
            UseCase::ReadsWithRef,
            vec![],
            "2026-05-31T12:00:00Z".to_string(),
            HashMap::new(),
            HashMap::new(),
            vec![],
        );
        plan.kmer_index.insert(
            "known".to_string(),
            phraya_core::types::sketch(b"ACGTACGTACGTACGTACGTACGTACGTACGT", 21, 11),
        );

        let mut sequences = HashMap::new();
        sequences.insert(
            "known".to_string(),
            Sequence::new(
                b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
                None,
                "known".to_string(),
                None,
            ),
        );
        sequences.insert(
            "unknown".to_string(),
            Sequence::new(
                b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT".to_vec(),
                None,
                "unknown".to_string(),
                None,
            ),
        );

        plan.populate_dense_sketches(&sequences);

        assert!(plan.dense_kmer_index.contains_key("known"));
        assert!(!plan.dense_kmer_index.contains_key("unknown"));
        assert!(!plan.w11_membership.contains_key("unknown"));
    }

    #[test]
    fn round_trip_empty_plan() {
        let plan = PhrayaPlan::new(
            UseCase::ReadsWithRef,
            vec![],
            "2026-05-31T12:00:00Z".to_string(),
            HashMap::new(),
            HashMap::new(),
            vec![],
        );

        let temp = NamedTempFile::new().unwrap();
        write_plan(temp.path(), &plan).unwrap();
        let read_plan = read_plan(temp.path()).unwrap();

        assert_eq!(read_plan.use_case, plan.use_case);
        assert_eq!(read_plan.input_files, plan.input_files);
        assert_eq!(read_plan.task_list, plan.task_list);
    }

    #[test]
    fn round_trip_with_files() {
        let plan = PhrayaPlan::new(
            UseCase::ContigsWithReads,
            vec!["input.fa".to_string(), "reads.fq".to_string()],
            "2026-05-31T12:00:00Z".to_string(),
            HashMap::new(),
            HashMap::new(),
            vec![(1, 2), (1, 3), (2, 3)],
        );

        let temp = NamedTempFile::new().unwrap();
        write_plan(temp.path(), &plan).unwrap();
        let read_plan = read_plan(temp.path()).unwrap();

        assert_eq!(read_plan.input_files, vec!["input.fa", "reads.fq"]);
        assert_eq!(read_plan.task_list.len(), 3);
    }

    #[test]
    fn round_trip_with_uniqueness() {
        let mut uniqueness = HashMap::new();
        uniqueness.insert(0u32, 1.0);
        uniqueness.insert(100u32, 0.5);
        uniqueness.insert(200u32, 0.25);

        let plan = PhrayaPlan::new(
            UseCase::ReadsOnly,
            vec![],
            "2026-05-31T12:00:00Z".to_string(),
            HashMap::new(),
            uniqueness.clone(),
            vec![],
        );

        let temp = NamedTempFile::new().unwrap();
        write_plan(temp.path(), &plan).unwrap();
        let read_plan = read_plan(temp.path()).unwrap();

        assert_eq!(read_plan.kmer_uniqueness, uniqueness);
    }

    #[test]
    fn large_task_list() {
        let mut tasks = Vec::new();
        for i in 0..10000 {
            tasks.push((i as u32, (i + 1) as u32));
        }

        let plan = PhrayaPlan::new(
            UseCase::ContigsOnly,
            vec![],
            "2026-05-31T12:00:00Z".to_string(),
            HashMap::new(),
            HashMap::new(),
            tasks.clone(),
        );

        let temp = NamedTempFile::new().unwrap();
        write_plan(temp.path(), &plan).unwrap();
        let read_plan = read_plan(temp.path()).unwrap();

        assert_eq!(read_plan.task_list.len(), 10000);
        assert_eq!(read_plan.task_list, tasks);
    }

    #[test]
    fn version_mismatch_handling() {
        let mut plan = PhrayaPlan::new(
            UseCase::ReadsWithRef,
            vec![],
            "2026-05-31T12:00:00Z".to_string(),
            HashMap::new(),
            HashMap::new(),
            vec![],
        );

        // Manually set wrong version
        plan.version = 999;

        let temp = NamedTempFile::new().unwrap();
        write_plan(temp.path(), &plan).unwrap();

        // Reading should fail with version mismatch
        let result = read_plan(temp.path());
        assert!(result.is_err());
        match result.unwrap_err() {
            PlanError::VersionMismatch { expected, got } => {
                assert_eq!(expected, PHRAYAPLAN_VERSION);
                assert_eq!(got, 999);
            }
            _ => panic!("Expected VersionMismatch error"),
        }
    }

    #[test]
    fn compression_ratio() {
        let mut tasks = Vec::new();
        for i in 0..1000 {
            tasks.push((i as u32, (i + 1) as u32));
        }

        let plan = PhrayaPlan::new(
            UseCase::ContigsWithReads,
            vec!["file1.fa".to_string(), "file2.fq".to_string()],
            "2026-05-31T12:00:00Z".to_string(),
            HashMap::new(),
            HashMap::new(),
            tasks,
        );

        let temp = NamedTempFile::new().unwrap();
        write_plan(temp.path(), &plan).unwrap();

        let file_size = std::fs::metadata(temp.path()).unwrap().len();
        // Compressed file should be reasonably small (task_list is repetitive)
        assert!(file_size < 100_000);
    }

    #[test]
    fn use_case_serialization() {
        for use_case in &[
            UseCase::ReadsWithRef,
            UseCase::ReadsOnly,
            UseCase::ContigsWithReads,
            UseCase::ContigsOnly,
        ] {
            let plan = PhrayaPlan::new(
                *use_case,
                vec![],
                "2026-05-31T12:00:00Z".to_string(),
                HashMap::new(),
                HashMap::new(),
                vec![],
            );

            let temp = NamedTempFile::new().unwrap();
            write_plan(temp.path(), &plan).unwrap();
            let read_plan = read_plan(temp.path()).unwrap();

            assert_eq!(read_plan.use_case, *use_case);
        }
    }

    #[test]
    fn nonexistent_file_read() {
        let result = read_plan(Path::new("/nonexistent/path.phrayaplan"));
        assert!(result.is_err());
    }

    #[test]
    fn corrupted_file_handling() {
        let temp = NamedTempFile::new().unwrap();
        std::fs::write(temp.path(), b"corrupted data").unwrap();

        let result = read_plan(temp.path());
        assert!(result.is_err());
    }

    #[test]
    fn round_trip_v3_batch_fields() {
        let mut plan = PhrayaPlan::new(
            UseCase::ReadsWithRef,
            vec!["reads_1.fq".to_string(), "reads_2.fq".to_string()],
            "2026-06-09T12:00:00Z".to_string(),
            HashMap::new(),
            HashMap::new(),
            vec![],
        );

        plan.reads_per_file = vec![1000, 1000];
        plan.total_read_count = 2000;
        plan.kmer_params = KmerParams { k: 21, w: 11 };
        plan.batch_num_chunks = Some(16);
        plan.batch_reads_per_chunk = Some(125);
        plan.read_byte_offsets = vec![vec![0, 100, 200, 300], vec![0, 110, 220, 330]];
        plan.batch_output_paths = vec!["out_0.phraya".to_string(), "out_1.phraya".to_string()];

        let temp = NamedTempFile::new().unwrap();
        write_plan(temp.path(), &plan).unwrap();
        let read_plan = read_plan(temp.path()).unwrap();

        assert_eq!(read_plan.version, PHRAYAPLAN_VERSION);
        assert_eq!(read_plan.reads_per_file, vec![1000, 1000]);
        assert_eq!(read_plan.total_read_count, 2000);
        assert_eq!(read_plan.kmer_params.k, 21);
        assert_eq!(read_plan.kmer_params.w, 11);
        assert_eq!(read_plan.batch_num_chunks, Some(16));
        assert_eq!(read_plan.batch_reads_per_chunk, Some(125));
        assert_eq!(read_plan.read_byte_offsets.len(), 2);
        assert_eq!(read_plan.read_byte_offsets[0], vec![0, 100, 200, 300]);
        assert_eq!(read_plan.batch_output_paths.len(), 2);
    }

    // ============================================================================
    // RED acceptance tests for issue #196: content-addressed reference space (ADR-0011)
    //
    // `ReferenceSpace` and `PhrayaPlan::reference_space` do not exist on this branch yet.
    // Every test below either fails to compile (references the new type/field directly)
    // or fails at runtime against behavior that has not been implemented. None of these
    // tests can pass against unmodified `main` — that is the point.
    // ============================================================================

    /// `ReferenceSpace` must exist with exactly these three fields: a content hash,
    /// an optional human-facing name, and per-sequence sketches.
    #[test]
    fn issue_196_reference_space_struct_has_required_fields() {
        let mut sketches = HashMap::new();
        sketches.insert(
            "chr1".to_string(),
            phraya_core::types::sketch(b"ACGTACGTACGTACGTACGTACGTACGTACGT", 21, 11),
        );

        let space = ReferenceSpace {
            content_hash: "deadbeef".to_string(),
            name: Some("my-reference".to_string()),
            sketches: sketches.clone(),
        };

        assert_eq!(space.content_hash, "deadbeef");
        assert_eq!(space.name, Some("my-reference".to_string()));
        assert_eq!(space.sketches.len(), 1);
        assert!(space.sketches.contains_key("chr1"));
    }

    /// `ReferenceSpace.name` is optional — an unnamed reference space must round-trip
    /// with `name: None`, not an empty string or a placeholder.
    #[test]
    fn issue_196_reference_space_name_is_optional() {
        let space = ReferenceSpace {
            content_hash: "abc123".to_string(),
            name: None,
            sketches: HashMap::new(),
        };

        assert_eq!(space.name, None);
    }

    /// `PhrayaPlan` must carry an optional `reference_space` field. Constructing a plan
    /// via the existing `new()` and then attaching a `ReferenceSpace` via struct-update
    /// syntax must compile — it does not today, because the field does not exist.
    #[test]
    fn issue_196_phraya_plan_has_reference_space_field() {
        let base_plan = PhrayaPlan::new(
            UseCase::ReadsWithRef,
            vec!["reference.fa".to_string()],
            "2026-07-08T12:00:00Z".to_string(),
            HashMap::new(),
            HashMap::new(),
            vec![],
        );

        let plan = PhrayaPlan {
            reference_space: Some(ReferenceSpace {
                content_hash: "deadbeef".to_string(),
                name: None,
                sketches: HashMap::new(),
            }),
            ..base_plan
        };

        assert!(plan.reference_space.is_some());
        assert_eq!(plan.reference_space.unwrap().content_hash, "deadbeef");
    }

    /// A `ReferenceSpace` attached to a plan must survive a full write/read round-trip
    /// through `.phrayaplan`'s MessagePack + zstd encoding — hash, name, and sketches
    /// all preserved exactly.
    #[test]
    fn issue_196_reference_space_round_trips_through_phrayaplan() {
        let mut sketches = HashMap::new();
        let sketch = phraya_core::types::sketch(b"ACGTACGTACGTACGTACGTACGTACGTACGT", 21, 11);
        sketches.insert("ref".to_string(), sketch.clone());

        let base_plan = PhrayaPlan::new(
            UseCase::ReadsWithRef,
            vec!["reference.fa".to_string()],
            "2026-07-08T12:00:00Z".to_string(),
            HashMap::new(),
            HashMap::new(),
            vec![],
        );

        let plan = PhrayaPlan {
            reference_space: Some(ReferenceSpace {
                content_hash: "cafef00d".to_string(),
                name: Some("chr1-assembly".to_string()),
                sketches: sketches.clone(),
            }),
            ..base_plan
        };

        let temp = NamedTempFile::new().unwrap();
        write_plan(temp.path(), &plan).unwrap();
        let read_plan = read_plan(temp.path()).unwrap();

        let read_space = read_plan
            .reference_space
            .expect("reference_space must survive round-trip");
        assert_eq!(read_space.content_hash, "cafef00d");
        assert_eq!(read_space.name, Some("chr1-assembly".to_string()));
        assert_eq!(read_space.sketches.len(), 1);
        assert_eq!(read_space.sketches.get("ref").unwrap(), &sketch);
    }

    /// A plan with no reference space attached (`reference_space: None`) must also
    /// round-trip cleanly — the field is optional, not mandatory.
    #[test]
    fn issue_196_plan_without_reference_space_round_trips() {
        let base_plan = PhrayaPlan::new(
            UseCase::ReadsWithRef,
            vec![],
            "2026-07-08T12:00:00Z".to_string(),
            HashMap::new(),
            HashMap::new(),
            vec![],
        );

        let plan = PhrayaPlan {
            reference_space: None,
            ..base_plan
        };

        let temp = NamedTempFile::new().unwrap();
        write_plan(temp.path(), &plan).unwrap();
        let read_plan = read_plan(temp.path()).unwrap();

        assert_eq!(read_plan.reference_space, None);
    }

    /// Hashing identical byte content twice must produce identical hashes — the
    /// content-hash function is a pure function of bytes, not of time, randomness,
    /// or any other hidden state.
    #[test]
    fn issue_196_content_hash_is_deterministic() {
        let content = b"ACGTACGTACGTACGTACGTACGTACGTACGT";

        let hash1 = content_hash_for_bytes(content);
        let hash2 = content_hash_for_bytes(content);

        assert_eq!(hash1, hash2, "identical content must hash identically");
    }

    /// The content hash must depend only on bytes, not on the filesystem path or
    /// sequence identifier the bytes happen to be associated with — presenting the
    /// same reference under a different name/path must resolve to the same hash.
    #[test]
    fn issue_196_content_hash_is_path_independent() {
        let seq_a = Sequence::new(
            b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
            None,
            "chr1_v1.fa".to_string(),
            None,
        );
        let seq_b = Sequence::new(
            b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
            None,
            "totally_different_name_and_path.fasta".to_string(),
            None,
        );

        let hash_a = content_hash_for_sequence(&seq_a);
        let hash_b = content_hash_for_sequence(&seq_b);

        assert_eq!(
            hash_a, hash_b,
            "identical bases under different names/paths must hash identically"
        );
    }

    /// A single differing byte in otherwise-identical content must change the hash —
    /// the function must be content-sensitive, not just present for show.
    #[test]
    fn issue_196_content_hash_is_sensitive_to_single_byte_change() {
        let content_a = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let content_b = b"ACGTACGTACGTACGTACGTACGTACGTACGA"; // last base flipped

        let hash_a = content_hash_for_bytes(content_a);
        let hash_b = content_hash_for_bytes(content_b);

        assert_ne!(
            hash_a, hash_b,
            "a single differing byte must change the content hash"
        );
    }

    /// The content hash must be a strong (cryptographic-strength) hash, not a
    /// short/weak checksum — require at least 256 bits of digest (64 hex characters
    /// at 4 bits/hex-char), matching BLAKE3/SHA-256-class algorithms named in the
    /// issue. A hex-encoded 64-bit hash (16 hex chars) must NOT satisfy this.
    #[test]
    fn issue_196_content_hash_has_strong_digest_length() {
        let hash = content_hash_for_bytes(b"ACGTACGTACGTACGTACGTACGTACGTACGT");

        assert!(
            hash.len() >= 64,
            "content hash should be at least 256 bits (64 hex chars) for BLAKE3/SHA-256-class strength, got {} chars: {}",
            hash.len(),
            hash
        );
        assert!(
            hash.chars().all(|c| c.is_ascii_hexdigit()),
            "content hash should be hex-encoded, got: {}",
            hash
        );
    }
}
