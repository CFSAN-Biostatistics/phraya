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

/// Fast 64-bit non-cryptographic hash for read content (ADR-0011).
/// Uses FNV-1a for speed and determinism. Fast hash is suitable for caching
/// by content within a single pipeline run; it is NOT a cryptographic identity.
pub fn read_content_hash(bytes: &[u8]) -> u64 {
    const FNV_OFFSET_BASIS: u64 = 0xcbf29ce484222325;
    const FNV_PRIME: u64 = 0x100000001b3;

    let mut hash = FNV_OFFSET_BASIS;
    for &byte in bytes {
        hash ^= byte as u64;
        hash = hash.wrapping_mul(FNV_PRIME);
    }
    hash
}

/// PhrayaPlan format version for forward compatibility
pub const PHRAYAPLAN_VERSION: u32 = 7;

/// Magic bytes at the start of a v7 plan file to distinguish from older single-frame format.
/// ASCII "PHR7" — chosen to be invalid zstd frame magic (0xFD2FB528) so old readers fail fast.
const V7_MAGIC: [u8; 4] = [b'P', b'H', b'R', b'7'];

/// Table of contents for a v7 chunk-addressable plan file.
/// Stored as the first frame after the 4-byte magic + 8-byte TOC-length prefix.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PlanToc {
    /// Format version (must be 7)
    pub version: u32,
    /// Reserved for future extension (set to 0)
    pub flags: u64,
    /// Number of read-sketch chunk frames
    pub num_chunks: u32,
    /// Byte offset (from file start) to the compressed shared frame
    pub shared_frame_offset: u64,
    /// Compressed byte length of the shared frame
    pub shared_frame_len: u64,
    /// (offset, compressed_len) for each chunk frame, in chunk order
    pub chunk_frame_offsets: Vec<(u64, u64)>,
}

/// The shared frame: metadata + reference palette. Loaded by every worker.
#[derive(Debug, Clone, Serialize, Deserialize)]
struct SharedFrame {
    pub version: u32,
    pub use_case: UseCase,
    pub reference_space: Vec<ReferenceSpace>,
    pub input_files: Vec<String>,
    pub timestamp: String,
    #[serde(serialize_with = "serialize_map_sorted")]
    pub kmer_index: HashMap<String, MinimizerSketch>,
    #[serde(serialize_with = "serialize_map_sorted")]
    pub kmer_uniqueness: HashMap<u32, f64>,
    pub task_list: Vec<(u32, u32)>,
    pub hotspot_intervals: Vec<(u32, u32)>,
    pub reads_per_file: Vec<usize>,
    pub total_read_count: usize,
    pub kmer_params: KmerParams,
    pub batch_num_chunks: Option<usize>,
    pub batch_reads_per_chunk: Option<usize>,
    pub batch_output_paths: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub insert_size_distribution: Option<InsertSizeDistribution>,
    #[serde(
        default,
        skip_serializing_if = "HashMap::is_empty",
        serialize_with = "serialize_map_sorted"
    )]
    pub dense_kmer_index: HashMap<String, MinimizerSketch>,
    #[serde(
        default,
        skip_serializing_if = "HashMap::is_empty",
        serialize_with = "serialize_map_sorted"
    )]
    pub w11_membership: HashMap<String, Vec<bool>>,
    #[serde(default, skip_serializing_if = "is_false")]
    pub sparse_mode: bool,
}

/// A single chunk frame: read sketches + byte offsets + mate info for one positional slice.
#[derive(Debug, Clone, Serialize, Deserialize)]
struct ChunkFrame {
    /// Read sketches for this chunk, keyed by content hash
    #[serde(serialize_with = "serialize_map_sorted")]
    pub read_sketches: HashMap<u64, MinimizerSketch>,
    /// Byte offsets for this chunk's reads, per input file
    pub read_byte_offsets: Vec<Vec<u64>>,
    /// Mate info for this chunk's reads
    #[serde(
        default,
        skip_serializing_if = "HashMap::is_empty",
        serialize_with = "serialize_map_sorted"
    )]
    pub mate_info: HashMap<String, phraya_core::types::MateInfo>,
}

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

fn default_read_sketches() -> HashMap<u64, MinimizerSketch> {
    HashMap::new()
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
    pub reference_space: Vec<ReferenceSpace>,
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
    /// Read sketches keyed by content hash, for reuse across pipeline stages
    /// (ADR-0011). Distinct from kmer_index, which is keyed by sequence ID
    /// and holds reference sketches. Empty by default.
    #[serde(default = "default_read_sketches")]
    pub read_sketches: HashMap<u64, MinimizerSketch>,
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
    /// Ordered list of read content hashes, in file-encounter (positional) order.
    /// Used by the v7 writer to partition read_sketches into chunks. Not serialized
    /// in the monolithic format — reconstructed on read from chunk frame ordering.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub read_hash_order: Vec<u64>,
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
            read_sketches: HashMap::new(),
            sparse_mode: false,
            reference_space: Vec::new(),
            read_hash_order: Vec::new(),
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

    /// Look up a reference space by its content hash.
    /// Returns None if no space with that hash exists in the palette.
    pub fn get_reference_space(&self, hash: &str) -> Option<&ReferenceSpace> {
        self.reference_space
            .iter()
            .find(|space| space.content_hash == hash)
    }

    /// Look up a stored read sketch by content hash
    pub fn get_read_sketch(&self, hash: u64) -> Option<&MinimizerSketch> {
        self.read_sketches.get(&hash)
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

/// Partition read_sketches into N chunk frames by positional order.
/// Uses `read_hash_order` to determine which hashes belong to which chunk.
/// If `read_hash_order` is empty, all sketches go into chunk 0.
fn partition_read_data(plan: &PhrayaPlan, num_chunks: usize) -> Vec<ChunkFrame> {
    if num_chunks == 0 {
        return vec![];
    }

    let total_reads = if plan.read_hash_order.is_empty() {
        plan.read_sketches.len()
    } else {
        plan.read_hash_order.len()
    };
    let chunk_size = if total_reads == 0 {
        0
    } else {
        (total_reads + num_chunks - 1) / num_chunks
    };

    let mut chunks: Vec<ChunkFrame> = (0..num_chunks)
        .map(|_| ChunkFrame {
            read_sketches: HashMap::new(),
            read_byte_offsets: Vec::new(),
            mate_info: HashMap::new(),
        })
        .collect();

    if plan.read_hash_order.is_empty() {
        // Fallback: no ordering info, dump all into chunk 0
        chunks[0].read_sketches = plan.read_sketches.clone();
        chunks[0].read_byte_offsets = plan.read_byte_offsets.clone();
        chunks[0].mate_info = plan.mate_info.clone();
    } else {
        // Partition by positional index
        for (pos, hash) in plan.read_hash_order.iter().enumerate() {
            let chunk_idx = if chunk_size == 0 {
                0
            } else {
                (pos / chunk_size).min(num_chunks - 1)
            };
            if let Some(sketch) = plan.read_sketches.get(hash) {
                chunks[chunk_idx].read_sketches.insert(*hash, sketch.clone());
            }
        }

        // Partition read_byte_offsets by positional range
        if !plan.read_byte_offsets.is_empty() {
            // read_byte_offsets[file_idx][read_idx] — partition inner vecs by chunk
            for chunk_idx in 0..num_chunks {
                let start = chunk_idx * chunk_size;
                let end = ((chunk_idx + 1) * chunk_size).min(total_reads);
                let mut chunk_offsets = Vec::new();
                for file_offsets in &plan.read_byte_offsets {
                    let slice_start = start.min(file_offsets.len());
                    let slice_end = end.min(file_offsets.len());
                    chunk_offsets.push(file_offsets[slice_start..slice_end].to_vec());
                }
                chunks[chunk_idx].read_byte_offsets = chunk_offsets;
            }
        }

        // Partition mate_info: we need to know which read names belong to which chunk.
        // mate_info is keyed by read name — we don't have a position→name map easily.
        // For now, put all mate_info in chunk 0 (conservative; a worker that doesn't own
        // a read's mate_info can still function — it just won't use insert-size info).
        // TODO: If mate_info partitioning becomes important, add a name→position map.
        if !plan.mate_info.is_empty() {
            chunks[0].mate_info = plan.mate_info.clone();
        }
    }

    chunks
}

/// Write PhrayaPlan to a v7 chunk-addressable file with seekable zstd frames.
///
/// File layout:
///   [4B magic "PHR7"] [8B LE toc_offset] [shared_zstd] [chunk_0_zstd] ... [chunk_N-1_zstd] [toc_msgpack]
///
/// TOC is at the end (after all data frames), so offsets are known before serializing it.
/// TOC is uncompressed msgpack. Shared and chunk frames are zstd-compressed.
/// All byte offsets in the TOC are absolute from file start.
pub fn write_plan(path: &Path, plan: &PhrayaPlan) -> Result<(), PlanError> {
    use std::io::Write;

    let num_chunks = plan.batch_num_chunks.unwrap_or(1).max(1);

    // Build shared frame
    let shared = SharedFrame {
        version: PHRAYAPLAN_VERSION,
        use_case: plan.use_case,
        reference_space: plan.reference_space.clone(),
        input_files: plan.input_files.clone(),
        timestamp: plan.timestamp.clone(),
        kmer_index: plan.kmer_index.clone(),
        kmer_uniqueness: plan.kmer_uniqueness.clone(),
        task_list: plan.task_list.clone(),
        hotspot_intervals: plan.hotspot_intervals.clone(),
        reads_per_file: plan.reads_per_file.clone(),
        total_read_count: plan.total_read_count,
        kmer_params: plan.kmer_params.clone(),
        batch_num_chunks: plan.batch_num_chunks,
        batch_reads_per_chunk: plan.batch_reads_per_chunk,
        batch_output_paths: plan.batch_output_paths.clone(),
        insert_size_distribution: plan.insert_size_distribution.clone(),
        dense_kmer_index: plan.dense_kmer_index.clone(),
        w11_membership: plan.w11_membership.clone(),
        sparse_mode: plan.sparse_mode,
    };

    let shared_bytes = rmp_serde::to_vec(&shared)
        .map_err(|e| PlanError::SerializationError(e.to_string()))?;
    let shared_compressed = zstd::encode_all(&shared_bytes[..], 3)
        .map_err(|e| PlanError::CompressionError(e.to_string()))?;

    // Build and compress chunk frames
    let chunk_frames = partition_read_data(plan, num_chunks);
    let mut chunk_compressed: Vec<Vec<u8>> = Vec::with_capacity(num_chunks);
    for frame in &chunk_frames {
        let bytes = rmp_serde::to_vec(frame)
            .map_err(|e| PlanError::SerializationError(e.to_string()))?;
        let compressed = zstd::encode_all(&bytes[..], 3)
            .map_err(|e| PlanError::CompressionError(e.to_string()))?;
        chunk_compressed.push(compressed);
    }

    // Compute layout: TOC goes at the end, so all data offsets are known first.
    // Header: [4B magic][8B toc_offset] = 12 bytes
    let header_len: u64 = 12;
    let shared_offset = header_len;
    let mut chunk_offsets: Vec<(u64, u64)> = Vec::with_capacity(num_chunks);
    let mut cursor = shared_offset + shared_compressed.len() as u64;
    for cc in &chunk_compressed {
        chunk_offsets.push((cursor, cc.len() as u64));
        cursor += cc.len() as u64;
    }
    let toc_offset = cursor;

    // Build TOC with final offsets (no fixpoint problem — TOC is written after all data)
    let toc = PlanToc {
        version: PHRAYAPLAN_VERSION,
        flags: 0,
        num_chunks: num_chunks as u32,
        shared_frame_offset: shared_offset,
        shared_frame_len: shared_compressed.len() as u64,
        chunk_frame_offsets: chunk_offsets,
    };
    let toc_bytes = rmp_serde::to_vec(&toc)
        .map_err(|e| PlanError::SerializationError(e.to_string()))?;

    // Write the file
    let file =
        std::fs::File::create(path).map_err(|e| PlanError::IoError(e.to_string()))?;
    let mut writer = std::io::BufWriter::new(file);

    writer
        .write_all(&V7_MAGIC)
        .map_err(|e| PlanError::IoError(e.to_string()))?;
    writer
        .write_all(&toc_offset.to_le_bytes())
        .map_err(|e| PlanError::IoError(e.to_string()))?;
    writer
        .write_all(&shared_compressed)
        .map_err(|e| PlanError::IoError(e.to_string()))?;
    for cc in &chunk_compressed {
        writer
            .write_all(cc)
            .map_err(|e| PlanError::IoError(e.to_string()))?;
    }
    // TOC at end — offsets are now finalized
    writer
        .write_all(&toc_bytes)
        .map_err(|e| PlanError::IoError(e.to_string()))?;
    writer
        .flush()
        .map_err(|e| PlanError::IoError(e.to_string()))?;

    Ok(())
}

/// Read the TOC from a v7 plan file.
/// Layout: [4B magic][8B toc_offset][data frames...][toc_msgpack at toc_offset]
fn read_toc(data: &[u8]) -> Result<PlanToc, PlanError> {
    if data.len() < 12 {
        return Err(PlanError::DecompressionError(
            "file too small for v7 header".to_string(),
        ));
    }
    if &data[0..4] != V7_MAGIC {
        return Err(PlanError::DecompressionError(
            "not a v7 plan file (bad magic)".to_string(),
        ));
    }
    let toc_offset = u64::from_le_bytes(data[4..12].try_into().unwrap()) as usize;
    if toc_offset > data.len() {
        return Err(PlanError::DecompressionError(
            "file truncated: TOC offset past EOF".to_string(),
        ));
    }
    let toc_raw = &data[toc_offset..];
    let toc: PlanToc = rmp_serde::from_slice(toc_raw)
        .map_err(|e| PlanError::SerializationError(e.to_string()))?;

    if toc.version != PHRAYAPLAN_VERSION {
        return Err(PlanError::VersionMismatch {
            expected: PHRAYAPLAN_VERSION,
            got: toc.version,
        });
    }

    Ok(toc)
}

/// Decompress and deserialize a frame at the given byte range.
fn read_frame<T: serde::de::DeserializeOwned>(data: &[u8], offset: u64, len: u64) -> Result<T, PlanError> {
    let start = offset as usize;
    let end = start + len as usize;
    if data.len() < end {
        return Err(PlanError::DecompressionError(
            "file truncated: frame extends past EOF".to_string(),
        ));
    }
    let decompressed = zstd::decode_all(&data[start..end])
        .map_err(|e| PlanError::DecompressionError(e.to_string()))?;
    rmp_serde::from_slice(&decompressed)
        .map_err(|e| PlanError::SerializationError(e.to_string()))
}

/// Reassemble a PhrayaPlan from a shared frame and chunk frames.
fn assemble_plan(shared: SharedFrame, chunks: Vec<ChunkFrame>) -> PhrayaPlan {
    let mut read_sketches = HashMap::new();
    let mut read_byte_offsets: Vec<Vec<u64>> = Vec::new();
    let mut mate_info = HashMap::new();
    let mut read_hash_order = Vec::new();

    for chunk in chunks {
        // Collect hashes in sorted order within each chunk for determinism
        let mut chunk_hashes: Vec<u64> = chunk.read_sketches.keys().copied().collect();
        chunk_hashes.sort();
        read_hash_order.extend_from_slice(&chunk_hashes);

        for (hash, sketch) in chunk.read_sketches {
            read_sketches.insert(hash, sketch);
        }
        // Merge byte offsets: extend each file's offset list
        if !chunk.read_byte_offsets.is_empty() {
            if read_byte_offsets.is_empty() {
                read_byte_offsets = chunk.read_byte_offsets;
            } else {
                for (i, offsets) in chunk.read_byte_offsets.into_iter().enumerate() {
                    if i < read_byte_offsets.len() {
                        read_byte_offsets[i].extend(offsets);
                    }
                }
            }
        }
        for (k, v) in chunk.mate_info {
            mate_info.insert(k, v);
        }
    }

    PhrayaPlan {
        version: shared.version,
        use_case: shared.use_case,
        reference_space: shared.reference_space,
        input_files: shared.input_files,
        timestamp: shared.timestamp,
        kmer_index: shared.kmer_index,
        kmer_uniqueness: shared.kmer_uniqueness,
        task_list: shared.task_list,
        read_sketches,
        hotspot_intervals: shared.hotspot_intervals,
        reads_per_file: shared.reads_per_file,
        total_read_count: shared.total_read_count,
        kmer_params: shared.kmer_params,
        batch_num_chunks: shared.batch_num_chunks,
        batch_reads_per_chunk: shared.batch_reads_per_chunk,
        read_byte_offsets,
        batch_output_paths: shared.batch_output_paths,
        insert_size_distribution: shared.insert_size_distribution,
        mate_info,
        dense_kmer_index: shared.dense_kmer_index,
        w11_membership: shared.w11_membership,
        sparse_mode: shared.sparse_mode,
        read_hash_order,
    }
}

/// Read PhrayaPlan from a v7 chunk-addressable file, loading all chunks.
/// This is the standard non-batch read path.
pub fn read_plan(path: &Path) -> Result<PhrayaPlan, PlanError> {
    let data = std::fs::read(path).map_err(|e| PlanError::IoError(e.to_string()))?;

    // Detect format: v7 starts with "PHR7" magic
    if data.len() >= 4 && &data[0..4] == V7_MAGIC {
        let toc = read_toc(&data)?;
        let shared: SharedFrame =
            read_frame(&data, toc.shared_frame_offset, toc.shared_frame_len)?;
        let mut chunks = Vec::with_capacity(toc.num_chunks as usize);
        for &(offset, len) in &toc.chunk_frame_offsets {
            let chunk: ChunkFrame = read_frame(&data, offset, len)?;
            chunks.push(chunk);
        }
        Ok(assemble_plan(shared, chunks))
    } else {
        // Not a v7 file — reject with version mismatch
        Err(PlanError::VersionMismatch {
            expected: PHRAYAPLAN_VERSION,
            got: 6,
        })
    }
}

/// Read PhrayaPlan from a v7 file, loading only the specified worker's chunk.
/// For batch mode: loads shared frame + one chunk frame, minimizing memory.
/// If the plan has only 1 chunk but worker_count > 1, loads the single chunk
/// and filters in-memory by positional range.
pub fn read_plan_worker(
    path: &Path,
    worker_id: usize,
    worker_count: usize,
) -> Result<PhrayaPlan, PlanError> {
    let data = std::fs::read(path).map_err(|e| PlanError::IoError(e.to_string()))?;

    if data.len() < 4 || &data[0..4] != V7_MAGIC {
        return Err(PlanError::VersionMismatch {
            expected: PHRAYAPLAN_VERSION,
            got: 6,
        });
    }

    let toc = read_toc(&data)?;
    let shared: SharedFrame =
        read_frame(&data, toc.shared_frame_offset, toc.shared_frame_len)?;

    let plan_chunks = toc.num_chunks as usize;

    if plan_chunks == worker_count && worker_id < plan_chunks {
        // Pre-split: load exactly our chunk
        let (offset, len) = toc.chunk_frame_offsets[worker_id];
        let chunk: ChunkFrame = read_frame(&data, offset, len)?;
        Ok(assemble_plan(shared, vec![chunk]))
    } else if plan_chunks == 1 && worker_count > 1 {
        // Fallback: load the single chunk, filter by positional range
        let (offset, len) = toc.chunk_frame_offsets[0];
        let mut chunk: ChunkFrame = read_frame(&data, offset, len)?;

        let total = chunk.read_sketches.len();
        let chunk_size = (total + worker_count - 1) / worker_count;
        let start = worker_id * chunk_size;
        let end = ((worker_id + 1) * chunk_size).min(total);

        // Sort keys for deterministic positional assignment
        let mut all_hashes: Vec<u64> = chunk.read_sketches.keys().copied().collect();
        all_hashes.sort();

        let my_hashes: std::collections::HashSet<u64> =
            all_hashes[start..end].iter().copied().collect();
        chunk.read_sketches.retain(|k, _| my_hashes.contains(k));

        // Filter byte offsets by range
        if !chunk.read_byte_offsets.is_empty() {
            for file_offsets in &mut chunk.read_byte_offsets {
                let slice_start = start.min(file_offsets.len());
                let slice_end = end.min(file_offsets.len());
                *file_offsets = file_offsets[slice_start..slice_end].to_vec();
            }
        }

        Ok(assemble_plan(shared, vec![chunk]))
    } else if worker_id < plan_chunks {
        // Mismatched chunk count — load our chunk by index
        let (offset, len) = toc.chunk_frame_offsets[worker_id];
        let chunk: ChunkFrame = read_frame(&data, offset, len)?;
        Ok(assemble_plan(shared, vec![chunk]))
    } else {
        Err(PlanError::IoError(format!(
            "worker_id {} exceeds plan chunk count {}",
            worker_id, plan_chunks
        )))
    }
}

/// Read just the TOC from a plan file (for inspection / plan-tasks without full load).
pub fn read_plan_toc(path: &Path) -> Result<PlanToc, PlanError> {
    let data = std::fs::read(path).map_err(|e| PlanError::IoError(e.to_string()))?;
    if data.len() < 4 || &data[0..4] != V7_MAGIC {
        return Err(PlanError::VersionMismatch {
            expected: PHRAYAPLAN_VERSION,
            got: 6,
        });
    }
    read_toc(&data)
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

    /// Generate a deterministic DNA sequence of given length, unique per seed.
    /// Uses a simple LCG to produce only ACGT characters.
    fn test_dna_sequence(seed: u64, len: usize) -> Vec<u8> {
        let bases = [b'A', b'C', b'G', b'T'];
        let mut state = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (0..len)
            .map(|_| {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                bases[((state >> 33) % 4) as usize]
            })
            .collect()
    }

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
        // Simulate a future-version v7 file: TOC at end with version=999
        let future_toc = PlanToc {
            version: 999,
            flags: 0,
            num_chunks: 1,
            shared_frame_offset: 12,
            shared_frame_len: 50,
            chunk_frame_offsets: vec![(62, 50)],
        };
        let toc_bytes = rmp_serde::to_vec(&future_toc).unwrap();

        let temp = NamedTempFile::new().unwrap();
        let mut file_bytes = Vec::new();
        file_bytes.extend_from_slice(&V7_MAGIC);
        // TOC offset = header (12) + padding (200) = 212
        let toc_offset: u64 = 12 + 200;
        file_bytes.extend_from_slice(&toc_offset.to_le_bytes());
        // Pad with dummy data frames (won't be read — version check fails first)
        file_bytes.extend(vec![0u8; 200]);
        // TOC at the end
        file_bytes.extend_from_slice(&toc_bytes);
        std::fs::write(temp.path(), &file_bytes).unwrap();

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
            reference_space: vec![ReferenceSpace {
                content_hash: "deadbeef".to_string(),
                name: None,
                sketches: HashMap::new(),
            }],
            ..base_plan
        };

        assert!(!plan.reference_space.is_empty());
        assert_eq!(plan.reference_space[0].content_hash, "deadbeef");
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
            reference_space: vec![ReferenceSpace {
                content_hash: "cafef00d".to_string(),
                name: Some("chr1-assembly".to_string()),
                sketches: sketches.clone(),
            }],
            ..base_plan
        };

        let temp = NamedTempFile::new().unwrap();
        write_plan(temp.path(), &plan).unwrap();
        let read_plan = read_plan(temp.path()).unwrap();

        let read_space = &read_plan.reference_space[0];
        assert_eq!(read_space.content_hash, "cafef00d");
        assert_eq!(read_space.name, Some("chr1-assembly".to_string()));
        assert_eq!(read_space.sketches.len(), 1);
        assert_eq!(read_space.sketches.get("ref").unwrap(), &sketch);
    }

    /// A plan with an empty reference palette must also
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
            reference_space: Vec::new(),
            ..base_plan
        };

        let temp = NamedTempFile::new().unwrap();
        write_plan(temp.path(), &plan).unwrap();
        let read_plan = read_plan(temp.path()).unwrap();

        assert!(read_plan.reference_space.is_empty());
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

    // ============================================================================
    // RED acceptance tests for issue #197 (plan-side slice): reference palette —
    // PhrayaPlan grows a single ReferenceSpace (#196) into a Vec<ReferenceSpace>
    // (ADR-0011). Deliberately scoped to plan.rs storage only; the align-side
    // `--reference` CLI/composability wiring is tracked as a follow-up (the
    // existing --reference flag/align modes need real design work to compose
    // with --worker/--ensure/traditional query-target modes, which is out of
    // scope for a test-immutable RED contract written without that design pass).
    // ============================================================================

    /// `PhrayaPlan.reference_space` must become a `Vec<ReferenceSpace>` (a palette),
    /// not `Option<ReferenceSpace>` (a single slot). Constructing a plan with two
    /// distinct reference spaces and reading both back must compile and succeed —
    /// it does not today, since the field is `Option`, which can hold at most one.
    #[test]
    fn issue_197_reference_space_is_a_vec_not_an_option() {
        let base_plan = PhrayaPlan::new(
            UseCase::ReadsWithRef,
            vec!["a.fa".to_string(), "b.fa".to_string()],
            "2026-07-13T00:00:00Z".to_string(),
            HashMap::new(),
            HashMap::new(),
            vec![],
        );

        let space_a = ReferenceSpace {
            content_hash: "hash_a".to_string(),
            name: Some("space-a".to_string()),
            sketches: HashMap::new(),
        };
        let space_b = ReferenceSpace {
            content_hash: "hash_b".to_string(),
            name: Some("space-b".to_string()),
            sketches: HashMap::new(),
        };

        let plan = PhrayaPlan {
            reference_space: vec![space_a.clone(), space_b.clone()],
            ..base_plan
        };

        assert_eq!(plan.reference_space.len(), 2);
        assert_eq!(plan.reference_space[0].content_hash, "hash_a");
        assert_eq!(plan.reference_space[1].content_hash, "hash_b");
    }

    /// An empty palette (no reference spaces at all) must be representable and
    /// round-trip cleanly — the empty-Vec default replaces `None` as the "no
    /// reference" case now that the field is a Vec.
    #[test]
    fn issue_197_empty_palette_round_trips() {
        let base_plan = PhrayaPlan::new(
            UseCase::ReadsOnly,
            vec![],
            "2026-07-13T00:00:00Z".to_string(),
            HashMap::new(),
            HashMap::new(),
            vec![],
        );

        let plan = PhrayaPlan {
            reference_space: vec![],
            ..base_plan
        };

        let temp = NamedTempFile::new().unwrap();
        write_plan(temp.path(), &plan).unwrap();
        let read_plan = read_plan(temp.path()).unwrap();

        assert!(read_plan.reference_space.is_empty());
    }

    /// A palette of N (N > 1) reference spaces must survive a full write/read
    /// round-trip through `.phrayaplan`'s MessagePack + zstd encoding — every
    /// space's hash, name, and sketches preserved, in order.
    #[test]
    fn issue_197_multi_space_palette_round_trips() {
        let base_plan = PhrayaPlan::new(
            UseCase::ContigsWithReads,
            vec!["a.fa".to_string(), "b.fa".to_string(), "c.fa".to_string()],
            "2026-07-13T00:00:00Z".to_string(),
            HashMap::new(),
            HashMap::new(),
            vec![],
        );

        let mut sketches_a = HashMap::new();
        sketches_a.insert(
            "seq_a".to_string(),
            phraya_core::types::sketch(b"ACGTACGTACGTACGTACGTACGTACGTACGT", 21, 11),
        );

        let spaces = vec![
            ReferenceSpace {
                content_hash: "hash_a".to_string(),
                name: Some("space-a".to_string()),
                sketches: sketches_a.clone(),
            },
            ReferenceSpace {
                content_hash: "hash_b".to_string(),
                name: None,
                sketches: HashMap::new(),
            },
            ReferenceSpace {
                content_hash: "hash_c".to_string(),
                name: Some("space-c".to_string()),
                sketches: HashMap::new(),
            },
        ];

        let plan = PhrayaPlan {
            reference_space: spaces.clone(),
            ..base_plan
        };

        let temp = NamedTempFile::new().unwrap();
        write_plan(temp.path(), &plan).unwrap();
        let read_plan = read_plan(temp.path()).unwrap();

        assert_eq!(read_plan.reference_space.len(), 3);
        assert_eq!(read_plan.reference_space, spaces);
    }

    /// `PhrayaPlan` must expose a lookup helper resolving a reference space by its
    /// content hash — the mechanism #197's align-side hit/miss resolution (a
    /// follow-up) will build on. Mirrors the existing `get_sketch`/`get_dense_sketch`
    /// accessor pattern.
    #[test]
    fn issue_197_get_reference_space_by_hash_looks_up_in_palette() {
        let base_plan = PhrayaPlan::new(
            UseCase::ContigsWithReads,
            vec!["a.fa".to_string(), "b.fa".to_string()],
            "2026-07-13T00:00:00Z".to_string(),
            HashMap::new(),
            HashMap::new(),
            vec![],
        );

        let space_a = ReferenceSpace {
            content_hash: "hash_a".to_string(),
            name: Some("space-a".to_string()),
            sketches: HashMap::new(),
        };
        let space_b = ReferenceSpace {
            content_hash: "hash_b".to_string(),
            name: Some("space-b".to_string()),
            sketches: HashMap::new(),
        };

        let plan = PhrayaPlan {
            reference_space: vec![space_a.clone(), space_b.clone()],
            ..base_plan
        };

        assert_eq!(
            plan.get_reference_space("hash_b"),
            Some(&space_b),
            "get_reference_space should find a space by its content hash"
        );
        assert_eq!(
            plan.get_reference_space("hash_nonexistent"),
            None,
            "get_reference_space should return None for an unknown hash"
        );
    }

    /// `phraya plan --reference X` (repeatable) must produce a palette with one
    /// `ReferenceSpace` per presented `--reference`, each with a real content hash
    /// computed from that file's actual bytes. This is the plan-side CLI wiring
    /// half of #197 — deliberately excludes the align-side --reference flag,
    /// which does not exist yet and needs its own design pass (see module doc).
    #[test]
    fn issue_197_plan_cli_accepts_repeatable_reference_and_stores_palette() {
        use std::process::Command;
        use tempfile::TempDir;

        fn manifest() -> std::path::PathBuf {
            // phraya-io's own manifest can't run the phraya binary; shell out via
            // the workspace's phraya-cli crate instead.
            let d = std::env::var("CARGO_MANIFEST_DIR").unwrap();
            std::path::Path::new(&d)
                .parent()
                .unwrap()
                .join("phraya-cli")
                .join("Cargo.toml")
        }

        let dir = TempDir::new().unwrap();
        let ref_a = dir.path().join("a.fa");
        let ref_b = dir.path().join("b.fa");
        let reads = dir.path().join("reads.fa");
        std::fs::write(&ref_a, ">a\nACGTACGTACGTACGTACGTACGTACGTACGT\n").unwrap();
        std::fs::write(&ref_b, ">b\nTGCATGCATGCATGCATGCATGCATGCATGCA\n").unwrap();
        std::fs::write(&reads, ">read0\nACGTACGTACGTACGTACGT\n").unwrap();
        let plan_path = dir.path().join("plan.phrayaplan");

        let out = Command::new("cargo")
            .arg("run")
            .arg("--manifest-path")
            .arg(manifest().to_str().unwrap())
            .arg("--")
            .arg("plan")
            .arg("--inputs")
            .arg(reads.to_str().unwrap())
            .arg("--reference")
            .arg(ref_a.to_str().unwrap())
            .arg("--reference")
            .arg(ref_b.to_str().unwrap())
            .arg("--output")
            .arg(plan_path.to_str().unwrap())
            .output()
            .expect("cargo run failed");

        assert!(
            out.status.success(),
            "phraya plan --reference A --reference B failed:\n{}",
            String::from_utf8_lossy(&out.stderr)
        );

        let plan = read_plan(&plan_path).unwrap();
        assert_eq!(
            plan.reference_space.len(),
            2,
            "plan --reference A --reference B should store a palette of 2 spaces, got {}",
            plan.reference_space.len()
        );

        let expected_hash_a = content_hash_for_bytes(b"ACGTACGTACGTACGTACGTACGTACGTACGT");
        let expected_hash_b = content_hash_for_bytes(b"TGCATGCATGCATGCATGCATGCATGCATGCA");
        assert!(
            plan.get_reference_space(&expected_hash_a).is_some(),
            "reference A's content hash should resolve to a stored space"
        );
        assert!(
            plan.get_reference_space(&expected_hash_b).is_some(),
            "reference B's content hash should resolve to a stored space"
        );
    }
    // ============================================================================
    // RED acceptance tests for issue #200: fat plan — read content hashing +
    // stored read sketches (ADR-0011)
    //
    // `read_content_hash`, `PhrayaPlan::read_sketches`, and
    // `PhrayaPlan::get_read_sketch` do not exist on this branch yet. Every test
    // below fails to compile or fails at runtime against unmodified main.
    // ============================================================================

    /// A fast, non-cryptographic 64-bit hash function for read content must exist,
    /// distinct from the strong (256-bit) reference content hash — the issue
    /// explicitly calls for a fast 64-bit hash (e.g. xxh3) for reads, not the
    /// cryptographic hash used for reference identity.
    #[test]
    fn issue_200_read_content_hash_is_64_bit() {
        let hash: u64 = read_content_hash(b"ACGTACGTACGTACGTACGTACGTACGTACGT");
        // A u64 return type is itself the strongest assertion here — this line
        // fails to compile if read_content_hash returns anything else (e.g. a
        // 256-bit hex String like the reference hash), and that's intentional:
        // the fast/strong hash distinction is the point of this test.
        let _: u64 = hash;
    }

    /// Hashing identical read bytes twice must produce identical hashes.
    #[test]
    fn issue_200_read_content_hash_is_deterministic() {
        let bases = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        assert_eq!(
            read_content_hash(bases),
            read_content_hash(bases),
            "identical read content must hash identically"
        );
    }

    /// A single differing byte must change the read hash — the function must be
    /// content-sensitive, not a constant or a length-only checksum.
    #[test]
    fn issue_200_read_content_hash_is_sensitive_to_content() {
        let a = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let b = b"ACGTACGTACGTACGTACGTACGTACGTACGA";
        assert_ne!(
            read_content_hash(a),
            read_content_hash(b),
            "a single differing byte must change the read content hash"
        );
    }

    /// `PhrayaPlan` must carry a `read_sketches: HashMap<u64, MinimizerSketch>`
    /// field, keyed by read content hash (distinct from `kmer_index`, which is
    /// keyed by sequence ID string and holds reference/target sketches).
    #[test]
    fn issue_200_phraya_plan_has_read_sketches_field() {
        let base_plan = PhrayaPlan::new(
            UseCase::ReadsWithRef,
            vec![],
            "2026-07-13T00:00:00Z".to_string(),
            HashMap::new(),
            HashMap::new(),
            vec![],
        );

        let read_bases = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let hash = read_content_hash(read_bases);
        let sketch = phraya_core::types::sketch(read_bases, 21, 11);

        let mut read_sketches = HashMap::new();
        read_sketches.insert(hash, sketch.clone());

        let plan = PhrayaPlan {
            read_sketches,
            ..base_plan
        };

        assert_eq!(plan.read_sketches.len(), 1);
        assert_eq!(plan.read_sketches.get(&hash), Some(&sketch));
    }

    /// `PhrayaPlan::get_read_sketch(hash)` looks up a stored read sketch by content
    /// hash, mirroring the existing `get_sketch(sequence_id)` accessor for
    /// reference sketches — returns `None` for an unstored hash.
    #[test]
    fn issue_200_get_read_sketch_looks_up_by_hash() {
        let base_plan = PhrayaPlan::new(
            UseCase::ReadsWithRef,
            vec![],
            "2026-07-13T00:00:00Z".to_string(),
            HashMap::new(),
            HashMap::new(),
            vec![],
        );

        let read_bases = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let hash = read_content_hash(read_bases);
        let sketch = phraya_core::types::sketch(read_bases, 21, 11);

        let mut read_sketches = HashMap::new();
        read_sketches.insert(hash, sketch.clone());

        let plan = PhrayaPlan {
            read_sketches,
            ..base_plan
        };

        assert_eq!(plan.get_read_sketch(hash), Some(&sketch));
        assert_eq!(plan.get_read_sketch(hash.wrapping_add(1)), None);
    }

    /// A plan with stored read sketches must round-trip through `.phrayaplan`'s
    /// write/read (MessagePack + zstd) unchanged — hash keys and sketch values
    /// both preserved exactly.
    #[test]
    fn issue_200_read_sketches_round_trip_through_phrayaplan() {
        let base_plan = PhrayaPlan::new(
            UseCase::ReadsWithRef,
            vec![],
            "2026-07-13T00:00:00Z".to_string(),
            HashMap::new(),
            HashMap::new(),
            vec![],
        );

        let read_a = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let read_b = b"TGCATGCATGCATGCATGCATGCATGCATGCA";
        let hash_a = read_content_hash(read_a);
        let hash_b = read_content_hash(read_b);
        let sketch_a = phraya_core::types::sketch(read_a, 21, 11);
        let sketch_b = phraya_core::types::sketch(read_b, 21, 11);

        let mut read_sketches = HashMap::new();
        read_sketches.insert(hash_a, sketch_a.clone());
        read_sketches.insert(hash_b, sketch_b.clone());

        let plan = PhrayaPlan {
            read_sketches,
            ..base_plan
        };

        let temp = NamedTempFile::new().unwrap();
        write_plan(temp.path(), &plan).unwrap();
        let read_plan = read_plan(temp.path()).unwrap();

        assert_eq!(read_plan.read_sketches.len(), 2);
        assert_eq!(read_plan.get_read_sketch(hash_a), Some(&sketch_a));
        assert_eq!(read_plan.get_read_sketch(hash_b), Some(&sketch_b));
    }

    /// A plan with no stored read sketches (the common case for plans that
    /// predate this feature, or reference-only plans) must still round-trip —
    /// the field is additive with a sensible empty default, not mandatory.
    #[test]
    fn issue_200_plan_without_read_sketches_round_trips() {
        let plan = PhrayaPlan::new(
            UseCase::ReadsWithRef,
            vec![],
            "2026-07-13T00:00:00Z".to_string(),
            HashMap::new(),
            HashMap::new(),
            vec![],
        );

        let temp = NamedTempFile::new().unwrap();
        write_plan(temp.path(), &plan).unwrap();
        let read_plan = read_plan(temp.path()).unwrap();

        assert!(read_plan.read_sketches.is_empty());
    }

    // ============================================================================
    // Issue #201: chunk-addressable plan v7 format tests
    // ============================================================================

    /// v7 round-trip: an empty plan (no reads, no chunks) writes and reads as v7.
    #[test]
    fn issue_201_v7_round_trip_empty_plan() {
        let plan = PhrayaPlan::new(
            UseCase::ReadsWithRef,
            vec!["ref.fa".to_string()],
            "2026-07-23T00:00:00Z".to_string(),
            HashMap::new(),
            HashMap::new(),
            vec![(1, 0)],
        );

        let temp = NamedTempFile::new().unwrap();
        write_plan(temp.path(), &plan).unwrap();

        // Verify magic bytes
        let raw = std::fs::read(temp.path()).unwrap();
        assert_eq!(&raw[0..4], b"PHR7", "v7 file must start with PHR7 magic");

        let loaded = read_plan(temp.path()).unwrap();
        assert_eq!(loaded.version, PHRAYAPLAN_VERSION);
        assert_eq!(loaded.use_case, UseCase::ReadsWithRef);
        assert_eq!(loaded.task_list, vec![(1, 0)]);
    }

    /// v7 round-trip with read sketches (single chunk, N=1 default).
    #[test]
    fn issue_201_v7_round_trip_with_read_sketches() {
        let mut plan = PhrayaPlan::new(
            UseCase::ReadsWithRef,
            vec!["ref.fa".to_string(), "reads.fq".to_string()],
            "2026-07-23T00:00:00Z".to_string(),
            HashMap::new(),
            HashMap::new(),
            vec![],
        );

        let read_a = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let read_b = b"TGCATGCATGCATGCATGCATGCATGCATGCA";
        let hash_a = read_content_hash(read_a);
        let hash_b = read_content_hash(read_b);
        let sketch_a = phraya_core::types::sketch(read_a, 21, 11);
        let sketch_b = phraya_core::types::sketch(read_b, 21, 11);

        plan.read_sketches.insert(hash_a, sketch_a.clone());
        plan.read_sketches.insert(hash_b, sketch_b.clone());
        plan.read_hash_order = vec![hash_a, hash_b];
        plan.total_read_count = 2;

        let temp = NamedTempFile::new().unwrap();
        write_plan(temp.path(), &plan).unwrap();
        let loaded = read_plan(temp.path()).unwrap();

        assert_eq!(loaded.read_sketches.len(), 2);
        assert_eq!(loaded.get_read_sketch(hash_a), Some(&sketch_a));
        assert_eq!(loaded.get_read_sketch(hash_b), Some(&sketch_b));
    }

    /// v7 with pre-split chunks (N=4): each chunk gets its subset.
    #[test]
    fn issue_201_v7_pre_split_chunks() {
        let mut plan = PhrayaPlan::new(
            UseCase::ReadsWithRef,
            vec!["ref.fa".to_string()],
            "2026-07-23T00:00:00Z".to_string(),
            HashMap::new(),
            HashMap::new(),
            vec![],
        );

        // Insert 8 read sketches
        let mut hashes = Vec::new();
        for i in 0u8..8 {
            let bases: Vec<u8> = test_dna_sequence(i as u64, 32);
            let hash = read_content_hash(&bases);
            let sketch = phraya_core::types::sketch(&bases, 21, 11);
            plan.read_sketches.insert(hash, sketch);
            hashes.push(hash);
        }
        plan.read_hash_order = hashes.clone();
        plan.total_read_count = 8;
        plan.batch_num_chunks = Some(4);

        let temp = NamedTempFile::new().unwrap();
        write_plan(temp.path(), &plan).unwrap();

        // Verify TOC has 4 chunks
        let toc = read_plan_toc(temp.path()).unwrap();
        assert_eq!(toc.num_chunks, 4);
        assert_eq!(toc.chunk_frame_offsets.len(), 4);

        // Full read should get all 8 sketches
        let loaded = read_plan(temp.path()).unwrap();
        assert_eq!(loaded.read_sketches.len(), 8);
    }

    /// Worker load: worker 0 of 4 gets only its chunk's sketches.
    #[test]
    fn issue_201_v7_worker_loads_only_own_chunk() {
        let mut plan = PhrayaPlan::new(
            UseCase::ReadsWithRef,
            vec!["ref.fa".to_string()],
            "2026-07-23T00:00:00Z".to_string(),
            HashMap::new(),
            HashMap::new(),
            vec![],
        );

        // 8 reads, 4 chunks → 2 reads per chunk
        let mut hashes = Vec::new();
        for i in 0u64..8 {
            let bases: Vec<u8> = test_dna_sequence(i, 36);
            let hash = read_content_hash(&bases);
            let sketch = phraya_core::types::sketch(&bases, 21, 11);
            plan.read_sketches.insert(hash, sketch);
            hashes.push(hash);
        }
        plan.read_hash_order = hashes.clone();
        plan.total_read_count = 8;
        plan.batch_num_chunks = Some(4);

        let temp = NamedTempFile::new().unwrap();
        write_plan(temp.path(), &plan).unwrap();

        // Worker 0 should get chunk 0 (first 2 reads)
        let worker_plan = read_plan_worker(temp.path(), 0, 4).unwrap();
        assert_eq!(
            worker_plan.read_sketches.len(),
            2,
            "Worker 0 of 4 should get 2 of 8 reads, got {}",
            worker_plan.read_sketches.len()
        );

        // Worker 3 should also get 2 reads
        let worker_plan = read_plan_worker(temp.path(), 3, 4).unwrap();
        assert_eq!(
            worker_plan.read_sketches.len(),
            2,
            "Worker 3 of 4 should get 2 of 8 reads, got {}",
            worker_plan.read_sketches.len()
        );

        // Shared data (metadata) should be present in all workers
        assert_eq!(worker_plan.version, PHRAYAPLAN_VERSION);
        assert_eq!(worker_plan.use_case, UseCase::ReadsWithRef);
    }

    /// Worker isolation: no overlap between workers' sketch sets.
    #[test]
    fn issue_201_v7_worker_chunks_are_disjoint() {
        let mut plan = PhrayaPlan::new(
            UseCase::ReadsWithRef,
            vec!["ref.fa".to_string()],
            "2026-07-23T00:00:00Z".to_string(),
            HashMap::new(),
            HashMap::new(),
            vec![],
        );

        let mut hashes = Vec::new();
        for i in 0u64..12 {
            let bases: Vec<u8> = test_dna_sequence(i + 100, 32);
            let hash = read_content_hash(&bases);
            let sketch = phraya_core::types::sketch(&bases, 21, 11);
            plan.read_sketches.insert(hash, sketch);
            hashes.push(hash);
        }
        plan.read_hash_order = hashes;
        plan.total_read_count = 12;
        plan.batch_num_chunks = Some(3);

        let temp = NamedTempFile::new().unwrap();
        write_plan(temp.path(), &plan).unwrap();

        let w0 = read_plan_worker(temp.path(), 0, 3).unwrap();
        let w1 = read_plan_worker(temp.path(), 1, 3).unwrap();
        let w2 = read_plan_worker(temp.path(), 2, 3).unwrap();

        // Each worker gets 4 reads
        assert_eq!(w0.read_sketches.len(), 4);
        assert_eq!(w1.read_sketches.len(), 4);
        assert_eq!(w2.read_sketches.len(), 4);

        // No overlap
        let keys_0: std::collections::HashSet<u64> =
            w0.read_sketches.keys().copied().collect();
        let keys_1: std::collections::HashSet<u64> =
            w1.read_sketches.keys().copied().collect();
        let keys_2: std::collections::HashSet<u64> =
            w2.read_sketches.keys().copied().collect();

        assert!(
            keys_0.is_disjoint(&keys_1),
            "Worker 0 and 1 must have disjoint sketch sets"
        );
        assert!(
            keys_1.is_disjoint(&keys_2),
            "Worker 1 and 2 must have disjoint sketch sets"
        );
        assert!(
            keys_0.is_disjoint(&keys_2),
            "Worker 0 and 2 must have disjoint sketch sets"
        );
    }

    /// Fallback: N=1 plan with worker_count > 1 filters in-memory.
    #[test]
    fn issue_201_v7_fallback_n1_with_multiple_workers() {
        let mut plan = PhrayaPlan::new(
            UseCase::ReadsWithRef,
            vec!["ref.fa".to_string()],
            "2026-07-23T00:00:00Z".to_string(),
            HashMap::new(),
            HashMap::new(),
            vec![],
        );

        let mut hashes = Vec::new();
        for i in 0u64..6 {
            let bases: Vec<u8> = test_dna_sequence(i + 200, 32);
            let hash = read_content_hash(&bases);
            let sketch = phraya_core::types::sketch(&bases, 21, 11);
            plan.read_sketches.insert(hash, sketch);
            hashes.push(hash);
        }
        plan.read_hash_order = hashes;
        plan.total_read_count = 6;
        // N=1 (default, no --chunks)

        let temp = NamedTempFile::new().unwrap();
        write_plan(temp.path(), &plan).unwrap();

        // Full load gets all 6
        let full = read_plan(temp.path()).unwrap();
        assert_eq!(full.read_sketches.len(), 6);

        // Worker 0/3 gets ~2, worker 1/3 gets ~2, worker 2/3 gets ~2
        let w0 = read_plan_worker(temp.path(), 0, 3).unwrap();
        let w1 = read_plan_worker(temp.path(), 1, 3).unwrap();
        let w2 = read_plan_worker(temp.path(), 2, 3).unwrap();

        let total_loaded = w0.read_sketches.len() + w1.read_sketches.len() + w2.read_sketches.len();
        assert_eq!(
            total_loaded, 6,
            "All 3 workers combined must cover all 6 reads, got {}",
            total_loaded
        );
    }

    /// v6 file rejection: a plan written with old format should be hard-rejected.
    #[test]
    fn issue_201_v7_rejects_non_v7_file() {
        // Write a raw zstd-compressed msgpack blob (simulating a v6 file)
        let fake_v6_plan = PhrayaPlan::new(
            UseCase::ReadsWithRef,
            vec![],
            "2026-07-23T00:00:00Z".to_string(),
            HashMap::new(),
            HashMap::new(),
            vec![],
        );
        let serialized = rmp_serde::to_vec(&fake_v6_plan).unwrap();
        let compressed = zstd::encode_all(&serialized[..], 3).unwrap();

        let temp = NamedTempFile::new().unwrap();
        std::fs::write(temp.path(), &compressed).unwrap();

        let result = read_plan(temp.path());
        assert!(result.is_err(), "v6 format must be rejected");
        match result.unwrap_err() {
            PlanError::VersionMismatch { expected, got } => {
                assert_eq!(expected, 7);
                assert_eq!(got, 6);
            }
            other => panic!("Expected VersionMismatch, got {:?}", other),
        }
    }

    /// TOC inspection: read_plan_toc returns correct metadata.
    #[test]
    fn issue_201_v7_toc_inspection() {
        let mut plan = PhrayaPlan::new(
            UseCase::ReadsWithRef,
            vec!["ref.fa".to_string()],
            "2026-07-23T00:00:00Z".to_string(),
            HashMap::new(),
            HashMap::new(),
            vec![],
        );
        plan.batch_num_chunks = Some(5);
        plan.total_read_count = 10;
        for i in 0u64..10 {
            let bases: Vec<u8> = test_dna_sequence(i + 300, 32);
            let hash = read_content_hash(&bases);
            plan.read_sketches
                .insert(hash, phraya_core::types::sketch(&bases, 21, 11));
            plan.read_hash_order.push(hash);
        }

        let temp = NamedTempFile::new().unwrap();
        write_plan(temp.path(), &plan).unwrap();

        let toc = read_plan_toc(temp.path()).unwrap();
        assert_eq!(toc.version, 7);
        assert_eq!(toc.flags, 0);
        assert_eq!(toc.num_chunks, 5);
        assert_eq!(toc.chunk_frame_offsets.len(), 5);
        // Shared frame must come before chunk frames
        assert!(toc.shared_frame_offset < toc.chunk_frame_offsets[0].0);
    }

    /// Byte determinism: writing the same plan twice produces identical bytes.
    #[test]
    fn issue_201_v7_byte_deterministic() {
        let mut plan = PhrayaPlan::new(
            UseCase::ReadsWithRef,
            vec!["ref.fa".to_string()],
            "2026-07-23T00:00:00Z".to_string(),
            HashMap::new(),
            HashMap::new(),
            vec![],
        );
        for i in 0u64..4 {
            let bases: Vec<u8> = test_dna_sequence(i + 200, 32);
            let hash = read_content_hash(&bases);
            plan.read_sketches
                .insert(hash, phraya_core::types::sketch(&bases, 21, 11));
            plan.read_hash_order.push(hash);
        }
        plan.batch_num_chunks = Some(2);

        let temp1 = NamedTempFile::new().unwrap();
        let temp2 = NamedTempFile::new().unwrap();
        write_plan(temp1.path(), &plan).unwrap();
        write_plan(temp2.path(), &plan).unwrap();

        let bytes1 = std::fs::read(temp1.path()).unwrap();
        let bytes2 = std::fs::read(temp2.path()).unwrap();
        assert_eq!(bytes1, bytes2, "Two writes of the same plan must be byte-identical");
    }
}
