use phraya_index::MinimimizerSketch;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;
use thiserror::Error;

/// PhrayaPlan format version for forward compatibility
pub const PHRAYAPLAN_VERSION: u32 = 1;

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

/// PhrayaPlan: read-only reference for alignment workers
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhrayaPlan {
    /// Format version
    pub version: u32,
    /// Detected use case
    pub use_case: UseCase,
    /// Input file paths
    pub input_files: Vec<String>,
    /// Timestamp (ISO8601)
    pub timestamp: String,
    /// K-mer sketches for all sequences
    pub kmer_index: Vec<MinimimizerSketch>,
    /// K-mer uniqueness: position → uniqueness score
    pub kmer_uniqueness: HashMap<u32, f64>,
    /// Task list: (query_id, target_id) pairs
    pub task_list: Vec<(u32, u32)>,
}

impl PhrayaPlan {
    /// Create a new plan
    pub fn new(
        use_case: UseCase,
        input_files: Vec<String>,
        timestamp: String,
        kmer_index: Vec<MinimimizerSketch>,
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

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    #[test]
    fn round_trip_empty_plan() {
        let plan = PhrayaPlan::new(
            UseCase::ReadsWithRef,
            vec![],
            "2026-05-31T12:00:00Z".to_string(),
            vec![],
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
            vec![],
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
            vec![],
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
            vec![],
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
            vec![],
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
            vec![],
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
                vec![],
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
}
