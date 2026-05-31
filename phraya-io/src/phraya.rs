use serde::{Deserialize, Serialize};
use phraya_core::types::{VariantObservation, CoverageTrack};
use thiserror::Error;

/// .phraya format version
pub const PHRAYA_VERSION: u32 = 1;

/// Phraya file format errors
#[derive(Debug, Error, Serialize, Deserialize)]
pub enum PhrayaError {
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

/// Header for .phraya file
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhrayaHeader {
    pub version: u32,
    pub reference_length: u32,
    pub sample_id: String,
    pub timestamp: String,
    pub observation_count: usize,
}

/// Phraya file: alignment results for a single sample
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhrayaFile {
    pub header: PhrayaHeader,
    pub observations: Vec<VariantObservation>,
    pub coverage_track: CoverageTrack,
}

impl PhrayaFile {
    /// Create a new Phraya file
    pub fn new(
        reference_length: u32,
        sample_id: String,
        timestamp: String,
        observations: Vec<VariantObservation>,
        coverage_track: CoverageTrack,
    ) -> Self {
        let observation_count = observations.len();
        let header = PhrayaHeader {
            version: PHRAYA_VERSION,
            reference_length,
            sample_id,
            timestamp,
            observation_count,
        };

        PhrayaFile {
            header,
            observations,
            coverage_track,
        }
    }
}

/// Write PhrayaFile to compressed binary format
pub fn write_phraya(path: &std::path::Path, file: &PhrayaFile) -> Result<(), PhrayaError> {
    let serialized = rmp_serde::to_vec(file)
        .map_err(|e| PhrayaError::SerializationError(e.to_string()))?;

    let compressed = zstd::encode_all(&serialized[..], 3)
        .map_err(|e| PhrayaError::CompressionError(e.to_string()))?;

    std::fs::write(path, compressed)
        .map_err(|e| PhrayaError::IoError(e.to_string()))?;

    Ok(())
}

/// Read PhrayaFile from compressed binary format
pub fn read_phraya(path: &std::path::Path) -> Result<PhrayaFile, PhrayaError> {
    let compressed = std::fs::read(path)
        .map_err(|e| PhrayaError::IoError(e.to_string()))?;

    let decompressed = zstd::decode_all(&compressed[..])
        .map_err(|e| PhrayaError::DecompressionError(e.to_string()))?;

    let file: PhrayaFile = rmp_serde::from_slice(&decompressed)
        .map_err(|e| PhrayaError::SerializationError(e.to_string()))?;

    if file.header.version != PHRAYA_VERSION {
        return Err(PhrayaError::VersionMismatch {
            expected: PHRAYA_VERSION,
            got: file.header.version,
        });
    }

    Ok(file)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;
    use std::collections::HashMap;

    #[test]
    fn round_trip_empty_observations() {
        let coverage = CoverageTrack::new(vec![10, 10, 5, 5]);
        let file = PhrayaFile::new(
            100,
            "sample1".to_string(),
            "2026-05-31T12:00:00Z".to_string(),
            vec![],
            coverage,
        );

        let temp = NamedTempFile::new().unwrap();
        write_phraya(temp.path(), &file).unwrap();
        let read_file = read_phraya(temp.path()).unwrap();

        assert_eq!(read_file.header.sample_id, "sample1");
        assert_eq!(read_file.observations.len(), 0);
    }

    #[test]
    fn round_trip_with_observations() {
        let mut alleles = HashMap::new();
        alleles.insert(b'A', 10);
        alleles.insert(b'T', 5);

        let obs = VariantObservation::new(
            100, b'A', alleles, 0.95, "10M".to_string(), 60, 2,
            vec![10, 12, 15, 18, 20], 35.5, "sample1:read42".to_string(),
        );

        let coverage = CoverageTrack::new(vec![10; 200]);
        let file = PhrayaFile::new(
            200,
            "sample1".to_string(),
            "2026-05-31T12:00:00Z".to_string(),
            vec![obs],
            coverage,
        );

        let temp = NamedTempFile::new().unwrap();
        write_phraya(temp.path(), &file).unwrap();
        let read_file = read_phraya(temp.path()).unwrap();

        assert_eq!(read_file.observations.len(), 1);
        assert_eq!(read_file.observations[0].position(), 100);
        assert_eq!(read_file.observations[0].mapq(), 60);
    }

    #[test]
    fn large_observation_list() {
        let mut observations = Vec::new();
        for i in 0..10000 {
            let mut alleles = HashMap::new();
            alleles.insert(b'A', (i % 100) as u32 + 1);

            let obs = VariantObservation::new(
                i as u32, b'A', alleles, 0.95, format!("{}M", 10 + (i % 5)),
                (i % 60) as u8, 0, vec![10], 35.0, format!("sample:read{}", i),
            );
            observations.push(obs);
        }

        let coverage = CoverageTrack::new(vec![10; 10000]);
        let file = PhrayaFile::new(
            10000,
            "large_sample".to_string(),
            "2026-05-31T12:00:00Z".to_string(),
            observations,
            coverage,
        );

        let temp = NamedTempFile::new().unwrap();
        write_phraya(temp.path(), &file).unwrap();
        let read_file = read_phraya(temp.path()).unwrap();

        assert_eq!(read_file.observations.len(), 10000);
        assert_eq!(read_file.header.observation_count, 10000);
    }

    #[test]
    fn variant_observation_fields_preserved() {
        let mut alleles = HashMap::new();
        alleles.insert(b'G', 15);
        alleles.insert(b'A', 3);

        let obs = VariantObservation::new(
            300, b'G',
            alleles.clone(),
            0.99,
            "25M".to_string(),
            60, 0,
            vec![15, 16, 18, 20],
            42.0,
            "sample4:read5".to_string(),
        );

        let coverage = CoverageTrack::new(vec![15; 300]);
        let file = PhrayaFile::new(
            300,
            "preserve_test".to_string(),
            "2026-05-31T12:00:00Z".to_string(),
            vec![obs],
            coverage,
        );

        let temp = NamedTempFile::new().unwrap();
        write_phraya(temp.path(), &file).unwrap();
        let read_file = read_phraya(temp.path()).unwrap();

        let read_obs = &read_file.observations[0];
        assert_eq!(read_obs.position(), 300);
        assert_eq!(read_obs.ref_base(), b'G');
        assert_eq!(read_obs.confidence(), 0.99);
        assert_eq!(read_obs.cigar(), "25M");
        assert_eq!(read_obs.mapq(), 60);
        assert_eq!(read_obs.edit_distance(), 0);
        assert_eq!(read_obs.avg_base_quality(), 42.0);
        assert_eq!(read_obs.provenance(), "sample4:read5");

        let all_alleles = read_obs.all_alleles();
        assert_eq!(all_alleles.get(&b'G'), Some(&15));
        assert_eq!(all_alleles.get(&b'A'), Some(&3));
    }

    #[test]
    fn version_mismatch() {
        let coverage = CoverageTrack::new(vec![10; 100]);
        let mut file = PhrayaFile::new(
            100,
            "version_test".to_string(),
            "2026-05-31T12:00:00Z".to_string(),
            vec![],
            coverage,
        );

        file.header.version = 999;

        let temp = NamedTempFile::new().unwrap();
        write_phraya(temp.path(), &file).unwrap();

        let result = read_phraya(temp.path());
        assert!(result.is_err());
        match result.unwrap_err() {
            PhrayaError::VersionMismatch { expected, got } => {
                assert_eq!(expected, PHRAYA_VERSION);
                assert_eq!(got, 999);
            }
            _ => panic!("Expected VersionMismatch"),
        }
    }

    #[test]
    fn compression_ratio() {
        let mut observations = Vec::new();
        for i in 0..5000 {
            let mut alleles = HashMap::new();
            alleles.insert(b'A', 50);

            let obs = VariantObservation::new(
                i as u32, b'A', alleles, 0.95, "50M".to_string(),
                60, 0, vec![50], 35.0, format!("sample:read{}", i),
            );
            observations.push(obs);
        }

        let coverage = CoverageTrack::new(vec![50; 5000]);
        let file = PhrayaFile::new(
            5000,
            "compression_test".to_string(),
            "2026-05-31T12:00:00Z".to_string(),
            observations,
            coverage,
        );

        let temp = NamedTempFile::new().unwrap();
        write_phraya(temp.path(), &file).unwrap();

        let file_size = std::fs::metadata(temp.path()).unwrap().len();
        // With repetitive data, file should be reasonably compressed
        assert!(file_size < 5_000_000); // Should be much smaller than 5MB for 5000 repetitive observations
    }

    #[test]
    fn nonexistent_file_read() {
        let result = read_phraya(std::path::Path::new("/nonexistent/path.phraya"));
        assert!(result.is_err());
    }

    #[test]
    fn corrupted_file() {
        let temp = NamedTempFile::new().unwrap();
        std::fs::write(temp.path(), b"corrupted data").unwrap();

        let result = read_phraya(temp.path());
        assert!(result.is_err());
    }

    #[test]
    fn header_consistency() {
        let coverage = CoverageTrack::new(vec![10; 100]);
        let file = PhrayaFile::new(
            100,
            "header_test".to_string(),
            "2026-05-31T12:00:00Z".to_string(),
            vec![],
            coverage,
        );

        let temp = NamedTempFile::new().unwrap();
        write_phraya(temp.path(), &file).unwrap();
        let read_file = read_phraya(temp.path()).unwrap();

        assert_eq!(read_file.header.version, PHRAYA_VERSION);
        assert_eq!(read_file.header.reference_length, 100);
        assert_eq!(read_file.header.sample_id, "header_test");
        assert_eq!(read_file.header.observation_count, 0);
    }
}
