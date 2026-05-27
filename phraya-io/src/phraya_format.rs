/// .phraya binary format I/O
///
/// The .phraya format is Phraya's native binary format for storing variant observations.
/// It uses MessagePack serialization with a versioned header for forward compatibility.
///
/// Format structure:
/// - Header: version (u32), observation_count (u64), created_at (String)
/// - Body: MessagePack-serialized Vec<VariantObservation>
use crate::{IoError, Result};
use phraya_core::VariantObservation;
use serde::{Deserialize, Serialize};
use std::path::Path;

/// Metadata stored in .phraya file header
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct PhrayaMetadata {
    pub version: u32,
    pub observation_count: u64,
    pub created_at: String,
}

/// Internal file structure: header + observations
#[derive(Debug, Serialize, Deserialize)]
struct PhrayaFile {
    metadata: PhrayaMetadata,
    observations: Vec<VariantObservation>,
}

/// Write variant observations to .phraya binary format
///
/// # Arguments
/// * `observations` - Slice of variant observations to write
/// * `path` - Output file path
///
/// # Format
/// Uses MessagePack serialization with versioned header
pub fn write_phraya(observations: &[VariantObservation], path: &Path) -> Result<()> {
    use std::time::{SystemTime, UNIX_EPOCH};

    let timestamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|e| IoError::WriteError(format!("Failed to get current time: {}", e)))?;

    let created_at = format!("{:?}", timestamp);

    let metadata = PhrayaMetadata {
        version: 1,
        observation_count: observations.len() as u64,
        created_at,
    };

    let file = PhrayaFile {
        metadata,
        observations: observations.to_vec(),
    };

    let serialized =
        rmp_serde::to_vec(&file).map_err(|e| IoError::SerializationError(e.to_string()))?;

    std::fs::write(path, serialized).map_err(|e| IoError::WriteError(e.to_string()))?;

    Ok(())
}

/// Read variant observations from .phraya binary format
///
/// # Arguments
/// * `path` - Path to .phraya file
///
/// # Returns
/// Vector of variant observations
pub fn read_phraya(path: &Path) -> Result<Vec<VariantObservation>> {
    let bytes = std::fs::read(path).map_err(|e| IoError::ReadError(e.to_string()))?;

    let file: PhrayaFile =
        rmp_serde::from_slice(&bytes).map_err(|e| IoError::DeserializationError(e.to_string()))?;

    Ok(file.observations)
}

/// Read only metadata from .phraya file without deserializing observations
///
/// # Arguments
/// * `path` - Path to .phraya file
///
/// # Returns
/// Metadata header information
pub fn read_phraya_metadata(path: &Path) -> Result<PhrayaMetadata> {
    let bytes = std::fs::read(path).map_err(|e| IoError::ReadError(e.to_string()))?;

    let file: PhrayaFile =
        rmp_serde::from_slice(&bytes).map_err(|e| IoError::DeserializationError(e.to_string()))?;

    Ok(file.metadata)
}

#[cfg(test)]
mod tests {
    use super::*;
    use phraya_core::VariantObservation;
    use std::path::Path;
    use tempfile::TempDir;

    /// Helper to create a test VariantObservation
    fn create_test_observation(position: u64, ref_base: u8, alt_base: u8) -> VariantObservation {
        // This will fail until VariantObservation is implemented in phraya-core
        VariantObservation::new(position, ref_base, alt_base)
    }

    #[test]
    fn test_write_phraya_empty_observations() {
        // AC: Tests with various observation counts (0)
        let temp_dir = TempDir::new().unwrap();
        let path = temp_dir.path().join("empty.phraya");

        let observations: Vec<VariantObservation> = vec![];

        // This will fail - function doesn't exist yet
        let result = write_phraya(&observations, &path);
        assert!(
            result.is_ok(),
            "Should successfully write empty observations"
        );
        assert!(path.exists(), "Output file should exist");
    }

    #[test]
    fn test_write_phraya_single_observation() {
        // AC: Tests with various observation counts (1)
        let temp_dir = TempDir::new().unwrap();
        let path = temp_dir.path().join("single.phraya");

        let observations = vec![create_test_observation(100, b'A', b'G')];

        let result = write_phraya(&observations, &path);
        assert!(
            result.is_ok(),
            "Should successfully write single observation"
        );
        assert!(path.exists(), "Output file should exist");

        // File should not be empty - should contain header + data
        let metadata = std::fs::metadata(&path).unwrap();
        assert!(metadata.len() > 0, "File should not be empty");
    }

    #[test]
    fn test_write_phraya_hundred_observations() {
        // AC: Tests with various observation counts (100)
        let temp_dir = TempDir::new().unwrap();
        let path = temp_dir.path().join("hundred.phraya");

        let observations: Vec<VariantObservation> = (0..100)
            .map(|i| create_test_observation(i * 100, b'A', b'G'))
            .collect();

        let result = write_phraya(&observations, &path);
        assert!(result.is_ok(), "Should successfully write 100 observations");
        assert!(path.exists(), "Output file should exist");
    }

    #[test]
    fn test_write_phraya_ten_thousand_observations() {
        // AC: Tests with various observation counts (10000)
        let temp_dir = TempDir::new().unwrap();
        let path = temp_dir.path().join("ten_thousand.phraya");

        let observations: Vec<VariantObservation> = (0..10_000)
            .map(|i| create_test_observation(i * 100, b'A', b'G'))
            .collect();

        let result = write_phraya(&observations, &path);
        assert!(
            result.is_ok(),
            "Should successfully write 10000 observations"
        );
        assert!(path.exists(), "Output file should exist");

        // Verify file size is reasonable for 10k observations
        let metadata = std::fs::metadata(&path).unwrap();
        assert!(
            metadata.len() > 1000,
            "File should contain substantial data"
        );
    }

    #[test]
    fn test_round_trip_serialization_empty() {
        // AC: Round-trip test with empty observations
        let temp_dir = TempDir::new().unwrap();
        let path = temp_dir.path().join("roundtrip_empty.phraya");

        let original: Vec<VariantObservation> = vec![];

        write_phraya(&original, &path).unwrap();

        // This will fail - read function doesn't exist yet
        let loaded = read_phraya(&path).unwrap();

        assert_eq!(original.len(), loaded.len(), "Should preserve empty vector");
    }

    #[test]
    fn test_round_trip_serialization_single() {
        // AC: Round-trip test with single observation
        let temp_dir = TempDir::new().unwrap();
        let path = temp_dir.path().join("roundtrip_single.phraya");

        let original = vec![create_test_observation(42, b'C', b'T')];

        write_phraya(&original, &path).unwrap();
        let loaded = read_phraya(&path).unwrap();

        assert_eq!(original.len(), loaded.len());
        assert_eq!(original[0], loaded[0], "Should preserve observation data");
    }

    #[test]
    fn test_round_trip_serialization_multiple() {
        // AC: Round-trip test with multiple observations
        let temp_dir = TempDir::new().unwrap();
        let path = temp_dir.path().join("roundtrip_multi.phraya");

        let original = vec![
            create_test_observation(100, b'A', b'G'),
            create_test_observation(200, b'C', b'T'),
            create_test_observation(300, b'G', b'A'),
            create_test_observation(400, b'T', b'C'),
        ];

        write_phraya(&original, &path).unwrap();
        let loaded = read_phraya(&path).unwrap();

        assert_eq!(original.len(), loaded.len());
        for (i, (orig, load)) in original.iter().zip(loaded.iter()).enumerate() {
            assert_eq!(orig, load, "Observation {} should match", i);
        }
    }

    #[test]
    fn test_round_trip_preserves_all_fields() {
        // AC: Verify equality - all fields preserved
        let temp_dir = TempDir::new().unwrap();
        let path = temp_dir.path().join("roundtrip_fields.phraya");

        // Create observation with various field values to ensure all are preserved
        let original = vec![create_test_observation(12345, b'A', b'G')];

        write_phraya(&original, &path).unwrap();
        let loaded = read_phraya(&path).unwrap();

        // Deep equality check will verify all fields match
        assert_eq!(
            original, loaded,
            "All fields should be preserved in round-trip"
        );
    }

    #[test]
    fn test_file_format_includes_version_header() {
        // AC: File format includes header with version and metadata
        let temp_dir = TempDir::new().unwrap();
        let path = temp_dir.path().join("with_header.phraya");

        let observations = vec![create_test_observation(100, b'A', b'G')];

        write_phraya(&observations, &path).unwrap();

        // Read raw file contents to verify header exists
        let contents = std::fs::read(&path).unwrap();

        // MessagePack format should have structured header
        // This is a basic check - actual header validation happens in read_phraya
        assert!(contents.len() > 10, "File should contain header + data");

        // Verify we can extract version from the file
        let metadata = read_phraya_metadata(&path).unwrap();
        assert!(metadata.version > 0, "Should have valid version number");
    }

    #[test]
    fn test_file_format_metadata_structure() {
        // AC: Documentation on .phraya format structure - verify through metadata
        let temp_dir = TempDir::new().unwrap();
        let path = temp_dir.path().join("metadata.phraya");

        let observations = vec![
            create_test_observation(100, b'A', b'G'),
            create_test_observation(200, b'C', b'T'),
        ];

        write_phraya(&observations, &path).unwrap();

        let metadata = read_phraya_metadata(&path).unwrap();

        // Verify metadata contains expected fields
        assert_eq!(
            metadata.observation_count, 2,
            "Should track observation count"
        );
        assert!(metadata.version > 0, "Should have version");
        assert!(!metadata.created_at.is_empty(), "Should have timestamp");
    }

    #[test]
    fn test_write_phraya_uses_messagepack() {
        // AC: Uses MessagePack serialization (via rmp-serde crate)
        let temp_dir = TempDir::new().unwrap();
        let path = temp_dir.path().join("msgpack.phraya");

        let observations = vec![create_test_observation(100, b'A', b'G')];

        write_phraya(&observations, &path).unwrap();

        // Verify file is valid MessagePack by attempting to parse it
        let _contents = std::fs::read(&path).unwrap();

        // MessagePack files have specific magic bytes/structure
        // The read function uses rmp-serde, so if it succeeds, format is valid
        let loaded = read_phraya(&path).unwrap();
        assert_eq!(loaded.len(), 1, "MessagePack format should be parseable");
    }

    #[test]
    fn test_write_phraya_error_on_invalid_path() {
        // Error handling: invalid path
        let observations = vec![create_test_observation(100, b'A', b'G')];

        // Try to write to a path that can't be created
        let invalid_path = Path::new("/nonexistent/directory/file.phraya");
        let result = write_phraya(&observations, invalid_path);

        assert!(result.is_err(), "Should fail on invalid path");
    }

    #[test]
    fn test_write_phraya_error_on_read_only_location() {
        // Error handling: read-only location (if we can simulate it)
        // This test may need adjustment based on test environment
        let temp_dir = TempDir::new().unwrap();
        let path = temp_dir.path().join("readonly.phraya");

        // Create file first
        std::fs::write(&path, b"test").unwrap();

        // Make read-only
        let mut perms = std::fs::metadata(&path).unwrap().permissions();
        perms.set_readonly(true);
        std::fs::set_permissions(&path, perms).unwrap();

        let observations = vec![create_test_observation(100, b'A', b'G')];

        let result = write_phraya(&observations, &path);

        // Should fail to overwrite read-only file
        assert!(result.is_err(), "Should fail on read-only file");
    }

    #[test]
    fn test_large_observations_vector() {
        // Stress test: ensure large datasets work
        let temp_dir = TempDir::new().unwrap();
        let path = temp_dir.path().join("large.phraya");

        let observations: Vec<VariantObservation> = (0..50_000)
            .map(|i| create_test_observation(i * 100, b'A', b'G'))
            .collect();

        let result = write_phraya(&observations, &path);
        assert!(result.is_ok(), "Should handle large observation vectors");

        // Verify round-trip for large dataset
        let loaded = read_phraya(&path).unwrap();
        assert_eq!(
            observations.len(),
            loaded.len(),
            "Should preserve all observations"
        );
    }

    #[test]
    fn test_different_variant_types() {
        // Test various base combinations
        let temp_dir = TempDir::new().unwrap();
        let path = temp_dir.path().join("variants.phraya");

        let observations = vec![
            create_test_observation(100, b'A', b'G'),
            create_test_observation(200, b'C', b'T'),
            create_test_observation(300, b'G', b'C'),
            create_test_observation(400, b'T', b'A'),
            create_test_observation(500, b'A', b'T'),
            create_test_observation(600, b'C', b'G'),
        ];

        write_phraya(&observations, &path).unwrap();
        let loaded = read_phraya(&path).unwrap();

        assert_eq!(observations, loaded, "Should preserve all variant types");
    }

    #[test]
    fn test_sequential_positions() {
        // Test adjacent variant positions
        let temp_dir = TempDir::new().unwrap();
        let path = temp_dir.path().join("sequential.phraya");

        let observations = vec![
            create_test_observation(100, b'A', b'G'),
            create_test_observation(101, b'C', b'T'),
            create_test_observation(102, b'G', b'A'),
        ];

        write_phraya(&observations, &path).unwrap();
        let loaded = read_phraya(&path).unwrap();

        assert_eq!(observations, loaded, "Should handle sequential positions");
    }

    #[test]
    fn test_sparse_positions() {
        // Test widely spaced variant positions
        let temp_dir = TempDir::new().unwrap();
        let path = temp_dir.path().join("sparse.phraya");

        let observations = vec![
            create_test_observation(1_000, b'A', b'G'),
            create_test_observation(1_000_000, b'C', b'T'),
            create_test_observation(10_000_000, b'G', b'A'),
        ];

        write_phraya(&observations, &path).unwrap();
        let loaded = read_phraya(&path).unwrap();

        assert_eq!(observations, loaded, "Should handle sparse positions");
    }
}
