use phraya_core::types::{CoverageTrack, VariantObservation};
use serde::{Deserialize, Serialize};
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
    let serialized =
        rmp_serde::to_vec(file).map_err(|e| PhrayaError::SerializationError(e.to_string()))?;

    let compressed = zstd::encode_all(&serialized[..], 3)
        .map_err(|e| PhrayaError::CompressionError(e.to_string()))?;

    std::fs::write(path, compressed).map_err(|e| PhrayaError::IoError(e.to_string()))?;

    Ok(())
}

/// Read PhrayaFile from compressed binary format
pub fn read_phraya(path: &std::path::Path) -> Result<PhrayaFile, PhrayaError> {
    let compressed = std::fs::read(path).map_err(|e| PhrayaError::IoError(e.to_string()))?;

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

/// Merge multiple .phraya files into a single file.
///
/// Combines observations from all inputs, grouping by position.
/// Preserves provenance for each observation.
/// Coverage tracks are summed element-wise.
/// Deduplicates identical observations.
/// Order-independent (merge is commutative).
pub fn merge_phraya_files(paths: &[&std::path::Path]) -> Result<PhrayaFile, PhrayaError> {
    if paths.is_empty() {
        return Err(PhrayaError::IoError("No files to merge".to_string()));
    }

    // Read all files
    let mut files = Vec::new();
    for path in paths {
        files.push(read_phraya(path)?);
    }

    // Verify all files have same reference length
    let ref_length = files[0].header.reference_length;
    for file in &files {
        if file.header.reference_length != ref_length {
            return Err(PhrayaError::IoError(
                "All files must have same reference length".to_string(),
            ));
        }
    }

    // Collect all observations and deduplicate
    let mut observations = Vec::new();
    let mut seen = std::collections::HashSet::new();

    for file in &files {
        for obs in &file.observations {
            // Create a key for deduplication: (position, alleles, provenance)
            let alleles_key = format!(
                "{:?}",
                obs.all_alleles()
                    .iter()
                    .map(|(k, v)| format!("{}:{}", *k as char, v))
                    .collect::<Vec<_>>()
            );
            let key = format!("{}:{}:{}", obs.position(), alleles_key, obs.provenance());

            if !seen.contains(&key) {
                observations.push(obs.clone());
                seen.insert(key);
            }
        }
    }

    // Sort observations by position
    observations.sort_by_key(|obs| obs.position());

    // Count observations per position (= reads supporting that variant after merge).
    // Update each observation's local_coverage so downstream filters see merged depth.
    // Also aggregate pair counts and insert stats so paired-end filters work post-merge.
    let mut obs_count: std::collections::HashMap<u32, u32> = std::collections::HashMap::new();
    let mut pair_totals: std::collections::HashMap<u32, (u32, u32)> = std::collections::HashMap::new();
    // (insert_size_sum, insert_size_count) keyed by position
    let mut insert_totals: std::collections::HashMap<u32, (i64, u32)> = std::collections::HashMap::new();
    // unmapped mate count keyed by position
    let mut unmapped_mate_totals: std::collections::HashMap<u32, u32> = std::collections::HashMap::new();
    for obs in &observations {
        *obs_count.entry(obs.position()).or_insert(0) += 1;
        let (total, proper) = pair_totals.entry(obs.position()).or_insert((0, 0));
        let (obs_total, obs_proper) = obs.pair_counts();
        *total += obs_total;
        *proper += obs_proper;
        let (isum, icount) = insert_totals.entry(obs.position()).or_insert((0, 0));
        let (obs_sum, obs_count_i) = obs.insert_stats();
        *isum += obs_sum;
        *icount += obs_count_i;
        *unmapped_mate_totals.entry(obs.position()).or_insert(0) += obs.unmapped_mate_count();
    }
    let observations: Vec<_> = observations
        .into_iter()
        .map(|obs| {
            let depth = *obs_count.get(&obs.position()).unwrap_or(&1);
            let (total_paired, proper_paired) = pair_totals.get(&obs.position()).copied().unwrap_or((0, 0));
            let (ins_sum, ins_count) = insert_totals.get(&obs.position()).copied().unwrap_or((0, 0));
            let unmapped_mates = unmapped_mate_totals.get(&obs.position()).copied().unwrap_or(0);
            phraya_core::types::VariantObservation::new(
                obs.position(),
                obs.ref_base(),
                obs.all_alleles().clone(),
                obs.confidence(),
                obs.cigar().to_string(),
                obs.mapq(),
                obs.edit_distance(),
                vec![depth],
                obs.avg_base_quality(),
                obs.provenance().to_string(),
            )
            .with_pair_counts(total_paired, proper_paired)
            .with_insert_stats(ins_sum, ins_count)
            .with_unmapped_mate_count(unmapped_mates)
        })
        .collect();

    // Merge coverage tracks
    let mut merged_coverage_vec = vec![0usize; ref_length as usize];
    for file in &files {
        let decompressed = file.coverage_track.decompress();
        for (i, &cov) in decompressed.iter().enumerate() {
            merged_coverage_vec[i] += cov as usize;
        }
    }

    let merged_coverage = CoverageTrack::new(merged_coverage_vec);

    // Create merged file header
    let merged_header = PhrayaHeader {
        version: PHRAYA_VERSION,
        reference_length: ref_length,
        sample_id: format!("merged_{}", files.len()),
        timestamp: chrono::Local::now().to_rfc3339(),
        observation_count: observations.len(),
    };

    Ok(PhrayaFile {
        header: merged_header,
        observations,
        coverage_track: merged_coverage,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;
    use tempfile::NamedTempFile;

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
            100,
            b'A',
            alleles,
            0.95,
            "10M".to_string(),
            60,
            2,
            vec![10, 12, 15, 18, 20],
            35.5,
            "sample1:read42".to_string(),
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
                i as u32,
                b'A',
                alleles,
                0.95,
                format!("{}M", 10 + (i % 5)),
                (i % 60) as u8,
                0,
                vec![10],
                35.0,
                format!("sample:read{}", i),
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
            300,
            b'G',
            alleles.clone(),
            0.99,
            "25M".to_string(),
            60,
            0,
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
                i as u32,
                b'A',
                alleles,
                0.95,
                "50M".to_string(),
                60,
                0,
                vec![50],
                35.0,
                format!("sample:read{}", i),
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

    #[test]
    fn merge_two_files_overlapping_positions() {
        let mut alleles1 = HashMap::new();
        alleles1.insert(b'A', 10);

        let obs1 = VariantObservation::new(
            50,
            b'A',
            alleles1,
            0.95,
            "10M".to_string(),
            60,
            0,
            vec![10],
            35.0,
            "sample1:read1".to_string(),
        );

        let mut alleles2 = HashMap::new();
        alleles2.insert(b'T', 5);

        let obs2 = VariantObservation::new(
            50,
            b'T',
            alleles2,
            0.90,
            "10M".to_string(),
            50,
            1,
            vec![5],
            30.0,
            "sample2:read1".to_string(),
        );

        let coverage1 = CoverageTrack::new(vec![10; 100]);
        let file1 = PhrayaFile::new(
            100,
            "sample1".to_string(),
            "2026-05-31T12:00:00Z".to_string(),
            vec![obs1],
            coverage1,
        );

        let coverage2 = CoverageTrack::new(vec![5; 100]);
        let file2 = PhrayaFile::new(
            100,
            "sample2".to_string(),
            "2026-05-31T12:00:00Z".to_string(),
            vec![obs2],
            coverage2,
        );

        let temp1 = NamedTempFile::new().unwrap();
        let temp2 = NamedTempFile::new().unwrap();
        write_phraya(temp1.path(), &file1).unwrap();
        write_phraya(temp2.path(), &file2).unwrap();

        let merged = merge_phraya_files(&[temp1.path(), temp2.path()]).unwrap();

        assert_eq!(merged.observations.len(), 2);
        assert_eq!(merged.header.reference_length, 100);
    }

    #[test]
    fn merge_two_files_disjoint_positions() {
        let mut alleles1 = HashMap::new();
        alleles1.insert(b'A', 10);

        let obs1 = VariantObservation::new(
            25,
            b'A',
            alleles1,
            0.95,
            "10M".to_string(),
            60,
            0,
            vec![10],
            35.0,
            "sample1:read1".to_string(),
        );

        let mut alleles2 = HashMap::new();
        alleles2.insert(b'T', 5);

        let obs2 = VariantObservation::new(
            75,
            b'T',
            alleles2,
            0.90,
            "10M".to_string(),
            50,
            1,
            vec![5],
            30.0,
            "sample2:read1".to_string(),
        );

        let coverage1 = CoverageTrack::new(vec![10; 100]);
        let file1 = PhrayaFile::new(
            100,
            "sample1".to_string(),
            "2026-05-31T12:00:00Z".to_string(),
            vec![obs1],
            coverage1,
        );

        let coverage2 = CoverageTrack::new(vec![5; 100]);
        let file2 = PhrayaFile::new(
            100,
            "sample2".to_string(),
            "2026-05-31T12:00:00Z".to_string(),
            vec![obs2],
            coverage2,
        );

        let temp1 = NamedTempFile::new().unwrap();
        let temp2 = NamedTempFile::new().unwrap();
        write_phraya(temp1.path(), &file1).unwrap();
        write_phraya(temp2.path(), &file2).unwrap();

        let merged = merge_phraya_files(&[temp1.path(), temp2.path()]).unwrap();

        assert_eq!(merged.observations.len(), 2);
        // Verify positions are sorted
        assert_eq!(merged.observations[0].position(), 25);
        assert_eq!(merged.observations[1].position(), 75);
    }

    #[test]
    fn merge_coverage_summing() {
        let coverage1 = CoverageTrack::new(vec![5, 5, 5, 5]);
        let file1 = PhrayaFile::new(
            4,
            "sample1".to_string(),
            "2026-05-31T12:00:00Z".to_string(),
            vec![],
            coverage1,
        );

        let coverage2 = CoverageTrack::new(vec![10, 10, 10, 10]);
        let file2 = PhrayaFile::new(
            4,
            "sample2".to_string(),
            "2026-05-31T12:00:00Z".to_string(),
            vec![],
            coverage2,
        );

        let temp1 = NamedTempFile::new().unwrap();
        let temp2 = NamedTempFile::new().unwrap();
        write_phraya(temp1.path(), &file1).unwrap();
        write_phraya(temp2.path(), &file2).unwrap();

        let merged = merge_phraya_files(&[temp1.path(), temp2.path()]).unwrap();

        let merged_coverage = merged.coverage_track.decompress();
        // Coverage should be summed (5 + 10 quantized → 15)
        for cov in merged_coverage {
            assert!(cov > 5); // Should be greater than individual values
        }
    }

    #[test]
    fn merge_commutativity() {
        let mut alleles = HashMap::new();
        alleles.insert(b'A', 10);

        let obs = VariantObservation::new(
            50,
            b'A',
            alleles,
            0.95,
            "10M".to_string(),
            60,
            0,
            vec![10],
            35.0,
            "sample:read".to_string(),
        );

        let coverage = CoverageTrack::new(vec![10; 100]);
        let file = PhrayaFile::new(
            100,
            "test".to_string(),
            "2026-05-31T12:00:00Z".to_string(),
            vec![obs],
            coverage,
        );

        let temp1 = NamedTempFile::new().unwrap();
        let temp2 = NamedTempFile::new().unwrap();
        write_phraya(temp1.path(), &file).unwrap();
        write_phraya(temp2.path(), &file).unwrap();

        let merged_12 = merge_phraya_files(&[temp1.path(), temp2.path()]).unwrap();
        let merged_21 = merge_phraya_files(&[temp2.path(), temp1.path()]).unwrap();

        // Observations should be identical (deduped)
        assert_eq!(merged_12.observations.len(), merged_21.observations.len());
    }

    #[test]
    fn merge_empty_files() {
        let coverage = CoverageTrack::new(vec![10; 100]);
        let file = PhrayaFile::new(
            100,
            "empty".to_string(),
            "2026-05-31T12:00:00Z".to_string(),
            vec![],
            coverage,
        );

        let temp1 = NamedTempFile::new().unwrap();
        let temp2 = NamedTempFile::new().unwrap();
        write_phraya(temp1.path(), &file).unwrap();
        write_phraya(temp2.path(), &file).unwrap();

        let merged = merge_phraya_files(&[temp1.path(), temp2.path()]).unwrap();

        assert_eq!(merged.observations.len(), 0);
    }

    #[test]
    fn merge_single_file() {
        let mut alleles = HashMap::new();
        alleles.insert(b'A', 10);

        let obs = VariantObservation::new(
            50,
            b'A',
            alleles,
            0.95,
            "10M".to_string(),
            60,
            0,
            vec![10],
            35.0,
            "sample:read".to_string(),
        );

        let coverage = CoverageTrack::new(vec![10; 100]);
        let file = PhrayaFile::new(
            100,
            "single".to_string(),
            "2026-05-31T12:00:00Z".to_string(),
            vec![obs],
            coverage,
        );

        let temp = NamedTempFile::new().unwrap();
        write_phraya(temp.path(), &file).unwrap();

        let merged = merge_phraya_files(&[temp.path()]).unwrap();

        assert_eq!(merged.observations.len(), 1);
        assert_eq!(merged.header.reference_length, 100);
    }

    #[test]
    fn merge_mismatched_reference_length_error() {
        let coverage1 = CoverageTrack::new(vec![10; 100]);
        let file1 = PhrayaFile::new(
            100,
            "sample1".to_string(),
            "2026-05-31T12:00:00Z".to_string(),
            vec![],
            coverage1,
        );

        let coverage2 = CoverageTrack::new(vec![5; 200]);
        let file2 = PhrayaFile::new(
            200,
            "sample2".to_string(),
            "2026-05-31T12:00:00Z".to_string(),
            vec![],
            coverage2,
        );

        let temp1 = NamedTempFile::new().unwrap();
        let temp2 = NamedTempFile::new().unwrap();
        write_phraya(temp1.path(), &file1).unwrap();
        write_phraya(temp2.path(), &file2).unwrap();

        let result = merge_phraya_files(&[temp1.path(), temp2.path()]);
        assert!(result.is_err());
    }
}
