use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use thiserror::Error;

/// Query index for multi-mapping analysis
pub type QueryIndex = HashMap<String, Vec<(u32, f64)>>;

/// Queries file format errors
#[derive(Debug, Error, Serialize, Deserialize)]
pub enum QueriesError {
    #[error("serialization error: {0}")]
    SerializationError(String),
    #[error("decompression error: {0}")]
    DecompressionError(String),
    #[error("compression error: {0}")]
    CompressionError(String),
    #[error("io error: {0}")]
    IoError(String),
}

/// Write query index to compressed binary format
///
/// # Arguments
/// * `path` - output file path
/// * `index` - QueryIndex: HashMap<query_id, Vec<(position, score)>>
///
/// Note: Filters positions to score_ratio >= 0.95 (hard-coded opinion), then drops any
/// query whose filtered list is empty. A query therefore appears in the file iff it placed
/// at least one alignment at score >= 0.95 — matching the documented contract that
/// downstream counters (e.g. the benchmark harness) rely on. Keeping empty-list keys
/// previously inflated the aligned-read count to 100% (issue #193).
pub fn write_queries(path: &std::path::Path, index: &QueryIndex) -> Result<(), QueriesError> {
    const SCORE_THRESHOLD: f64 = 0.95;

    // Filter index to keep only high-confidence alignments, dropping reads left with none
    // (issue #193: empty-list keys otherwise inflate the aligned-read count to 100%).
    //
    // Serialize through a BTreeMap keyed by query id so output is deterministic run-to-run:
    // the in-memory QueryIndex is a HashMap whose iteration order is randomized per process,
    // which would otherwise make the `.phraya.queries` bytes vary between identical runs.
    // Within each query, positions are sorted by (position, score) for the same reason (the
    // alignment order depends on anchor-extension order, which we do not want leaking into
    // the serialized form).
    let filtered_index: std::collections::BTreeMap<String, Vec<(u32, f64)>> = index
        .iter()
        .filter_map(|(query_id, alignments)| {
            let mut filtered_alignments: Vec<(u32, f64)> = alignments
                .iter()
                .filter(|(_, score)| *score >= SCORE_THRESHOLD)
                .copied()
                .collect();
            if filtered_alignments.is_empty() {
                None
            } else {
                filtered_alignments.sort_by(|a, b| {
                    a.0.cmp(&b.0)
                        .then(b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal))
                });
                Some((query_id.clone(), filtered_alignments))
            }
        })
        .collect();

    let serialized = rmp_serde::to_vec(&filtered_index)
        .map_err(|e| QueriesError::SerializationError(e.to_string()))?;

    let compressed = zstd::encode_all(&serialized[..], 3)
        .map_err(|e| QueriesError::CompressionError(e.to_string()))?;

    std::fs::write(path, compressed).map_err(|e| QueriesError::IoError(e.to_string()))?;

    Ok(())
}

/// Read query index from compressed binary format
///
/// # Arguments
/// * `path` - input file path
///
/// # Returns
/// QueryIndex: HashMap<query_id, Vec<(position, score)>>
pub fn read_queries(path: &std::path::Path) -> Result<QueryIndex, QueriesError> {
    let compressed = std::fs::read(path).map_err(|e| QueriesError::IoError(e.to_string()))?;

    let decompressed = zstd::decode_all(&compressed[..])
        .map_err(|e| QueriesError::DecompressionError(e.to_string()))?;

    let index: QueryIndex = rmp_serde::from_slice(&decompressed)
        .map_err(|e| QueriesError::SerializationError(e.to_string()))?;

    Ok(index)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    #[test]
    fn round_trip_empty_index() {
        let index: QueryIndex = HashMap::new();

        let temp = NamedTempFile::new().unwrap();
        write_queries(temp.path(), &index).unwrap();
        let read_index = read_queries(temp.path()).unwrap();

        assert_eq!(read_index.len(), 0);
    }

    #[test]
    fn round_trip_single_query() {
        let mut index = HashMap::new();
        index.insert(
            "query1".to_string(),
            vec![(100u32, 0.98f64), (200u32, 0.95f64)],
        );

        let temp = NamedTempFile::new().unwrap();
        write_queries(temp.path(), &index).unwrap();
        let read_index = read_queries(temp.path()).unwrap();

        assert_eq!(read_index.len(), 1);
        let positions = &read_index["query1"];
        assert_eq!(positions.len(), 2);
        assert_eq!(positions[0], (100, 0.98));
        assert_eq!(positions[1], (200, 0.95));
    }

    /// The serialized `.queries` bytes must not depend on `HashMap` iteration order:
    /// two indices holding the same logical content, inserted in different orders, must
    /// serialize to identical bytes. This locks in the deterministic-output guarantee so a
    /// content hash can be used as a benchmark/regression gate.
    #[test]
    fn write_queries_is_order_independent() {
        // Enough distinct keys that HashMap iteration order almost certainly varies between
        // the two maps below (and between runs).
        let keys: Vec<String> = (0..64).map(|i| format!("read_{i}")).collect();

        let mut forward = HashMap::new();
        for (i, k) in keys.iter().enumerate() {
            forward.insert(k.clone(), vec![(i as u32 * 7, 0.99f64), (i as u32 * 7 + 3, 0.96f64)]);
        }
        let mut reverse = HashMap::new();
        for (i, k) in keys.iter().enumerate().rev() {
            // Same key -> same value, inserted in the opposite order and with the positions
            // provided in the opposite (score) order to also exercise the per-query sort.
            reverse.insert(k.clone(), vec![(i as u32 * 7 + 3, 0.96f64), (i as u32 * 7, 0.99f64)]);
        }

        let ta = NamedTempFile::new().unwrap();
        let tb = NamedTempFile::new().unwrap();
        write_queries(ta.path(), &forward).unwrap();
        write_queries(tb.path(), &reverse).unwrap();

        let bytes_a = std::fs::read(ta.path()).unwrap();
        let bytes_b = std::fs::read(tb.path()).unwrap();
        assert_eq!(
            bytes_a, bytes_b,
            "identical content in different insertion/position orders must serialize identically"
        );
    }

    #[test]
    fn round_trip_multiple_queries() {
        let mut index = HashMap::new();
        index.insert(
            "query1".to_string(),
            vec![(100u32, 0.98f64), (200u32, 0.95f64)],
        );
        index.insert("query2".to_string(), vec![(50u32, 0.99f64)]);
        index.insert(
            "query3".to_string(),
            vec![(150u32, 0.96f64), (250u32, 0.95f64), (350u32, 0.97f64)],
        );

        let temp = NamedTempFile::new().unwrap();
        write_queries(temp.path(), &index).unwrap();
        let read_index = read_queries(temp.path()).unwrap();

        assert_eq!(read_index.len(), 3);
        assert_eq!(read_index["query1"].len(), 2);
        assert_eq!(read_index["query2"].len(), 1);
        assert_eq!(read_index["query3"].len(), 3);
    }

    #[test]
    fn large_query_index() {
        let mut index = HashMap::new();

        for q in 0..10000 {
            let query_id = format!("query_{}", q);
            let mut positions = Vec::new();
            for p in 0..(q % 10 + 1) {
                positions.push(((q * 100 + p) as u32, 0.95f64 + (p as f64) * 0.001));
            }
            index.insert(query_id, positions);
        }

        let temp = NamedTempFile::new().unwrap();
        write_queries(temp.path(), &index).unwrap();
        let read_index = read_queries(temp.path()).unwrap();

        assert_eq!(read_index.len(), 10000);

        // Spot check a few entries
        assert_eq!(read_index["query_0"].len(), 1);
        assert_eq!(read_index["query_99"].len(), 10);
    }

    #[test]
    fn query_positions_preserved() {
        let mut index = HashMap::new();
        let positions = vec![
            (100u32, 0.99f64),
            (200u32, 0.98f64),
            (300u32, 0.97f64),
            (400u32, 0.96f64),
            (500u32, 0.95f64),
        ];
        index.insert("test_query".to_string(), positions.clone());

        let temp = NamedTempFile::new().unwrap();
        write_queries(temp.path(), &index).unwrap();
        let read_index = read_queries(temp.path()).unwrap();

        let read_positions = &read_index["test_query"];
        assert_eq!(read_positions, &positions);
    }

    #[test]
    fn empty_query_alignment_list_is_dropped() {
        // A read whose only alignments are sub-threshold (or which placed nothing) leaves an
        // empty filtered list; that key must NOT be written (issue #193 — it inflated the
        // aligned-read count). Only reads with >=1 placement appear in the file.
        let mut index = HashMap::new();
        index.insert("query_no_hits".to_string(), vec![]);
        index.insert("query_subthreshold".to_string(), vec![(10u32, 0.40f64)]);

        let temp = NamedTempFile::new().unwrap();
        write_queries(temp.path(), &index).unwrap();
        let read_index = read_queries(temp.path()).unwrap();

        assert_eq!(read_index.len(), 0, "empty and all-sub-threshold reads must be dropped");
    }

    #[test]
    fn subthreshold_positions_dropped_but_read_kept_if_any_placement() {
        // A read with one placement and one sub-threshold hit keeps the key, with only the
        // >=0.95 position retained.
        let mut index = HashMap::new();
        index.insert(
            "mixed".to_string(),
            vec![(100u32, 0.98f64), (200u32, 0.40f64)],
        );

        let temp = NamedTempFile::new().unwrap();
        write_queries(temp.path(), &index).unwrap();
        let read_index = read_queries(temp.path()).unwrap();

        assert_eq!(read_index.len(), 1);
        assert_eq!(read_index["mixed"], vec![(100u32, 0.98f64)]);
    }

    #[test]
    fn nonexistent_file_read() {
        let result = read_queries(std::path::Path::new("/nonexistent/path.queries"));
        assert!(result.is_err());
    }

    #[test]
    fn corrupted_file() {
        let temp = NamedTempFile::new().unwrap();
        std::fs::write(temp.path(), b"corrupted data").unwrap();

        let result = read_queries(temp.path());
        assert!(result.is_err());
    }

    #[test]
    fn compression_efficiency() {
        let mut index = HashMap::new();

        // Create repetitive data that compresses well
        for q in 0..1000 {
            let query_id = format!("query_{:04}", q);
            let mut positions = Vec::new();
            for p in 0..100 {
                positions.push(((p * 10) as u32, 0.95f64));
            }
            index.insert(query_id, positions);
        }

        let temp = NamedTempFile::new().unwrap();
        write_queries(temp.path(), &index).unwrap();

        let file_size = std::fs::metadata(temp.path()).unwrap().len();
        // With repetitive data, file should compress well
        assert!(file_size < 1_000_000); // Should be much smaller than 1MB
    }

    #[test]
    fn high_score_threshold_values() {
        let mut index = HashMap::new();
        let mut positions = Vec::new();

        // Test various score values around 0.95 threshold
        for i in 0..100 {
            let score = 0.95f64 + (i as f64) * 0.0001;
            positions.push(((i * 10) as u32, score));
        }

        index.insert("high_scores".to_string(), positions);

        let temp = NamedTempFile::new().unwrap();
        write_queries(temp.path(), &index).unwrap();
        let read_index = read_queries(temp.path()).unwrap();

        let read_positions = &read_index["high_scores"];
        assert_eq!(read_positions.len(), 100);

        // Verify precision is maintained
        for (i, (pos, score)) in read_positions.iter().enumerate() {
            assert_eq!(*pos, (i as u32 * 10));
            assert!((score - (0.95f64 + (i as f64) * 0.0001)).abs() < 1e-10);
        }
    }
}
