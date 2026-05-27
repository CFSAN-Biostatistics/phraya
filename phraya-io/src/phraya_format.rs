// Tests for mergeable .phraya binary format (Issue #19)
//
// This module tests the .phraya format's merge capabilities:
// - Metadata storage (sources, evidence hash, alignment params)
// - merge() function for combining multiple .phraya files
// - Overlapping observation handling with provenance
// - Merge correctness, idempotence, and commutativity
// - Performance requirements

#[cfg(test)]
mod tests {
    use super::*;

    // Test 1: .phraya format includes required metadata
    #[test]
    fn test_phraya_format_includes_metadata() {
        // Given: a variant observation with alignment parameters
        let observation = VariantObservation {
            position: 12345,
            reference_allele: "A".to_string(),
            alternate_allele: "G".to_string(),
            confidence: 0.95,
            provenance: None,
        };

        let metadata = PhrayaMetadata {
            input_sources: vec!["sample1.fasta".to_string()],
            evidence_layer_hash: "abc123def456".to_string(),
            alignment_params: AlignmentParams {
                strategy: "balanced".to_string(),
                min_mapq: 30,
            },
        };

        // When: writing to .phraya format
        let mut buffer: Vec<u8> = Vec::new();
        let result = PhrayaWriter::new(&mut buffer)
            .with_metadata(metadata.clone())
            .write_observation(&observation);

        assert!(result.is_ok(), "Failed to write .phraya format");

        // Then: reading back should recover the metadata
        let reader = PhrayaReader::new(&buffer[..]).unwrap();
        let recovered_metadata = reader.metadata();

        assert_eq!(recovered_metadata.input_sources, vec!["sample1.fasta"]);
        assert_eq!(recovered_metadata.evidence_layer_hash, "abc123def456");
        assert_eq!(recovered_metadata.alignment_params.strategy, "balanced");
        assert_eq!(recovered_metadata.alignment_params.min_mapq, 30);
    }

    // Test 2: merge() function combines multiple .phraya files
    #[test]
    fn test_merge_combines_multiple_files() {
        // Given: three separate .phraya files
        let file_a = create_phraya_file(
            vec![VariantObservation {
                position: 100,
                reference_allele: "A".to_string(),
                alternate_allele: "G".to_string(),
                confidence: 0.95,
                provenance: None,
            }],
            "sample_a.fasta",
        );

        let file_b = create_phraya_file(
            vec![VariantObservation {
                position: 200,
                reference_allele: "C".to_string(),
                alternate_allele: "T".to_string(),
                confidence: 0.92,
                provenance: None,
            }],
            "sample_b.fasta",
        );

        let file_c = create_phraya_file(
            vec![VariantObservation {
                position: 300,
                reference_allele: "G".to_string(),
                alternate_allele: "A".to_string(),
                confidence: 0.88,
                provenance: None,
            }],
            "sample_c.fasta",
        );

        // When: merging the files
        let merged = merge(&[&file_a[..], &file_b[..], &file_c[..]]).unwrap();

        // Then: merged file contains all observations
        let reader = PhrayaReader::new(&merged[..]).unwrap();
        let observations: Vec<_> = reader.observations().collect();

        assert_eq!(observations.len(), 3);
        assert!(observations.iter().any(|o| o.position == 100));
        assert!(observations.iter().any(|o| o.position == 200));
        assert!(observations.iter().any(|o| o.position == 300));
    }

    // Test 3: merge handles overlapping observations with provenance
    #[test]
    fn test_merge_handles_overlapping_observations() {
        // Given: two files with observations at the same position
        let file_a = create_phraya_file(
            vec![VariantObservation {
                position: 100,
                reference_allele: "A".to_string(),
                alternate_allele: "G".to_string(),
                confidence: 0.95,
                provenance: None,
            }],
            "sample_a.fasta",
        );

        let file_b = create_phraya_file(
            vec![VariantObservation {
                position: 100,
                reference_allele: "A".to_string(),
                alternate_allele: "G".to_string(),
                confidence: 0.88,
                provenance: None,
            }],
            "sample_b.fasta",
        );

        // When: merging files with overlapping observations
        let merged = merge(&[&file_a[..], &file_b[..]]).unwrap();

        // Then: both observations are kept with provenance markers
        let reader = PhrayaReader::new(&merged[..]).unwrap();
        let observations: Vec<_> = reader.observations().collect();

        // Should have 2 observations at position 100
        let pos_100_obs: Vec<_> = observations.iter().filter(|o| o.position == 100).collect();
        assert_eq!(
            pos_100_obs.len(),
            2,
            "Both overlapping observations should be kept"
        );

        // Each should have provenance information
        assert!(pos_100_obs[0].provenance.is_some());
        assert!(pos_100_obs[1].provenance.is_some());

        // Provenances should be different
        let provenances: Vec<_> = pos_100_obs
            .iter()
            .map(|o| &o.provenance.as_ref().unwrap().source)
            .collect();
        assert!(provenances.contains(&&"sample_a.fasta".to_string()));
        assert!(provenances.contains(&&"sample_b.fasta".to_string()));
    }

    // Test 4: merge correctness - align all-vs-all vs separate then merge
    #[test]
    fn test_merge_correctness_vs_monolithic_alignment() {
        // Given: three samples
        let samples = vec!["sampleA.fasta", "sampleB.fasta", "sampleC.fasta"];
        let reference = "reference.fasta";

        // When: aligning all three together (monolithic)
        let monolithic_result = align_all_samples(&samples, reference);

        // And: aligning each separately then merging
        let separate_a = align_single_sample("sampleA.fasta", reference);
        let separate_b = align_single_sample("sampleB.fasta", reference);
        let separate_c = align_single_sample("sampleC.fasta", reference);
        let merged_result = merge(&[&separate_a[..], &separate_b[..], &separate_c[..]]).unwrap();

        // Then: variant observations should be identical (modulo provenance)
        let mono_obs = extract_observations(&monolithic_result);
        let merged_obs = extract_observations(&merged_result);

        assert_eq!(
            mono_obs.len(),
            merged_obs.len(),
            "Monolithic and merged approaches should find same variants"
        );

        for obs in &mono_obs {
            let matching = merged_obs.iter().find(|m| {
                m.position == obs.position
                    && m.reference_allele == obs.reference_allele
                    && m.alternate_allele == obs.alternate_allele
            });
            assert!(
                matching.is_some(),
                "Observation at {} should exist in merged result",
                obs.position
            );

            // Confidence should match within tolerance
            let merged_conf = matching.unwrap().confidence;
            assert!(
                (obs.confidence - merged_conf).abs() < 0.01,
                "Confidence mismatch: monolithic={}, merged={}",
                obs.confidence,
                merged_conf
            );
        }
    }

    // Test 5: idempotence - merge(A, A) == A
    #[test]
    fn test_merge_idempotence() {
        // Given: a .phraya file with multiple observations
        let file_a = create_phraya_file(
            vec![
                VariantObservation {
                    position: 100,
                    reference_allele: "A".to_string(),
                    alternate_allele: "G".to_string(),
                    confidence: 0.95,
                    provenance: None,
                },
                VariantObservation {
                    position: 200,
                    reference_allele: "C".to_string(),
                    alternate_allele: "T".to_string(),
                    confidence: 0.92,
                    provenance: None,
                },
            ],
            "sample_a.fasta",
        );

        // When: merging A with itself
        let merged = merge(&[&file_a[..], &file_a[..]]).unwrap();

        // Then: result should be equivalent to A (no duplication)
        let original_obs = extract_observations(&file_a);
        let merged_obs = extract_observations(&merged);

        // Should have same number of observations (deduplication)
        assert_eq!(
            original_obs.len(),
            merged_obs.len(),
            "merge(A, A) should have same observation count as A"
        );

        // Each observation in A should appear exactly once in merged result
        for obs in &original_obs {
            let matches: Vec<_> = merged_obs
                .iter()
                .filter(|m| {
                    m.position == obs.position
                        && m.reference_allele == obs.reference_allele
                        && m.alternate_allele == obs.alternate_allele
                })
                .collect();
            assert_eq!(
                matches.len(),
                1,
                "Observation at {} should appear exactly once after merge(A, A)",
                obs.position
            );
        }
    }

    // Test 6: commutativity - merge(A, B) == merge(B, A)
    #[test]
    fn test_merge_commutativity() {
        // Given: two different .phraya files
        let file_a = create_phraya_file(
            vec![
                VariantObservation {
                    position: 100,
                    reference_allele: "A".to_string(),
                    alternate_allele: "G".to_string(),
                    confidence: 0.95,
                    provenance: None,
                },
                VariantObservation {
                    position: 300,
                    reference_allele: "T".to_string(),
                    alternate_allele: "C".to_string(),
                    confidence: 0.90,
                    provenance: None,
                },
            ],
            "sample_a.fasta",
        );

        let file_b = create_phraya_file(
            vec![
                VariantObservation {
                    position: 200,
                    reference_allele: "C".to_string(),
                    alternate_allele: "T".to_string(),
                    confidence: 0.92,
                    provenance: None,
                },
                VariantObservation {
                    position: 400,
                    reference_allele: "G".to_string(),
                    alternate_allele: "A".to_string(),
                    confidence: 0.88,
                    provenance: None,
                },
            ],
            "sample_b.fasta",
        );

        // When: merging A then B vs B then A
        let ab = merge(&[&file_a[..], &file_b[..]]).unwrap();
        let ba = merge(&[&file_b[..], &file_a[..]]).unwrap();

        // Then: results should be equivalent
        let ab_obs = extract_observations(&ab);
        let ba_obs = extract_observations(&ba);

        assert_eq!(
            ab_obs.len(),
            ba_obs.len(),
            "merge(A,B) and merge(B,A) should have same observation count"
        );

        // Sort by position for comparison
        let mut ab_sorted = ab_obs.clone();
        ab_sorted.sort_by_key(|o| o.position);

        let mut ba_sorted = ba_obs.clone();
        ba_sorted.sort_by_key(|o| o.position);

        for (ab_obs, ba_obs) in ab_sorted.iter().zip(ba_sorted.iter()) {
            assert_eq!(ab_obs.position, ba_obs.position);
            assert_eq!(ab_obs.reference_allele, ba_obs.reference_allele);
            assert_eq!(ab_obs.alternate_allele, ba_obs.alternate_allele);
            assert!((ab_obs.confidence - ba_obs.confidence).abs() < 0.001);
        }
    }

    // Test 7: performance - merge 100 files with 1000 variants each in <5 seconds
    #[test]
    fn test_merge_performance() {
        use std::time::Instant;

        // Given: 100 .phraya files, each with 1000 variants
        let mut files = Vec::new();
        for file_idx in 0..100 {
            let mut observations = Vec::new();
            for var_idx in 0..1000 {
                observations.push(VariantObservation {
                    position: (file_idx * 10000) + var_idx,
                    reference_allele: "A".to_string(),
                    alternate_allele: "G".to_string(),
                    confidence: 0.90,
                    provenance: None,
                });
            }
            let file = create_phraya_file(observations, &format!("sample_{}.fasta", file_idx));
            files.push(file);
        }

        // When: merging all files
        let start = Instant::now();
        let file_refs: Vec<&[u8]> = files.iter().map(|f| f.as_slice()).collect();
        let merged = merge(&file_refs).unwrap();
        let elapsed = start.elapsed();

        // Then: should complete in less than 5 seconds
        assert!(
            elapsed.as_secs() < 5,
            "Merging 100 files with 1000 variants each took {:?}, expected <5s",
            elapsed
        );

        // Verify correctness: should have 100,000 observations
        let obs = extract_observations(&merged);
        assert_eq!(
            obs.len(),
            100_000,
            "Merged result should contain all 100,000 observations"
        );
    }

    // Test 8: merge preserves all metadata from constituent files
    #[test]
    fn test_merge_preserves_metadata() {
        // Given: files with different alignment parameters
        let file_a = create_phraya_file_with_params(
            vec![VariantObservation {
                position: 100,
                reference_allele: "A".to_string(),
                alternate_allele: "G".to_string(),
                confidence: 0.95,
                provenance: None,
            }],
            "sample_a.fasta",
            AlignmentParams {
                strategy: "exact".to_string(),
                min_mapq: 30,
            },
        );

        let file_b = create_phraya_file_with_params(
            vec![VariantObservation {
                position: 200,
                reference_allele: "C".to_string(),
                alternate_allele: "T".to_string(),
                confidence: 0.92,
                provenance: None,
            }],
            "sample_b.fasta",
            AlignmentParams {
                strategy: "balanced".to_string(),
                min_mapq: 20,
            },
        );

        // When: merging files with different parameters
        let merged = merge(&[&file_a[..], &file_b[..]]).unwrap();

        // Then: merged file should track all source parameters
        let reader = PhrayaReader::new(&merged[..]).unwrap();
        let metadata = reader.metadata();

        assert_eq!(metadata.input_sources.len(), 2);
        assert!(
            metadata
                .input_sources
                .contains(&"sample_a.fasta".to_string())
        );
        assert!(
            metadata
                .input_sources
                .contains(&"sample_b.fasta".to_string())
        );

        // Alignment params should be preserved per-source
        let source_params = reader.source_alignment_params();
        assert_eq!(source_params.len(), 2);
        assert!(source_params.iter().any(|p| p.strategy == "exact"));
        assert!(source_params.iter().any(|p| p.strategy == "balanced"));
    }

    // Test 9: empty file merge behavior
    #[test]
    fn test_merge_empty_files() {
        // Given: an empty .phraya file and a file with observations
        let empty_file = create_phraya_file(vec![], "empty.fasta");
        let non_empty = create_phraya_file(
            vec![VariantObservation {
                position: 100,
                reference_allele: "A".to_string(),
                alternate_allele: "G".to_string(),
                confidence: 0.95,
                provenance: None,
            }],
            "sample.fasta",
        );

        // When: merging empty with non-empty
        let merged = merge(&[&empty_file[..], &non_empty[..]]).unwrap();

        // Then: result should match non-empty file
        let merged_obs = extract_observations(&merged);
        let non_empty_obs = extract_observations(&non_empty);

        assert_eq!(merged_obs.len(), non_empty_obs.len());
        assert_eq!(merged_obs[0].position, 100);
    }

    // Test 10: merge with conflicting evidence layer hashes
    #[test]
    fn test_merge_conflicting_evidence_hashes() {
        // Given: files with different evidence layer hashes (incompatible)
        let file_a = create_phraya_file_with_hash(
            vec![VariantObservation {
                position: 100,
                reference_allele: "A".to_string(),
                alternate_allele: "G".to_string(),
                confidence: 0.95,
                provenance: None,
            }],
            "sample_a.fasta",
            "hash_v1",
        );

        let file_b = create_phraya_file_with_hash(
            vec![VariantObservation {
                position: 200,
                reference_allele: "C".to_string(),
                alternate_allele: "T".to_string(),
                confidence: 0.92,
                provenance: None,
            }],
            "sample_b.fasta",
            "hash_v2", // Different hash - incompatible format
        );

        // When/Then: merging should fail with clear error
        let result = merge(&[&file_a[..], &file_b[..]]);
        assert!(result.is_err());

        let error_msg = result.unwrap_err().to_string();
        assert!(
            error_msg.contains("incompatible evidence layer"),
            "Error should mention incompatible evidence layers"
        );
    }

    // Helper functions (these will fail to compile - that's expected for RED tests)

    fn create_phraya_file(observations: Vec<VariantObservation>, source: &str) -> Vec<u8> {
        let metadata = PhrayaMetadata {
            input_sources: vec![source.to_string()],
            evidence_layer_hash: "default_hash".to_string(),
            alignment_params: AlignmentParams {
                strategy: "balanced".to_string(),
                min_mapq: 30,
            },
        };

        let mut buffer = Vec::new();
        let mut writer = PhrayaWriter::new(&mut buffer).with_metadata(metadata);

        for obs in observations {
            writer.write_observation(&obs).unwrap();
        }
        writer.finish().unwrap();

        buffer
    }

    fn create_phraya_file_with_params(
        observations: Vec<VariantObservation>,
        source: &str,
        params: AlignmentParams,
    ) -> Vec<u8> {
        let metadata = PhrayaMetadata {
            input_sources: vec![source.to_string()],
            evidence_layer_hash: "default_hash".to_string(),
            alignment_params: params,
        };

        let mut buffer = Vec::new();
        let mut writer = PhrayaWriter::new(&mut buffer).with_metadata(metadata);

        for obs in observations {
            writer.write_observation(&obs).unwrap();
        }
        writer.finish().unwrap();

        buffer
    }

    fn create_phraya_file_with_hash(
        observations: Vec<VariantObservation>,
        source: &str,
        hash: &str,
    ) -> Vec<u8> {
        let metadata = PhrayaMetadata {
            input_sources: vec![source.to_string()],
            evidence_layer_hash: hash.to_string(),
            alignment_params: AlignmentParams {
                strategy: "balanced".to_string(),
                min_mapq: 30,
            },
        };

        let mut buffer = Vec::new();
        let mut writer = PhrayaWriter::new(&mut buffer).with_metadata(metadata);

        for obs in observations {
            writer.write_observation(&obs).unwrap();
        }
        writer.finish().unwrap();

        buffer
    }

    fn align_all_samples(_samples: &[&str], _reference: &str) -> Vec<u8> {
        // Placeholder that will fail - implementation agent will define this
        unimplemented!("align_all_samples - to be implemented")
    }

    fn align_single_sample(_sample: &str, _reference: &str) -> Vec<u8> {
        // Placeholder that will fail - implementation agent will define this
        unimplemented!("align_single_sample - to be implemented")
    }

    fn extract_observations(phraya_bytes: &[u8]) -> Vec<VariantObservation> {
        let reader = PhrayaReader::new(phraya_bytes).unwrap();
        reader.observations().collect()
    }
}

// Stub types that should fail to compile (implementations belong in main code)
// These stubs are here only to show what the tests expect

#[allow(dead_code)]
#[derive(Clone)]
struct VariantObservation {
    position: u64,
    reference_allele: String,
    alternate_allele: String,
    confidence: f64,
    provenance: Option<Provenance>,
}

#[allow(dead_code)]
#[derive(Clone)]
struct Provenance {
    source: String,
    alignment_params: AlignmentParams,
}

#[allow(dead_code)]
#[derive(Clone)]
struct PhrayaMetadata {
    input_sources: Vec<String>,
    evidence_layer_hash: String,
    alignment_params: AlignmentParams,
}

#[allow(dead_code)]
#[derive(Clone)]
struct AlignmentParams {
    strategy: String,
    min_mapq: u32,
}

#[allow(dead_code)]
struct PhrayaWriter<W> {
    writer: W,
}

#[allow(dead_code)]
impl<W> PhrayaWriter<W> {
    fn new(_writer: W) -> Self {
        unimplemented!()
    }

    fn with_metadata(self, _metadata: PhrayaMetadata) -> Self {
        unimplemented!()
    }

    fn write_observation(&mut self, _obs: &VariantObservation) -> Result<(), std::io::Error> {
        unimplemented!()
    }

    fn finish(self) -> Result<(), std::io::Error> {
        unimplemented!()
    }
}

#[allow(dead_code)]
struct PhrayaReader<R> {
    reader: R,
}

#[allow(dead_code)]
impl<R> PhrayaReader<R> {
    fn new(_reader: R) -> Result<Self, std::io::Error> {
        unimplemented!()
    }

    fn metadata(&self) -> &PhrayaMetadata {
        unimplemented!()
    }

    fn observations(&self) -> impl Iterator<Item = VariantObservation> {
        std::iter::empty()
    }

    fn source_alignment_params(&self) -> Vec<AlignmentParams> {
        unimplemented!()
    }
}

#[allow(dead_code)]
fn merge(_files: &[&[u8]]) -> Result<Vec<u8>, MergeError> {
    unimplemented!()
}

#[allow(dead_code)]
#[derive(Debug)]
struct MergeError {
    message: String,
}

impl std::fmt::Display for MergeError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.message)
    }
}

impl std::error::Error for MergeError {}
