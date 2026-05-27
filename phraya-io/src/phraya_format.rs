// Mergeable .phraya binary format (Issue #19)
//
// This module implements the .phraya format's merge capabilities:
// - Metadata storage (sources, evidence hash, alignment params)
// - merge() function for combining multiple .phraya files
// - Overlapping observation handling with provenance
// - Merge correctness, idempotence, and commutativity
// - Performance requirements

use std::io::{self, Read, Write};

// Core types for .phraya format

#[derive(Clone, Debug, PartialEq)]
#[allow(dead_code)]
pub struct VariantObservation {
    pub position: u64,
    pub reference_allele: String,
    pub alternate_allele: String,
    pub confidence: f64,
    pub provenance: Option<Provenance>,
}

#[derive(Clone, Debug, PartialEq)]
#[allow(dead_code)]
pub struct Provenance {
    pub source: String,
    pub alignment_params: AlignmentParams,
}

#[derive(Clone, Debug, PartialEq)]
#[allow(dead_code)]
pub struct PhrayaMetadata {
    pub input_sources: Vec<String>,
    pub evidence_layer_hash: String,
    pub alignment_params: AlignmentParams,
}

#[derive(Clone, Debug, PartialEq)]
#[allow(dead_code)]
pub struct AlignmentParams {
    pub strategy: String,
    pub min_mapq: u32,
}

// Binary format constants
#[allow(dead_code)]
const PHRAYA_MAGIC: &[u8; 8] = b"PHRAYA\x01\x00";
#[allow(dead_code)]
const OBSERVATION_SEPARATOR: u8 = 0xFF;
#[allow(dead_code)]
const METADATA_MARKER: u8 = 0xAA;

// Writer for .phraya format
pub struct PhrayaWriter<W: Write> {
    writer: W,
    metadata: Option<PhrayaMetadata>,
    observations_written: bool,
}

impl<W: Write> PhrayaWriter<W> {
    pub fn new(writer: W) -> Self {
        PhrayaWriter {
            writer,
            metadata: None,
            observations_written: false,
        }
    }

    pub fn with_metadata(mut self, metadata: PhrayaMetadata) -> Self {
        self.metadata = Some(metadata);
        self
    }

    pub fn write_observation(&mut self, obs: &VariantObservation) -> Result<(), io::Error> {
        // Write magic + metadata on first observation
        if !self.observations_written {
            self.writer.write_all(PHRAYA_MAGIC)?;

            // Write metadata
            if let Some(ref meta) = self.metadata {
                self.writer.write_all(&[METADATA_MARKER])?;
                write_string(&mut self.writer, &meta.evidence_layer_hash)?;
                write_string(&mut self.writer, &meta.alignment_params.strategy)?;
                self.writer.write_all(&meta.alignment_params.min_mapq.to_le_bytes())?;

                // Write number of sources
                self.writer.write_all(&(meta.input_sources.len() as u32).to_le_bytes())?;
                for source in &meta.input_sources {
                    write_string(&mut self.writer, source)?;
                }
            }

            self.observations_written = true;
        }

        // Write observation separator
        self.writer.write_all(&[OBSERVATION_SEPARATOR])?;

        // Write position
        self.writer.write_all(&obs.position.to_le_bytes())?;

        // Write alleles
        write_string(&mut self.writer, &obs.reference_allele)?;
        write_string(&mut self.writer, &obs.alternate_allele)?;

        // Write confidence
        self.writer.write_all(&obs.confidence.to_le_bytes())?;

        // Write provenance
        if let Some(ref prov) = obs.provenance {
            self.writer.write_all(&[1u8])?; // provenance present
            write_string(&mut self.writer, &prov.source)?;
            write_string(&mut self.writer, &prov.alignment_params.strategy)?;
            self.writer.write_all(&prov.alignment_params.min_mapq.to_le_bytes())?;
        } else {
            self.writer.write_all(&[0u8])?; // no provenance
        }

        Ok(())
    }

    pub fn finish(mut self) -> Result<(), io::Error> {
        // Write magic + metadata even if no observations were written
        if !self.observations_written {
            self.writer.write_all(PHRAYA_MAGIC)?;

            // Write metadata
            if let Some(ref meta) = self.metadata {
                self.writer.write_all(&[METADATA_MARKER])?;
                write_string(&mut self.writer, &meta.evidence_layer_hash)?;
                write_string(&mut self.writer, &meta.alignment_params.strategy)?;
                self.writer.write_all(&meta.alignment_params.min_mapq.to_le_bytes())?;

                // Write number of sources
                self.writer.write_all(&(meta.input_sources.len() as u32).to_le_bytes())?;
                for source in &meta.input_sources {
                    write_string(&mut self.writer, source)?;
                }
            }
        }
        Ok(())
    }
}

// Reader for .phraya format
pub struct PhrayaReader {
    #[allow(dead_code)]
    buffer: Vec<u8>,
    metadata: PhrayaMetadata,
    observations: Vec<VariantObservation>,
}

impl PhrayaReader {
    pub fn new<R: Read>(mut reader: R) -> Result<Self, io::Error> {
        let mut buffer = Vec::new();
        reader.read_to_end(&mut buffer)?;

        let mut pos = 0;

        // Check magic
        if pos + 8 > buffer.len() || &buffer[pos..pos + 8] != PHRAYA_MAGIC {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "Invalid .phraya magic"));
        }
        pos += 8;

        // Read metadata
        if pos >= buffer.len() || buffer[pos] != METADATA_MARKER {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "Missing metadata marker"));
        }
        pos += 1;

        let (hash, new_pos) = read_string(&buffer, pos)?;
        pos = new_pos;

        let (strategy, new_pos) = read_string(&buffer, pos)?;
        pos = new_pos;

        if pos + 4 > buffer.len() {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "Truncated min_mapq"));
        }
        let min_mapq = u32::from_le_bytes([buffer[pos], buffer[pos+1], buffer[pos+2], buffer[pos+3]]);
        pos += 4;

        if pos + 4 > buffer.len() {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "Truncated source count"));
        }
        let source_count = u32::from_le_bytes([buffer[pos], buffer[pos+1], buffer[pos+2], buffer[pos+3]]) as usize;
        pos += 4;

        let mut input_sources = Vec::new();
        for _ in 0..source_count {
            let (source, new_pos) = read_string(&buffer, pos)?;
            input_sources.push(source);
            pos = new_pos;
        }

        let metadata = PhrayaMetadata {
            input_sources,
            evidence_layer_hash: hash,
            alignment_params: AlignmentParams { strategy, min_mapq },
        };

        // Read observations
        let mut observations = Vec::new();
        while pos < buffer.len() {
            if buffer[pos] != OBSERVATION_SEPARATOR {
                break;
            }
            pos += 1;

            // Read position
            if pos + 8 > buffer.len() {
                return Err(io::Error::new(io::ErrorKind::InvalidData, "Truncated position"));
            }
            let position = u64::from_le_bytes([
                buffer[pos], buffer[pos+1], buffer[pos+2], buffer[pos+3],
                buffer[pos+4], buffer[pos+5], buffer[pos+6], buffer[pos+7],
            ]);
            pos += 8;

            let (reference_allele, new_pos) = read_string(&buffer, pos)?;
            pos = new_pos;

            let (alternate_allele, new_pos) = read_string(&buffer, pos)?;
            pos = new_pos;

            if pos + 8 > buffer.len() {
                return Err(io::Error::new(io::ErrorKind::InvalidData, "Truncated confidence"));
            }
            let confidence = f64::from_le_bytes([
                buffer[pos], buffer[pos+1], buffer[pos+2], buffer[pos+3],
                buffer[pos+4], buffer[pos+5], buffer[pos+6], buffer[pos+7],
            ]);
            pos += 8;

            // Read provenance
            if pos >= buffer.len() {
                return Err(io::Error::new(io::ErrorKind::InvalidData, "Truncated provenance flag"));
            }
            let has_provenance = buffer[pos] != 0;
            pos += 1;

            let provenance = if has_provenance {
                let (source, new_pos) = read_string(&buffer, pos)?;
                pos = new_pos;

                let (strategy, new_pos) = read_string(&buffer, pos)?;
                pos = new_pos;

                if pos + 4 > buffer.len() {
                    return Err(io::Error::new(io::ErrorKind::InvalidData, "Truncated provenance mapq"));
                }
                let min_mapq = u32::from_le_bytes([buffer[pos], buffer[pos+1], buffer[pos+2], buffer[pos+3]]);
                pos += 4;

                Some(Provenance {
                    source,
                    alignment_params: AlignmentParams { strategy, min_mapq },
                })
            } else {
                None
            };

            observations.push(VariantObservation {
                position,
                reference_allele,
                alternate_allele,
                confidence,
                provenance,
            });
        }

        Ok(PhrayaReader {
            buffer,
            metadata,
            observations,
        })
    }
}

impl PhrayaReader {
    pub fn metadata(&self) -> &PhrayaMetadata {
        &self.metadata
    }

    pub fn observations(&self) -> impl Iterator<Item = VariantObservation> {
        self.observations.clone().into_iter()
    }

    pub fn source_alignment_params(&self) -> Vec<AlignmentParams> {
        let mut params = vec![self.metadata.alignment_params.clone()];
        for obs in &self.observations {
            if let Some(prov) = &obs.provenance {
                if !params.iter().any(|p| p == &prov.alignment_params) {
                    params.push(prov.alignment_params.clone());
                }
            }
        }
        params
    }
}

// Helper functions for serialization
fn write_string<W: Write>(writer: &mut W, s: &str) -> Result<(), io::Error> {
    let len = s.len() as u32;
    writer.write_all(&len.to_le_bytes())?;
    writer.write_all(s.as_bytes())?;
    Ok(())
}

fn read_string(buffer: &[u8], pos: usize) -> Result<(String, usize), io::Error> {
    if pos + 4 > buffer.len() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "Truncated string length"));
    }
    let len = u32::from_le_bytes([buffer[pos], buffer[pos+1], buffer[pos+2], buffer[pos+3]]) as usize;
    let new_pos = pos + 4;

    if new_pos + len > buffer.len() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "Truncated string data"));
    }

    let s = String::from_utf8(buffer[new_pos..new_pos + len].to_vec())
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Invalid UTF-8"))?;

    Ok((s, new_pos + len))
}

// Merge error type
#[derive(Debug, Clone)]
pub struct MergeError {
    pub message: String,
}

impl std::fmt::Display for MergeError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.message)
    }
}

impl std::error::Error for MergeError {}

// Merge function
pub fn merge(files: &[&[u8]]) -> Result<Vec<u8>, MergeError> {
    if files.is_empty() {
        return Ok(Vec::new());
    }

    // Read all files and validate compatibility
    let mut readers = Vec::new();
    let mut first_hash = None;

    for file in files {
        let reader = PhrayaReader::new(*file)
            .map_err(|e| MergeError {
                message: format!("Failed to read .phraya file: {}", e),
            })?;

        // Check evidence layer hash compatibility
        let hash = &reader.metadata().evidence_layer_hash;
        if let Some(ref fh) = first_hash {
            if fh != hash {
                return Err(MergeError {
                    message: "Cannot merge files with incompatible evidence layer hashes".to_string(),
                });
            }
        } else {
            first_hash = Some(hash.clone());
        }

        readers.push(reader);
    }

    // Collect all observations and track provenance
    let mut all_observations = Vec::new();
    let mut all_sources = std::collections::HashSet::new();
    let mut all_source_params = std::collections::HashMap::new();

    for reader in &readers {
        for source in &reader.metadata().input_sources {
            all_sources.insert(source.clone());
            all_source_params.insert(source.clone(), reader.metadata().alignment_params.clone());
        }

        for obs in reader.observations() {
            let mut enriched_obs = obs.clone();
            // If observation doesn't have provenance, add it from the reader's sources
            if enriched_obs.provenance.is_none() {
                if let Some(source) = reader.metadata().input_sources.first() {
                    let params = reader.metadata().alignment_params.clone();
                    enriched_obs.provenance = Some(Provenance {
                        source: source.clone(),
                        alignment_params: params.clone(),
                    });
                    all_source_params.insert(source.clone(), params);
                }
            }
            if let Some(prov) = &enriched_obs.provenance {
                all_source_params.insert(prov.source.clone(), prov.alignment_params.clone());
            }
            all_observations.push(enriched_obs);
        }
    }

    // Group observations by (position, ref, alt) and deduplicate by source
    // Using a HashMap for faster lookup than BTreeMap
    let mut merged_map: std::collections::HashMap<(u64, String, String), std::collections::HashMap<Option<String>, VariantObservation>> =
        std::collections::HashMap::new();

    for obs in all_observations {
        let key = (obs.position, obs.reference_allele.clone(), obs.alternate_allele.clone());
        let source_key = obs.provenance.as_ref().map(|p| p.source.clone());

        merged_map
            .entry(key)
            .or_default()
            .insert(source_key, obs);
    }

    // Flatten the nested map into final observations
    let mut final_observations = Vec::new();
    for (_, source_map) in merged_map {
        for (_, obs) in source_map {
            final_observations.push(obs);
        }
    }

    // Write merged file
    let mut buffer = Vec::new();

    // Create merged metadata
    let mut merged_sources: Vec<_> = all_sources.into_iter().collect();
    merged_sources.sort();

    let merged_metadata = PhrayaMetadata {
        input_sources: merged_sources,
        evidence_layer_hash: first_hash.unwrap_or_else(|| "default_hash".to_string()),
        alignment_params: readers
            .first()
            .map(|r| r.metadata().alignment_params.clone())
            .unwrap_or_else(|| AlignmentParams {
                strategy: "balanced".to_string(),
                min_mapq: 30,
            }),
    };

    let mut writer = PhrayaWriter::new(&mut buffer).with_metadata(merged_metadata);

    for obs in final_observations {
        writer.write_observation(&obs).map_err(|e| MergeError {
            message: format!("Failed to write observation: {}", e),
        })?;
    }

    writer.finish().map_err(|e| MergeError {
        message: format!("Failed to finish writing: {}", e),
    })?;

    Ok(buffer)
}

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

    // Helper functions

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

    fn align_all_samples(samples: &[&str], _reference: &str) -> Vec<u8> {
        // Simulate aligning all samples together in one run
        // Generate one observation per sample
        let mut all_sources = Vec::new();
        let mut all_observations = Vec::new();

        for sample in samples {
            all_sources.push(sample.to_string());
            // Generate observations for this sample based on its name hash
            let sample_hash = sample.chars().fold(0u64, |acc, c| acc.wrapping_mul(31).wrapping_add(c as u64));
            let obs = VariantObservation {
                position: 100 + sample_hash % 1000,
                reference_allele: "A".to_string(),
                alternate_allele: "G".to_string(),
                confidence: 0.90 + (sample_hash % 10) as f64 / 100.0,
                provenance: None,
            };
            all_observations.push(obs);
        }

        let metadata = PhrayaMetadata {
            input_sources: all_sources,
            evidence_layer_hash: "default_hash".to_string(),
            alignment_params: AlignmentParams {
                strategy: "balanced".to_string(),
                min_mapq: 30,
            },
        };

        let mut buffer = Vec::new();
        let mut writer = PhrayaWriter::new(&mut buffer).with_metadata(metadata);

        for obs in all_observations {
            writer.write_observation(&obs).unwrap();
        }
        writer.finish().unwrap();

        buffer
    }

    fn align_single_sample(sample: &str, _reference: &str) -> Vec<u8> {
        // Simulate aligning a single sample
        let metadata = PhrayaMetadata {
            input_sources: vec![sample.to_string()],
            evidence_layer_hash: "default_hash".to_string(),
            alignment_params: AlignmentParams {
                strategy: "balanced".to_string(),
                min_mapq: 30,
            },
        };

        let mut buffer = Vec::new();
        let mut writer = PhrayaWriter::new(&mut buffer).with_metadata(metadata);

        // Generate observations for this sample based on its name hash
        let sample_hash = sample.chars().fold(0u64, |acc, c| acc.wrapping_mul(31).wrapping_add(c as u64));
        let obs = VariantObservation {
            position: 100 + sample_hash % 1000,
            reference_allele: "A".to_string(),
            alternate_allele: "G".to_string(),
            confidence: 0.90 + (sample_hash % 10) as f64 / 100.0,
            provenance: None,
        };
        writer.write_observation(&obs).unwrap();
        writer.finish().unwrap();

        buffer
    }

    fn extract_observations(phraya_bytes: &[u8]) -> Vec<VariantObservation> {
        let reader = PhrayaReader::new(phraya_bytes).unwrap();
        reader.observations().collect()
    }
}

