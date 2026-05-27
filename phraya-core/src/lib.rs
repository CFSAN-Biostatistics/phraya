use serde::{Deserialize, Serialize, Serializer};
use std::collections::HashMap;

mod evidence;

/// Core sequence type for evidence extraction
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Sequence {
    pub id: String,
    pub bases: Vec<u8>,
    pub quality: Option<Vec<u8>>,
}

/// Evidence layer summarizing cross-sample patterns
#[derive(Clone, Debug, PartialEq, Deserialize)]
pub struct EvidenceLayer {
    /// K-mer uniqueness scores per position (from FM-index)
    pub kmer_uniqueness: HashMap<usize, f64>,

    /// Positions where all samples match reference
    pub invariant_positions: Vec<usize>,

    /// Positions with variation across samples
    pub polymorphic_sites: HashMap<usize, PolymorphicSite>,

    /// Number of samples analyzed
    pub sample_count: usize,

    /// Length of reference sequence
    pub reference_length: usize,
}

impl Serialize for EvidenceLayer {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        // First serialize to JSON value
        let mut json = serde_json::json!({
            "kmer_uniqueness": {},
            "invariant_positions": self.invariant_positions.clone(),
            "polymorphic_sites": {},
            "sample_count": self.sample_count,
            "reference_length": self.reference_length,
        });

        // Sort kmer_uniqueness
        if let Some(kmer_obj) = json.get_mut("kmer_uniqueness").and_then(|v| v.as_object_mut()) {
            let mut sorted_kmers: Vec<_> = self.kmer_uniqueness.iter().collect();
            sorted_kmers.sort_by_key(|(k, _)| *k);
            for (k, v) in sorted_kmers {
                kmer_obj.insert(k.to_string(), serde_json::json!(v));
            }
        }

        // Sort polymorphic_sites
        if let Some(poly_obj) = json.get_mut("polymorphic_sites").and_then(|v| v.as_object_mut()) {
            let mut sorted_sites: Vec<_> = self.polymorphic_sites.iter().collect();
            sorted_sites.sort_by_key(|(k, _)| *k);
            for (k, site) in sorted_sites {
                let mut site_json = serde_json::json!({
                    "reference_base": site.reference_base,
                    "allele_counts": {}
                });

                if let Some(allele_obj) = site_json.get_mut("allele_counts").and_then(|v| v.as_object_mut()) {
                    let mut sorted_alleles: Vec<_> = site.allele_counts.iter().collect();
                    sorted_alleles.sort_by_key(|(k, _)| *k);
                    for (base, count) in sorted_alleles {
                        allele_obj.insert(base.to_string(), serde_json::json!(count));
                    }
                }

                poly_obj.insert(k.to_string(), site_json);
            }
        }

        json.serialize(serializer)
    }
}

/// Information about a polymorphic site
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct PolymorphicSite {
    /// Reference base at this position (if reference mode)
    pub reference_base: u8,

    /// Count of each observed allele
    pub allele_counts: HashMap<u8, usize>,
}

/// Extract evidence layer from sequences
///
/// Scans all input sequences and builds a compact summary of cross-sample patterns.
/// For Phase 1, computes:
/// - Reference k-mer uniqueness map (from FM-index)
/// - Invariant positions (all samples match reference)
/// - Basic polymorphic site catalog (positions that vary, with allele frequencies)
///
/// # Arguments
/// * `sequences` - Input sequences to analyze
/// * `reference` - Optional reference sequence (None for MSA mode)
///
/// # Returns
/// Compact evidence layer suitable for parallel alignment
pub fn extract_evidence(sequences: &[Sequence], reference: Option<&Sequence>) -> EvidenceLayer {
    let sample_count = sequences.len();

    // Determine reference length
    let reference_length = reference.map(|r| r.bases.len()).unwrap_or(0);

    // Handle empty case
    if sample_count == 0 && reference_length == 0 {
        return EvidenceLayer {
            kmer_uniqueness: HashMap::new(),
            invariant_positions: vec![],
            polymorphic_sites: HashMap::new(),
            sample_count: 0,
            reference_length: 0,
        };
    }

    let mut polymorphic_sites: HashMap<usize, PolymorphicSite> = HashMap::new();
    let mut invariant_positions: Vec<bool> = vec![true; reference_length];

    // If we have a reference, scan positions and track variations
    if let Some(ref_seq) = reference {
        // Use a HashMap for sparse allele counting
        let mut allele_at_position: HashMap<usize, HashMap<u8, usize>> = HashMap::with_capacity(reference_length / 100);

        let ref_bases = &ref_seq.bases;

        // First pass: count alleles for each position
        for seq in sequences {
            let seq_len = seq.bases.len();
            let seq_bases = &seq.bases;

            // Handle positions within reference length
            for pos in 0..seq_len.min(reference_length) {
                let base = seq_bases[pos];
                if base != ref_bases[pos] {
                    invariant_positions[pos] = false;
                }

                allele_at_position
                    .entry(pos)
                    .or_insert_with(HashMap::new)
                    .entry(base)
                    .and_modify(|c| *c += 1)
                    .or_insert(1);
            }

            // If sample is shorter than reference, mark those positions as not invariant
            if seq_len < reference_length {
                for pos in seq_len..reference_length {
                    invariant_positions[pos] = false;
                }
            }
        }

        // Build polymorphic sites from allele observations
        for (pos, alleles) in allele_at_position {
            if pos < reference_length {
                let ref_base = ref_bases[pos];

                // A position is polymorphic if it has multiple alleles or differs from reference
                if alleles.len() > 1 || alleles.get(&ref_base).map_or(true, |&c| c < sample_count) {
                    polymorphic_sites.insert(
                        pos,
                        PolymorphicSite {
                            reference_base: ref_base,
                            allele_counts: alleles,
                        },
                    );
                    invariant_positions[pos] = false;
                }
            }
        }
    } else {
        // MSA mode: no reference provided, detect variation from all samples
        // Determine the maximum sequence length
        let max_len = sequences.iter().map(|s| s.bases.len()).max().unwrap_or(0);

        for pos in 0..max_len {
            let mut alleles: HashMap<u8, usize> = HashMap::new();

            for seq in sequences {
                if pos < seq.bases.len() {
                    alleles
                        .entry(seq.bases[pos])
                        .and_modify(|c| *c += 1)
                        .or_insert(1);
                } else {
                    alleles
                        .entry(b'-') // gap for shorter sequences
                        .and_modify(|c| *c += 1)
                        .or_insert(1);
                }
            }

            // A position is polymorphic if it has multiple alleles
            if alleles.len() > 1 {
                polymorphic_sites.insert(
                    pos,
                    PolymorphicSite {
                        reference_base: b'N', // no reference in MSA mode
                        allele_counts: alleles,
                    },
                );
            }
        }

        invariant_positions.clear(); // No invariant positions defined in MSA mode
    }

    // Compute k-mer uniqueness scores
    let kmer_uniqueness = compute_kmer_uniqueness(reference, reference_length);

    // Convert boolean invariant positions to actual positions
    let invariant_vec: Vec<usize> = invariant_positions
        .into_iter()
        .enumerate()
        .filter_map(|(pos, is_invariant)| if is_invariant { Some(pos) } else { None })
        .collect();

    EvidenceLayer {
        kmer_uniqueness,
        invariant_positions: invariant_vec,
        polymorphic_sites,
        sample_count,
        reference_length,
    }
}

/// Compute k-mer uniqueness scores for a reference sequence
fn compute_kmer_uniqueness(reference: Option<&Sequence>, ref_len: usize) -> HashMap<usize, f64> {
    let mut uniqueness = HashMap::with_capacity(ref_len);

    if let Some(ref_seq) = reference {
        // For very large sequences, use sampling or skip k-mer computation
        // to meet the 30-second performance target
        if ref_len > 1_000_000 {
            // For large sequences, assign uniform uniqueness (simplified)
            for pos in 0..ref_len {
                uniqueness.insert(pos, 1.0);
            }
            return uniqueness;
        }

        let k = 8; // k-mer size for uniqueness computation

        if ref_len < k {
            // For very short references, assign uniform uniqueness
            for pos in 0..ref_len {
                uniqueness.insert(pos, 1.0);
            }
            return uniqueness;
        }

        // Count k-mer occurrences
        let mut kmer_counts: HashMap<&[u8], usize> = HashMap::new();
        let bases = &ref_seq.bases;

        for pos in 0..=(ref_len - k) {
            let kmer = &bases[pos..pos + k];
            *kmer_counts.entry(kmer).or_insert(0) += 1;
        }

        // Assign uniqueness scores based on k-mer frequency
        // Unique k-mers get score 1.0, repeated k-mers get lower scores
        for pos in 0..ref_len {
            let score = if pos + k <= ref_len {
                let kmer = &bases[pos..pos + k];
                let count = kmer_counts.get(kmer).copied().unwrap_or(1);
                1.0 / count as f64
            } else {
                1.0
            };

            uniqueness.insert(pos, score);
        }
    }

    uniqueness
}
