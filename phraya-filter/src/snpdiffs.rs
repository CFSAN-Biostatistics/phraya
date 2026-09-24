//! CSP2-compatible `.snpdiffs` output formatter.
//!
//! Produces a text serialization of `.phraya` alignment data in the exact format
//! that `compileMUMmer.py` emits and `csp2_common.py`'s `parseSNPDiffs()` consumes.
//!
//! This is a pure serialization/interpretation layer — no filtering decisions are made
//! here; all variant observations that reach this function are emitted verbatim.

use phraya_core::cigar::CigarStats;
use phraya_core::types::{Strand, VariantObservation, VariantType};
use phraya_io::plan::content_hash_for_bytes;
use phraya_io::SequenceParser;
use std::collections::HashMap;
use std::path::Path;

/// Metadata about a FASTA file, computed for the snpdiffs `#` header.
#[derive(Debug, Clone)]
struct FastaMeta {
    contig_count: usize,
    total_bases: usize,
    n50: usize,
    n90: usize,
    l50: usize,
    l90: usize,
    sha256: String,
    /// (contig_name, length) pairs for per-contig lookups (e.g. Dist_to_Query_End).
    contig_lengths: Vec<(String, usize)>,
}

/// Compute FASTA metadata: contig count, total bases, N50/L90, L50/L90, SHA-256.
///
/// SHA-256 is computed over the **raw file bytes** (including FASTA headers and
/// whitespace), matching `compileMUMmer.py`'s `fasta_info()` which does
/// `hashlib.sha256(file.read()).hexdigest()` on the raw file.
fn compute_fasta_meta(path: &Path) -> Result<FastaMeta, String> {
    // Read raw file bytes for SHA256 (matching compileMUMmer.py)
    let raw_bytes = std::fs::read(path)
        .map_err(|e| format!("failed to read FASTA {}: {}", path.display(), e))?;
    let sha256 = content_hash_for_bytes(&raw_bytes);

    // Parse sequences for contig stats
    let mut parser = SequenceParser::from_path(path).map_err(|e| e.to_string())?;
    let mut contig_lengths: Vec<(String, usize)> = Vec::new();

    while let Some(result) = parser.next() {
        let seq = result.map_err(|e| e.to_string())?;
        let name = seq.id().to_string();
        let len = seq.len();
        contig_lengths.push((name, len));
    }

    let total: usize = contig_lengths.iter().map(|(_, l)| *l).sum();

    // N50/L50 and N90/L90 — sort descending, walk until cumulative >= fraction
    let mut sorted_desc: Vec<usize> = contig_lengths.iter().map(|(_, l)| *l).collect();
    sorted_desc.sort_unstable_by(|a, b| b.cmp(a));

    let (n50, l50) = n_stats(&sorted_desc, total, 0.5);
    let (n90, l90) = n_stats(&sorted_desc, total, 0.9);

    Ok(FastaMeta {
        contig_count: contig_lengths.len(),
        total_bases: total,
        n50,
        n90,
        l50,
        l90,
        sha256,
        contig_lengths,
    })
}

/// Look up a contig's length by name from FASTA metadata.
/// If not found and the FASTA has a single contig, returns that length.
/// Returns None if not found and multiple contigs exist.
fn lookup_contig_length(meta: &FastaMeta, name: &str) -> Option<usize> {
    if let Some(&(_, len)) = meta.contig_lengths.iter().find(|(n, _)| n == name) {
        return Some(len);
    }
    // Single-contig FASTA: use total_bases (which == the single contig length)
    if meta.contig_count == 1 {
        return Some(meta.total_bases);
    }
    None
}

/// Compute N-x and L-x statistics for a given fraction (e.g. 0.5 for N50/L50).
fn n_stats(sorted_desc_lengths: &[usize], total_bases: usize, fraction: f64) -> (usize, usize) {
    if total_bases == 0 || sorted_desc_lengths.is_empty() {
        return (0, 0);
    }
    let threshold = (total_bases as f64 * fraction).ceil() as usize;
    let mut cumsum = 0usize;
    for (i, &len) in sorted_desc_lengths.iter().enumerate() {
        cumsum += len;
        if cumsum >= threshold {
            return (len, i + 1);
        }
    }
    (0, 0)
}

/// Walk a CIGAR string forward to recover the alignment's reference start position.
///
/// Given a `VariantObservation`'s `position` (absolute reference coordinate, 0-indexed),
/// `query_position` (absolute query coordinate, 0-indexed), and `cigar` (full alignment
/// CIGAR), this computes `target_start` — where the alignment begins on the reference.
///
/// The algorithm walks the CIGAR tracking cumulative reference (`r_cum`) and query (`q_cum`)
/// consumption from the alignment start. At the CIGAR operation containing the variant:
/// - SNPs (`X` op): the variant is at some offset within the X run; both ref and query
///   advance in lockstep, so `r_cum_at_variant = r_cum + offset` where `offset = query_pos - q_cum`
/// - Deletions (`I` op, Phraya = target-extra): query_position is fixed; `r_cum` is read
///   at the I op boundary
/// - Insertions (`D` op, Phraya = query-extra): reference position is fixed; `r_cum` is
///   read at the D op boundary
///
/// `query_position` is 0-indexed from the start of the query sequence, and the alignment
/// always starts at query position 0 (CIGAR describes the full query from its beginning).
/// So `q_cum` tracks query position directly.
fn compute_alignment_start(cigar: &str, position: u32, query_position: u32) -> u32 {
    let mut q_cum = 0u32;
    let mut r_cum = 0u32;

    for (count, op) in phraya_core::cigar::parse_ops(cigar) {
        let count = count as u32;
        match op {
            'M' | 'X' => {
                // Both ref and query consumed. For X ops, check if the variant falls here.
                if op == 'X'
                    && query_position >= q_cum
                    && query_position < q_cum + count
                {
                    let offset = query_position - q_cum;
                    return position.saturating_sub(r_cum + offset);
                }
                q_cum += count;
                r_cum += count;
            }
            'I' => {
                // Reference-only (deletion in query). query_position is fixed at q_cum.
                if q_cum == query_position {
                    return position.saturating_sub(r_cum);
                }
                r_cum += count;
            }
            'D' => {
                // Query-only (insertion in query). query_position is fixed at q_cum.
                if q_cum == query_position {
                    return position.saturating_sub(r_cum);
                }
                q_cum += count;
            }
            _ => {}
        }
    }

    // Fallback: if we couldn't locate the variant in the CIGAR, return 0.
    0
}

/// Get the alternate allele with the highest count from `all_alleles`.
///
/// For SNPs: the alternate base (non-ref) with the most supporting observations.
/// For insertions: the inserted base(s) sorted by count desc.
/// For deletions: `.` (no query allele).
fn best_alt_allele(
    ref_base: u8,
    all_alleles: &HashMap<u8, u32>,
    variant_type: VariantType,
) -> String {
    match variant_type {
        VariantType::Snp => {
            let mut best_base: Option<(u8, u32)> = None;
            for (&allele, &count) in all_alleles {
                if allele != ref_base && allele != b'.' {
                    match best_base {
                        None => best_base = Some((allele, count)),
                        Some((_, best_count)) if count > best_count => {
                            best_base = Some((allele, count));
                        }
                        _ => {}
                    }
                }
            }
            match best_base {
                Some((base, _)) => (base as char).to_string(),
                None => ".".to_string(),
            }
        }
        VariantType::Insertion => {
            let mut bases: Vec<(u8, u32)> = all_alleles
                .iter()
                .filter(|(&b, _)| b != b'.')
                .map(|(&b, &c)| (b, c))
                .collect();
            bases.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
            bases.iter().map(|(b, _)| *b as char).collect()
        }
        VariantType::Deletion => ".".to_string(),
    }
}

/// Format the reference base for snpdiffs output.
///
/// For SNPs: the reference base character.
/// For insertions: `.` (reference has no base at an insertion).
/// For deletions: the first deleted base (stored as `ref_base` on the observation).
fn format_ref_base(ref_base: u8, variant_type: VariantType) -> String {
    match variant_type {
        VariantType::Snp => (ref_base as char).to_string(),
        VariantType::Insertion => ".".to_string(),
        VariantType::Deletion => (ref_base as char).to_string(),
    }
}

/// Compute the alignment-block identity from a CIGAR and edit distance.
///
/// Uses `CigarStats::identity(edit_distance) * 100.0` — the raw alignment identity
/// (fraction of query-aligned bases that are matches), NOT the kmer-weighted `confidence`
/// field. This matches MUMmer's `%ID` which is uniform within an alignment block.
///
/// ASSUMPTION (to confirm against CSP2): `Perc_Iden` in snpdiffs variant rows is the
/// block-level identity `100 × (1 - edit_distance / query_aligned_len)`, not the
/// kmer-uniqueness-penalized `confidence`.
fn obs_identity(obs: &VariantObservation) -> f64 {
    let stats = CigarStats::parse(obs.cigar());
    stats.identity(obs.edit_distance()) * 100.0
}

/// Compute the median of a mutable slice, or 0.0 if empty.
fn median_value(v: &mut [f64]) -> f64 {
    if v.is_empty() {
        0.0
    } else {
        v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let mid = v.len() / 2;
        if v.len() % 2 == 1 {
            v[mid]
        } else {
            (v[mid - 1] + v[mid]) / 2.0
        }
    }
}

/// Derive the query contig name from an observation's provenance.
///
/// Provenance is typically `"contig_name:read_id"` or just `"contig_name"`.
/// Falls back to `query_id` if provenance is empty, or `"query"` if both are absent.
fn derive_query_contig(obs: &VariantObservation, query_id: Option<&str>) -> String {
    let prov = obs.provenance();
    if let Some(contig) = prov.split(':').next() {
        if !contig.is_empty() {
            return contig.to_string();
        }
    }
    query_id.map(str::to_string).unwrap_or_else(|| "query".to_string())
}

/// Format a FASTA metadata section as tab-separated key:value pairs.
/// Matches `compileMUMmer.py`'s `query_string`/`reference_string` construction:
/// `Query_ID`, `Query_Assembly`, `Query_Contig_Count`, `Query_Assembly_Bases`,
/// `Query_N50`, `Query_N90`, `Query_L50`, `Query_L90`, `Query_SHA256`.
fn format_fasta_meta_section(
    meta: &Option<FastaMeta>,
    assembly_id: &str,
    assembly_path: &str,
    is_query: bool,
) -> String {
    let prefix = if is_query { "Query" } else { "Reference" };

    match meta {
        Some(m) => format!(
            "{p}_ID:{id}\t{p}_Assembly:{asm}\t{p}_Contig_Count:{cc}\t{p}_Assembly_Bases:{ab}\t{p}_N50:{n50}\t{p}_N90:{n90}\t{p}_L50:{l50}\t{p}_L90:{l90}\t{p}_SHA256:{sha}",
            p = prefix,
            id = assembly_id,
            asm = assembly_path,
            cc = m.contig_count,
            ab = m.total_bases,
            n50 = m.n50,
            n90 = m.n90,
            l50 = m.l50,
            l90 = m.l90,
            sha = m.sha256,
        ),
        None => format!(
            "{p}_ID:{id}\t{p}_Assembly:NA\t{p}_Contig_Count:NA\t{p}_Assembly_Bases:NA\t{p}_N50:NA\t{p}_N90:NA\t{p}_L50:NA\t{p}_L90:NA\t{p}_SHA256:NA",
            p = prefix,
            id = assembly_id,
        ),
    }
}

/// Summary metrics for the `#` header.
struct SummaryMetrics {
    snp_count: usize,
    indel_count: usize,
    median_identity: f64,
    median_aln_len: f64,
    percent_ref_aligned: f64,
    percent_query_aligned: f64,
}

fn compute_summary(
    observations: &[VariantObservation],
    covered_positions: Option<u32>,
    reference_length: u32,
    query_length: Option<u32>,
) -> SummaryMetrics {
    let snp_count = observations
        .iter()
        .filter(|obs| obs.variant_type() == VariantType::Snp)
        .count();
    let indel_count = observations
        .iter()
        .filter(|obs| {
            matches!(
                obs.variant_type(),
                VariantType::Insertion | VariantType::Deletion
            )
        })
        .count();

    let mut identities: Vec<f64> = observations.iter().map(obs_identity).collect();
    let mut aln_lens: Vec<f64> = observations
        .iter()
        .map(|obs| CigarStats::parse(obs.cigar()).query_aligned_len() as f64)
        .collect();

    let median_identity = median_value(&mut identities);
    let median_aln_len = median_value(&mut aln_lens);

    let percent_ref_aligned = if reference_length > 0 {
        let covered = covered_positions.unwrap_or(0) as f64;
        100.0 * covered / reference_length as f64
    } else {
        0.0
    };

    let percent_query_aligned = match query_length {
        Some(qlen) if qlen > 0 => {
            let covered = covered_positions.unwrap_or(0) as f64;
            100.0 * covered / qlen as f64
        }
        _ => 0.0,
    };

    SummaryMetrics {
        snp_count,
        indel_count,
        median_identity,
        median_aln_len,
        percent_ref_aligned,
        percent_query_aligned,
    }
}

/// Format the `#` header line (tab-separated key:value pairs).
///
/// Matches `compileMUMmer.py`'s header: `#` + `\t` + query meta + `\t` +
/// reference meta + `\t` + summary metrics. All fields are `key:value`.
fn format_header(
    query_meta: &Option<FastaMeta>,
    query_id: &str,
    query_fasta_path: &str,
    ref_meta: &Option<FastaMeta>,
    ref_id: &str,
    ref_fasta_path: &str,
    summary: &SummaryMetrics,
) -> String {
    let query_str = format_fasta_meta_section(query_meta, query_id, query_fasta_path, true);
    let ref_str = format_fasta_meta_section(ref_meta, ref_id, ref_fasta_path, false);

    let summary_str = format!(
        "SNPs:{}\tReference_Percent_Aligned:{:.2}\tQuery_Percent_Aligned:{:.2}\tMedian_Percent_Identity:{:.2}\tMedian_Alignment_Length:{:.2}\tKmer_Similarity:NA\tShared_Kmers:NA\tReference_Unique_Kmers:NA\tQuery_Unique_Kmers:NA\tReference_Breakpoints:0\tQuery_Breakpoints:0\tReference_Relocations:0\tQuery_Relocations:0\tReference_Translocations:0\tQuery_Translocations:0\tReference_Inversions:0\tQuery_Inversions:0\tReference_Insertions:0\tQuery_Insertions:0\tReference_Tandem:0\tQuery_Tandem:0\tIndels:{}\tInvalid:0\tgSNPs:NA\tgIndels:NA",
        summary.snp_count,
        summary.percent_ref_aligned,
        summary.percent_query_aligned,
        summary.median_identity,
        summary.median_aln_len,
        summary.indel_count,
    );

    format!("#\t{}\t{}\t{}", query_str, ref_str, summary_str)
}

/// Format a `##` BED row (one per reference contig).
///
/// Single tab-separated line: `##\t<11 fields>`
/// Columns: Ref_Contig, Ref_Start, Ref_End, Ref_Length, Ref_Aligned,
///          Query_Contig, Query_Start, Query_End, Query_Length, Query_Aligned, Perc_Iden
///
/// For the single-contig-vs-reference use case, one row covers the full contig span.
/// `Ref_Aligned` = covered positions (from CoverageTrack header).
/// `Perc_Iden` = median of per-variant identity scores (or 100.0 if no variants but covered).
fn format_bed_row(
    ref_contig: &str,
    ref_length: u32,
    covered_positions: Option<u32>,
    query_contig: &str,
    query_length: Option<u32>,
    observations: &[VariantObservation],
) -> String {
    let ref_aligned = covered_positions.unwrap_or(0);

    // For single contig vs reference, query_aligned ≈ ref_aligned.
    let query_aligned = ref_aligned;

    // Overall Perc_Iden: median identity across all variant observations.
    // Uses raw CigarStats identity (not kmer-weighted confidence).
    let perc_iden = if !observations.is_empty() {
        let mut idens: Vec<f64> = observations.iter().map(obs_identity).collect();
        median_value(&mut idens)
    } else if ref_aligned > 0 {
        100.0
    } else {
        0.0
    };

    let query_end = query_length.unwrap_or(0);
    let query_len = query_length.unwrap_or(0);

    // MUMmer convention: Ref_Start is 0-indexed, Ref_End is end position
    // (for full contig span: 0 to ref_length).
    format!(
        "##\t{}\t0\t{}\t{}\t{}\t{}\t0\t{}\t{}\t{}\t{:.1}",
        ref_contig,
        ref_length,       // Ref_End (full span end)
        ref_length,       // Ref_Length
        ref_aligned,      // Ref_Aligned
        query_contig,
        query_end,        // Query_End (full span end)
        query_len,        // Query_Length
        query_aligned,    // Query_Aligned
        perc_iden,        // Perc_Iden
    )
}

/// Format a single variant observation as one tab-separated line.
///
/// Columns (21, tab-separated, matching `csp2_common.py`'s `snp_columns`):
/// Ref_Contig, Start_Ref, Ref_Pos, Query_Contig, Start_Query, Query_Pos,
/// Ref_Loc, Query_Loc, Ref_Start, Ref_End, Query_Start, Query_End,
/// Ref_Base, Query_Base, Dist_to_Ref_End, Dist_to_Query_End,
/// Ref_Aligned, Query_Aligned, Query_Direction, Perc_Iden, Cat
///
/// NOTE: `Query_Contig` and `Query_Length` are derived **per-observation** from
/// `obs.provenance()` (which carries the query contig name), then looked up in
/// `query_meta`. This is essential for multi-contig query FASTAs where different
/// observations originate from different contigs.
fn format_variant_row(
    obs: &VariantObservation,
    ref_contig: &str,
    query_meta: &Option<FastaMeta>,
    query_id_fallback: Option<&str>,
    ref_length: u32,
) -> String {
    let position = obs.position();
    let ref_pos_1indexed = position + 1; // 1-indexed
    let start_ref = position; // 0-indexed (Ref_Pos - 1)
    let query_position = obs.query_position();
    let query_pos_1indexed = query_position + 1; // 1-indexed
    let start_query = query_position; // 0-indexed (Query_Pos - 1)

    let cigar_str = obs.cigar();
    let cigar_stats = CigarStats::parse(cigar_str);
    let ref_aligned = cigar_stats.target_aligned_len();
    let query_aligned = cigar_stats.query_aligned_len();

    // Compute alignment block start from CIGAR + positions.
    let alignment_start = compute_alignment_start(cigar_str, position, query_position);
    let ref_start_block = alignment_start;
    // Ref_End: alignment start + aligned reference length (1-indexed end convention)
    let ref_end_block = alignment_start.saturating_add(ref_aligned);
    let query_start_block = 0u32; // Query starts at position 0
    let query_end_block = query_aligned;

    let ref_base = format_ref_base(obs.ref_base(), obs.variant_type());
    let query_base = best_alt_allele(obs.ref_base(), obs.all_alleles(), obs.variant_type());

    // Dist_to_Ref_End = min(Ref_Pos_1indexed, Ref_Length - Ref_Pos_1indexed)
    // Matches compileMUMmer.py: min([x, y]) where x=Ref_Pos, y=Ref_Length-Ref_Pos
    let dist_to_ref_end = ref_pos_1indexed.min(ref_length.saturating_sub(ref_pos_1indexed));

    // Derive query contig per-observation from provenance, then look up its length.
    let query_contig = derive_query_contig(obs, query_id_fallback);
    let query_length = query_meta.as_ref()
        .and_then(|meta| lookup_contig_length(meta, &query_contig))
        .map(|l| l as u32);

    // Dist_to_Query_End = min(Query_Pos_1indexed, Query_Length - Query_Pos_1indexed)
    let dist_to_query_end = match query_length {
        Some(qlen) => query_pos_1indexed.min(qlen.saturating_sub(query_pos_1indexed)),
        None => 0,
    };

    // csp2_common.py checks for '-1' in reverse orientation — use 1/-1, not F/R
    let direction = match obs.strand() {
        Strand::Forward => "1",
        Strand::Reverse => "-1",
    };

    let cat = match obs.variant_type() {
        VariantType::Snp => "SNP",
        VariantType::Insertion | VariantType::Deletion => "Indel",
    };

    // Perc_Iden: raw CigarStats identity (not kmer-weighted confidence)
    let perc_iden = cigar_stats.identity(obs.edit_distance()) * 100.0;

    let ref_loc = format!("{}/{}", ref_contig, ref_pos_1indexed);
    let query_loc = format!("{}/{}", query_contig, query_pos_1indexed);

    // All 21 fields, tab-separated, single line
    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1}\t{}",
        ref_contig,                    // 1. Ref_Contig
        start_ref,                     // 2. Start_Ref (0-indexed)
        ref_pos_1indexed,              // 3. Ref_Pos (1-indexed)
        query_contig,                  // 4. Query_Contig (per-observation)
        start_query,                   // 5. Start_Query (0-indexed)
        query_pos_1indexed,            // 6. Query_Pos (1-indexed)
        ref_loc,                       // 7. Ref_Loc
        query_loc,                     // 8. Query_Loc
        ref_start_block,               // 9. Ref_Start (alignment block start)
        ref_end_block,                 // 10. Ref_End
        query_start_block,             // 11. Query_Start
        query_end_block,               // 12. Query_End
        ref_base,                      // 13. Ref_Base
        query_base,                    // 14. Query_Base
        dist_to_ref_end,               // 15. Dist_to_Ref_End
        dist_to_query_end,             // 16. Dist_to_Query_End
        ref_aligned,                   // 17. Ref_Aligned
        query_aligned,                 // 18. Query_Aligned
        direction,                     // 19. Query_Direction (1 or -1)
        perc_iden,                     // 20. Perc_Iden
        cat,                           // 21. Cat (SNP or Indel)
    )
}

/// Error type for snpdiffs formatting.
#[derive(Debug)]
pub struct SnpdiffsError {
    pub message: String,
}

impl std::fmt::Display for SnpdiffsError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.message)
    }
}

impl std::error::Error for SnpdiffsError {}

/// Format a `PhrayaFile` as CSP2-compatible snpdiffs text.
///
/// # Arguments
/// * `phraya_file` — the full `.phraya` file (header + observations + coverage track)
/// * `observations` — filtered variant observations to emit (may be a subset of the file's full set)
/// * `reference_id` — reference sequence label for header and BED rows (falls back to `header.sample_id`)
/// * `query_id` — query sequence label for header (falls back to first observation's provenance contig)
/// * `reference_fasta` — path to reference FASTA for N50/L90/contig count/SHA256 metadata
/// * `query_fasta` — path to query FASTA for same, AND for per-contig length lookup
///
/// # Output
/// Three sections, newline-separated:
/// 1. `#` header line (tab-separated key:value)
/// 2. `##` BED row(s) (one per reference contig, tab-separated)
/// 3. Variant data rows (one tab-separated line per variant)
pub fn format_snpdiffs(
    phraya_file: &phraya_io::phraya::PhrayaFile,
    observations: &[VariantObservation],
    reference_id: Option<&str>,
    query_id: Option<&str>,
    reference_fasta: Option<&Path>,
    query_fasta: Option<&Path>,
) -> Result<String, SnpdiffsError> {
    let ref_length = phraya_file.header.reference_length;
    let covered = phraya_file.header.covered_positions;

    // IDs for header labels (always present, matching compileMUMmer.py's Query_ID/Reference_ID)
    let ref_label = reference_id
        .unwrap_or(&phraya_file.header.sample_id)
        .to_string();
    let ref_fasta_path = reference_fasta
        .map(|p| p.to_string_lossy().to_string())
        .unwrap_or_default();

    // For the header and BED row, use query_id if provided, else first obs provenance
    let query_label = query_id
        .map(str::to_string)
        .or_else(|| {
            observations.first().map(|obs| {
                let prov = obs.provenance();
                prov.split(':').next().unwrap_or(prov).to_string()
            })
        })
        .unwrap_or_else(|| "query".to_string());

    let query_fasta_path = query_fasta
        .map(|p| p.to_string_lossy().to_string())
        .unwrap_or_default();

    // FASTA metadata
    let ref_meta = match reference_fasta {
        Some(path) => Some(compute_fasta_meta(path).map_err(|e| SnpdiffsError {
            message: e,
        })?),
        None => None,
    };
    let query_meta = match query_fasta {
        Some(path) => Some(compute_fasta_meta(path).map_err(|e| SnpdiffsError {
            message: e,
        })?),
        None => None,
    };

    // For the BED row, derive query_length from the first observation's contig
    // (or query_label if no observations). This is an aggregate value.
    let query_length_for_bed = match &query_meta {
        Some(meta) => {
            let bed_query_contig = if observations.is_empty() {
                query_label.clone()
            } else {
                derive_query_contig(observations.first().unwrap(), Some(&query_label))
            };
            lookup_contig_length(meta, &bed_query_contig).map(|l| l as u32)
        }
        None => None,
    };

    // For the header's Query_Percent_Aligned, use the same query_length
    let summary = compute_summary(observations, covered, ref_length, query_length_for_bed);

    // Build output
    let mut output = String::new();

    // 1. Header
    output.push_str(&format_header(
        &query_meta, &query_label, &query_fasta_path,
        &ref_meta, &ref_label, &ref_fasta_path,
        &summary,
    ));
    output.push('\n');
    output.push_str(&format_bed_row(
        &ref_label,
        ref_length,
        covered,
        &query_label,
        query_length_for_bed,
        observations,
    ));
    output.push('\n');

    // 3. Variant rows (one tab-separated line per variant)
    // Each row derives its own Query_Contig from provenance, then looks up
    // the per-contig length from query_meta.
    let query_id_fallback: Option<&str> = query_id;
    for obs in observations {
        output.push_str(&format_variant_row(
            obs,
            &ref_label,
            &query_meta,
            query_id_fallback,
            ref_length,
        ));
        output.push('\n');
    }

    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;
    use phraya_core::types::CoverageTrack;
    use std::io::Write;

    /// Helper: build a VariantObservation with the correct number of args.
    /// `edit_distance` is the key parameter for snpdiffs identity computation;
    /// `confidence` is set to a dummy value since snpdiffs uses CigarStats identity.
    fn make_obs(
        position: u32,
        ref_base: u8,
        alleles: HashMap<u8, u32>,
        confidence: f64,
        edit_distance: u32,
        cigar: &str,
        query_position: u32,
        variant_type: VariantType,
        strand: Strand,
    ) -> VariantObservation {
        VariantObservation::new(
            position,
            ref_base,
            alleles,
            confidence,
            cigar.to_string(),
            60,            // mapq
            edit_distance, // edit_distance
            vec![10],      // local_coverage
            35.0,          // avg_base_quality
            "query:read1".to_string(), // provenance
        )
        .with_variant_type(variant_type)
        .with_strand(strand)
        .with_query_position(query_position)
        .with_coverage_window_offset(50)
    }

    // ── compute_alignment_start ──────────────────────────────────────

    #[test]
    fn test_compute_alignment_start_snp_in_x_run() {
        // CIGAR: 5M1X14M, SNP at first X (position=10, query_position=5)
        // 5M: q_cum=5, r_cum=5; 1X: q_cum=5, query_pos=5 → offset=0, r_cum=5
        // target_start = 10 - (5 + 0) = 5
        let start = compute_alignment_start("5M1X14M", 10, 5);
        assert_eq!(start, 5);
    }

    #[test]
    fn test_compute_alignment_start_snp_in_x_run_offset() {
        // CIGAR: 5M3X12M, SNP at second X (position=11, query_position=6)
        // 5M: q_cum=5, r_cum=5; 3X: q_cum=5, query_pos=6 → offset=1, r_cum=5
        // target_start = 11 - (5 + 1) = 5
        let start = compute_alignment_start("5M3X12M", 11, 6);
        assert_eq!(start, 5);
    }

    #[test]
    fn test_compute_alignment_start_deletion() {
        // CIGAR: 5M2I3M (Phraya I = target-extra = deletion in query)
        // Deletion at position=10, query_position=5
        // 5M: q_cum=5, r_cum=5; 2I: q_cum=5==5 → target_start = 10 - 5 = 5
        let start = compute_alignment_start("5M2I3M", 10, 5);
        assert_eq!(start, 5);
    }

    #[test]
    fn test_compute_alignment_start_insertion() {
        // CIGAR: 5M2D3M (Phraya D = query-extra = insertion in query)
        // Insertion at position=10, query_position=5
        // 5M: q_cum=5, r_cum=5; 2D: q_cum=5==5 → target_start = 10 - 5 = 5
        let start = compute_alignment_start("5M2D3M", 10, 5);
        assert_eq!(start, 5);
    }

    #[test]
    fn test_compute_alignment_start_multiple_ops() {
        // CIGAR: 3M1I2M1D2X3M
        // Trace:
        // 3M: q_cum=3, r_cum=3
        // 1I: q_cum stays 3, check 3==7? No. r_cum=4
        // 2M: q_cum=5, r_cum=6
        // 1D: q_cum starts at 5, check 5==7? No. q_cum=6
        // 2X: q_cum=6, query_pos=7 is in [6, 8). offset=1, r_cum=6
        // target_start = 20 - (6 + 1) = 13
        let start = compute_alignment_start("3M1I2M1D2X3M", 20, 7);
        assert_eq!(start, 13);
    }

    #[test]
    fn test_compute_alignment_start_at_sequence_start() {
        // SNP at position 0, query_position 0
        // CIGAR: 1X4M
        // At the X op: q_cum=0, r_cum=0, offset=0
        // target_start = 0 - (0 + 0) = 0
        let start = compute_alignment_start("1X4M", 0, 0);
        assert_eq!(start, 0);
    }

    // ── best_alt_allele ─────────────────────────────────────────────

    #[test]
    fn test_best_alt_allele_snp() {
        let mut alleles = HashMap::new();
        alleles.insert(b'A', 10);
        alleles.insert(b'T', 3);
        let result = best_alt_allele(b'A', &alleles, VariantType::Snp);
        assert_eq!(result, "T");
    }

    #[test]
    fn test_best_alt_allele_snp_highest_count() {
        let mut alleles = HashMap::new();
        alleles.insert(b'A', 10);
        alleles.insert(b'C', 20);
        alleles.insert(b'G', 5);
        let result = best_alt_allele(b'A', &alleles, VariantType::Snp);
        assert_eq!(result, "C");
    }

    #[test]
    fn test_best_alt_allele_insertion() {
        let mut alleles = HashMap::new();
        alleles.insert(b'C', 1);
        alleles.insert(b'G', 1);
        let result = best_alt_allele(b'.', &alleles, VariantType::Insertion);
        assert_eq!(result, "CG");
    }

    #[test]
    fn test_best_alt_allele_deletion() {
        let alleles = HashMap::new();
        let result = best_alt_allele(b'A', &alleles, VariantType::Deletion);
        assert_eq!(result, ".");
    }

    // ── obs_identity ────────────────────────────────────────────────

    #[test]
    fn test_obs_identity_perfect() {
        let mut alleles = HashMap::new();
        alleles.insert(b'T', 5);
        let obs = make_obs(10, b'A', alleles, 0.95, 0, "100M", 5,
            VariantType::Snp, Strand::Forward);
        assert_eq!(obs_identity(&obs), 100.0);
    }

    #[test]
    fn test_obs_identity_with_mismatch() {
        let mut alleles = HashMap::new();
        alleles.insert(b'T', 5);
        // CIGAR: 99M1X, 1 mismatch in 100bp → identity = 1 - 1/100 = 0.99 → 99.0%
        let obs = make_obs(10, b'A', alleles, 0.95, 1, "99M1X", 5,
            VariantType::Snp, Strand::Forward);
        assert!((obs_identity(&obs) - 99.0).abs() < 0.01);
    }

    // ── format_variant_row — single line, tab-separated ─────────────

    #[test]
    fn test_format_variant_row_snp() {
        // CIGAR: 10M1X39M, SNP at first X: position=60, query_position=10
        // 10M: q=10, r=10; 1X: q_cum=10==10 → offset=0, r_cum=10
        // target_start = 60 - (10 + 0) = 50
        // CigarStats: matches=49, mismatches=1
        // target_aligned=50, query_aligned=50
        // edit_distance=1, identity = 1 - 1/50 = 0.98 → 98.0%
        let mut alleles = HashMap::new();
        alleles.insert(b'T', 5);
        let obs = make_obs(
            60, b'A', alleles, 0.95, 1, "10M1X39M", 10,
            VariantType::Snp, Strand::Forward,
        );

        // No query FASTA meta → query_length = None → Dist_to_Query_End = 0
        let query_meta: Option<FastaMeta> = None;
        let row = format_variant_row(&obs, "ref_ctg", &query_meta, None, 200);
        let fields: Vec<&str> = row.split('\t').collect();
        assert_eq!(fields.len(), 21);

        assert_eq!(fields[0], "ref_ctg");       // Ref_Contig
        assert_eq!(fields[1], "60");             // Start_Ref (0-indexed)
        assert_eq!(fields[2], "61");             // Ref_Pos (1-indexed)
        assert_eq!(fields[3], "query");          // Query_Contig (from provenance "query:read1")
        assert_eq!(fields[4], "10");             // Start_Query (0-indexed)
        assert_eq!(fields[5], "11");             // Query_Pos (1-indexed)
        assert_eq!(fields[6], "ref_ctg/61");     // Ref_Loc
        assert_eq!(fields[7], "query/11");       // Query_Loc
        assert_eq!(fields[8], "50");             // Ref_Start (alignment block start)
        assert_eq!(fields[9], "100");            // Ref_End = 50 + 50
        assert_eq!(fields[10], "0");             // Query_Start
        assert_eq!(fields[11], "50");            // Query_End = 50

        assert_eq!(fields[12], "A");             // Ref_Base
        assert_eq!(fields[13], "T");             // Query_Base

        // Dist_to_Ref_End = min(61, 200-61) = min(61, 139) = 61
        assert_eq!(fields[14], "61");
        // Dist_to_Query_End = 0 (no query_meta → query_length = None)
        assert_eq!(fields[15], "0");

        assert_eq!(fields[16], "50");           // Ref_Aligned
        assert_eq!(fields[17], "50");           // Query_Aligned
        assert_eq!(fields[18], "1");            // Query_Direction (forward = 1)
        assert_eq!(fields[19], "98.0");         // Perc_Iden = 100 * (1 - 1/50)
        assert_eq!(fields[20], "SNP");          // Cat
    }

    #[test]
    fn test_format_variant_row_snp_with_query_meta() {
        // Same as above but with query FASTA metadata so query_length is resolved.
        // Provenance is "query:read1" → contig name "query" → looked up in FASTA.
        let mut alleles = HashMap::new();
        alleles.insert(b'T', 5);
        let obs = make_obs(
            60, b'A', alleles, 0.95, 1, "10M1X39M", 10,
            VariantType::Snp, Strand::Forward,
        );

        let query_meta = Some(FastaMeta {
            contig_count: 1,
            total_bases: 200,
            n50: 200,
            n90: 200,
            l50: 1,
            l90: 1,
            sha256: "abc".to_string(),
            contig_lengths: vec![("query".to_string(), 200)],
        });

        let row = format_variant_row(&obs, "ref_ctg", &query_meta, None, 200);
        let fields: Vec<&str> = row.split('\t').collect();

        // Dist_to_Query_End = min(11, 200-11) = min(11, 189) = 11
        assert_eq!(fields[15], "11");
    }

    #[test]
    fn test_format_variant_row_deletion() {
        // CIGAR: 50M2I3M (Phraya I = target-extra = deletion)
        // Deletion at position=50, query_position=50
        // 50M: q=50, r=50; 2I: q_cum=50==50 → target_start = 50-50 = 0
        // CigarStats: matches=53, mismatches=0, target_extra=2
        // target_aligned = 53+0+2 = 55, query_aligned = 53+0+0 = 53
        // edit_distance=0 → identity = 1.0 → 100.0%
        let mut alleles = HashMap::new();
        alleles.insert(b'.', 1);
        let obs = make_obs(
            50, b'T', alleles, 0.90, 0, "50M2I3M", 50,
            VariantType::Deletion, Strand::Reverse,
        );

        let query_meta: Option<FastaMeta> = None;
        let row = format_variant_row(&obs, "ref", &query_meta, None, 200);
        let fields: Vec<&str> = row.split('\t').collect();
        assert_eq!(fields.len(), 21);

        assert_eq!(fields[8], "0");      // Ref_Start = 0
        assert_eq!(fields[9], "55");     // Ref_End = 0 + 55
        assert_eq!(fields[11], "53");    // Query_End

        assert_eq!(fields[12], "T");     // Ref_Base (deleted base)
        assert_eq!(fields[13], ".");     // Query_Base (deletion → .)

        assert_eq!(fields[16], "55");    // Ref_Aligned
        assert_eq!(fields[17], "53");    // Query_Aligned
        assert_eq!(fields[18], "-1");    // Query_Direction (reverse = -1)
        assert_eq!(fields[19], "100.0"); // Perc_Iden = identity(0) = 100.0%
        assert_eq!(fields[20], "Indel");
    }

    #[test]
    fn test_format_variant_row_insertion() {
        // CIGAR: 50M1D3M (Phraya D = query-extra = insertion)
        // Insertion at position=50, query_position=50
        // 50M: q=50, r=50; 1D: q_cum=50==50 → target_start = 50-50 = 0
        // CigarStats: matches=53, mismatches=0, query_extra=1
        // target_aligned = 53+0+0 = 53, query_aligned = 53+0+1 = 54
        // edit_distance=0 → identity = 1.0 → 100.0%
        let mut alleles = HashMap::new();
        alleles.insert(b'G', 1);
        let obs = make_obs(
            50, b'.', alleles, 0.92, 0, "50M1D3M", 50,
            VariantType::Insertion, Strand::Forward,
        );

        let query_meta: Option<FastaMeta> = None;
        let row = format_variant_row(&obs, "ref", &query_meta, None, 200);
        let fields: Vec<&str> = row.split('\t').collect();
        assert_eq!(fields.len(), 21);

        assert_eq!(fields[8], "0");      // Ref_Start
        assert_eq!(fields[9], "53");     // Ref_End = 0 + 53
        assert_eq!(fields[11], "54");    // Query_End

        assert_eq!(fields[12], ".");     // Ref_Base (insertion → .)
        assert_eq!(fields[13], "G");     // Query_Base (inserted base)

        assert_eq!(fields[16], "53");    // Ref_Aligned
        assert_eq!(fields[17], "54");    // Query_Aligned
        assert_eq!(fields[18], "1");     // Query_Direction (forward)
        assert_eq!(fields[19], "100.0"); // Perc_Iden = identity(0) = 100.0%
        assert_eq!(fields[20], "Indel");
    }

    #[test]
    fn test_format_variant_row_reverse_strand() {
        // CIGAR: 10M1X29M, 1 mismatch → identity = 1 - 1/40 = 0.975 → 97.5%
        let mut alleles = HashMap::new();
        alleles.insert(b'G', 5);
        let obs = make_obs(
            50, b'A', alleles, 0.98, 1, "10M1X29M", 10,
            VariantType::Snp, Strand::Reverse,
        );

        let query_meta: Option<FastaMeta> = None;
        let row = format_variant_row(&obs, "ref", &query_meta, None, 100);
        let fields: Vec<&str> = row.split('\t').collect();
        assert_eq!(fields[18], "-1");   // Reverse direction
        assert_eq!(fields[19], "97.5"); // Perc_Iden = 100 * (1 - 1/40)
        assert_eq!(fields[20], "SNP");
    }

    #[test]
    fn test_format_variant_row_distance_to_edge() {
        // CIGAR: 1X9M, 1 mismatch → identity = 1 - 1/10 = 0.9 → 90.0%
        let mut alleles = HashMap::new();
        alleles.insert(b'T', 5);
        let obs = make_obs(
            2, b'A', alleles, 0.95, 1, "1X9M", 0,
            VariantType::Snp, Strand::Forward,
        );

        let query_meta = Some(FastaMeta {
            contig_count: 1,
            total_bases: 10,
            n50: 10,
            n90: 10,
            l50: 1,
            l90: 1,
            sha256: "abc".to_string(),
            contig_lengths: vec![("query".to_string(), 10)],
        });

        let row = format_variant_row(&obs, "ref", &query_meta, None, 100);
        let fields: Vec<&str> = row.split('\t').collect();

        // Ref_Pos_1indexed = 3, Ref_Length = 100
        // Dist_to_Ref_End = min(3, 100-3) = min(3, 97) = 3
        assert_eq!(fields[14], "3");
        // Query_Pos_1indexed = 1, Query_Length = 10 (looked up from "query" contig)
        // Dist_to_Query_End = min(1, 10-1) = min(1, 9) = 1
        assert_eq!(fields[15], "1");
    }

    #[test]
    fn test_format_variant_row_distance_at_end() {
        // Variant near end of reference
        // CIGAR: 2X8M, 2 mismatches → identity = 1 - 2/10 = 0.8 → 80.0%
        let mut alleles = HashMap::new();
        alleles.insert(b'T', 5);
        let obs = make_obs(
            97, b'A', alleles, 0.95, 2, "2X8M", 0,
            VariantType::Snp, Strand::Forward,
        );

        let query_meta = Some(FastaMeta {
            contig_count: 1,
            total_bases: 10,
            n50: 10,
            n90: 10,
            l50: 1,
            l90: 1,
            sha256: "abc".to_string(),
            contig_lengths: vec![("query".to_string(), 10)],
        });

        let row = format_variant_row(&obs, "ref", &query_meta, None, 100);
        let fields: Vec<&str> = row.split('\t').collect();

        // Ref_Pos_1indexed = 98, Ref_Length = 100
        // Dist_to_Ref_End = min(98, 100-98) = min(98, 2) = 2
        assert_eq!(fields[14], "2");
        // Query_Pos_1indexed = 1, Query_Length = 10
        // Dist_to_Query_End = min(1, 10-1) = min(1, 9) = 1
        assert_eq!(fields[15], "1");
    }

    #[test]
    fn test_format_variant_row_multi_contig_query() {
        // Verify that query contig is derived per-observation from provenance,
        // not from a single query_label. Two observations with different provenance
        // should get different Query_Contig and Dist_to_Query_End values.
        let obs_a = VariantObservation::new(
            50,
            b'A',
            {
                let mut m = HashMap::new();
                m.insert(b'T', 10);
                m
            },
            0.95,
            "100M".to_string(),
            60,
            0,     // edit_distance
            vec![10],
            35.0,
            "contig_a:read1".to_string(),
        )
        .with_variant_type(VariantType::Snp)
        .with_strand(Strand::Forward)
        .with_query_position(50)
        .with_coverage_window_offset(50);

        let obs_b = VariantObservation::new(
            50,
            b'G',
            {
                let mut m = HashMap::new();
                m.insert(b'C', 10);
                m
            },
            0.95,
            "100M".to_string(),
            60,
            0,     // edit_distance
            vec![10],
            35.0,
            "contig_b:read2".to_string(),
        )
        .with_variant_type(VariantType::Snp)
        .with_strand(Strand::Forward)
        .with_query_position(50)
        .with_coverage_window_offset(50);

        // FASTA with two contigs of different lengths
        let query_meta = Some(FastaMeta {
            contig_count: 2,
            total_bases: 200,
            n50: 100,
            n90: 100,
            l50: 1,
            l90: 2,
            sha256: "abc".to_string(),
            contig_lengths: vec![
                ("contig_a".to_string(), 100),
                ("contig_b".to_string(), 80),
            ],
        });

        let row_a = format_variant_row(&obs_a, "ref", &query_meta, Some("user_label"), 200);
        let row_b = format_variant_row(&obs_b, "ref", &query_meta, Some("user_label"), 200);

        let fields_a: Vec<&str> = row_a.split('\t').collect();
        let fields_b: Vec<&str> = row_b.split('\t').collect();

        // Each observation gets its own query contig name from provenance
        assert_eq!(fields_a[3], "contig_a");
        assert_eq!(fields_b[3], "contig_b");

        // Dist_to_Query_End differs because contig lengths differ
        // Ref_Pos = 51, Query_Pos = 51
        // contig_a (len=100): min(51, 100-51) = min(51, 49) = 49
        // contig_b (len=80): min(51, 80-51) = min(51, 29) = 29
        assert_eq!(fields_a[15], "49");
        assert_eq!(fields_b[15], "29");

        // Query_Loc also differs
        assert_eq!(fields_a[7], "contig_a/51");
        assert_eq!(fields_b[7], "contig_b/51");
    }

    // ── format_bed_row ──────────────────────────────────────────────

    #[test]
    fn test_format_bed_row() {
        // CIGAR: 100M, edit_distance=0 → identity = 100.0%
        let obs = make_obs(
            50, b'A', {
                let mut m = HashMap::new();
                m.insert(b'T', 5);
                m
            }, 0.95, 0, "100M", 50,
            VariantType::Snp, Strand::Forward,
        );

        let bed = format_bed_row("ref", 1000, Some(500), "query", Some(800), &[obs]);
        let fields: Vec<&str> = bed.split('\t').collect();
        assert_eq!(fields.len(), 12); // ## + 11 fields
        assert_eq!(fields[0], "##");
        assert_eq!(fields[1], "ref");
        assert_eq!(fields[2], "0");      // Ref_Start
        assert_eq!(fields[3], "1000");   // Ref_End
        assert_eq!(fields[4], "1000");   // Ref_Length
        assert_eq!(fields[5], "500");    // Ref_Aligned (covered_positions)
        assert_eq!(fields[6], "query");
        assert_eq!(fields[7], "0");      // Query_Start
        assert_eq!(fields[8], "800");     // Query_End
        assert_eq!(fields[9], "800");     // Query_Length
        assert_eq!(fields[10], "500");   // Query_Aligned
        assert_eq!(fields[11], "100.0");  // Perc_Iden = identity(0) = 100.0%
    }

    #[test]
    fn test_format_bed_row_no_coverage() {
        let bed = format_bed_row("ref", 1000, Some(0), "query", Some(800), &[]);
        let fields: Vec<&str> = bed.split('\t').collect();
        assert_eq!(fields[0], "##");
        assert_eq!(fields[5], "0");     // Ref_Aligned
        assert_eq!(fields[11], "0.0");  // Perc_Iden (no coverage, no obs)
    }

    // ── format_header ───────────────────────────────────────────────

    #[test]
    fn test_format_header_with_fasta() {
        let query_meta = Some(FastaMeta {
            contig_count: 3,
            total_bases: 500,
            n50: 200,
            n90: 100,
            l50: 2,
            l90: 3,
            sha256: "def456".to_string(),
            contig_lengths: vec![("q1".to_string(), 200), ("q2".to_string(), 200), ("q3".to_string(), 100)],
        });
        let ref_meta = Some(FastaMeta {
            contig_count: 1,
            total_bases: 1000,
            n50: 1000,
            n90: 1000,
            l50: 1,
            l90: 1,
            sha256: "abc123".to_string(),
            contig_lengths: vec![("r1".to_string(), 1000)],
        });

        let summary = SummaryMetrics {
            snp_count: 10,
            indel_count: 3,
            median_identity: 98.5,
            median_aln_len: 95.0,
            percent_ref_aligned: 95.0,
            percent_query_aligned: 90.0,
        };

        let header = format_header(
            &query_meta, "query1", "query.fasta",
            &ref_meta, "ref1", "ref.fasta",
            &summary,
        );

        let fields: Vec<&str> = header.split('\t').collect();
        assert_eq!(fields[0], "#");
        assert!(fields.iter().any(|f| f == &"Query_ID:query1"));
        assert!(fields.iter().any(|f| f == &"Query_Assembly:query.fasta"));
        assert!(fields.iter().any(|f| f == &"Query_Contig_Count:3"));
        assert!(fields.iter().any(|f| f == &"Query_Assembly_Bases:500"));
        assert!(fields.iter().any(|f| f == &"Query_N50:200"));
        assert!(fields.iter().any(|f| f == &"Query_N90:100"));
        assert!(fields.iter().any(|f| f == &"Query_L50:2"));
        assert!(fields.iter().any(|f| f == &"Query_L90:3"));
        assert!(fields.iter().any(|f| f == &"Query_SHA256:def456"));
        assert!(fields.iter().any(|f| f == &"Reference_ID:ref1"));
        assert!(fields.iter().any(|f| f == &"Reference_Assembly:ref.fasta"));
        assert!(fields.iter().any(|f| f == &"Reference_Contig_Count:1"));
        assert!(fields.iter().any(|f| f == &"Reference_Assembly_Bases:1000"));
        assert!(fields.iter().any(|f| f == &"Reference_SHA256:abc123"));
        assert!(fields.iter().any(|f| f.starts_with("SNPs:10")));
        assert!(fields.iter().any(|f| f.starts_with("Indels:3")));
        assert!(fields.iter().any(|f| f.starts_with("Median_Percent_Identity:98.50")));
    }

    #[test]
    fn test_format_header_no_fasta() {
        let summary = SummaryMetrics {
            snp_count: 0,
            indel_count: 0,
            median_identity: 0.0,
            median_aln_len: 0.0,
            percent_ref_aligned: 0.0,
            percent_query_aligned: 0.0,
        };

        let header = format_header(&None, "qry", "", &None, "ref", "", &summary);
        let fields: Vec<&str> = header.split('\t').collect();
        assert_eq!(fields[0], "#");
        assert!(fields.iter().any(|f| f == &"Query_ID:qry"));
        assert!(fields.iter().any(|f| f == &"Reference_ID:ref"));
        assert!(fields.iter().any(|f| f == &"Query_Contig_Count:NA"));
        assert!(fields.iter().any(|f| f == &"Reference_SHA256:NA"));
        assert!(fields.iter().any(|f| f.starts_with("Kmer_Similarity:NA")));
        assert!(fields.iter().any(|f| f.starts_with("gSNPs:NA")));
    }

    // ── compute_fasta_meta ──────────────────────────────────────────

    #[test]
    fn test_compute_fasta_meta() {
        let dir = tempfile::TempDir::new().unwrap();
        let fasta_path = dir.path().join("test.fasta");
        {
            let mut f = std::fs::File::create(&fasta_path).unwrap();
            writeln!(f, ">contig1").unwrap();
            writeln!(f, "ACGTACGTAC").unwrap(); // 10 bp
            writeln!(f, ">contig2").unwrap();
            writeln!(f, "TTTTGGGGCC").unwrap(); // 10 bp
            writeln!(f, ">contig3").unwrap();
            writeln!(f, "AAAAAAAAAA").unwrap(); // 10 bp
        }

        let meta = compute_fasta_meta(&fasta_path).unwrap();
        assert_eq!(meta.contig_count, 3);
        assert_eq!(meta.total_bases, 30);
        // Sorted desc: [10, 10, 10]. Half=15. After second: 20 >= 15
        assert_eq!(meta.n50, 10);
        assert_eq!(meta.l50, 2);
        assert_eq!(meta.n90, 10);
        assert_eq!(meta.l90, 3);
        assert!(!meta.sha256.is_empty());
        assert_eq!(meta.contig_lengths.len(), 3);
        assert_eq!(meta.contig_lengths[0], ("contig1".to_string(), 10));
    }

    #[test]
    fn test_compute_fasta_meta_different_lengths() {
        let dir = tempfile::TempDir::new().unwrap();
        let fasta_path = dir.path().join("test.fasta");
        {
            let mut f = std::fs::File::create(&fasta_path).unwrap();
            writeln!(f, ">big").unwrap();
            writeln!(f, "{}", "A".repeat(100)).unwrap(); // 100 bp
            writeln!(f, ">small1").unwrap();
            writeln!(f, "{}", "C".repeat(10)).unwrap(); // 10 bp
            writeln!(f, ">small2").unwrap();
            writeln!(f, "{}", "G".repeat(10)).unwrap(); // 10 bp
        }

        let meta = compute_fasta_meta(&fasta_path).unwrap();
        assert_eq!(meta.contig_count, 3);
        assert_eq!(meta.total_bases, 120);
        // Sorted desc: [100, 10, 10]. Half=60. After first: 100 >= 60
        assert_eq!(meta.n50, 100);
        assert_eq!(meta.l50, 1);
        // 90% threshold = ceil(120*0.9) = 108. After first: 100 < 108, after second: 110 >= 108
        assert_eq!(meta.n90, 10);
        assert_eq!(meta.l90, 2);
    }

    // ── lookup_contig_length ────────────────────────────────────────

    #[test]
    fn test_lookup_contig_length_single_contig() {
        let meta = FastaMeta {
            contig_count: 1,
            total_bases: 200,
            n50: 200,
            n90: 200,
            l50: 1,
            l90: 1,
            sha256: "abc".to_string(),
            contig_lengths: vec![("contig1".to_string(), 200)],
        };
        // Name matches
        assert_eq!(lookup_contig_length(&meta, "contig1"), Some(200));
        // Name doesn't match but single contig → fallback to total_bases
        assert_eq!(lookup_contig_length(&meta, "other"), Some(200));
    }

    #[test]
    fn test_lookup_contig_length_multi_contig() {
        let meta = FastaMeta {
            contig_count: 3,
            total_bases: 300,
            n50: 100,
            n90: 100,
            l50: 1,
            l90: 3,
            sha256: "abc".to_string(),
            contig_lengths: vec![
                ("big".to_string(), 100),
                ("mid".to_string(), 100),
                ("small".to_string(), 100),
            ],
        };
        assert_eq!(lookup_contig_length(&meta, "big"), Some(100));
        assert_eq!(lookup_contig_length(&meta, "mid"), Some(100));
        assert_eq!(lookup_contig_length(&meta, "nonexistent"), None);
    }

    // ── n_stats ─────────────────────────────────────────────────────

    #[test]
    fn test_n_stats_uniform() {
        let lengths = vec![10, 10, 10, 10];
        let (n50, l50) = n_stats(&lengths, 40, 0.5);
        assert_eq!(n50, 10);
        assert_eq!(l50, 2);
    }

    #[test]
    fn test_n_stats_skewed() {
        let lengths = vec![100, 50, 25];
        let (n50, l50) = n_stats(&lengths, 175, 0.5);
        // After first: 100 >= 88 → N50=100, L50=1
        assert_eq!(n50, 100);
        assert_eq!(l50, 1);
        let (n90, l90) = n_stats(&lengths, 175, 0.9);
        // 90% threshold = ceil(175*0.9) = 158. After third: 175 >= 158
        assert_eq!(n90, 25);
        assert_eq!(l90, 3);
    }

    // ── derive_query_contig ─────────────────────────────────────────

    #[test]
    fn test_derive_query_contig_from_provenance() {
        let obs = make_obs(
            10, b'A', { let mut m = HashMap::new(); m.insert(b'T', 1); m },
            0.95, 0, "50M", 5, VariantType::Snp, Strand::Forward,
        );
        assert_eq!(derive_query_contig(&obs, None), "query");
    }

    #[test]
    fn test_derive_query_contig_with_fallback() {
        let obs = VariantObservation::new(
            10, b'A', { let mut m = HashMap::new(); m.insert(b'T', 1); m },
            0.95, "50M".to_string(), 60, 0, vec![10], 35.0,
            "".to_string(), // empty provenance
        )
        .with_variant_type(VariantType::Snp)
        .with_strand(Strand::Forward)
        .with_query_position(5)
        .with_coverage_window_offset(50);

        assert_eq!(derive_query_contig(&obs, Some("fallback_id")), "fallback_id");
    }

    #[test]
    fn test_derive_query_contig_default() {
        let obs = VariantObservation::new(
            10, b'A', { let mut m = HashMap::new(); m.insert(b'T', 1); m },
            0.95, "50M".to_string(), 60, 0, vec![10], 35.0,
            "".to_string(),
        )
        .with_variant_type(VariantType::Snp)
        .with_strand(Strand::Forward)
        .with_query_position(5)
        .with_coverage_window_offset(50);

        assert_eq!(derive_query_contig(&obs, None), "query");
    }

    // ── Integration: format_snpdiffs end-to-end ─────────────────────

    #[test]
    fn test_format_snpdiffs_end_to_end() {
        // CIGAR: 100M, edit_distance=0 → identity = 100.0%
        let obs = VariantObservation::new(
            50,
            b'A',
            {
                let mut m = HashMap::new();
                m.insert(b'T', 10);
                m
            },
            0.95,
            "100M".to_string(),
            60,
            0,      // edit_distance
            vec![10],
            35.0,
            "query:read1".to_string(),
        )
        .with_variant_type(VariantType::Snp)
        .with_strand(Strand::Forward)
        .with_query_position(50)
        .with_coverage_window_offset(50);

        let mut coverage_vec = vec![0u32; 200];
        for i in 0..100 {
            coverage_vec[i] = 1;
        }
        let coverage = CoverageTrack::new(&coverage_vec);
        let phraya_file = phraya_io::phraya::PhrayaFile::new(
            200,
            "ref_contig".to_string(),
            "2024-01-01T00:00:00Z".to_string(),
            vec![obs.clone()],
            coverage,
        )
        .with_coverage_breadth(100, 0);

        let output = format_snpdiffs(
            &phraya_file,
            &[obs],
            Some("my_ref"),
            Some("my_query"),
            None,
            None,
        )
        .unwrap();

        let lines: Vec<&str> = output.lines().collect();

        // Line 0: Header
        let header_fields: Vec<&str> = lines[0].split('\t').collect();
        assert!(header_fields.iter().any(|f| f == &"Query_ID:my_query"));
        assert!(header_fields.iter().any(|f| f == &"Reference_ID:my_ref"));
        assert!(header_fields.iter().any(|f| f == &"Query_Contig_Count:NA"));
        assert!(header_fields.iter().any(|f| f == &"Reference_SHA256:NA"));
        assert!(header_fields.iter().any(|f| f.starts_with("SNPs:1")));
        assert!(header_fields.iter().any(|f| f.starts_with("Indels:0")));
        assert!(header_fields.iter().any(|f| f == &"Kmer_Similarity:NA"));

        // Line 1: BED row
        let bed_fields: Vec<&str> = lines[1].split('\t').collect();
        assert_eq!(bed_fields.len(), 12);
        assert_eq!(bed_fields[0], "##");
        assert_eq!(bed_fields[1], "my_ref");
        assert_eq!(bed_fields[2], "0");      // Ref_Start
        assert_eq!(bed_fields[3], "200");    // Ref_End
        assert_eq!(bed_fields[4], "200");    // Ref_Length
        assert_eq!(bed_fields[5], "100");    // Ref_Aligned
        assert_eq!(bed_fields[6], "my_query");
        assert_eq!(bed_fields[7], "0");      // Query_Start
        assert_eq!(bed_fields[8], "0");      // Query_End (no query_length → 0)
        assert_eq!(bed_fields[9], "0");      // Query_Length
        assert_eq!(bed_fields[10], "100");   // Query_Aligned
        assert_eq!(bed_fields[11], "100.0"); // Perc_Iden = identity(0) = 100.0%

        // Line 2: Variant row (single line, 21 tab-separated fields)
        let var_fields: Vec<&str> = lines[2].split('\t').collect();
        assert_eq!(var_fields.len(), 21);
        assert_eq!(var_fields[0], "my_ref");     // Ref_Contig
        assert_eq!(var_fields[1], "50");          // Start_Ref
        assert_eq!(var_fields[2], "51");          // Ref_Pos
        assert_eq!(var_fields[3], "query");        // Query_Contig (from provenance, not header label)
        assert_eq!(var_fields[4], "50");          // Start_Query
        assert_eq!(var_fields[5], "51");          // Query_Pos
        assert_eq!(var_fields[6], "my_ref/51");   // Ref_Loc
        assert_eq!(var_fields[7], "query/51");     // Query_Loc (from provenance contig)
        assert_eq!(var_fields[12], "A");          // Ref_Base
        assert_eq!(var_fields[13], "T");          // Query_Base
        assert_eq!(var_fields[18], "1");           // Query_Direction (forward = 1)
        assert_eq!(var_fields[19], "100.0");      // Perc_Iden = identity(0) = 100.0%
        assert_eq!(var_fields[20], "SNP");        // Cat
    }

    #[test]
    fn test_format_snpdiffs_empty_observations() {
        let coverage = CoverageTrack::uncovered(200);
        let phraya_file = phraya_io::phraya::PhrayaFile::new(
            200,
            "ref".to_string(),
            "2024-01-01T00:00:00Z".to_string(),
            vec![],
            coverage,
        )
        .with_coverage_breadth(0, 0);

        let output = format_snpdiffs(
            &phraya_file,
            &[],
            Some("ref"),
            Some("query"),
            None,
            None,
        )
        .unwrap();

        let lines: Vec<&str> = output.lines().collect();

        // Header + BED row only (no variants)
        assert_eq!(lines.len(), 2, "should have header + 1 BED line only");
        assert!(lines[0].starts_with("#\t"));
        assert_eq!(lines[1], "##\tref\t0\t200\t200\t0\tquery\t0\t0\t0\t0\t0.0");
    }

    #[test]
    fn test_format_snpdiffs_query_id_from_provenance() {
        let obs = VariantObservation::new(
            50,
            b'A',
            {
                let mut m = HashMap::new();
                m.insert(b'T', 5);
                m
            },
            0.95,
            "100M".to_string(),
            60,
            0,      // edit_distance
            vec![10],
            35.0,
            "assembly42:contig1".to_string(),
        )
        .with_variant_type(VariantType::Snp)
        .with_strand(Strand::Forward)
        .with_query_position(50)
        .with_coverage_window_offset(50);

        let coverage = CoverageTrack::uncovered(200);
        let phraya_file = phraya_io::phraya::PhrayaFile::new(
            200,
            "ref_contig".to_string(),
            "2024-01-01T00:00:00Z".to_string(),
            vec![obs.clone()],
            coverage,
        )
        .with_coverage_breadth(100, 0);

        let output = format_snpdiffs(
            &phraya_file,
            &[obs],
            Some("ref"),
            None, // No query_id → falls back to provenance "assembly42"
            None,
            None,
        )
        .unwrap();

        let lines: Vec<&str> = output.lines().collect();
        let header_fields: Vec<&str> = lines[0].split('\t').collect();
        // Provenance is "assembly42:contig1", split on ':' gives "assembly42"
        assert!(header_fields.iter().any(|f| f == &"Query_ID:assembly42"));

        let bed_fields: Vec<&str> = lines[1].split('\t').collect();
        assert_eq!(bed_fields[6], "assembly42");
    }

    #[test]
    fn test_format_snpdiffs_with_fasta() {
        use std::io::Write;

        // Create reference FASTA
        let dir = tempfile::TempDir::new().unwrap();
        let ref_fasta = dir.path().join("ref.fasta");
        {
            let mut f = std::fs::File::create(&ref_fasta).unwrap();
            writeln!(f, ">ref_contig").unwrap();
            writeln!(f, "{}", "A".repeat(200)).unwrap();
        }

        // Create query FASTA
        let query_fasta = dir.path().join("query.fasta");
        {
            let mut f = std::fs::File::create(&query_fasta).unwrap();
            writeln!(f, ">query_contig").unwrap();
            writeln!(f, "{}", "T".repeat(200)).unwrap();
        }

        // edit_distance=0, CIGAR 100M → identity = 100.0%
        // Provenance: "query_contig:read1" → contig name "query_contig"
        let obs = VariantObservation::new(
            50,
            b'A',
            {
                let mut m = HashMap::new();
                m.insert(b'T', 10);
                m
            },
            0.95,
            "100M".to_string(),
            60,
            0,      // edit_distance
            vec![10],
            35.0,
            "query_contig:read1".to_string(),
        )
        .with_variant_type(VariantType::Snp)
        .with_strand(Strand::Forward)
        .with_query_position(50)
        .with_coverage_window_offset(50);

        let mut coverage_vec = vec![0u32; 200];
        for i in 0..100 {
            coverage_vec[i] = 1;
        }
        let coverage = CoverageTrack::new(&coverage_vec);
        let phraya_file = phraya_io::phraya::PhrayaFile::new(
            200,
            "ref_contig".to_string(),
            "2024-01-01T00:00:00Z".to_string(),
            vec![obs.clone()],
            coverage,
        )
        .with_coverage_breadth(100, 0);

        let output = format_snpdiffs(
            &phraya_file,
            &[obs],
            Some("my_ref"),
            Some("my_query"),
            Some(&ref_fasta),
            Some(&query_fasta),
        )
        .unwrap();

        let lines: Vec<&str> = output.lines().collect();

        // Header should have FASTA metadata now
        let header_fields: Vec<&str> = lines[0].split('\t').collect();
        assert!(header_fields.iter().any(|f| f == &"Query_ID:my_query"));
        assert!(header_fields.iter().any(|f| f == &"Reference_ID:my_ref"));
        assert!(header_fields.iter().any(|f| f == &"Query_Contig_Count:1"));
        assert!(header_fields.iter().any(|f| f == &"Query_Assembly_Bases:200"));
        assert!(header_fields.iter().any(|f| f == &"Reference_Contig_Count:1"));
        assert!(header_fields.iter().any(|f| f == &"Reference_Assembly_Bases:200"));
        // SHA256 should be a real hash, not NA
        assert!(header_fields.iter().any(|f| f.starts_with("Query_SHA256:") && f != &"Query_SHA256:NA"));
        assert!(header_fields.iter().any(|f| f.starts_with("Reference_SHA256:") && f != &"Reference_SHA256:NA"));

        // BED row should have query length now (single-contig FASTA → lookup works)
        let bed_fields: Vec<&str> = lines[1].split('\t').collect();
        assert_eq!(bed_fields[8], "200");   // Query_End
        assert_eq!(bed_fields[9], "200");   // Query_Length
        assert_eq!(bed_fields[10], "100");  // Query_Aligned

        // Variant row should use per-observation query contig from provenance
        let var_fields: Vec<&str> = lines[2].split('\t').collect();
        assert_eq!(var_fields[3], "query_contig");  // Query_Contig from provenance
        assert_eq!(var_fields[7], "query_contig/51"); // Query_Loc
    }

    #[test]
    fn test_parseSNPDiffs_compatibility() {
        // Verify that the output can be parsed by the same logic as csp2_common.py's
        // parseSNPDiffs: header starts with "#\t", BED rows start with "##\t",
        // variant rows have no prefix and are tab-separated.
        // CIGAR: 49M1X, edit_distance=1 → identity = 1 - 1/50 = 0.98 → 98.0%
        let obs = VariantObservation::new(
            100,
            b'A',
            {
                let mut m = HashMap::new();
                m.insert(b'T', 5);
                m
            },
            0.98,
            "49M1X".to_string(),
            30,
            1,      // edit_distance
            vec![10],
            35.0,
            "qry:read1".to_string(),
        )
        .with_variant_type(VariantType::Snp)
        .with_strand(Strand::Forward)
        .with_query_position(30)
        .with_coverage_window_offset(30);

        let coverage = CoverageTrack::uncovered(200);
        let phraya_file = phraya_io::phraya::PhrayaFile::new(
            200,
            "ref_seq".to_string(),
            "2024-01-01T00:00:00Z".to_string(),
            vec![obs.clone()],
            coverage,
        )
        .with_coverage_breadth(100, 0);

        let output = format_snpdiffs(
            &phraya_file,
            &[obs],
            Some("ref_seq"),
            Some("qry"),
            None,
            None,
        )
        .unwrap();

        // Simulate csp2_common.py's parseSNPDiffs logic
        let mut bed_rows = Vec::new();
        let mut snp_rows = Vec::new();

        for line in output.lines() {
            if line.starts_with("#\t") {
                // Header — skip during data parsing
            } else if line.starts_with("##\t") {
                bed_rows.push(line.strip_prefix("##\t").unwrap().split('\t').collect::<Vec<_>>());
            } else {
                snp_rows.push(line.split('\t').collect::<Vec<_>>());
            }
        }

        // BED row should have 11 columns (matching csp2_common.py bed_columns)
        assert_eq!(bed_rows.len(), 1);
        assert_eq!(bed_rows[0].len(), 11);

        // Variant row should have 21 columns (matching csp2_common.py snp_columns)
        assert_eq!(snp_rows.len(), 1);
        assert_eq!(snp_rows[0].len(), 21);

        // Verify key fields
        assert_eq!(snp_rows[0][0], "ref_seq");   // Ref_Contig
        assert_eq!(snp_rows[0][19], "98.0");      // Perc_Iden = 100 * (1 - 1/50)
        assert_eq!(snp_rows[0][20], "SNP");       // Cat
        assert_eq!(snp_rows[0][18], "1");          // Query_Direction (forward = 1, not 'F')

        // Verify header parsing (csp2_common.py fetchHeaders: split on ':', expect 1 colon per field)
        let header_line = output.lines().next().unwrap();
        let header_fields: Vec<&str> = header_line.strip_prefix("#\t").unwrap().split('\t').collect();
        for field in &header_fields {
            // Each field must be key:value (exactly 1 colon, no colons in key or value)
            assert_eq!(
                field.matches(':').count(), 1,
                "field '{}' has != 1 colon", field
            );
        }
    }

    // ── Per-observation query contig in multi-contig FASTA ──────────

    #[test]
    fn test_format_snpdiffs_multicontig_query_per_variant() {
        use std::io::Write;

        // Create query FASTA with two contigs of different lengths
        let dir = tempfile::TempDir::new().unwrap();
        let query_fasta = dir.path().join("query.fasta");
        {
            let mut f = std::fs::File::create(&query_fasta).unwrap();
            writeln!(f, ">contigA").unwrap();
            writeln!(f, "{}", "A".repeat(100)).unwrap();
            writeln!(f, ">contigB").unwrap();
            writeln!(f, "{}", "T".repeat(80)).unwrap();
        }

        let ref_fasta = dir.path().join("ref.fasta");
        {
            let mut f = std::fs::File::create(&ref_fasta).unwrap();
            writeln!(f, ">ref").unwrap();
            writeln!(f, "{}", "N".repeat(200)).unwrap();
        }

        // Two observations from different query contigs
        let obs_a = VariantObservation::new(
            50, b'A', { let mut m = HashMap::new(); m.insert(b'T', 10); m },
            0.95, "100M".to_string(), 60, 0, vec![10], 35.0,
            "contigA:read1".to_string(),
        )
        .with_variant_type(VariantType::Snp)
        .with_strand(Strand::Forward)
        .with_query_position(50)
        .with_coverage_window_offset(50);

        let obs_b = VariantObservation::new(
            50, b'G', { let mut m = HashMap::new(); m.insert(b'C', 10); m },
            0.95, "100M".to_string(), 60, 0, vec![10], 35.0,
            "contigB:read2".to_string(),
        )
        .with_variant_type(VariantType::Snp)
        .with_strand(Strand::Forward)
        .with_query_position(50)
        .with_coverage_window_offset(50);

        let coverage = CoverageTrack::new(&vec![1u32; 200]);
        let phraya_file = phraya_io::phraya::PhrayaFile::new(
            200,
            "ref".to_string(),
            "2024-01-01T00:00:00Z".to_string(),
            vec![obs_a.clone(), obs_b.clone()],
            coverage,
        )
        .with_coverage_breadth(100, 0);

        let output = format_snpdiffs(
            &phraya_file,
            &[obs_a.clone(), obs_b.clone()],
            Some("ref"),
            None, // No query_id → falls back to first provenance "contigA"
            Some(&ref_fasta),
            Some(&query_fasta),
        )
        .unwrap();

        let lines: Vec<&str> = output.lines().collect();

        // Header uses "contigA" as query_id (from first obs provenance)
        let header_fields: Vec<&str> = lines[0].split('\t').collect();
        assert!(header_fields.iter().any(|f| f == &"Query_ID:contigA"));

        // Two variant rows
        let var_lines: Vec<&str> = lines.iter().copied().filter(|l| !l.starts_with("#\t")).collect();
        assert_eq!(var_lines.len(), 3); // 1 BED + 2 variants

        let fields_a: Vec<&str> = var_lines[1].split('\t').collect();
        let fields_b: Vec<&str> = var_lines[2].split('\t').collect();

        // Each variant row has its own query contig from provenance
        assert_eq!(fields_a[3], "contigA");
        assert_eq!(fields_b[3], "contigB");

        // Dist_to_Query_End uses per-contig lengths:
        // contigA (len=100): min(51, 100-51) = min(51, 49) = 49
        // contigB (len=80): min(51, 80-51) = min(51, 29) = 29
        assert_eq!(fields_a[15], "49");
        assert_eq!(fields_b[15], "29");

        // Query_Contig field in BED row uses first obs contig
        let bed_fields: Vec<&str> = var_lines[0].split('\t').collect();
        assert_eq!(bed_fields[6], "contigA");
    }
}
