//! CSP2-compatible `.snpdiffs` output (the `dnadiff` replacement, ADR-0015).
//!
//! Produces the text format `compileMUMmer.py` emits and `csp2_common.py`'s
//! `parseSNPDiffs()` consumes: a `#` header of `key:value` fields, `##` BED rows, and one
//! 21-column row per variant. It is computed from the alignment blocks `phraya align
//! --strategy assembly` stores in each `.phraya` file:
//!
//! - blocks are selected one-to-one across all inputs ([`crate::blocks`], MUMmer's
//!   `delta-filter -1`), so every reference and query base is explained at most once;
//! - one BED row per kept block, plus `.` rows for uncovered reference and query regions;
//! - variant rows only for observations inside a kept block, with that block's coordinates
//!   and identity, query positions in forward query coordinates;
//! - header coverage is kept bases over assembly bases, medians are over kept blocks.

use crate::blocks::{select_one_to_one, KeptBlock, SpaceBlocks};
use phraya_core::types::{AlignmentBlock, Strand, VariantObservation, VariantType};
use phraya_io::phraya::PhrayaFile;
use phraya_io::plan::content_hash_for_bytes;
use phraya_io::SequenceParser;
use std::collections::{HashMap, HashSet};
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
    /// (contig_name, length) pairs, in file order.
    contig_lengths: Vec<(String, usize)>,
}

/// Compute FASTA metadata: contig count, total bases, N50/N90, L50/L90, SHA-256.
///
/// SHA-256 is computed over the **raw file bytes** (including FASTA headers and
/// whitespace), matching `compileMUMmer.py`'s `fasta_info()` which does
/// `hashlib.sha256(file.read()).hexdigest()` on the raw file.
fn compute_fasta_meta(path: &Path) -> Result<FastaMeta, String> {
    let raw_bytes = std::fs::read(path)
        .map_err(|e| format!("failed to read FASTA {}: {}", path.display(), e))?;
    let sha256 = content_hash_for_bytes(&raw_bytes);

    let mut parser = SequenceParser::from_path(path).map_err(|e| e.to_string())?;
    let mut contig_lengths: Vec<(String, usize)> = Vec::new();
    while let Some(result) = parser.next() {
        let seq = result.map_err(|e| e.to_string())?;
        contig_lengths.push((seq.id().to_string(), seq.len()));
    }

    let total: usize = contig_lengths.iter().map(|(_, l)| *l).sum();
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

/// Median of `v`, or `None` if empty.
fn median(mut v: Vec<f64>) -> Option<f64> {
    if v.is_empty() {
        return None;
    }
    v.sort_by(f64::total_cmp);
    let mid = v.len() / 2;
    Some(if v.len() % 2 == 1 { v[mid] } else { (v[mid - 1] + v[mid]) / 2.0 })
}

/// `{:.2}`, or `NA` when the value is undefined (as `compileMUMmer.py` writes it).
fn two_decimals_or_na(v: Option<f64>) -> String {
    v.map_or_else(|| "NA".to_string(), |x| format!("{x:.2}"))
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
struct Summary {
    snps: usize,
    indels: usize,
    percent_ref_aligned: Option<f64>,
    percent_query_aligned: Option<f64>,
    median_identity: Option<f64>,
    median_alignment_length: Option<f64>,
}

/// Format the `#` header line: `#`, query metadata, reference metadata, summary — all
/// tab-separated `key:value` fields in `compileMUMmer.py`'s order. Structural and k-mer
/// fields are not computed and stay `0` / `NA`.
fn format_header(
    query_meta: &Option<FastaMeta>,
    query_id: &str,
    query_fasta_path: &str,
    ref_meta: &Option<FastaMeta>,
    ref_id: &str,
    ref_fasta_path: &str,
    s: &Summary,
) -> String {
    let query_str = format_fasta_meta_section(query_meta, query_id, query_fasta_path, true);
    let ref_str = format_fasta_meta_section(ref_meta, ref_id, ref_fasta_path, false);
    let summary_str = format!(
        "SNPs:{}\tReference_Percent_Aligned:{}\tQuery_Percent_Aligned:{}\tMedian_Percent_Identity:{}\tMedian_Alignment_Length:{}\tKmer_Similarity:NA\tShared_Kmers:NA\tReference_Unique_Kmers:NA\tQuery_Unique_Kmers:NA\tReference_Breakpoints:0\tQuery_Breakpoints:0\tReference_Relocations:0\tQuery_Relocations:0\tReference_Translocations:0\tQuery_Translocations:0\tReference_Inversions:0\tQuery_Inversions:0\tReference_Insertions:0\tQuery_Insertions:0\tReference_Tandem:0\tQuery_Tandem:0\tIndels:{}\tInvalid:0\tgSNPs:NA\tgIndels:NA",
        s.snps,
        two_decimals_or_na(s.percent_ref_aligned),
        two_decimals_or_na(s.percent_query_aligned),
        two_decimals_or_na(s.median_identity),
        two_decimals_or_na(s.median_alignment_length),
        s.indels,
    );
    format!("#\t{}\t{}\t{}", query_str, ref_str, summary_str)
}

/// `##` BED row for a kept block. Columns: Ref_Contig, Ref_Start, Ref_End, Ref_Length,
/// Ref_Aligned, Query_Contig, Query_Start, Query_End, Query_Length, Query_Aligned,
/// Perc_Iden; starts 0-based, ends exclusive, as `compileMUMmer.py` writes `.1coords`.
fn format_block_row(kept: &KeptBlock, ref_contig: &str, ref_len: u32) -> String {
    let b = &kept.block;
    format!(
        "##\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}",
        ref_contig,
        b.t_start,
        b.t_end,
        ref_len,
        b.t_end - b.t_start,
        b.query_id,
        b.q_start,
        b.q_end,
        b.query_len,
        b.q_end - b.q_start,
        kept.identity * 100.0,
    )
}

/// One tab-separated variant row (21 columns, `csp2_common.py`'s `snp_columns`):
/// Ref_Contig, Start_Ref, Ref_Pos, Query_Contig, Start_Query, Query_Pos, Ref_Loc, Query_Loc,
/// Ref_Start, Ref_End, Query_Start, Query_End, Ref_Base, Query_Base, Dist_to_Ref_End,
/// Dist_to_Query_End, Ref_Aligned, Query_Aligned, Query_Direction, Perc_Iden, Cat.
///
/// Alignment columns come from the kept block containing the variant; `Query_Pos` is the
/// forward query coordinate. `Query_Base` is in reference orientation (the complement of the
/// query FASTA base on `-1` rows), as `show-snps` reports it.
fn format_variant_row(obs: &VariantObservation, kept: &KeptBlock, ref_contig: &str, ref_len: u32) -> String {
    let b: &AlignmentBlock = &kept.block;
    let ref_pos = obs.position() + 1;
    let query_pos = b.forward_query_position(obs.query_position()) + 1;
    let direction = match b.strand {
        Strand::Forward => "1",
        Strand::Reverse => "-1",
    };
    let cat = match obs.variant_type() {
        VariantType::Snp => "SNP",
        VariantType::Insertion | VariantType::Deletion => "Indel",
    };
    format!(
        "{rc}\t{}\t{ref_pos}\t{qc}\t{}\t{query_pos}\t{rc}/{ref_pos}\t{qc}/{query_pos}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{direction}\t{:.2}\t{cat}",
        ref_pos - 1,
        query_pos - 1,
        b.t_start,
        b.t_end,
        b.q_start,
        b.q_end,
        format_ref_base(obs.ref_base(), obs.variant_type()),
        best_alt_allele(obs.ref_base(), obs.all_alleles(), obs.variant_type()),
        ref_pos.min(ref_len.saturating_sub(ref_pos)),
        query_pos.min(b.query_len.saturating_sub(query_pos)),
        b.t_end - b.t_start,
        b.q_end - b.q_start,
        kept.identity * 100.0,
        rc = ref_contig,
        qc = b.query_id,
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

/// Format `.phraya` files from `phraya align --reference ... --strategy assembly` (one per
/// reference space, `header.sample_id` = reference contig) as one CSP2-compatible snpdiffs
/// block. Each file is paired with the observations to emit (filtering may have dropped
/// some); an observation is reported only if it lies in a block kept by one-to-one selection.
///
/// * `reference_id` / `query_id` — assembly labels for the header (`Reference_ID`,
///   `Query_ID`); default to the FASTA file stems.
/// * `reference_fasta` / `query_fasta` — assemblies for header metadata, coverage
///   denominators and `.` rows. Without the query FASTA, `Query_Percent_Aligned` is `NA` and
///   no uncovered-query rows are written.
pub fn format_snpdiffs(
    files: &[(&PhrayaFile, &[VariantObservation])],
    reference_id: Option<&str>,
    query_id: Option<&str>,
    reference_fasta: Option<&Path>,
    query_fasta: Option<&Path>,
) -> Result<String, SnpdiffsError> {
    if files.is_empty() {
        return Err(SnpdiffsError { message: "snpdiffs needs at least one .phraya file".to_string() });
    }
    if let Some((f, _)) = files.iter().find(|(f, obs)| f.alignments.is_empty() && !obs.is_empty()) {
        return Err(SnpdiffsError {
            message: format!(
                "{} has variant observations but no alignment blocks; snpdiffs is computed from \
                 assembly alignment blocks: re-run `phraya align --reference ... --strategy assembly`",
                f.header.sample_id
            ),
        });
    }

    let read_meta = |path: Option<&Path>| -> Result<Option<FastaMeta>, SnpdiffsError> {
        path.map(|p| compute_fasta_meta(p).map_err(|message| SnpdiffsError { message })).transpose()
    };
    let ref_meta = read_meta(reference_fasta)?;
    let query_meta = read_meta(query_fasta)?;
    let path_str = |p: Option<&Path>| p.map(|p| p.to_string_lossy().to_string()).unwrap_or_default();
    let label = |id: Option<&str>, path: Option<&Path>, default: &str| {
        id.map(str::to_string)
            .or_else(|| path.and_then(|p| p.file_stem()).map(|s| s.to_string_lossy().to_string()))
            .unwrap_or_else(|| default.to_string())
    };
    let ref_label = label(reference_id, reference_fasta, "reference");
    let query_label = label(query_id, query_fasta, "query");

    let spaces: Vec<SpaceBlocks<'_>> = files
        .iter()
        .map(|(f, _)| SpaceBlocks {
            ref_contig: &f.header.sample_id,
            ref_len: f.header.reference_length,
            blocks: &f.alignments,
        })
        .collect();
    let selection = select_one_to_one(&spaces);

    // Variant rows: observations inside a kept block, one per (space, position, kind).
    let mut variant_rows = Vec::new();
    let mut seen = HashSet::new();
    let (mut snps, mut indels) = (0usize, 0usize);
    for (space, (file, observations)) in files.iter().enumerate() {
        for obs in observations.iter() {
            let Some(kept) = selection.block_for(space, obs) else { continue };
            let kind = match obs.variant_type() {
                VariantType::Snp => 0u8,
                VariantType::Insertion => 1,
                VariantType::Deletion => 2,
            };
            if !seen.insert((space, obs.position(), kind)) {
                continue;
            }
            if kind == 0 {
                snps += 1;
            } else {
                indels += 1;
            }
            variant_rows.push(format_variant_row(obs, kept, &file.header.sample_id, file.header.reference_length));
        }
    }

    let ref_total = match &ref_meta {
        Some(m) => m.total_bases as f64,
        None => files.iter().map(|(f, _)| f.header.reference_length as f64).sum(),
    };
    let summary = Summary {
        snps,
        indels,
        percent_ref_aligned: (ref_total > 0.0).then(|| 100.0 * selection.reference_bases_covered() as f64 / ref_total),
        percent_query_aligned: query_meta
            .as_ref()
            .filter(|m| m.total_bases > 0)
            .map(|m| 100.0 * selection.query_bases_covered() as f64 / m.total_bases as f64),
        median_identity: median(selection.kept.iter().map(|k| k.identity * 100.0).collect()),
        median_alignment_length: median(selection.kept.iter().map(|k| (k.block.t_end - k.block.t_start) as f64).collect()),
    };

    let mut out = format_header(
        &query_meta,
        &query_label,
        &path_str(query_fasta),
        &ref_meta,
        &ref_label,
        &path_str(reference_fasta),
        &summary,
    );
    out.push('\n');
    for kept in &selection.kept {
        let (file, _) = files[kept.space];
        out.push_str(&format_block_row(kept, &file.header.sample_id, file.header.reference_length));
        out.push('\n');
    }

    // Uncovered regions, as compileMUMmer.py writes them (consumers drop these rows).
    let space_of: HashMap<&str, usize> =
        files.iter().enumerate().map(|(i, (f, _))| (f.header.sample_id.as_str(), i)).collect();
    let ref_contigs: Vec<(String, u32)> = match &ref_meta {
        Some(m) => m.contig_lengths.iter().map(|(n, l)| (n.clone(), *l as u32)).collect(),
        None => files.iter().map(|(f, _)| (f.header.sample_id.clone(), f.header.reference_length)).collect(),
    };
    for (contig, len) in &ref_contigs {
        let gaps = match space_of.get(contig.as_str()) {
            Some(&space) => selection.uncovered_reference(space),
            None if *len > 0 => vec![(0, *len)],
            None => Vec::new(),
        };
        for (s, e) in gaps {
            out.push_str(&format!("##\t{contig}\t{s}\t{e}\t.\t.\t.\t.\t.\t.\t.\t.\n"));
        }
    }
    if let Some(m) = &query_meta {
        for (contig, len) in &m.contig_lengths {
            for (s, e) in selection.uncovered_query(contig, *len as u32) {
                out.push_str(&format!("##\t.\t.\t.\t.\t.\t{contig}\t{s}\t{e}\t.\t.\t.\n"));
            }
        }
    }

    for row in variant_rows {
        out.push_str(&row);
        out.push('\n');
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use phraya_core::types::CoverageTrack;
    use std::io::Write;
    use std::path::PathBuf;

    fn write_fasta(dir: &Path, name: &str, records: &[(&str, usize)]) -> PathBuf {
        let path = dir.join(name);
        let mut f = std::fs::File::create(&path).unwrap();
        for (id, len) in records {
            writeln!(f, ">{id}").unwrap();
            writeln!(f, "{}", "ACGT".repeat(len / 4 + 1)[..*len].to_string()).unwrap();
        }
        path
    }

    /// A block with mismatches at the given column offsets and optional extra CIGAR.
    fn block(query: &str, query_len: u32, q_start: u32, t_start: u32, cigar: &str, strand: Strand) -> AlignmentBlock {
        let ops = phraya_core::cigar::parse_ops(cigar);
        let t_len: u32 = ops.iter().filter(|(_, op)| "MXI".contains(*op)).map(|(n, _)| *n as u32).sum();
        let q_len: u32 = ops.iter().filter(|(_, op)| "MXD".contains(*op)).map(|(n, _)| *n as u32).sum();
        let edit_distance: u32 = ops.iter().filter(|(_, op)| *op != 'M').map(|(n, _)| *n as u32).sum();
        AlignmentBlock {
            query_id: query.into(),
            query_len,
            q_start,
            q_end: q_start + q_len,
            t_start,
            t_end: t_start + t_len,
            strand,
            edit_distance,
            cigar: cigar.into(),
        }
    }

    fn variant(kind: VariantType, position: u32, ref_base: u8, alt: u8, query: &str, query_position: u32, strand: Strand) -> VariantObservation {
        VariantObservation::new(position, ref_base, HashMap::from([(alt, 1)]), 1.0, "unused".into(), 60, 1, vec![1], 60.0, query.into())
            .with_variant_type(kind)
            .with_query_position(query_position)
            .with_strand(strand)
    }

    fn space(contig: &str, len: u32, observations: Vec<VariantObservation>, blocks: Vec<AlignmentBlock>) -> PhrayaFile {
        PhrayaFile::new(len, contig.into(), "t".into(), observations, CoverageTrack::uncovered(len)).with_alignments(blocks)
    }

    struct Parsed {
        header: HashMap<String, String>,
        bed: Vec<String>,
        variants: Vec<String>,
    }

    fn parse(text: &str) -> Parsed {
        let mut header = HashMap::new();
        let (mut bed, mut variants) = (Vec::new(), Vec::new());
        for line in text.lines() {
            if line.starts_with("##\t") {
                bed.push(line.to_string());
            } else if let Some(rest) = line.strip_prefix("#\t") {
                for f in rest.split('\t') {
                    let (k, v) = f.split_once(':').unwrap();
                    header.insert(k.to_string(), v.to_string());
                }
            } else {
                variants.push(line.to_string());
            }
        }
        Parsed { header, bed, variants }
    }

    fn render(files: &[PhrayaFile], ref_fa: Option<&Path>, qry_fa: Option<&Path>) -> Result<String, SnpdiffsError> {
        let pairs: Vec<(&PhrayaFile, &[VariantObservation])> = files.iter().map(|f| (f, f.observations.as_slice())).collect();
        format_snpdiffs(&pairs, Some("REF"), Some("QRY"), ref_fa, qry_fa)
    }

    #[test]
    fn kept_block_drives_bed_rows_coverage_and_variant_fields() {
        let dir = tempfile::TempDir::new().unwrap();
        let ref_fa = write_fasta(dir.path(), "ref.fa", &[("r0", 1_000), ("r1", 500)]);
        let qry_fa = write_fasta(dir.path(), "qry.fa", &[("q0", 1_200)]);
        let snps = vec![
            variant(VariantType::Snp, 200, b'A', b'C', "q0", 300, Strand::Forward),
            variant(VariantType::Snp, 700, b'G', b'T', "q0", 800, Strand::Forward),
        ];
        let files = [
            space("r0", 1_000, snps, vec![block("q0", 1_200, 100, 0, "200M1X499M1X299M", Strand::Forward)]),
            space("r1", 500, vec![], vec![]),
        ];
        let p = parse(&render(&files, Some(&ref_fa), Some(&qry_fa)).unwrap());

        assert_eq!(p.header["Reference_ID"], "REF");
        assert_eq!(p.header["Query_ID"], "QRY");
        assert_eq!((p.header["SNPs"].as_str(), p.header["Indels"].as_str()), ("2", "0"));
        assert_eq!(p.header["Reference_Percent_Aligned"], "66.67", "1000 of 1500 reference bases");
        assert_eq!(p.header["Query_Percent_Aligned"], "83.33", "1000 of 1200 query bases");
        assert_eq!(p.header["Median_Percent_Identity"], "99.80");
        assert_eq!(p.header["Median_Alignment_Length"], "1000.00");
        assert_eq!(
            p.bed,
            vec![
                "##\tr0\t0\t1000\t1000\t1000\tq0\t100\t1100\t1200\t1000\t99.80",
                "##\tr1\t0\t500\t.\t.\t.\t.\t.\t.\t.\t.",
                "##\t.\t.\t.\t.\t.\tq0\t0\t100\t.\t.\t.",
                "##\t.\t.\t.\t.\t.\tq0\t1100\t1200\t.\t.\t.",
            ]
        );
        assert_eq!(
            p.variants[0],
            "r0\t200\t201\tq0\t300\t301\tr0/201\tq0/301\t0\t1000\t100\t1100\tA\tC\t201\t301\t1000\t1000\t1\t99.80\tSNP"
        );
        assert_eq!(p.variants.len(), 2);
    }

    #[test]
    fn observations_outside_kept_blocks_are_not_reported() {
        let observations = vec![
            variant(VariantType::Snp, 200, b'A', b'C', "q0", 200, Strand::Forward),
            variant(VariantType::Snp, 800, b'A', b'C', "q0", 800, Strand::Forward), // past the block
            variant(VariantType::Snp, 300, b'A', b'C', "other", 300, Strand::Forward), // other query
        ];
        let files = [space("r0", 1_000, observations, vec![block("q0", 1_000, 0, 0, "200M1X299M", Strand::Forward)])];
        let p = parse(&render(&files, None, None).unwrap());
        assert_eq!(p.header["SNPs"], "1");
        assert_eq!(p.variants.len(), 1);
        assert!(p.variants[0].starts_with("r0\t200\t201\t"));
    }

    #[test]
    fn reverse_strand_rows_report_forward_query_positions() {
        // Strand-local offset 100 in the reverse-complemented 1000 bp query is forward 899.
        let obs = variant(VariantType::Snp, 100, b'A', b'C', "q0", 100, Strand::Reverse);
        let files = [space("r0", 1_000, vec![obs], vec![block("q0", 1_000, 0, 0, "100M1X899M", Strand::Reverse)])];
        let p = parse(&render(&files, None, None).unwrap());
        let row: Vec<&str> = p.variants[0].split('\t').collect();
        assert_eq!((row[4], row[5], row[7]), ("899", "900", "q0/900"));
        assert_eq!((row[10], row[11], row[18]), ("0", "1000", "-1"));
        assert_eq!(row[15], "100", "Dist_to_Query_End = min(900, 1000 - 900)");
    }

    #[test]
    fn indels_are_counted_and_formatted() {
        // 2 bp deletion from the query at reference 500-501.
        let del = variant(VariantType::Deletion, 500, b'G', b'.', "q0", 500, Strand::Forward);
        let files = [space("r0", 1_000, vec![del], vec![block("q0", 998, 0, 0, "500M2I498M", Strand::Forward)])];
        let p = parse(&render(&files, None, None).unwrap());
        assert_eq!((p.header["SNPs"].as_str(), p.header["Indels"].as_str()), ("0", "1"));
        let row: Vec<&str> = p.variants[0].split('\t').collect();
        assert_eq!((row[12], row[13], row[20]), ("G", ".", "Indel"));
        assert_eq!((row[16], row[17]), ("1000", "998"));
    }

    #[test]
    fn a_position_is_reported_once() {
        let obs = variant(VariantType::Snp, 200, b'A', b'C', "q0", 200, Strand::Forward);
        let files = [space("r0", 1_000, vec![obs.clone(), obs], vec![block("q0", 1_000, 0, 0, "200M1X299M", Strand::Forward)])];
        let p = parse(&render(&files, None, None).unwrap());
        assert_eq!((p.header["SNPs"].as_str(), p.variants.len()), ("1", 1));
    }

    #[test]
    fn observations_without_alignment_blocks_are_rejected() {
        let obs = variant(VariantType::Snp, 200, b'A', b'C', "q0", 200, Strand::Forward);
        let err = render(&[space("r0", 1_000, vec![obs], vec![])], None, None).unwrap_err();
        assert!(err.message.contains("--strategy assembly"), "{}", err.message);
    }

    /// An observation whose stored query offset is not the offset the kept block's CIGAR pairs
    /// with its reference position must be rejected. This is the §6 failure mode: an off-diagonal
    /// local-repeat block emits an observation at a reference position that the kept full contig
    /// block pairs with a *different* query offset, and the emitter must not adopt it. The spurious
    /// observation is sorted before the true one to check the before-`seen.insert` ordering.
    #[test]
    fn off_diagonal_observation_is_rejected() {
        // Kept block: 100M1X99M, forward. ref 100 (0-based) is the X column, paired with q 100.
        // Off-diagonal observations carry ref 100 at q 130 — not what this block pairs.
        let observations = vec![
            // Spurious: ref 50 / q 130. Block pairs ref 50 with q 50 -> reject.
            variant(VariantType::Snp, 50, b'G', b'A', "q0", 130, Strand::Forward),
            // Spurious: ref 100 / q 130. True offset at ref 100 is q 100 -> reject.
            variant(VariantType::Snp, 100, b'G', b'A', "q0", 130, Strand::Forward),
            // True: ref 100 / q 100. Must survive (and not be deduped by the spurious above).
            variant(VariantType::Snp, 100, b'G', b'A', "q0", 100, Strand::Forward),
        ];
        let files = [space("r0", 200, observations, vec![block("q0", 200, 0, 0, "100M1X99M", Strand::Forward)])];
        let p = parse(&render(&files, None, None).unwrap());
        assert_eq!((p.header["SNPs"].as_str(), p.variants.len()), ("1", 1));
        let row: Vec<&str> = p.variants[0].split('\t').collect();
        assert_eq!((row[0], row[1], row[4], row[5]), ("r0", "100", "100", "101"));
    }

    #[test]
    fn indel_observations_are_accepted_at_consumer_coordinates() {
        // Block "100M2I50M2D48M", t_start=0, q_start=0, forward:
        //   M run t=0..100 q=0..100
        //   I run t=100..102 q=100 (deletion of ref 100,101; obs carries q=100)
        //   M run t=102..152 q=100..150
        //   D run t=152 q=150 (insertion anchored at ref 152; obs carries q=150)
        //   M run t=152..200 q=152..200
        let b = block("q0", 200, 0, 0, "100M2I50M2D48M", Strand::Forward);
        let files = [space("r0", 300, vec![
            // Deletion: first deleted ref base t=100, q=100.
            variant(VariantType::Deletion, 100, b'G', b'.', "q0", 100, Strand::Forward),
            // Deletion off-diagonal (q=130): not the I run's offset -> reject.
            variant(VariantType::Deletion, 100, b'G', b'.', "q0", 130, Strand::Forward),
            // Insertion at anchor ref 152, q=150.
            variant(VariantType::Insertion, 152, b'.', b'A', "q0", 150, Strand::Forward),
        ], vec![b])];
        let p = parse(&render(&files, None, None).unwrap());
        assert_eq!((p.header["SNPs"].as_str(), p.header["Indels"].as_str()), ("0", "2"));
    }

    #[test]
    fn nothing_aligned_is_an_empty_comparison() {
        let dir = tempfile::TempDir::new().unwrap();
        let qry_fa = write_fasta(dir.path(), "qry.fa", &[("q0", 300)]);
        let p = parse(&render(&[space("r0", 1_000, vec![], vec![])], None, Some(&qry_fa)).unwrap());
        assert_eq!(p.header["SNPs"], "0");
        assert_eq!(p.header["Reference_Percent_Aligned"], "0.00");
        assert_eq!(p.header["Query_Percent_Aligned"], "0.00");
        assert_eq!(p.header["Median_Percent_Identity"], "NA");
        assert_eq!(p.bed, vec!["##\tr0\t0\t1000\t.\t.\t.\t.\t.\t.\t.\t.", "##\t.\t.\t.\t.\t.\tq0\t0\t300\t.\t.\t."]);
        assert!(p.variants.is_empty());
    }

    #[test]
    fn query_coverage_is_na_without_the_query_fasta() {
        let files = [space("r0", 1_000, vec![], vec![block("q0", 1_000, 0, 0, "500M", Strand::Forward)])];
        let p = parse(&render(&files, None, None).unwrap());
        assert_eq!(p.header["Query_Percent_Aligned"], "NA");
        assert_eq!(p.header["Reference_Percent_Aligned"], "50.00", "falls back to the space lengths");
        assert_eq!(p.header["Query_Assembly_Bases"], "NA");
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
}
