use crate::VariantObservation;
use phraya_core::types::VariantType;
use std::collections::BTreeMap;

/// Generate VCF 4.2 header
pub fn vcf_header(reference_name: &str, reference_length: u32) -> String {
    let mut header = String::new();
    header.push_str("##fileformat=VCFv4.2\n");
    header.push_str(&format!(
        "##contig=<ID={},length={}>\n",
        reference_name, reference_length
    ));
    header.push_str("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">\n");
    header.push_str("##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"Mapping quality\">\n");
    header.push_str("##INFO=<ID=CIGAR,Number=1,Type=String,Description=\"CIGAR string\">\n");
    header.push_str("##INFO=<ID=ED,Number=1,Type=Integer,Description=\"Edit distance\">\n");
    header.push_str("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    header.push_str("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n");
    header.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n");
    header
}

/// Convert VariantObservations to VCF records
pub fn format_vcf(
    observations: impl Iterator<Item = VariantObservation>,
    reference_name: &str,
    reference_length: u32,
) -> String {
    let mut output = vcf_header(reference_name, reference_length);

    // Group observations by position for multi-allelic handling
    let mut by_position: BTreeMap<u32, Vec<VariantObservation>> = BTreeMap::new();
    for obs in observations {
        by_position
            .entry(obs.position())
            .or_insert_with(Vec::new)
            .push(obs);
    }

    // Convert grouped observations to VCF records
    for (position, obs_list) in by_position {
        let primary = &obs_list[0];
        let chrom = reference_name;
        let variant_type = primary.variant_type();

        // Handle indel-specific positioning and REF/ALT encoding
        let (vcf_pos, ref_seq, alt_seq) = match variant_type {
            VariantType::Snp => {
                // SNPs use the position directly
                let pos = position + 1; // VCF is 1-indexed
                let ref_base = (primary.ref_base() as char).to_string();

                // Collect all unique ALT alleles
                let mut alt_bases = std::collections::HashSet::new();
                for obs in &obs_list {
                    for (&allele, _) in obs.all_alleles().iter() {
                        if allele != obs.ref_base() {
                            alt_bases.insert(allele as char);
                        }
                    }
                }

                let alt = if alt_bases.is_empty() {
                    ".".to_string()
                } else {
                    let mut alts: Vec<char> = alt_bases.into_iter().collect();
                    alts.sort();
                    alts.iter().collect::<String>()
                };

                (pos, ref_base, alt)
            }
            VariantType::Deletion => {
                // For deletion: ref_base contains the deleted bases, alt is "."
                // In VCF format, we encode as ref=deleted_bases, alt="."
                let pos = position + 1; // VCF is 1-indexed
                let ref_base = (primary.ref_base() as char).to_string();
                let alt = ".".to_string();
                (pos, ref_base, alt)
            }
            VariantType::Insertion => {
                // For insertion: ref_base is ".", alt contains inserted bases
                // In VCF format, we encode as ref=".", alt=inserted_bases
                let pos = position + 1; // VCF is 1-indexed
                let ref_base = ".".to_string();

                // Collect inserted bases from alleles
                let mut inserted = String::new();
                for (&allele, _) in primary.all_alleles().iter() {
                    if allele != b'.' {
                        inserted.push(allele as char);
                    }
                }

                let alt = if inserted.is_empty() {
                    ".".to_string()
                } else {
                    inserted
                };

                (pos, ref_base, alt)
            }
        };

        let id = ".";
        let qual = format!("{:.0}", primary.confidence() * 100.0);
        let filter = "PASS";

        // Build INFO field
        let coverage = primary.coverage_at_variant().unwrap_or(0);
        let mapq = primary.mapq() as u32;
        let cigar = primary.cigar();
        let edit_dist = primary.edit_distance();

        let info = format!(
            "DP={};MQ={};CIGAR={};ED={}",
            coverage, mapq, cigar, edit_dist
        );

        let format = "GT:GQ";
        let sample = format!("0/1:{}", (primary.confidence() * 100.0) as u32);

        output.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            chrom, vcf_pos, id, ref_seq, alt_seq, qual, filter, info, format, sample
        ));
    }

    output
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    fn create_observation(
        position: u32,
        ref_base: u8,
        alleles: HashMap<u8, u32>,
        mapq: u8,
    ) -> VariantObservation {
        VariantObservation::new(
            position,
            ref_base,
            alleles,
            0.95,
            "10M".to_string(),
            mapq,
            0,
            vec![10],
            35.0,
            "sample:read".to_string(),
        )
    }

    #[test]
    fn vcf_header_format() {
        let header = vcf_header("chr1", 1000);

        assert!(header.contains("##fileformat=VCFv4.2"));
        assert!(header.contains("##contig=<ID=chr1,length=1000>"));
        assert!(header.contains("##INFO=<ID=DP"));
        assert!(header.contains("##INFO=<ID=MQ"));
        assert!(header.contains("##INFO=<ID=CIGAR"));
        assert!(header.contains("##INFO=<ID=ED"));
        assert!(header.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"));
    }

    #[test]
    fn single_snp() {
        let mut alleles = HashMap::new();
        alleles.insert(b'A', 10);
        alleles.insert(b'T', 5);

        let obs = create_observation(99, b'A', alleles, 60);
        let vcf = format_vcf(std::iter::once(obs), "chr1", 1000);

        // Check for the record (1-indexed position)
        let record_line = vcf.lines().find(|l| l.starts_with("chr1\t100\t"));
        assert!(record_line.is_some());

        let record = record_line.unwrap();
        // Check basic structure: should have both REF A and ALT T
        assert!(record.contains("\tA\t")); // REF field
        assert!(record.contains("\tT\t")); // ALT field
    }

    #[test]
    fn multi_allelic() {
        let mut alleles = HashMap::new();
        alleles.insert(b'A', 10);
        alleles.insert(b'C', 3);
        alleles.insert(b'G', 2);

        let obs = create_observation(49, b'A', alleles, 60);
        let vcf = format_vcf(std::iter::once(obs), "chr1", 1000);

        let record_line = vcf.lines().find(|l| l.starts_with("chr1\t50\t"));
        assert!(record_line.is_some());

        let record = record_line.unwrap();
        // ALT should contain both C and G (may be in any order)
        let has_c = record.contains("C");
        let has_g = record.contains("G");
        assert!(has_c && has_g);
    }

    #[test]
    fn no_variants() {
        let vcf = format_vcf(std::iter::empty(), "chr1", 1000);

        // Should have header but no records
        let lines: Vec<&str> = vcf.lines().collect();
        assert!(lines.len() > 0);
        assert!(lines.last().unwrap().starts_with("#CHROM"));
    }

    #[test]
    fn vcf_position_1indexed() {
        let mut alleles = HashMap::new();
        alleles.insert(b'A', 10);

        let obs = create_observation(0, b'A', alleles, 60); // 0-indexed
        let vcf = format_vcf(std::iter::once(obs), "chr1", 1000);

        let record_line = vcf.lines().find(|l| l.starts_with("chr1\t1\t"));
        assert!(record_line.is_some()); // Should be position 1 in VCF (1-indexed)
    }

    #[test]
    fn info_field_present() {
        let mut alleles = HashMap::new();
        alleles.insert(b'A', 10);

        let obs = create_observation(99, b'A', alleles, 50);
        let vcf = format_vcf(std::iter::once(obs), "chr1", 1000);

        let record_line = vcf.lines().find(|l| l.starts_with("chr1\t100\t"));
        let record = record_line.unwrap();

        // Check INFO field contains expected data
        assert!(record.contains("DP="));
        assert!(record.contains("MQ="));
        assert!(record.contains("CIGAR="));
        assert!(record.contains("ED="));
    }

    #[test]
    fn quality_encoding() {
        let mut alleles = HashMap::new();
        alleles.insert(b'A', 10);

        // Create observation with 0.80 confidence -> 80 QUAL
        let obs = VariantObservation::new(
            50,
            b'A',
            alleles,
            0.80,
            "10M".to_string(),
            60,
            0,
            vec![10],
            35.0,
            "sample:read".to_string(),
        );

        let vcf = format_vcf(std::iter::once(obs), "chr1", 1000);
        let record_line = vcf.lines().find(|l| l.starts_with("chr1\t51\t"));

        // Quality should be encoded as confidence * 100
        assert!(record_line.unwrap().contains("\t80\t"));
    }

    #[test]
    fn format_field() {
        let mut alleles = HashMap::new();
        alleles.insert(b'A', 10);

        let obs = create_observation(99, b'A', alleles, 60);
        let vcf = format_vcf(std::iter::once(obs), "chr1", 1000);

        let record_line = vcf.lines().find(|l| l.starts_with("chr1\t100\t"));
        let record = record_line.unwrap();

        // FORMAT and SAMPLE fields should be present
        let parts: Vec<&str> = record.split('\t').collect();
        assert!(parts.len() >= 10);
        assert_eq!(parts[8], "GT:GQ"); // FORMAT field
        assert!(parts[9].contains("0/1")); // SAMPLE field
    }

    #[test]
    fn multiple_positions() {
        let mut alleles1 = HashMap::new();
        alleles1.insert(b'A', 10);
        let obs1 = create_observation(49, b'A', alleles1, 60);

        let mut alleles2 = HashMap::new();
        alleles2.insert(b'C', 5);
        let obs2 = create_observation(99, b'C', alleles2, 50);

        let vcf = format_vcf(vec![obs1, obs2].into_iter(), "chr1", 1000);

        let record_count = vcf.lines().filter(|l| l.starts_with("chr1")).count();
        assert_eq!(record_count, 2);
    }

    #[test]
    fn deletion_variant_encoding() {
        let mut alleles = HashMap::new();
        alleles.insert(b'A', 10);

        let obs = create_observation(99, b'A', alleles, 60).with_variant_type(VariantType::Deletion);
        let vcf = format_vcf(std::iter::once(obs), "chr1", 1000);

        let record_line = vcf.lines().find(|l| l.starts_with("chr1\t100\t")).unwrap();
        let parts: Vec<&str> = record_line.split('\t').collect();
        // REF carries the deleted base(s), ALT is "."
        assert_eq!(parts[3], "A");
        assert_eq!(parts[4], ".");
    }

    #[test]
    fn insertion_variant_encoding() {
        let mut alleles = HashMap::new();
        alleles.insert(b'T', 10);

        let obs = create_observation(99, b'.', alleles, 60).with_variant_type(VariantType::Insertion);
        let vcf = format_vcf(std::iter::once(obs), "chr1", 1000);

        let record_line = vcf.lines().find(|l| l.starts_with("chr1\t100\t")).unwrap();
        let parts: Vec<&str> = record_line.split('\t').collect();
        // REF is ".", ALT carries the inserted base(s)
        assert_eq!(parts[3], ".");
        assert_eq!(parts[4], "T");
    }

    #[test]
    fn insertion_with_all_dot_alleles_yields_dot_alt() {
        let mut alleles = HashMap::new();
        alleles.insert(b'.', 10);

        let obs = create_observation(99, b'.', alleles, 60).with_variant_type(VariantType::Insertion);
        let vcf = format_vcf(std::iter::once(obs), "chr1", 1000);

        let record_line = vcf.lines().find(|l| l.starts_with("chr1\t100\t")).unwrap();
        let parts: Vec<&str> = record_line.split('\t').collect();
        assert_eq!(parts[3], ".");
        assert_eq!(parts[4], ".");
    }
}
