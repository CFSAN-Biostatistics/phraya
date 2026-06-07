use clap::{Parser, Subcommand};
use log::info;
use phraya_align::executor::{align_task_with_config, AlignConfig, Strategy};
use phraya_core::types::{
    compute_kmer_uniqueness, detect_hotspot_intervals, select_centroid, sketch_sequence_default,
    CoverageTrack, MinimizerSketch, Sequence,
};
use phraya_filter::{vcf, FilterBuilder, FilterPreset};
use phraya_io::{
    phraya,
    plan::{self, write_plan, PhrayaPlan, UseCase},
    queries,
    SequenceParser,
};
use std::collections::HashMap;
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "phraya")]
#[command(about = "Phraya: pairwise sequence aligner for bacterial genomics")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Plan alignment tasks from input sequences
    Plan {
        /// Input files (FASTA/FASTQ/BAM/CRAM)
        #[arg(long, value_name = "FILE", required = true)]
        inputs: Vec<PathBuf>,

        /// Reference file (optional, auto-detects use case)
        #[arg(long, value_name = "FILE")]
        reference: Option<PathBuf>,

        /// Output plan file
        #[arg(long, value_name = "FILE", required = true)]
        output: PathBuf,
    },
    /// Extract task list from a .phrayaplan file and output as TSV
    PlanTasks {
        /// Path to the .phrayaplan file
        #[arg(value_name = "PLAN_FILE")]
        plan_file: PathBuf,
    },
    /// Merge multiple .phraya files into a single file
    Merge {
        /// Input .phraya files to merge
        #[arg(value_name = "FILE", required = true)]
        inputs: Vec<PathBuf>,

        /// Output merged .phraya file
        #[arg(long, value_name = "FILE", required = true)]
        output: PathBuf,
    },
    /// Align a query sequence against a target using a plan file
    Align {
        /// Plan file (.phrayaplan)
        #[arg(value_name = "PLAN_FILE")]
        plan_file: PathBuf,

        /// Query sequence ID
        #[arg(value_name = "QUERY_ID")]
        query_id: String,

        /// Target sequence ID
        #[arg(value_name = "TARGET_ID")]
        target_id: String,

        /// Output .phraya file
        #[arg(long, value_name = "FILE", required = true)]
        output: PathBuf,

        /// Alignment strategy: fast (±150bp window), balanced (±50bp, default), exact (±25bp)
        #[arg(long, value_name = "STRATEGY", default_value = "balanced")]
        strategy: String,
    },
    /// Filter .phraya file by thresholds and output in specified format
    Filter {
        /// Input .phraya file
        #[arg(value_name = "INPUT")]
        input: PathBuf,

        /// Minimum coverage threshold
        #[arg(long, value_name = "N")]
        min_coverage: Option<u32>,

        /// Maximum coverage threshold
        #[arg(long, value_name = "N")]
        max_coverage: Option<u32>,

        /// Minimum MAPQ threshold
        #[arg(long, value_name = "N")]
        min_mapq: Option<u8>,

        /// Maximum MAPQ threshold
        #[arg(long, value_name = "N")]
        max_mapq: Option<u8>,

        /// Output format (vcf, tsv, phraya)
        #[arg(long, value_name = "FORMAT", default_value = "vcf")]
        format: String,

        /// Output file (required for phraya format)
        #[arg(long, value_name = "FILE")]
        output: Option<PathBuf>,

        /// Named filter preset (conservative or sensitive). Individual threshold flags override preset values.
        #[arg(long, value_name = "PRESET")]
        preset: Option<String>,
    },
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::builder().is_test(false).try_init().ok();

    let cli = Cli::parse();

    match cli.command {
        Commands::Plan {
            inputs,
            reference,
            output,
        } => {
            run_plan(&inputs, reference.as_deref(), &output)?;
        }
        Commands::PlanTasks { plan_file } => {
            plan_tasks(&plan_file)?;
        }
        Commands::Align {
            plan_file,
            query_id,
            target_id,
            output,
            strategy,
        } => {
            let strat = match strategy.as_str() {
                "fast" => Strategy::Fast,
                "balanced" => Strategy::Balanced,
                "exact" => Strategy::Exact,
                other => return Err(format!("unknown strategy: {other}; expected fast, balanced, or exact").into()),
            };
            run_align(&plan_file, &query_id, &target_id, &output, AlignConfig::new(strat))?;
        }
        Commands::Merge { inputs, output } => {
            run_merge(&inputs, &output)?;
        }
        Commands::Filter {
            input,
            min_coverage,
            max_coverage,
            min_mapq,
            max_mapq,
            format,
            output,
            preset,
        } => {
            run_filter(
                &input,
                min_coverage,
                max_coverage,
                min_mapq,
                max_mapq,
                &format,
                output.as_deref(),
                preset.as_deref(),
            )?;
        }
    }

    Ok(())
}

fn run_align(
    plan_path: &std::path::Path,
    query_id: &str,
    target_id: &str,
    output_path: &std::path::Path,
    config: AlignConfig,
) -> Result<(), Box<dyn std::error::Error>> {
    let plan = plan::read_plan(plan_path)?;

    // Read all sequences from every input file listed in the plan.
    let mut seqs: HashMap<String, Sequence> = HashMap::new();
    for file_path in &plan.input_files {
        let mut parser = SequenceParser::from_path(file_path)?;
        while let Some(result) = parser.next() {
            let seq = result?;
            seqs.insert(seq.id().to_string(), seq);
        }
    }

    let query = seqs
        .get(query_id)
        .ok_or_else(|| format!("unknown query_id: {query_id}"))?;
    let target = seqs
        .get(target_id)
        .ok_or_else(|| format!("unknown target_id: {target_id}"))?;

    eprintln!("Aligning {query_id} to {target_id}");

    let result = align_task_with_config(query, target, &plan, &config)
        .ok_or_else(|| format!("alignment failed for {query_id} vs {target_id}"))?;

    // Build .phraya file
    let coverage = CoverageTrack::new(result.coverage_track.iter().map(|&v| v as usize).collect());
    let phraya_file = phraya::PhrayaFile::new(
        target.len() as u32,
        query_id.to_string(),
        chrono::Utc::now().to_rfc3339(),
        result.variants,
        coverage,
    );
    phraya::write_phraya(output_path, &phraya_file)?;

    // Build .phraya.queries sidecar
    let mut index = queries::QueryIndex::new();
    index.insert(query_id.to_string(), result.query_positions);
    let queries_path = {
        let mut p = output_path.as_os_str().to_owned();
        p.push(".queries");
        std::path::PathBuf::from(p)
    };
    queries::write_queries(&queries_path, &index)?;

    Ok(())
}

fn run_plan(
    input_paths: &[PathBuf],
    reference_path: Option<&std::path::Path>,
    output_path: &PathBuf,
) -> Result<(), Box<dyn std::error::Error>> {
    // Read all input sequences and track which file they came from
    let mut all_sequences: Vec<(Sequence, String)> = Vec::new();
    let mut sequence_to_file_index: Vec<usize> = Vec::new(); // Track which input file each sequence came from
    let mut input_file_list = Vec::new();
    let mut ref_seq_index: Option<usize> = None;

    // First, read reference if provided
    if let Some(ref_path) = reference_path {
        let ref_path_str = ref_path.to_string_lossy().to_string();
        input_file_list.push(ref_path_str.clone());

        let mut parser = SequenceParser::from_path(ref_path)?;
        while let Some(seq_result) = parser.next() {
            let seq = seq_result?;
            ref_seq_index = Some(all_sequences.len());
            all_sequences.push((seq, ref_path_str.clone()));
            sequence_to_file_index.push(0); // File index 0 for reference
            break; // Take only the first sequence as reference
        }
    }

    // Read input sequences
    for (file_idx, input_path) in input_paths.iter().enumerate() {
        let input_path_str = input_path.to_string_lossy().to_string();
        input_file_list.push(input_path_str.clone());

        let mut parser = SequenceParser::from_path(input_path)?;
        while let Some(seq_result) = parser.next() {
            let seq = seq_result?;
            all_sequences.push((seq, input_path_str.clone()));
            // File indices: 1+ for input files (offset by 1 if there's a reference)
            let file_index = if ref_seq_index.is_some() { file_idx + 1 } else { file_idx };
            sequence_to_file_index.push(file_index);
        }
    }

    if all_sequences.is_empty() {
        return Err("No sequences found in input files".into());
    }

    // Compute sketches for all sequences (ordered Vec for task generation)
    let sketches: Vec<MinimizerSketch> = all_sequences
        .iter()
        .map(|(seq, _)| sketch_sequence_default(seq))
        .collect();

    // Detect use case
    let use_case = detect_use_case(
        reference_path.is_some(),
        input_paths.len(),
        &sketches,
        &all_sequences,
    );

    eprintln!(
        "Detected Case {}: {:?}",
        use_case_number(&use_case),
        format!("{:?}", use_case)
    );

    // Compute k-mer uniqueness
    let kmer_uniqueness = compute_kmer_uniqueness(&sketches);

    // Detect hotspot intervals from k-mer uniqueness (threshold 0.5)
    let hotspot_intervals = detect_hotspot_intervals(&kmer_uniqueness, 0.5);

    // Generate task list based on use case
    let task_list = generate_task_list(
        &use_case,
        input_paths.len(),
        reference_path.is_some(),
        &sketches,
        &sequence_to_file_index,
        ref_seq_index,
    );

    // Build HashMap kmer_index keyed by sequence ID for reuse during alignment
    let kmer_index: HashMap<String, MinimizerSketch> = all_sequences
        .iter()
        .zip(sketches.iter())
        .map(|((seq, _), sketch)| (seq.id().to_string(), sketch.clone()))
        .collect();

    // Create and write plan
    let mut plan = PhrayaPlan::new(
        use_case,
        input_file_list,
        chrono::Utc::now().to_rfc3339_opts(chrono::SecondsFormat::Secs, true),
        kmer_index,
        kmer_uniqueness,
        task_list,
    );
    plan.hotspot_intervals = hotspot_intervals;

    write_plan(output_path, &plan)?;

    info!("Plan written to {:?}", output_path);
    Ok(())
}

fn plan_tasks(plan_file: &PathBuf) -> Result<(), Box<dyn std::error::Error>> {
    // Read the plan file
    let plan = plan::read_plan(plan_file)?;

    // Output TSV header
    println!("query_id\ttarget_id");

    // Output each task as a TSV line
    for (query_id, target_id) in &plan.task_list {
        println!("{}\t{}", query_id, target_id);
    }

    Ok(())
}

/// Detect the use case from input characteristics
fn detect_use_case(
    has_reference: bool,
    num_input_files: usize,
    _sketches: &[MinimizerSketch],
    all_sequences: &[(Sequence, String)],
) -> UseCase {
    // Count sequences from input files only (exclude reference if provided)
    let num_input_sequences = all_sequences.len() - if has_reference { 1 } else { 0 };

    if has_reference {
        // Reference is stored first; check whether the remaining inputs are contigs (≥5kb).
        let inputs_are_contigs = all_sequences
            .iter()
            .skip(1) // skip reference
            .all(|(seq, _)| seq.len() >= 5000);

        if inputs_are_contigs && num_input_sequences > 0 {
            UseCase::ContigsOnly // Case 4: contigs + reference
        } else {
            UseCase::ReadsWithRef // Case 2: reads + reference
        }
    } else if num_input_files > 1 || (num_input_files == 1 && num_input_sequences > 1) {
        // Case 3 or 4: We have multiple sequences without reference
        // Case 3: contigs + reads (multiple sequences in input files)
        // Case 4: contigs only (but we still treat as contigs for simplicity)

        // Simple heuristic: if we have exactly one input file with multiple sequences,
        // it's likely contigs. If multiple input files, likely contigs + reads.
        // For now, we classify as ContigsWithReads if we have multiple files,
        // otherwise ContigsOnly
        if num_input_files > 1 {
            UseCase::ContigsWithReads
        } else {
            UseCase::ContigsOnly
        }
    } else {
        // Reads only, no reference — not supported
        UseCase::ReadsOnly
    }
}

/// Generate task list based on use case and input characteristics
fn generate_task_list(
    use_case: &UseCase,
    _num_input_files: usize,
    has_reference: bool,
    sketches: &[MinimizerSketch],
    sequence_to_file_index: &[usize],
    ref_seq_index: Option<usize>,
) -> Vec<(u32, u32)> {
    match use_case {
        UseCase::ReadsWithRef => {
            // Case 2: N reads + 1 reference
            // Tasks: (query_id, target_id) where target is always 0 (reference)
            let mut tasks = Vec::new();
            let num_reads = sketches.len() - 1; // All sequences except reference
            for read_id in 1..=num_reads {
                tasks.push((read_id as u32, 0));
            }
            tasks
        }
        UseCase::ContigsWithReads => {
            // Case 3: M contigs + N reads, select centroid from contigs
            // Contigs are in the first input file (or first two if there's a reference)
            let mut tasks = Vec::new();

            // Determine contig indices (all sequences from the first input file)
            let first_input_file_index = if ref_seq_index.is_some() { 1 } else { 0 };
            let contig_indices: Vec<usize> = sequence_to_file_index
                .iter()
                .enumerate()
                .filter(|(_, &file_idx)| file_idx == first_input_file_index)
                .map(|(seq_idx, _)| seq_idx)
                .collect();

            if contig_indices.is_empty() {
                return tasks; // No contigs found
            }

            // Select centroid from contigs only
            let contig_sketches: Vec<_> = contig_indices
                .iter()
                .map(|&i| sketches[i].clone())
                .collect();
            let centroid_offset = select_centroid(&contig_sketches)
                .unwrap_or(0);
            let centroid_idx = contig_indices[centroid_offset];

            eprintln!("Case 3: Selected centroid contig at index {}", centroid_idx);

            // Generate tasks: all other contigs and reads to centroid
            for (seq_idx, _) in sketches.iter().enumerate() {
                if seq_idx != centroid_idx {
                    tasks.push((seq_idx as u32, centroid_idx as u32));
                }
            }

            tasks
        }
        UseCase::ContigsOnly => {
            let mut tasks = Vec::new();
            if has_reference {
                // Case 4 with reference: reference is at index 0, contigs at 1..M
                for contig_id in 1..sketches.len() {
                    tasks.push((contig_id as u32, 0));
                }
            } else {
                // Case 4 MSA: all-vs-all pairs
                for i in 0..sketches.len() {
                    for j in (i + 1)..sketches.len() {
                        tasks.push((i as u32, j as u32));
                    }
                }
            }
            tasks
        }
        UseCase::ReadsOnly => {
            // Reads only, no reference — not supported
            // For now, return empty task list
            Vec::new()
        }
    }
}

fn use_case_number(use_case: &UseCase) -> u32 {
    match use_case {
        UseCase::ReadsWithRef => 2,
        UseCase::ReadsOnly => 1,
        UseCase::ContigsWithReads => 3,
        UseCase::ContigsOnly => 4,
    }
}

fn run_merge(
    input_paths: &[PathBuf],
    output_path: &PathBuf,
) -> Result<(), Box<dyn std::error::Error>> {
    if input_paths.is_empty() {
        return Err("No input files specified".into());
    }

    eprintln!("Merging {} samples...", input_paths.len());

    // Convert PathBuf references to &Path for the merge function
    let paths: Vec<&std::path::Path> = input_paths.iter().map(|p| p.as_path()).collect();

    // Merge the files
    let merged = phraya::merge_phraya_files(&paths)?;

    // Write the merged file
    phraya::write_phraya(output_path, &merged)?;

    eprintln!("Merged file written to {:?}", output_path);

    Ok(())
}

fn run_filter(
    input_path: &PathBuf,
    min_coverage: Option<u32>,
    max_coverage: Option<u32>,
    min_mapq: Option<u8>,
    max_mapq: Option<u8>,
    format: &str,
    output_path: Option<&std::path::Path>,
    preset: Option<&str>,
) -> Result<(), Box<dyn std::error::Error>> {
    // Validate format
    if !["vcf", "tsv", "phraya"].contains(&format) {
        return Err(format!(
            "Invalid format '{}'. Must be one of: vcf, tsv, phraya",
            format
        )
        .into());
    }

    // Read the .phraya file
    let phraya_file = phraya::read_phraya(input_path)?;

    let initial_count = phraya_file.observations.len();

    // Start from preset defaults if specified, then apply explicit overrides on top.
    let mut filter_builder = match preset {
        Some("conservative") => FilterPreset::Conservative.builder(),
        Some("sensitive") => FilterPreset::Sensitive.builder(),
        Some(other) => {
            return Err(format!(
                "Unknown preset '{}'. Valid presets: conservative, sensitive",
                other
            )
            .into())
        }
        None => FilterBuilder::new(),
    };
    if let Some(min_cov) = min_coverage {
        filter_builder = filter_builder.min_coverage(min_cov);
    }
    if let Some(max_cov) = max_coverage {
        filter_builder = filter_builder.max_coverage(max_cov);
    }
    if let Some(min_mq) = min_mapq {
        filter_builder = filter_builder.min_mapq(min_mq);
    }
    if let Some(max_mq) = max_mapq {
        filter_builder = filter_builder.max_mapq(max_mq);
    }

    let filter = filter_builder.build();

    // Apply filter to observations
    let filtered_observations: Vec<_> = filter.filter(&phraya_file.observations).cloned().collect();

    let final_count = filtered_observations.len();
    eprintln!("Filtered {} → {} observations", initial_count, final_count);

    // Output in specified format
    match format {
        "vcf" => {
            let vcf_output = vcf::format_vcf(
                filtered_observations.into_iter(),
                &phraya_file.header.sample_id,
                phraya_file.header.reference_length,
            );
            println!("{}", vcf_output);
        }
        "tsv" => {
            output_tsv(&filtered_observations)?;
        }
        "phraya" => {
            if let Some(out_path) = output_path {
                let filtered_file = phraya::PhrayaFile::new(
                    phraya_file.header.reference_length,
                    phraya_file.header.sample_id.clone(),
                    phraya_file.header.timestamp.clone(),
                    filtered_observations,
                    phraya_file.coverage_track.clone(),
                );
                phraya::write_phraya(out_path, &filtered_file)?;
            } else {
                return Err("--output is required when format is 'phraya'".into());
            }
        }
        _ => return Err(format!("Unsupported format: {}", format).into()),
    }

    Ok(())
}

fn output_tsv(
    observations: &[phraya_core::types::VariantObservation],
) -> Result<(), Box<dyn std::error::Error>> {
    // Output TSV header
    println!("position\tref_base\tall_alleles\tmapq\tconfidence\tcigar\tedit_distance\tcoverage\tavg_base_quality\tprovenance");

    // Output each observation as a TSV line
    for obs in observations {
        let coverage = if !obs.local_coverage().is_empty() {
            obs.local_coverage()[0]
        } else {
            0
        };

        // Format all_alleles as key:value pairs
        let alleles_str = obs
            .all_alleles()
            .iter()
            .map(|(&base, &count)| format!("{}:{}", base as char, count))
            .collect::<Vec<_>>()
            .join(",");

        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            obs.position(),
            obs.ref_base() as char,
            alleles_str,
            obs.mapq(),
            obs.confidence(),
            obs.cigar(),
            obs.edit_distance(),
            coverage,
            obs.avg_base_quality(),
            obs.provenance(),
        );
    }

    Ok(())
}
