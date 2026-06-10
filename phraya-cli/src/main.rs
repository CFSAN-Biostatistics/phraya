use clap::{Parser, Subcommand};
use log::info;
use phraya_align::executor::{align_task_with_config, AlignConfig, Strategy};
use phraya_core::types::{
    compute_kmer_uniqueness, detect_hotspot_intervals, select_centroid, sketch_sequence_default,
    CoverageTrack, MinimizerSketch, Sequence,
};
use phraya_filter::{vcf, FilterBuilder, FilterPreset};
use phraya_io::{
    bam_cram::{BamCramParser, ParsedReads},
    phraya,
    plan::{self, write_plan, InsertSizeDistribution, PhrayaPlan, UseCase},
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

        /// Batch mode: divide reads into N chunks
        #[arg(long, value_name = "N", conflicts_with = "batch_by")]
        batch_to: Option<usize>,

        /// Batch mode: X reads per chunk
        #[arg(long, value_name = "N", conflicts_with = "batch_to")]
        batch_by: Option<usize>,

        /// Batch output pattern with {worker} placeholder (required if batch_to or batch_by specified)
        #[arg(long, value_name = "PATTERN")]
        batch_output_pattern: Option<String>,
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

        /// Query sequence ID (not used in batch mode)
        #[arg(value_name = "QUERY_ID")]
        query_id: Option<String>,

        /// Target sequence ID (not used in batch mode)
        #[arg(value_name = "TARGET_ID")]
        target_id: Option<String>,

        /// Output .phraya file (not used in batch --worker mode, taken from plan)
        #[arg(long, value_name = "FILE")]
        output: Option<PathBuf>,

        /// Alignment strategy: fast (±150bp window), balanced (±50bp, default), exact (±25bp)
        #[arg(long, value_name = "STRATEGY", default_value = "balanced")]
        strategy: String,

        /// Batch mode: worker ID (0-indexed)
        #[arg(long, value_name = "N")]
        worker: Option<usize>,

        /// Batch mode: process all missing chunks sequentially
        #[arg(long)]
        ensure: bool,

        /// Number of parallel threads for --ensure mode (default: auto-detect)
        #[arg(long, value_name = "N")]
        threads: Option<usize>,
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

        /// Minimum k-mer uniqueness threshold (0.0-1.0)
        #[arg(long, value_name = "F")]
        min_kmer_uniqueness: Option<f64>,

        /// Output format (vcf, tsv, phraya)
        #[arg(long, value_name = "FORMAT", default_value = "vcf")]
        format: String,

        /// Output file (required for phraya format)
        #[arg(long, value_name = "FILE")]
        output: Option<PathBuf>,

        /// Named filter preset (conservative or sensitive). Individual threshold flags override preset values.
        #[arg(long, value_name = "PRESET")]
        preset: Option<String>,

        /// Exclude discordant pairs (insert size beyond mean ± 3σ)
        #[arg(long)]
        exclude_discordant_pairs: bool,

        /// Sigma threshold for discordant pair detection (default: 3.0)
        #[arg(long, value_name = "F", default_value = "3.0")]
        discordant_sigma: f64,

        /// Require proper pairs (SAM flag 0x2)
        #[arg(long)]
        require_proper_pairs: bool,

        /// Minimum insert size (absolute value)
        #[arg(long, value_name = "N")]
        min_insert_size: Option<i32>,

        /// Maximum insert size (absolute value)
        #[arg(long, value_name = "N")]
        max_insert_size: Option<i32>,

        /// Require both mates mapped
        #[arg(long)]
        require_both_mates_mapped: bool,
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
            batch_to,
            batch_by,
            batch_output_pattern,
        } => {
            // Validate batch flags
            if (batch_to.is_some() || batch_by.is_some()) && batch_output_pattern.is_none() {
                return Err("--batch-output-pattern required when --batch-to or --batch-by specified".into());
            }
            run_plan(&inputs, reference.as_deref(), &output, batch_to, batch_by, batch_output_pattern.as_deref())?;
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
            worker,
            ensure,
            threads,
        } => {
            let strat = match strategy.as_str() {
                "fast" => Strategy::Fast,
                "balanced" => Strategy::Balanced,
                "exact" => Strategy::Exact,
                other => return Err(format!("unknown strategy: {other}; expected fast, balanced, or exact").into()),
            };

            if ensure {
                run_align_ensure(&plan_file, AlignConfig::new(strat), threads)?;
            } else if let Some(worker_id) = worker {
                if output.is_some() {
                    return Err("--output not allowed with --worker (output path from plan)".into());
                }
                run_align_worker(&plan_file, worker_id, AlignConfig::new(strat))?;
            } else {
                // Traditional mode
                let query = query_id.ok_or("QUERY_ID required in traditional mode")?;
                let target = target_id.ok_or("TARGET_ID required in traditional mode")?;
                let out = output.ok_or("--output required in traditional mode")?;
                run_align(&plan_file, &query, &target, &out, AlignConfig::new(strat))?;
            }
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
            min_kmer_uniqueness,
            format,
            output,
            preset,
            exclude_discordant_pairs,
            discordant_sigma,
            require_proper_pairs,
            min_insert_size,
            max_insert_size,
            require_both_mates_mapped,
        } => {
            run_filter(
                &input,
                min_coverage,
                max_coverage,
                min_mapq,
                max_mapq,
                min_kmer_uniqueness,
                &format,
                output.as_deref(),
                preset.as_deref(),
                exclude_discordant_pairs,
                discordant_sigma,
                require_proper_pairs,
                min_insert_size,
                max_insert_size,
                require_both_mates_mapped,
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

/// Align specific worker chunk in batch mode
fn run_align_worker(
    plan_path: &std::path::Path,
    worker_id: usize,
    config: AlignConfig,
) -> Result<(), Box<dyn std::error::Error>> {
    let plan = plan::read_plan(plan_path)?;
    run_align_worker_with_plan(worker_id, &plan, config)
}

/// Align worker chunk with borrowed plan (for parallel execution)
fn run_align_worker_with_plan(
    worker_id: usize,
    plan: &PhrayaPlan,
    config: AlignConfig,
) -> Result<(), Box<dyn std::error::Error>> {

    // Validate worker ID
    let num_chunks = plan.batch_num_chunks.ok_or("Plan has no batch configuration")?;
    if worker_id >= num_chunks {
        return Err(format!("Worker {} out of range for {} chunks", worker_id, num_chunks).into());
    }

    // Get output path from plan
    let output_path = plan.batch_output_paths.get(worker_id)
        .ok_or(format!("No output path for worker {}", worker_id))?;

    // Calculate chunk boundaries (exclude reference from read pool if present)
    // Batch mode only processes query reads, not the reference
    let ref_count = if !plan.input_files.is_empty() && plan.reads_per_file.len() > 1 {
        plan.reads_per_file[0] // Assume first file is reference
    } else {
        0
    };
    let query_read_count = plan.total_read_count.saturating_sub(ref_count);

    let chunk_size = if let Some(reads_per) = plan.batch_reads_per_chunk {
        reads_per
    } else {
        (query_read_count + num_chunks - 1) / num_chunks
    };
    let start_idx = worker_id * chunk_size;
    let end_idx = std::cmp::min(start_idx + chunk_size, query_read_count);

    eprintln!("Worker {} processing reads [{}, {})", worker_id, start_idx, end_idx);

    // Read target sequence (reference or centroid)
    // For now, assume target is first sequence in first file
    let mut target_seq: Option<Sequence> = None;
    let mut parser = SequenceParser::from_path(&plan.input_files[0])?;
    if let Some(result) = parser.next() {
        target_seq = Some(result?);
    }
    let target = target_seq.ok_or("No target sequence found")?;

    // Extract and align reads in this chunk
    let mut all_variants = Vec::new();
    let mut all_query_positions = HashMap::new();
    let mut coverage_track = vec![0u32; target.len()];

    for global_idx in start_idx..end_idx {
        // Map global_idx to (file_idx, local_idx)
        let (file_idx, local_idx) = map_global_to_local(&plan, global_idx)?;

        // Extract sequence at byte offset
        let seq = extract_sequence_at_offset(&plan, file_idx, local_idx)?;
        let query_id = seq.id().to_string();

        // Align
        if let Some(result) = align_task_with_config(&seq, &target, &plan, &config) {
            all_variants.extend(result.variants);
            all_query_positions.insert(query_id, result.query_positions);
            for (i, &cov) in result.coverage_track.iter().enumerate() {
                coverage_track[i] += cov;
            }
        }
    }

    // Write output
    let coverage = CoverageTrack::new(coverage_track.iter().map(|&v| v as usize).collect());
    let phraya_file = phraya::PhrayaFile::new(
        target.len() as u32,
        format!("worker_{}", worker_id),
        chrono::Utc::now().to_rfc3339(),
        all_variants,
        coverage,
    );
    phraya::write_phraya(std::path::Path::new(output_path), &phraya_file)?;

    // Write queries sidecar
    let mut index = queries::QueryIndex::new();
    for (qid, positions) in all_query_positions {
        index.insert(qid, positions);
    }
    let queries_path = format!("{}.queries", output_path);
    queries::write_queries(std::path::Path::new(&queries_path), &index)?;

    eprintln!("Worker {} complete: {} observations written to {}", worker_id, phraya_file.observations.len(), output_path);
    Ok(())
}

/// Map global read index to (file_idx, local_idx), skipping reference (file 0)
fn map_global_to_local(plan: &PhrayaPlan, global_idx: usize) -> Result<(usize, usize), Box<dyn std::error::Error>> {
    let mut cumulative = 0;
    // Start from file 1 (skip reference at file 0)
    for (file_idx, &count) in plan.reads_per_file.iter().enumerate().skip(1) {
        if global_idx < cumulative + count {
            return Ok((file_idx, global_idx - cumulative));
        }
        cumulative += count;
    }
    Err(format!("Global index {} out of bounds", global_idx).into())
}

/// Process all missing chunks in batch mode
fn run_align_ensure(
    plan_path: &std::path::Path,
    config: AlignConfig,
    threads: Option<usize>,
) -> Result<(), Box<dyn std::error::Error>> {
    use rayon::prelude::*;
    use std::sync::Arc;

    // Configure thread pool if requested
    if let Some(n) = threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .map_err(|e| format!("Failed to configure thread pool: {}", e))?;
    }

    // Load plan once, share via Arc
    let plan = Arc::new(plan::read_plan(plan_path)?);

    let _num_chunks = plan.batch_num_chunks.ok_or("Plan has no batch configuration")?;
    let mut missing_chunks = Vec::new();

    // Find missing outputs
    for (worker_id, output_path) in plan.batch_output_paths.iter().enumerate() {
        if !std::path::Path::new(output_path).exists() {
            missing_chunks.push(worker_id);
        }
    }

    if missing_chunks.is_empty() {
        eprintln!("All chunks already complete");
        return Ok(());
    }

    eprintln!("Processing {} missing chunks in parallel (threads: {})",
        missing_chunks.len(), rayon::current_num_threads());

    // Parallel execution with error collection
    let results: Vec<_> = missing_chunks
        .par_iter()
        .map(|&worker_id| {
            let result = run_align_worker_with_plan(worker_id, &*plan, config);
            (worker_id, result.map_err(|e| e.to_string()))
        })
        .collect();

    // Report failures
    let failures: Vec<_> = results.iter()
        .filter_map(|(id, r)| r.as_ref().err().map(|e| (*id, e.clone())))
        .collect();

    if !failures.is_empty() {
        for (id, err) in &failures {
            eprintln!("Worker {} failed: {}", id, err);
        }
        return Err(format!("{} of {} workers failed",
            failures.len(), missing_chunks.len()).into());
    }

    eprintln!("Ensure mode complete");
    Ok(())
}

/// Extract sequence at specific byte offset
fn extract_sequence_at_offset(
    plan: &PhrayaPlan,
    file_idx: usize,
    local_idx: usize,
) -> Result<Sequence, Box<dyn std::error::Error>> {
    use std::fs::File;
    use std::io::{BufRead, BufReader, Seek, SeekFrom};

    let file_path = &plan.input_files[file_idx];
    let offset = plan.read_byte_offsets[file_idx][local_idx];

    let file = File::open(file_path)?;
    let mut reader = BufReader::new(file);
    reader.seek(SeekFrom::Start(offset))?;

    // Parse one record manually (FASTA or FASTQ)
    let mut line = String::new();
    reader.read_line(&mut line)?;

    if line.starts_with('>') {
        // FASTA format
        let id = line[1..].trim().to_string();
        let mut sequence = String::new();
        line.clear();
        while reader.read_line(&mut line)? > 0 {
            if line.starts_with('>') {
                break;
            }
            sequence.push_str(line.trim());
            line.clear();
        }
        Ok(Sequence::new(sequence.into_bytes(), None, id, None))
    } else if line.starts_with('@') {
        // FASTQ format
        let id = line[1..].trim().to_string();
        line.clear();
        reader.read_line(&mut line)?;
        let sequence = line.trim().to_string();
        // Skip + and quality lines
        line.clear();
        reader.read_line(&mut line)?; // +
        line.clear();
        let mut quality_line = String::new();
        reader.read_line(&mut quality_line)?;
        let quality = quality_line.trim().as_bytes().to_vec();
        Ok(Sequence::new(sequence.into_bytes(), Some(quality), id, None))
    } else {
        Err(format!("Unknown format at offset {} in {}", offset, file_path).into())
    }
}

/// Infer insert size distribution from BAM files in inputs
fn infer_insert_size_from_inputs(
    input_paths: &[PathBuf],
) -> Result<Option<InsertSizeDistribution>, Box<dyn std::error::Error>> {
    // Find first BAM file
    let bam_path = input_paths.iter()
        .find(|p| p.to_string_lossy().ends_with(".bam"));

    if let Some(bam_path) = bam_path {
        infer_insert_size_from_bam(bam_path)
    } else {
        Ok(None)
    }
}

/// Infer insert size distribution from BAM proper pairs
fn infer_insert_size_from_bam(
    bam_path: &std::path::Path
) -> Result<Option<InsertSizeDistribution>, Box<dyn std::error::Error>> {
    use noodles_bam::io::Reader;
    use std::fs::File;
    use std::io::BufReader;

    let file = File::open(bam_path)?;
    let reader = BufReader::new(file);
    let mut bam_reader = Reader::new(reader);
    let _header = bam_reader.read_header()?;

    let mut tlens = Vec::new();
    const MAX_SAMPLES: usize = 10000; // Sample first 10k proper pairs

    for result in bam_reader.records() {
        if tlens.len() >= MAX_SAMPLES {
            break;
        }

        let record = result?;
        let flags = record.flags();

        // Only sample proper pairs for insert size estimation
        if flags.is_properly_segmented() {
            let tlen = record.template_length();
            if tlen != 0 {
                tlens.push(tlen);
            }
        }
    }

    Ok(InsertSizeDistribution::from_bam_proper_pairs(&tlens))
}

/// Helper to parse sequences from any format (FASTA/FASTQ/BAM/CRAM)
fn parse_sequences_from_file(
    path: &std::path::Path,
) -> Result<(Vec<Sequence>, HashMap<String, phraya_core::types::MateInfo>), Box<dyn std::error::Error>> {
    let path_str = path.to_string_lossy();

    if path_str.ends_with(".bam") {
        let parsed = BamCramParser::from_bam_path(path)?;
        Ok((parsed.sequences, parsed.mate_info))
    } else if path_str.ends_with(".cram") {
        let parsed = BamCramParser::from_cram_path(path)?;
        Ok((parsed.sequences, parsed.mate_info))
    } else {
        // FASTA/FASTQ - no mate info
        let mut sequences = Vec::new();
        let mut parser = SequenceParser::from_path(path)?;
        while let Some(seq_result) = parser.next() {
            sequences.push(seq_result?);
        }
        Ok((sequences, HashMap::new()))
    }
}

fn run_plan(
    input_paths: &[PathBuf],
    reference_path: Option<&std::path::Path>,
    output_path: &PathBuf,
    batch_to: Option<usize>,
    batch_by: Option<usize>,
    batch_output_pattern: Option<&str>,
) -> Result<(), Box<dyn std::error::Error>> {
    // Read all input sequences and track which file they came from
    let mut all_sequences: Vec<(Sequence, String)> = Vec::new();
    let mut sequence_to_file_index: Vec<usize> = Vec::new(); // Track which input file each sequence came from
    let mut input_file_list = Vec::new();
    let mut ref_seq_index: Option<usize> = None;
    let mut all_mate_info: HashMap<String, phraya_core::types::MateInfo> = HashMap::new();

    // First, read reference if provided
    if let Some(ref_path) = reference_path {
        let ref_path_str = ref_path.to_string_lossy().to_string();
        input_file_list.push(ref_path_str.clone());

        let (sequences, mate_info) = parse_sequences_from_file(ref_path)?;
        if let Some(seq) = sequences.into_iter().next() {
            ref_seq_index = Some(all_sequences.len());
            all_sequences.push((seq, ref_path_str.clone()));
            sequence_to_file_index.push(0); // File index 0 for reference
        }
        all_mate_info.extend(mate_info);
    }

    // Read input sequences
    for (file_idx, input_path) in input_paths.iter().enumerate() {
        let input_path_str = input_path.to_string_lossy().to_string();
        input_file_list.push(input_path_str.clone());

        let (sequences, mate_info) = parse_sequences_from_file(input_path)?;
        for seq in sequences {
            all_sequences.push((seq, input_path_str.clone()));
            // File indices: 1+ for input files (offset by 1 if there's a reference)
            let file_index = if ref_seq_index.is_some() { file_idx + 1 } else { file_idx };
            sequence_to_file_index.push(file_index);
        }
        all_mate_info.extend(mate_info);
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

    // Build byte offset index and read counts for batch mode
    let (read_byte_offsets, reads_per_file, total_read_count) =
        build_byte_offset_index(reference_path, input_paths)?;

    // Infer insert size distribution from BAM files if present
    let insert_size_distribution = infer_insert_size_from_inputs(input_paths)?;
    if let Some(ref dist) = insert_size_distribution {
        eprintln!("Inferred insert size: mean={}, std_dev={}, n={}",
            dist.mean, dist.std_dev, dist.sample_size);
    }

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
    plan.read_byte_offsets = read_byte_offsets;
    plan.reads_per_file = reads_per_file;
    plan.total_read_count = total_read_count;
    plan.kmer_params = phraya_io::plan::KmerParams { k: 21, w: 11 };
    plan.insert_size_distribution = insert_size_distribution;
    plan.mate_info = all_mate_info;

    // Handle batch mode configuration
    if let Some(batch_pattern) = batch_output_pattern {
        configure_batch_mode(&mut plan, batch_to, batch_by, batch_pattern)?;
    }

    write_plan(output_path, &plan)?;

    info!("Plan written to {:?}", output_path);
    Ok(())
}

/// Build byte offset index for all input files
fn build_byte_offset_index(
    reference_path: Option<&std::path::Path>,
    input_paths: &[PathBuf],
) -> Result<(Vec<Vec<u64>>, Vec<usize>, usize), Box<dyn std::error::Error>> {
    let mut all_offsets = Vec::new();
    let mut all_counts = Vec::new();
    let mut total = 0;

    // Process reference first if present
    if let Some(ref_path) = reference_path {
        let (offsets, count) = index_file_offsets(ref_path)?;
        all_offsets.push(offsets);
        all_counts.push(count);
        total += count;
    }

    // Process input files
    for input_path in input_paths {
        let (offsets, count) = index_file_offsets(input_path)?;
        all_offsets.push(offsets);
        all_counts.push(count);
        total += count;
    }

    Ok((all_offsets, all_counts, total))
}

/// Configure batch mode for the plan
fn configure_batch_mode(
    plan: &mut PhrayaPlan,
    batch_to: Option<usize>,
    batch_by: Option<usize>,
    batch_pattern: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    use std::collections::HashSet;

    // Determine num_chunks based on flags
    let num_chunks = if let Some(n) = batch_to {
        plan.batch_num_chunks = Some(n);
        if let Some(reads_per) = batch_by {
            plan.batch_reads_per_chunk = Some(reads_per);
            // Validate coverage
            let expected_total = n * reads_per;
            if expected_total < plan.total_read_count {
                eprintln!(
                    "Warning: under-provisioned batch config. {} chunks × {} reads/chunk = {} < {} total reads",
                    n, reads_per, expected_total, plan.total_read_count
                );
            } else if expected_total > plan.total_read_count {
                eprintln!(
                    "Warning: over-provisioned batch config. {} chunks × {} reads/chunk = {} > {} total reads",
                    n, reads_per, expected_total, plan.total_read_count
                );
            }
        }
        n
    } else if let Some(reads_per) = batch_by {
        plan.batch_reads_per_chunk = Some(reads_per);
        let n = (plan.total_read_count + reads_per - 1) / reads_per; // Ceiling division
        plan.batch_num_chunks = Some(n);
        n
    } else {
        return Err("Either --batch-to or --batch-by must be specified".into());
    };

    // Expand output pattern
    let mut expanded_paths = Vec::new();
    for worker_id in 0..num_chunks {
        let path = batch_pattern.replace("{worker}", &worker_id.to_string());
        expanded_paths.push(path);
    }

    // Check for collisions
    let unique_paths: HashSet<_> = expanded_paths.iter().collect();
    if unique_paths.len() != expanded_paths.len() {
        return Err("Batch output paths contain collisions (non-unique paths)".into());
    }

    plan.batch_output_paths = expanded_paths;
    Ok(())
}

/// Index byte offsets for a single file
fn index_file_offsets(path: &std::path::Path) -> Result<(Vec<u64>, usize), Box<dyn std::error::Error>> {
    use std::fs::File;
    use std::io::{BufRead, BufReader, Seek, SeekFrom};

    let file = File::open(path)?;
    let mut reader = BufReader::new(file);
    let mut offsets = Vec::new();
    let mut position: u64 = 0;
    let mut line = String::new();

    // Detect format
    reader.read_line(&mut line)?;
    let is_fasta = line.starts_with('>');
    let is_fastq = line.starts_with('@');

    reader.seek(SeekFrom::Start(0))?;
    position = 0;
    line.clear();

    if is_fasta {
        // FASTA: record starts with '>'
        loop {
            let start_pos = position;
            let bytes_read = reader.read_line(&mut line)?;
            if bytes_read == 0 {
                break;
            }
            if line.starts_with('>') {
                offsets.push(start_pos);
            }
            position += bytes_read as u64;
            line.clear();
        }
    } else if is_fastq {
        // FASTQ: 4-line records, starts with '@'
        loop {
            let start_pos = position;
            // Read header
            let bytes_read = reader.read_line(&mut line)?;
            if bytes_read == 0 {
                break;
            }
            if line.starts_with('@') {
                offsets.push(start_pos);
                position += bytes_read as u64;
                line.clear();
                // Skip sequence, +, quality lines
                for _ in 0..3 {
                    let n = reader.read_line(&mut line)?;
                    position += n as u64;
                    line.clear();
                }
            } else {
                position += bytes_read as u64;
                line.clear();
            }
        }
    } else {
        // BAM/CRAM: not text-based, can't build simple byte offset index
        // For now, return empty offsets (batch mode won't work for BAM/CRAM inputs)
        eprintln!("Warning: byte offset indexing not supported for BAM/CRAM files");
        return Ok((Vec::new(), 0));
    }

    let count = offsets.len();
    Ok((offsets, count))
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

    // Auto-detect plan file
    if input_paths.len() == 1 && input_paths[0].extension().and_then(|s| s.to_str()) == Some("phrayaplan") {
        eprintln!("Detected plan file, merging batch outputs...");
        let plan = plan::read_plan(&input_paths[0])?;

        if plan.batch_output_paths.is_empty() {
            return Err("Plan has no batch outputs to merge".into());
        }

        // Verify all outputs exist
        for (idx, output_path_str) in plan.batch_output_paths.iter().enumerate() {
            let path = std::path::Path::new(output_path_str);
            if !path.exists() {
                return Err(format!(
                    "Missing batch output {} (worker {}). Run with --ensure to complete.",
                    output_path_str, idx
                ).into());
            }
        }

        eprintln!("Merging {} batch outputs...", plan.batch_output_paths.len());
        let paths: Vec<&std::path::Path> = plan.batch_output_paths.iter()
            .map(|s| std::path::Path::new(s.as_str()))
            .collect();
        let merged = phraya::merge_phraya_files(&paths)?;
        phraya::write_phraya(output_path, &merged)?;
        eprintln!("Merged file written to {:?}", output_path);
        return Ok(());
    }

    // Traditional mode
    eprintln!("Merging {} samples...", input_paths.len());
    let paths: Vec<&std::path::Path> = input_paths.iter().map(|p| p.as_path()).collect();
    let merged = phraya::merge_phraya_files(&paths)?;
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
    min_kmer_uniqueness: Option<f64>,
    format: &str,
    output_path: Option<&std::path::Path>,
    preset: Option<&str>,
    exclude_discordant_pairs: bool,
    discordant_sigma: f64,
    require_proper_pairs: bool,
    min_insert_size: Option<i32>,
    max_insert_size: Option<i32>,
    require_both_mates_mapped: bool,
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
    if let Some(min_km) = min_kmer_uniqueness {
        filter_builder = filter_builder.min_kmer_uniqueness(min_km);
    }

    // Paired-end filters
    if exclude_discordant_pairs {
        filter_builder = filter_builder
            .exclude_discordant_pairs(true)
            .discordant_sigma_threshold(discordant_sigma);
    }
    if require_proper_pairs {
        filter_builder = filter_builder.require_proper_pairs(true);
    }
    if let Some(min) = min_insert_size {
        filter_builder = filter_builder.min_insert_size(min);
    }
    if let Some(max) = max_insert_size {
        filter_builder = filter_builder.max_insert_size(max);
    }
    if require_both_mates_mapped {
        filter_builder = filter_builder.require_both_mates_mapped(true);
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
