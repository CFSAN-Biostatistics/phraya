use clap::Parser;
use log::info;
use phraya_core::types::Sequence;
use phraya_index::{compute_kmer_uniqueness, sketch_sequence_default};
use phraya_io::{plan::PhrayaPlan, plan::UseCase, plan::write_plan, SequenceParser};
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "phraya")]
#[command(about = "Phraya: pairwise sequence aligner for bacterial genomics")]
enum Cli {
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
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::builder()
        .is_test(false)
        .try_init()
        .ok();

    let cli = Cli::parse();

    match cli {
        Cli::Plan {
            inputs,
            reference,
            output,
        } => {
            run_plan(&inputs, reference.as_deref(), &output)?;
        }
    }

    Ok(())
}

fn run_plan(
    input_paths: &[PathBuf],
    reference_path: Option<&std::path::Path>,
    output_path: &PathBuf,
) -> Result<(), Box<dyn std::error::Error>> {
    // Read all input sequences
    let mut all_sequences: Vec<(Sequence, String)> = Vec::new();
    let mut input_file_list = Vec::new();

    // First, read reference if provided
    let mut ref_sequence: Option<(Sequence, String)> = None;
    if let Some(ref_path) = reference_path {
        let ref_path_str = ref_path.to_string_lossy().to_string();
        input_file_list.push(ref_path_str.clone());

        let mut parser = SequenceParser::from_path(ref_path)?;
        while let Some(seq_result) = parser.next() {
            let seq = seq_result?;
            ref_sequence = Some((seq, ref_path_str.clone()));
            break; // Take only the first sequence as reference
        }

        if let Some((ref_seq, _)) = &ref_sequence {
            all_sequences.push((ref_seq.clone(), ref_path_str));
        }
    }

    // Read input sequences
    for input_path in input_paths {
        let input_path_str = input_path.to_string_lossy().to_string();
        input_file_list.push(input_path_str.clone());

        let mut parser = SequenceParser::from_path(input_path)?;
        while let Some(seq_result) = parser.next() {
            let seq = seq_result?;
            all_sequences.push((seq, input_path_str.clone()));
        }
    }

    if all_sequences.is_empty() {
        return Err("No sequences found in input files".into());
    }

    // Compute sketches for all sequences
    let sketches: Vec<_> = all_sequences
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

    eprintln!("Detected Case {}: {:?}", use_case_number(&use_case), format!("{:?}", use_case));

    // Compute k-mer uniqueness
    let kmer_uniqueness = compute_kmer_uniqueness(&sketches);

    // Generate task list based on use case
    let task_list = generate_task_list(&use_case, input_paths.len(), reference_path.is_some(), &sketches);

    // Create and write plan
    let plan = PhrayaPlan::new(
        use_case,
        input_file_list,
        chrono::Utc::now().to_rfc3339_opts(chrono::SecondsFormat::Secs, true),
        sketches,
        kmer_uniqueness,
        task_list,
    );

    write_plan(output_path, &plan)?;

    info!("Plan written to {:?}", output_path);
    Ok(())
}

/// Detect the use case from input characteristics
fn detect_use_case(
    has_reference: bool,
    num_input_files: usize,
    _sketches: &[phraya_index::MinimimizerSketch],
    all_sequences: &[(Sequence, String)],
) -> UseCase {
    // Count sequences from input files only (exclude reference if provided)
    let num_input_sequences = all_sequences.len() - if has_reference { 1 } else { 0 };

    if has_reference {
        // Case 2: Reads + reference
        UseCase::ReadsWithRef
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
        // Case 1: Reads only (should not happen in MVP)
        UseCase::ReadsOnly
    }
}

/// Generate task list based on use case and input characteristics
fn generate_task_list(
    use_case: &UseCase,
    _num_input_files: usize,
    _has_reference: bool,
    sketches: &[phraya_index::MinimimizerSketch],
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
            // Case 3: M contigs + N reads, select centroid as reference
            // For simplicity, we generate all pairwise tasks
            let mut tasks = Vec::new();
            // Generate all pairwise alignments (reads vs contigs and reads vs reads)
            for i in 0..sketches.len() {
                for j in 0..sketches.len() {
                    if i != j {
                        tasks.push((i as u32, j as u32));
                    }
                }
            }
            tasks
        }
        UseCase::ContigsOnly => {
            // Case 4: M contigs only
            // Generate all pairwise alignments
            let mut tasks = Vec::new();
            for i in 0..sketches.len() {
                for j in (i + 1)..sketches.len() {
                    tasks.push((i as u32, j as u32));
                }
            }
            tasks
        }
        UseCase::ReadsOnly => {
            // Case 1: N reads only (not in MVP)
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
