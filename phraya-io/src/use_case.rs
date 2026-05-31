use crate::SequenceParser;
/// Use case detection module for classifying alignment workflows.
///
/// Detects workflow type from input files (reads vs contigs) and reference presence.
/// Case 2: N reads + reference
/// Case 3: M contigs + N reads, no reference (with centroid selection)
/// Case 4: M contigs ± reference (or contig MSA)
/// Case 1: reads only, no reference (deferred to Phase 2)
use phraya_index::{select_centroid, sketch_sequence_default};
use std::path::Path;
use thiserror::Error;

/// Detected input classification
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputType {
    /// Short sequence (<5kb, likely from sequencing reads)
    Read,
    /// Long sequence (≥5kb, likely an assembly contig)
    Contig,
}

/// Detected use case from inputs and reference
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum UseCase {
    /// Case 2: N reads + reference → N alignment tasks
    Case2 {
        /// Number of alignment tasks (typically equals number of input sequences)
        n_tasks: usize,
    },
    /// Case 3: M contigs + N reads, no reference → centroid + M+N tasks
    Case3 {
        /// Index of the centroid contig in inputs
        centroid_idx: usize,
        /// Total number of tasks: 1 (centroid ref) + (M-1) (contigs) + N (reads)
        n_tasks: usize,
    },
    /// Case 4: M contigs ± reference
    Case4 {
        /// Number of alignment tasks
        n_tasks: usize,
        /// True if MSA mode (contigs only, no reference)
        is_msa: bool,
    },
}

/// Errors during use case detection
#[derive(Debug, Clone, Error)]
pub enum UseCaseError {
    /// Case 1 (reads only, no reference) deferred to Phase 2
    #[error("Case 1 (reads only, no reference) deferred to Phase 2")]
    Case1Deferred,

    /// Invalid input specification
    #[error("invalid input: {0}")]
    InvalidInput(String),

    /// I/O errors during file parsing
    #[error("io error: {0}")]
    IoError(String),

    /// No sequences found in input
    #[error("no sequences found in input")]
    NoSequences,
}

/// Detect input type (read or contig) from the first sequence in a file.
///
/// Returns InputType::Read if first sequence is <5kb, InputType::Contig if ≥5kb.
///
/// # Arguments
/// * `path` - Path to FASTA/FASTQ file
///
/// # Returns
/// InputType or UseCaseError on I/O failure or empty file
pub fn classify_input(path: &Path) -> Result<InputType, UseCaseError> {
    let mut parser =
        SequenceParser::from_path(path).map_err(|e| UseCaseError::IoError(e.to_string()))?;

    match parser.next() {
        Some(Ok(seq)) => {
            if seq.bases().len() >= 5000 {
                Ok(InputType::Contig)
            } else {
                Ok(InputType::Read)
            }
        }
        Some(Err(e)) => Err(UseCaseError::IoError(e.to_string())),
        None => Err(UseCaseError::NoSequences),
    }
}

/// Detect use case from inputs and optional reference.
///
/// Classifies each input file as READ (<5kb) or CONTIG (≥5kb) based on the first sequence,
/// then determines which workflow case applies:
///
/// - Case 2: has reference, inputs are reads → N tasks
/// - Case 3: no reference, mixed contigs + reads → select centroid, M+N tasks
/// - Case 4: contigs only → M tasks (with ref) or M×(M-1)/2 (without ref, MSA)
/// - Case 1: reads only, no reference → ERROR (deferred)
///
/// # Arguments
/// * `inputs` - Slice of paths to input FASTA/FASTQ files (can be &[Path] or &[PathBuf])
/// * `reference` - Optional path to reference sequence
///
/// # Returns
/// (UseCase, Vec<InputType>) on success, UseCaseError on failure
pub fn detect_use_case<P: AsRef<Path>>(
    inputs: &[P],
    reference: Option<&Path>,
) -> Result<(UseCase, Vec<InputType>), UseCaseError> {
    // Validate inputs not empty
    if inputs.is_empty() {
        return Err(UseCaseError::InvalidInput(
            "no input files provided".to_string(),
        ));
    }

    // Classify all inputs
    let input_types: Result<Vec<_>, _> = inputs
        .iter()
        .map(|input| classify_input(input.as_ref()))
        .collect();
    let input_types = input_types?;

    // Count reads and contigs
    let num_reads = input_types
        .iter()
        .filter(|&&t| t == InputType::Read)
        .count();
    let num_contigs = input_types.len() - num_reads;

    // Determine use case
    let use_case = match (reference, num_contigs, num_reads) {
        // Case 2: has reference, inputs can be reads (reads with reference) or mixed
        (Some(_), _, n_reads) if n_reads > 0 => {
            // If reference present and there are reads, treat as Case 2
            // Count only reads as tasks
            UseCase::Case2 { n_tasks: num_reads }
        }
        // Case 3: no reference, mixed contigs + reads
        (None, n_ctg, n_rd) if n_ctg > 0 && n_rd > 0 => {
            let centroid_idx = select_centroid_index(inputs, &input_types)?;
            // Total tasks: M contigs + N reads
            UseCase::Case3 {
                centroid_idx,
                n_tasks: num_contigs + num_reads,
            }
        }
        // Case 4: contigs only (with or without reference)
        (_, n_ctg, 0) if n_ctg > 0 => {
            if reference.is_some() {
                // With reference: M tasks
                UseCase::Case4 {
                    n_tasks: num_contigs,
                    is_msa: false,
                }
            } else {
                // Without reference (MSA mode): M×(M-1)/2 tasks
                let n_tasks = num_contigs * (num_contigs - 1) / 2;
                UseCase::Case4 {
                    n_tasks,
                    is_msa: true,
                }
            }
        }
        // Case 1: reads only, no reference (deferred)
        (None, 0, _) if num_reads > 0 => {
            return Err(UseCaseError::Case1Deferred);
        }
        // All other combinations are invalid
        _ => {
            return Err(UseCaseError::InvalidInput(
                "invalid combination of inputs and reference".to_string(),
            ));
        }
    };

    Ok((use_case, input_types))
}

/// Select the centroid contig index for Case 3.
///
/// For Case 3, we need to identify which contig is the centroid (will serve as reference coordinate space).
/// This requires sketching contigs and computing Jaccard similarity to find the median.
fn select_centroid_index<P: AsRef<Path>>(
    inputs: &[P],
    input_types: &[InputType],
) -> Result<usize, UseCaseError> {
    // Sketch all contigs
    let mut contig_sketches = Vec::new();
    let mut contig_indices = Vec::new();

    for (idx, input) in inputs.iter().enumerate() {
        if input_types[idx] == InputType::Contig {
            let mut parser = SequenceParser::from_path(input.as_ref())
                .map_err(|e| UseCaseError::IoError(e.to_string()))?;

            if let Some(Ok(seq)) = parser.next() {
                let sketch = sketch_sequence_default(&seq);
                contig_sketches.push(sketch);
                contig_indices.push(idx);
            } else {
                return Err(UseCaseError::NoSequences);
            }
        }
    }

    // Select centroid from contig sketches
    let centroid_sketch_idx = select_centroid(&contig_sketches)
        .ok_or_else(|| UseCaseError::InvalidInput("no contigs to select from".to_string()))?;

    // Map back to original input index
    Ok(contig_indices[centroid_sketch_idx])
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use std::path::PathBuf;
    use tempfile::NamedTempFile;

    // ===== Test utilities =====

    fn write_fasta(content: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        write!(f, "{}", content).unwrap();
        f.flush().unwrap();
        f
    }

    fn write_fastq(content: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        write!(f, "{}", content).unwrap();
        f.flush().unwrap();
        f
    }

    // ===== classify_input tests =====

    #[test]
    fn test_issue_67_classify_short_sequence_as_read() {
        // Sequence < 5kb should be classified as Read
        let fasta = ">read1\nACGT\n";
        let f = write_fasta(fasta);

        let result = classify_input(f.path()).unwrap();
        assert_eq!(result, InputType::Read);
    }

    #[test]
    fn test_issue_67_classify_long_sequence_as_contig() {
        // Sequence >= 5kb should be classified as Contig
        let long_seq = "A".repeat(5000);
        let fasta = format!(">contig1\n{}\n", long_seq);
        let f = write_fasta(&fasta);

        let result = classify_input(f.path()).unwrap();
        assert_eq!(result, InputType::Contig);
    }

    #[test]
    fn test_issue_67_classify_exactly_5kb_as_contig() {
        // Boundary: exactly 5kb should be Contig
        let seq = "A".repeat(5000);
        let fasta = format!(">boundary\n{}\n", seq);
        let f = write_fasta(&fasta);

        let result = classify_input(f.path()).unwrap();
        assert_eq!(result, InputType::Contig);
    }

    #[test]
    fn test_issue_67_classify_4999bp_as_read() {
        // Just below 5kb should be Read
        let seq = "A".repeat(4999);
        let fasta = format!(">just_under\n{}\n", seq);
        let f = write_fasta(&fasta);

        let result = classify_input(f.path()).unwrap();
        assert_eq!(result, InputType::Read);
    }

    #[test]
    fn test_issue_67_classify_reads_from_fastq() {
        // FASTQ format short reads should be classified as Read
        let fastq = "@read1\nACGT\n+\nIIII\n";
        let f = write_fastq(fastq);

        let result = classify_input(f.path()).unwrap();
        assert_eq!(result, InputType::Read);
    }

    #[test]
    fn test_issue_67_classify_wrapped_long_sequence() {
        // Long sequence wrapped across lines should be Contig
        let seq_part = "A".repeat(1250); // 4*1250 = 5000
        let fasta = format!(
            ">contig1\n{}\n{}\n{}\n{}\n",
            seq_part, seq_part, seq_part, seq_part
        );
        let f = write_fasta(&fasta);

        let result = classify_input(f.path()).unwrap();
        assert_eq!(result, InputType::Contig);
    }

    #[test]
    fn test_issue_67_classify_empty_file_error() {
        // Empty file should return NoSequences error
        let f = NamedTempFile::new().unwrap();

        let result = classify_input(f.path());
        assert!(matches!(result, Err(UseCaseError::NoSequences)));
    }

    // ===== detect_use_case Case 2 tests =====

    #[test]
    fn test_issue_67_case2_single_read_with_reference() {
        // Case 2: 1 read + reference = 1 task
        let read = ">read1\nACGT\n";
        let reference = ">ref\nACGTACGT\n";

        let read_f = write_fasta(read);
        let ref_f = write_fasta(reference);

        let (use_case, types) =
            detect_use_case(&[read_f.path().to_path_buf()], Some(ref_f.path())).unwrap();

        assert_eq!(use_case, UseCase::Case2 { n_tasks: 1 });
        assert_eq!(types, vec![InputType::Read]);
    }

    #[test]
    fn test_issue_67_case2_multiple_reads_with_reference() {
        // Case 2: N reads + reference = N tasks
        let read1 = ">read1\nACGT\n";
        let read2 = ">read2\nTGCA\n";
        let reference = ">ref\nACGTACGTACGT\n";

        let read1_f = write_fasta(read1);
        let read2_f = write_fasta(read2);
        let ref_f = write_fasta(reference);

        let (use_case, types) = detect_use_case(
            &[read1_f.path().to_path_buf(), read2_f.path().to_path_buf()],
            Some(ref_f.path()),
        )
        .unwrap();

        assert_eq!(use_case, UseCase::Case2 { n_tasks: 2 });
        assert_eq!(types, vec![InputType::Read, InputType::Read]);
    }

    #[test]
    fn test_issue_67_case2_reads_from_fasta() {
        // Case 2 with FASTA short reads
        let fasta_reads = ">read_1\nACGT\n";
        let reference = ">ref\nACGTACGT\n";

        let read_f = write_fasta(fasta_reads);
        let ref_f = write_fasta(reference);

        let (use_case, _) =
            detect_use_case(&[read_f.path().to_path_buf()], Some(ref_f.path())).unwrap();

        assert_eq!(use_case, UseCase::Case2 { n_tasks: 1 });
    }

    // ===== detect_use_case Case 3 tests =====

    #[test]
    fn test_issue_67_case3_contigs_and_reads_no_reference() {
        // Case 3: M contigs + N reads, no reference
        // 1 contig + 1 read = 2 tasks, centroid should be the single contig
        let contig = format!(">contig1\n{}\n", "A".repeat(5000));
        let read = ">read1\nACGT\n";

        let contig_f = write_fasta(&contig);
        let read_f = write_fasta(read);

        let (use_case, types) = detect_use_case(
            &[contig_f.path().to_path_buf(), read_f.path().to_path_buf()],
            None,
        )
        .unwrap();

        assert_eq!(types, vec![InputType::Contig, InputType::Read]);
        if let UseCase::Case3 {
            centroid_idx,
            n_tasks,
        } = use_case
        {
            assert_eq!(centroid_idx, 0); // First input is the contig
            assert_eq!(n_tasks, 2); // 1 contig + 1 read
        } else {
            panic!("expected Case3");
        }
    }

    #[test]
    fn test_issue_67_case3_multiple_contigs_and_reads_selects_centroid() {
        // Case 3: Multiple contigs + reads, should select centroid
        let contig1 = format!(">contig1\n{}\n", "A".repeat(5000));
        let contig2 = format!(">contig2\n{}\n", "A".repeat(5000)); // Identical to contig1
        let read = ">read1\nACGT\n";

        let contig1_f = write_fasta(&contig1);
        let contig2_f = write_fasta(&contig2);
        let read_f = write_fasta(read);

        let (use_case, types) = detect_use_case(
            &[
                contig1_f.path().to_path_buf(),
                contig2_f.path().to_path_buf(),
                read_f.path().to_path_buf(),
            ],
            None,
        )
        .unwrap();

        assert_eq!(
            types,
            vec![InputType::Contig, InputType::Contig, InputType::Read]
        );
        if let UseCase::Case3 {
            centroid_idx,
            n_tasks,
        } = use_case
        {
            // Centroid should be one of the first two indices (both are contigs)
            assert!(centroid_idx == 0 || centroid_idx == 1);
            assert_eq!(n_tasks, 3); // 2 contigs + 1 read
        } else {
            panic!("expected Case3");
        }
    }

    #[test]
    fn test_issue_67_case3_diverse_contigs() {
        // Case 3 with diverse contigs (different sequences)
        let contig1 = format!(">contig1\n{}\n", "A".repeat(5000));
        let contig2 = format!(">contig2\n{}\n", "T".repeat(5000));
        let read = ">read1\nACGT\n";

        let contig1_f = write_fasta(&contig1);
        let contig2_f = write_fasta(&contig2);
        let read_f = write_fasta(read);

        let (use_case, _) = detect_use_case(
            &[
                contig1_f.path().to_path_buf(),
                contig2_f.path().to_path_buf(),
                read_f.path().to_path_buf(),
            ],
            None,
        )
        .unwrap();

        if let UseCase::Case3 {
            centroid_idx: _,
            n_tasks,
        } = use_case
        {
            assert_eq!(n_tasks, 3);
        } else {
            panic!("expected Case3");
        }
    }

    // ===== detect_use_case Case 4 tests =====

    #[test]
    fn test_issue_67_case4_single_contig_with_reference() {
        // Case 4: 1 contig + reference = 1 task
        let contig = format!(">contig1\n{}\n", "A".repeat(5000));
        let reference = format!(">ref\n{}\n", "A".repeat(5000));

        let contig_f = write_fasta(&contig);
        let ref_f = write_fasta(&reference);

        let (use_case, types) =
            detect_use_case(&[contig_f.path().to_path_buf()], Some(ref_f.path())).unwrap();

        assert_eq!(
            use_case,
            UseCase::Case4 {
                n_tasks: 1,
                is_msa: false
            }
        );
        assert_eq!(types, vec![InputType::Contig]);
    }

    #[test]
    fn test_issue_67_case4_multiple_contigs_with_reference() {
        // Case 4: M contigs + reference = M tasks
        let contig1 = format!(">contig1\n{}\n", "A".repeat(5000));
        let contig2 = format!(">contig2\n{}\n", "T".repeat(5000));
        let reference = format!(">ref\n{}\n", "A".repeat(5000));

        let contig1_f = write_fasta(&contig1);
        let contig2_f = write_fasta(&contig2);
        let ref_f = write_fasta(&reference);

        let (use_case, types) = detect_use_case(
            &[
                contig1_f.path().to_path_buf(),
                contig2_f.path().to_path_buf(),
            ],
            Some(ref_f.path()),
        )
        .unwrap();

        assert_eq!(
            use_case,
            UseCase::Case4 {
                n_tasks: 2,
                is_msa: false
            }
        );
        assert_eq!(types, vec![InputType::Contig, InputType::Contig]);
    }

    #[test]
    fn test_issue_67_case4_contigs_only_msa_mode() {
        // Case 4: M contigs, no reference = MSA mode, M×(M-1)/2 tasks
        let contig1 = format!(">contig1\n{}\n", "A".repeat(5000));
        let contig2 = format!(">contig2\n{}\n", "T".repeat(5000));
        let contig3 = format!(">contig3\n{}\n", "G".repeat(5000));

        let contig1_f = write_fasta(&contig1);
        let contig2_f = write_fasta(&contig2);
        let contig3_f = write_fasta(&contig3);

        let (use_case, types) = detect_use_case(
            &[
                contig1_f.path().to_path_buf(),
                contig2_f.path().to_path_buf(),
                contig3_f.path().to_path_buf(),
            ],
            None,
        )
        .unwrap();

        // 3 contigs, no ref → 3×2/2 = 3 tasks
        assert_eq!(
            use_case,
            UseCase::Case4 {
                n_tasks: 3,
                is_msa: true
            }
        );
        assert_eq!(
            types,
            vec![InputType::Contig, InputType::Contig, InputType::Contig]
        );
    }

    #[test]
    fn test_issue_67_case4_two_contigs_msa() {
        // Case 4: 2 contigs only = MSA, 2×1/2 = 1 task
        let contig1 = format!(">contig1\n{}\n", "A".repeat(5000));
        let contig2 = format!(">contig2\n{}\n", "T".repeat(5000));

        let contig1_f = write_fasta(&contig1);
        let contig2_f = write_fasta(&contig2);

        let (use_case, _) = detect_use_case(
            &[
                contig1_f.path().to_path_buf(),
                contig2_f.path().to_path_buf(),
            ],
            None,
        )
        .unwrap();

        assert_eq!(
            use_case,
            UseCase::Case4 {
                n_tasks: 1,
                is_msa: true
            }
        );
    }

    // ===== Error cases =====

    #[test]
    fn test_issue_67_case1_error_reads_only_no_reference() {
        // Case 1 should return error: deferred to Phase 2
        let read = ">read1\nACGT\n";
        let read_f = write_fasta(read);

        let result = detect_use_case(&[read_f.path().to_path_buf()], None);

        assert!(matches!(result, Err(UseCaseError::Case1Deferred)));
    }

    #[test]
    fn test_issue_67_multiple_reads_only_error() {
        // Multiple reads with no reference → Case 1 error
        let read1 = ">read1\nACGT\n";
        let read2 = ">read2\nTGCA\n";
        let read1_f = write_fasta(read1);
        let read2_f = write_fasta(read2);

        let result = detect_use_case(
            &[read1_f.path().to_path_buf(), read2_f.path().to_path_buf()],
            None,
        );

        assert!(matches!(result, Err(UseCaseError::Case1Deferred)));
    }

    #[test]
    fn test_issue_67_no_inputs_error() {
        // Empty inputs should error
        let empty: &[PathBuf] = &[];
        let result = detect_use_case(empty, None);
        assert!(matches!(result, Err(UseCaseError::InvalidInput(_))));
    }

    #[test]
    fn test_issue_67_invalid_file_io_error() {
        // Non-existent file should return IoError
        let result = detect_use_case(&[PathBuf::from("/nonexistent/path/file.fa")], None);
        assert!(matches!(result, Err(UseCaseError::IoError(_))));
    }

    #[test]
    fn test_issue_67_malformed_fasta_error() {
        // Invalid FASTA (no > or @) should error
        let mut f = NamedTempFile::new().unwrap();
        write!(f, "ACGTACGTACGT\n").unwrap();
        f.flush().unwrap();

        let result = detect_use_case(&[f.path().to_path_buf()], None);
        assert!(matches!(result, Err(UseCaseError::IoError(_))));
    }

    // ===== Edge cases =====

    #[test]
    fn test_issue_67_read_and_contig_with_reference() {
        // If reference present and there are reads, should be Case 2
        let read = ">read1\nACGT\n";
        let contig = format!(">contig1\n{}\n", "A".repeat(5000));
        let reference = ">ref\nACGT\n";

        let read_f = write_fasta(read);
        let contig_f = write_fasta(&contig);
        let ref_f = write_fasta(reference);

        let (use_case, _) = detect_use_case(
            &[read_f.path().to_path_buf(), contig_f.path().to_path_buf()],
            Some(ref_f.path()),
        )
        .unwrap();

        // With reference and reads, should be Case 2
        if let UseCase::Case2 { n_tasks } = use_case {
            assert_eq!(n_tasks, 1); // Only the read counts
        } else {
            panic!("expected Case2");
        }
    }

    #[test]
    fn test_issue_67_classification_consistency() {
        // Input types should match the detected use case
        let contig = format!(">contig1\n{}\n", "A".repeat(5000));
        let read = ">read1\nACGT\n";

        let contig_f = write_fasta(&contig);
        let read_f = write_fasta(read);

        let (_, types) = detect_use_case(
            &[contig_f.path().to_path_buf(), read_f.path().to_path_buf()],
            None,
        )
        .unwrap();

        assert_eq!(types.len(), 2);
        assert_eq!(types[0], InputType::Contig);
        assert_eq!(types[1], InputType::Read);
    }

    #[test]
    fn test_issue_67_fastq_reads_with_reference() {
        // FASTQ reads should be classified as Read
        let fastq = "@read1\nACGTACGT\n+\nIIIIIIII\n";
        let reference = ">ref\nACGT\n";

        let read_f = write_fastq(fastq);
        let ref_f = write_fasta(reference);

        let (use_case, types) =
            detect_use_case(&[read_f.path().to_path_buf()], Some(ref_f.path())).unwrap();

        assert_eq!(use_case, UseCase::Case2 { n_tasks: 1 });
        assert_eq!(types, vec![InputType::Read]);
    }
}
