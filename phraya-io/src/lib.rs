// FASTA/FASTQ parser placeholder
// This function signature is expected by the tests but not yet implemented (TDD RED phase)

use phraya_core::Sequence;
use std::path::Path;

/// Parse sequences from a FASTA or FASTQ file (auto-detecting format).
/// Returns an iterator of Sequence objects for memory-efficient streaming.
///
/// This is a placeholder signature - implementation will be added in the GREEN phase.
pub fn parse_sequences<P: AsRef<Path>>(
    _path: P,
) -> Result<Box<dyn Iterator<Item = Sequence>>, phraya_core::ParseError> {
    // Placeholder that will fail all tests
    Err(phraya_core::ParseError::InvalidFormat(
        "Not yet implemented".to_string(),
    ))
}

#[cfg(test)]
mod fasta_fastq_tests;
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn placeholder_fails_as_expected() {
        // This test documents that parse_sequences is not yet implemented
        let result = parse_sequences("nonexistent.fa");
        assert!(result.is_err());
    }
}
