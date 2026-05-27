/// FM-Index for fast substring searching in DNA sequences
#[derive(Clone, Debug)]
pub struct FmIndex {
    reference: Vec<u8>,
    // Placeholder fields for actual FM-index implementation
    // In a full implementation, this would contain:
    // - BWT (Burrows-Wheeler Transform)
    // - Suffix array or SA samples
    // - Occurrence count tables
}

impl FmIndex {
    /// Create a new FM-index from a reference sequence
    pub fn new(reference: &[u8]) -> Self {
        FmIndex {
            reference: reference.to_vec(),
        }
    }

    /// Get reference sequence
    pub fn reference(&self) -> &[u8] {
        &self.reference
    }

    /// Count occurrences of a pattern (placeholder)
    pub fn count_occurrences(&self, _pattern: &[u8]) -> usize {
        // Placeholder: returns 1 for any pattern
        1
    }
}

pub fn add(left: u64, right: u64) -> u64 {
    left + right
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }
}
