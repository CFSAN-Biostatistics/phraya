//! CIGAR parsing and per-operation statistics.
//!
//! This codebase's CIGAR alphabet is `{M, X, I, D}` and its meaning is **not** the SAM
//! convention: `M` real match, `X` real mismatch, `D` query-extra (insertion in the query),
//! `I` target-extra (deletion in the query). See `phraya_align::executor` (the sole
//! producer) for the WFA traceback that emits these ops and the assertion oracle that
//! guarantees no other op character ever appears.

/// Parse a CIGAR string into ordered `(count, op)` pairs.
///
/// A bare operation letter with no preceding digits (malformed input) is treated as count
/// 1, matching the historical behavior of the parser this replaces.
pub fn parse_ops(cigar: &str) -> Vec<(usize, char)> {
    let mut ops = Vec::new();
    let mut count_str = String::new();
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            count_str.push(ch);
        } else {
            let count: usize = count_str.parse().unwrap_or(1);
            count_str.clear();
            ops.push((count, ch));
        }
    }
    ops
}

/// Per-operation totals for a CIGAR string.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct CigarStats {
    /// `M` — real matches.
    pub matches: u32,
    /// `X` — real mismatches.
    pub mismatches: u32,
    /// `D` — query-extra (insertion in the query relative to the target).
    pub query_extra: u32,
    /// `I` — target-extra (deletion in the query relative to the target).
    pub target_extra: u32,
}

impl CigarStats {
    /// Accumulate op totals from a CIGAR string. Unknown op characters are parsed (so
    /// `parse_ops`'s count consumption stays correct) but contribute to no field.
    pub fn parse(cigar: &str) -> Self {
        let mut stats = CigarStats::default();
        for (count, op) in parse_ops(cigar) {
            let count = count as u32;
            match op {
                'M' => stats.matches += count,
                'X' => stats.mismatches += count,
                'D' => stats.query_extra += count,
                'I' => stats.target_extra += count,
                _ => {}
            }
        }
        stats
    }

    /// Bases of the query consumed by the alignment: `M + X + D`.
    pub fn query_aligned_len(&self) -> u32 {
        self.matches + self.mismatches + self.query_extra
    }

    /// Bases of the target spanned by the alignment: `M + X + I`.
    pub fn target_aligned_len(&self) -> u32 {
        self.matches + self.mismatches + self.target_extra
    }

    /// Alignment columns: `M + X + I + D`.
    pub fn alignment_columns(&self) -> u32 {
        self.matches + self.mismatches + self.query_extra + self.target_extra
    }

    /// Fraction of alignment columns that are exact matches (BLAST/MUMmer convention:
    /// `M / (M+X+I+D)`). `0.0` when there are no columns.
    pub fn match_fraction(&self) -> f64 {
        let columns = self.alignment_columns();
        if columns == 0 {
            0.0
        } else {
            self.matches as f64 / columns as f64
        }
    }

    /// Absolute normalized identity: `1 - edit_distance / query_aligned_len` — the same
    /// quantity ADR-0011 stores per placement in `.phraya.queries` /
    /// `cross_space.phraya.queries`. Clamped to `[0.0, 1.0]`; `0.0` when
    /// `query_aligned_len() == 0`.
    pub fn identity(&self, edit_distance: u32) -> f64 {
        let qlen = self.query_aligned_len();
        if qlen == 0 {
            0.0
        } else {
            (1.0 - edit_distance as f64 / qlen as f64).clamp(0.0, 1.0)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_ops_reproduces_legacy_parser() {
        assert_eq!(parse_ops("10M2X3D5I"), vec![(10, 'M'), (2, 'X'), (3, 'D'), (5, 'I')]);
        // Bare op letter (no digits) parses as count 1 — legacy `unwrap_or(1)` behavior.
        assert_eq!(parse_ops("M"), vec![(1, 'M')]);
        assert_eq!(parse_ops(""), vec![]);
    }

    #[test]
    fn cigar_stats_op_accounting() {
        let stats = CigarStats::parse("10M2X3D5I");
        assert_eq!(stats.matches, 10);
        assert_eq!(stats.mismatches, 2);
        assert_eq!(stats.query_extra, 3);
        assert_eq!(stats.target_extra, 5);
        assert_eq!(stats.query_aligned_len(), 15); // M+X+D
        assert_eq!(stats.target_aligned_len(), 17); // M+X+I
        assert_eq!(stats.alignment_columns(), 20); // M+X+I+D
        assert_eq!(stats.match_fraction(), 0.5);
        assert_eq!(stats.identity(0), 1.0);
    }

    #[test]
    fn cigar_stats_empty_guards_against_nan() {
        let stats = CigarStats::default();
        assert_eq!(stats.alignment_columns(), 0);
        assert_eq!(stats.match_fraction(), 0.0);
        assert_eq!(stats.identity(0), 0.0);
    }

    #[test]
    fn cigar_stats_identity_clamped() {
        // edit_distance > query_aligned_len would otherwise go negative.
        let stats = CigarStats::parse("5M5X");
        assert_eq!(stats.identity(100), 0.0);
    }
}
