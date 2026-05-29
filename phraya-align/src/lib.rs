pub mod wfa_simd;
pub mod wfa_simd_dispatch;
pub mod wfa_simd_safety;

/// Seed anchor position for WFA extension (legacy, kept for compatibility).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SeedAnchor {
    pub query_pos: usize,
    pub target_pos: usize,
}

/// Alignment result from WFA.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alignment {
    pub cigar: String,
    pub edit_distance: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub target_start: usize,
    pub target_end: usize,
}

/// Result type for WFA operations (legacy, kept for compatibility).
pub type WfaResult = Result<Alignment, WfaError>;

/// Errors that can occur during WFA alignment (legacy, kept for compatibility).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WfaError {
    InvalidInput(String),
    AlignmentFailed(String),
}

/// Runtime-dispatched WFA extension using multiversion.
///
/// Automatically selects SSE4.2 (if available) or naive implementation
/// based on runtime CPU feature detection via the multiversion crate.
///
/// # Arguments
/// * `query` - Query sequence
/// * `target` - Target sequence
/// * `seed_pos` - Seed position (starting point for alignment, typically 0 for full alignment)
///
/// # Returns
/// Alignment struct with CIGAR string, edit distance, and alignment positions.
#[cfg(target_arch = "x86_64")]
pub fn wfa_extend(query: &[u8], target: &[u8], seed_pos: usize) -> Alignment {
    // Use runtime dispatch: try SSE4.2 first, fall back to naive
    if is_x86_feature_detected!("sse4.2") {
        wfa_simd::wfa_extend_simd_impl(query, target, seed_pos)
    } else {
        wfa_simd::wfa_extend_naive_impl(query, target, seed_pos)
    }
}

#[cfg(not(target_arch = "x86_64"))]
pub fn wfa_extend(query: &[u8], target: &[u8], seed_pos: usize) -> Alignment {
    // Non-x86 platforms use naive implementation
    wfa_simd::wfa_extend_naive_impl(query, target, seed_pos)
}
