pub mod wfa_simd;
pub mod wfa_simd_dispatch;
pub mod wfa_simd_safety;

/// Seed anchor position for WFA extension.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SeedAnchor {
    pub query_pos: usize,
    pub target_pos: usize,
}

/// Alignment result from WFA.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alignment {
    pub cigar: String,
    pub score: i32,
}

/// Result type for WFA operations.
pub type WfaResult = Result<Alignment, WfaError>;

/// Errors that can occur during WFA alignment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WfaError {
    InvalidInput(String),
    AlignmentFailed(String),
}

// Core WFA implementation - naive baseline
/// Naive scalar WFA extension implementation.
///
/// This is the baseline implementation that all optimized versions match against.
/// Uses simple dynamic programming without SIMD optimization.
pub fn wfa_extend_naive(query: &[u8], target: &[u8], seed: SeedAnchor) -> WfaResult {
    wfa_simd::wfa_extend_naive_impl(query, target, seed)
}

/// SSE4.2-accelerated WFA extension.
///
/// Uses SSE4.2 intrinsics for diagonal fill operations.
/// Falls back to naive implementation if SSE4.2 is not available.
pub fn wfa_extend_simd(query: &[u8], target: &[u8], seed: SeedAnchor) -> WfaResult {
    #[cfg(target_arch = "x86_64")]
    {
        wfa_simd::wfa_extend_simd_impl(query, target, seed)
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        // On non-x86 platforms, fall back to naive
        wfa_simd::wfa_extend_naive_impl(query, target, seed)
    }
}

/// Runtime-dispatched WFA extension using multiversion.
///
/// Automatically selects SSE4.2 (if available) or naive implementation
/// based on runtime CPU feature detection via the multiversion crate.
#[cfg(target_arch = "x86_64")]
pub fn wfa_extend(query: &[u8], target: &[u8], seed: SeedAnchor) -> WfaResult {
    // Use runtime dispatch: try SSE4.2 first, fall back to naive
    if is_x86_feature_detected!("sse4.2") {
        wfa_simd::wfa_extend_simd_impl(query, target, seed)
    } else {
        wfa_simd::wfa_extend_naive_impl(query, target, seed)
    }
}

#[cfg(not(target_arch = "x86_64"))]
pub fn wfa_extend(query: &[u8], target: &[u8], seed: SeedAnchor) -> WfaResult {
    // Non-x86 platforms use naive implementation
    wfa_simd::wfa_extend_naive_impl(query, target, seed)
}
