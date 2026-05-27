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

// These function signatures are expected by the tests.
// They do not exist yet - tests will fail (RED phase).

/// Naive WFA extension (baseline scalar implementation).
/// This function does not exist yet - tests should fail.
pub fn wfa_extend_naive(_query: &[u8], _target: &[u8], _seed: SeedAnchor) -> WfaResult {
    unimplemented!("wfa_extend_naive not yet implemented - from issue #5")
}

/// SSE4.2-accelerated WFA extension.
/// This function does not exist yet - tests should fail.
pub fn wfa_extend_simd(_query: &[u8], _target: &[u8], _seed: SeedAnchor) -> WfaResult {
    unimplemented!("wfa_extend_simd not yet implemented")
}

/// Runtime-dispatched WFA extension (uses multiversion for SSE4.2 vs naive).
/// This function does not exist yet - tests should fail.
pub fn wfa_extend(_query: &[u8], _target: &[u8], _seed: SeedAnchor) -> WfaResult {
    unimplemented!("wfa_extend not yet implemented")
}
