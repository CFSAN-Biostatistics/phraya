use serde::{Deserialize, Serialize};

/// Represents a biological sequence (DNA, RNA, or protein)
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Sequence {
    /// Sequence identifier (first token in header)
    pub id: Option<String>,
    /// Full header description (everything after the first space)
    pub description: Option<String>,
    /// The sequence data as bytes
    pub data: Vec<u8>,
    /// Quality scores (FASTQ only), same length as data
    pub quality: Option<Vec<u8>>,
    /// Optional pairing information for paired-end reads
    pub pairing_info: Option<PairingInfo>,
}

/// Information about paired-end reads
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct PairingInfo {
    /// R1 or R2
    pub read_number: u8,
    /// Optional mate pair identifier
    pub mate_id: Option<String>,
}
