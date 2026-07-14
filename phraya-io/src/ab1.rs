/// AB1 (ABIF) parser for Sanger capillary sequencing traces.
///
/// ABIF format structure:
/// - Header: "ABIF" magic (4 bytes)
/// - Version: u32 big-endian (usually 101)
/// - Number of elements: u32 big-endian
/// - Number of elements again: u32 big-endian (repeated)
/// - Directory: array of 28-byte tag entries
/// - Data section: variable-length records
///
/// Each tag entry is 24 bytes:
/// - Tag name: 4 ASCII characters (e.g., "PBAS", "PCON")
/// - Tag number: u32 big-endian
/// - Element type: u16 big-endian (4=byte, 5=char, 7=int, 8=float, etc.)
/// - Element size: u16 big-endian (bytes per element)
/// - Number of elements: u32 big-endian
/// - Data size: u32 big-endian (total bytes)
/// - Data offset or value: u32 big-endian (file offset if data > 4 bytes, else inline)

use phraya_core::types::{ParseError, Sequence};
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;

const AB1_MAGIC: &[u8; 4] = b"ABIF";
const AB1_VERSION: u32 = 101;
const TAG_SIZE: usize = 24;
const PBAS_TAG: &[u8; 4] = b"PBAS";
const PCON_TAG: &[u8; 4] = b"PCON";

/// Parse an AB1 file and extract the basecalled sequence.
///
/// Returns a Sequence with:
/// - bases: from PBAS tag
/// - quality_scores: from PCON tag (if present)
/// - id: from filename or "ab1_sequence"
pub fn parse_ab1_file<P: AsRef<Path>>(path: P) -> Result<Sequence, ParseError> {
    let mut file =
        std::fs::File::open(path.as_ref()).map_err(|e| ParseError::InvalidFormat(e.to_string()))?;

    let mut header = [0u8; 16];
    file.read_exact(&mut header)
        .map_err(|e| ParseError::InvalidFormat(format!("failed to read AB1 header: {}", e)))?;

    // Verify magic bytes
    if &header[0..4] != AB1_MAGIC {
        return Err(ParseError::InvalidFormat(format!(
            "invalid AB1 magic bytes: expected ABIF, got {:?}",
            String::from_utf8_lossy(&header[0..4])
        )));
    }

    // Verify version (should be 101)
    let version = u32::from_be_bytes([header[4], header[5], header[6], header[7]]);
    if version != AB1_VERSION {
        return Err(ParseError::InvalidFormat(format!(
            "unsupported AB1 version: {}, expected {}",
            version, AB1_VERSION
        )));
    }

    // Read number of directory elements
    let num_elements = u32::from_be_bytes([header[8], header[9], header[10], header[11]]) as usize;

    // The second num_elements field (header[12..16]) is a duplicate and should match

    // Parse directory and find PBAS and PCON tags
    let mut pbas_offset = None;
    let mut pbas_len = None;
    let mut pcon_offset = None;
    let mut pcon_len = None;

    for i in 0..num_elements {
        let mut tag_entry = [0u8; TAG_SIZE];
        file.read_exact(&mut tag_entry).map_err(|e| {
            ParseError::InvalidFormat(format!("failed to read tag entry {}: {}", i, e))
        })?;

        let tag_name = &tag_entry[0..4];
        // let tag_number = u32::from_be_bytes([tag_entry[4], tag_entry[5], tag_entry[6], tag_entry[7]]);
        // let element_type = u16::from_be_bytes([tag_entry[8], tag_entry[9]]);
        // let element_size = u16::from_be_bytes([tag_entry[10], tag_entry[11]]);
        let _num_elements_in_tag =
            u32::from_be_bytes([tag_entry[12], tag_entry[13], tag_entry[14], tag_entry[15]]);
        let data_size = u32::from_be_bytes([tag_entry[16], tag_entry[17], tag_entry[18], tag_entry[19]]);
        let data_offset_or_value =
            u32::from_be_bytes([tag_entry[20], tag_entry[21], tag_entry[22], tag_entry[23]]);

        if tag_name == PBAS_TAG {
            pbas_offset = Some(data_offset_or_value as usize);
            pbas_len = Some(data_size as usize);
        } else if tag_name == PCON_TAG {
            pcon_offset = Some(data_offset_or_value as usize);
            pcon_len = Some(data_size as usize);
        }
    }

    // Extract PBAS (required)
    let pbas_offset = pbas_offset.ok_or_else(|| {
        ParseError::InvalidFormat("missing required PBAS tag (basecalls)".to_string())
    })?;
    let pbas_len = pbas_len.unwrap_or(0);

    if pbas_len == 0 {
        return Err(ParseError::InvalidFormat(
            "PBAS tag is empty (no bases)".to_string(),
        ));
    }

    // Read PBAS data
    file.seek(SeekFrom::Start(pbas_offset as u64))
        .map_err(|e| ParseError::InvalidFormat(format!("failed to seek to PBAS: {}", e)))?;

    let mut bases = vec![0u8; pbas_len];
    file.read_exact(&mut bases)
        .map_err(|e| ParseError::InvalidFormat(format!("failed to read PBAS data: {}", e)))?;

    // Extract PCON (optional)
    let quality_scores = if let (Some(pcon_offset), Some(pcon_len)) = (pcon_offset, pcon_len) {
        if pcon_len > 0 {
            file.seek(SeekFrom::Start(pcon_offset as u64))
                .map_err(|e| ParseError::InvalidFormat(format!("failed to seek to PCON: {}", e)))?;

            let mut quality = vec![0u8; pcon_len];
            file.read_exact(&mut quality).map_err(|e| {
                ParseError::InvalidFormat(format!("failed to read PCON data: {}", e))
            })?;

            // Verify quality length matches bases length
            if quality.len() != bases.len() {
                return Err(ParseError::InvalidFormat(format!(
                    "PCON quality length ({}) does not match PBAS bases length ({})",
                    quality.len(),
                    bases.len()
                )));
            }

            Some(quality)
        } else {
            None
        }
    } else {
        None
    };

    // Extract ID from filename (or use a default)
    let id = path
        .as_ref()
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("ab1_sequence")
        .to_string();

    Ok(Sequence::new(bases, quality_scores, id, None))
}
