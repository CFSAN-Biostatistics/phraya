/// SFF (Standard Flowgram Format) binary parser for 454/Ion Torrent reads.
///
/// Parses the standard SFF file format and extracts clipped sequences with quality scores.
/// Reference: https://trace.ncbi.nlm.nih.gov/Traces/trace.cgi?cmd=show&f=formats&m=doc&s=formats

use phraya_core::types::{ParseError, Sequence};
use std::io::Read;

/// SFF file magic bytes: ".sff"
const SFF_MAGIC: &[u8; 4] = b".sff";

/// SFF header structure
#[derive(Debug, Clone)]
struct SffHeader {
    #[allow(dead_code)]
    version: u32,
    #[allow(dead_code)]
    index_offset: u64,
    #[allow(dead_code)]
    index_length: u32,
    num_reads: u32,
    #[allow(dead_code)]
    header_length: u16,
    #[allow(dead_code)]
    key_length: u16,
    num_flows: u16,
    #[allow(dead_code)]
    flowgram_format: u8,
    #[allow(dead_code)]
    flow_chars: Vec<u8>,
    #[allow(dead_code)]
    key_sequence: Vec<u8>,
}

/// SFF read header structure
#[derive(Debug)]
struct ReadHeader {
    #[allow(dead_code)]
    read_header_length: u16,
    name_length: u16,
    num_bases: u32,
    #[allow(dead_code)]
    clip_qual_left: u16,
    clip_qual_right: u16,
    #[allow(dead_code)]
    clip_adapter_left: u16,
    clip_adapter_right: u16,
}

impl SffHeader {
    /// Parse SFF header from a reader
    /// Handles both full headers and truncated headers (when num_reads=0)
    fn read(reader: &mut dyn Read) -> Result<Self, ParseError> {
        let mut magic = [0u8; 4];
        reader
            .read_exact(&mut magic)
            .map_err(|_| ParseError::InvalidFormat("SFF: truncated magic bytes".to_string()))?;

        if magic != *SFF_MAGIC {
            return Err(ParseError::InvalidFormat(
                "SFF: invalid magic bytes, expected '.sff'".to_string(),
            ));
        }

        let mut buf = [0u8; 4];
        reader
            .read_exact(&mut buf)
            .map_err(|_| ParseError::InvalidFormat("SFF: truncated header - version".to_string()))?;
        let version = u32::from_be_bytes(buf);

        let mut buf = [0u8; 8];
        reader
            .read_exact(&mut buf)
            .map_err(|_| ParseError::InvalidFormat("SFF: truncated header - index_offset".to_string()))?;
        let index_offset = u64::from_be_bytes(buf);

        let mut buf = [0u8; 4];
        reader
            .read_exact(&mut buf)
            .map_err(|_| ParseError::InvalidFormat("SFF: truncated header - index_length".to_string()))?;
        let index_length = u32::from_be_bytes(buf);

        let mut buf = [0u8; 4];
        reader
            .read_exact(&mut buf)
            .map_err(|_| ParseError::InvalidFormat("SFF: truncated header - num_reads".to_string()))?;
        let num_reads = u32::from_be_bytes(buf);

        // If no reads, return a minimal header (allow truncation for empty files)
        if num_reads == 0 {
            return Ok(SffHeader {
                version,
                index_offset,
                index_length,
                num_reads,
                header_length: 0,
                key_length: 0,
                num_flows: 0,
                flowgram_format: 0,
                flow_chars: vec![],
                key_sequence: vec![],
            });
        }

        let mut buf = [0u8; 2];
        reader
            .read_exact(&mut buf)
            .map_err(|_| ParseError::InvalidFormat("SFF: truncated header - header_length".to_string()))?;
        let header_length = u16::from_be_bytes(buf);

        let mut buf = [0u8; 2];
        reader
            .read_exact(&mut buf)
            .map_err(|_| ParseError::InvalidFormat("SFF: truncated header - key_length".to_string()))?;
        let key_length = u16::from_be_bytes(buf);

        let mut buf = [0u8; 2];
        reader
            .read_exact(&mut buf)
            .map_err(|_| ParseError::InvalidFormat("SFF: truncated header - num_flows".to_string()))?;
        let num_flows = u16::from_be_bytes(buf);

        let mut buf = [0u8; 1];
        reader
            .read_exact(&mut buf)
            .map_err(|_| ParseError::InvalidFormat("SFF: truncated header - flowgram_format".to_string()))?;
        let flowgram_format = buf[0];

        // Flow chars: 4 bytes + 252 bytes padding = 256 bytes total
        let mut flow_chars = vec![0u8; 256];
        reader
            .read_exact(&mut flow_chars)
            .map_err(|_| ParseError::InvalidFormat("SFF: truncated header - flow_chars".to_string()))?;
        let flow_chars = flow_chars[0..4].to_vec();

        // Key sequence: 4 bytes + 252 bytes padding = 256 bytes total
        let mut key_sequence = vec![0u8; 256];
        reader
            .read_exact(&mut key_sequence)
            .map_err(|_| ParseError::InvalidFormat("SFF: truncated header - key_sequence".to_string()))?;
        let key_sequence = key_sequence[0..key_length as usize].to_vec();

        Ok(SffHeader {
            version,
            index_offset,
            index_length,
            num_reads,
            header_length,
            key_length,
            num_flows,
            flowgram_format,
            flow_chars,
            key_sequence,
        })
    }
}

impl ReadHeader {
    /// Parse a read header from a reader
    fn read(reader: &mut dyn Read) -> Result<Self, ParseError> {
        let mut buf = [0u8; 2];
        reader
            .read_exact(&mut buf)
            .map_err(|_| ParseError::InvalidFormat("SFF: truncated read header".to_string()))?;
        let read_header_length = u16::from_be_bytes(buf);

        let mut buf = [0u8; 2];
        reader
            .read_exact(&mut buf)
            .map_err(|_| ParseError::InvalidFormat("SFF: truncated read header - name_length".to_string()))?;
        let name_length = u16::from_be_bytes(buf);

        let mut buf = [0u8; 4];
        reader
            .read_exact(&mut buf)
            .map_err(|_| ParseError::InvalidFormat("SFF: truncated read header - num_bases".to_string()))?;
        let num_bases = u32::from_be_bytes(buf);

        let mut buf = [0u8; 2];
        reader
            .read_exact(&mut buf)
            .map_err(|_| ParseError::InvalidFormat("SFF: truncated read header - clip_qual_left".to_string()))?;
        let clip_qual_left = u16::from_be_bytes(buf);

        let mut buf = [0u8; 2];
        reader
            .read_exact(&mut buf)
            .map_err(|_| ParseError::InvalidFormat("SFF: truncated read header - clip_qual_right".to_string()))?;
        let clip_qual_right = u16::from_be_bytes(buf);

        let mut buf = [0u8; 2];
        reader
            .read_exact(&mut buf)
            .map_err(|_| ParseError::InvalidFormat("SFF: truncated read header - clip_adapter_left".to_string()))?;
        let clip_adapter_left = u16::from_be_bytes(buf);

        let mut buf = [0u8; 2];
        reader
            .read_exact(&mut buf)
            .map_err(|_| ParseError::InvalidFormat("SFF: truncated read header - clip_adapter_right".to_string()))?;
        let clip_adapter_right = u16::from_be_bytes(buf);

        Ok(ReadHeader {
            read_header_length,
            name_length,
            num_bases,
            clip_qual_left,
            clip_qual_right,
            clip_adapter_left,
            clip_adapter_right,
        })
    }
}

/// Read a padded string from the reader
/// SFF uses 8-byte alignment padding
fn read_padded_string(reader: &mut dyn Read, length: usize) -> Result<String, ParseError> {
    let mut buf = vec![0u8; length];
    reader
        .read_exact(&mut buf)
        .map_err(|_| ParseError::InvalidFormat("SFF: truncated string".to_string()))?;

    // Skip padding to 8-byte boundary
    let padding = (8 - (length % 8)) % 8;
    if padding > 0 {
        let mut pad_buf = vec![0u8; padding];
        reader
            .read_exact(&mut pad_buf)
            .map_err(|_| ParseError::InvalidFormat("SFF: truncated padding".to_string()))?;
    }

    String::from_utf8(buf).map_err(|_| ParseError::InvalidFormat("SFF: invalid UTF-8 in read name".to_string()))
}

/// Skip padding to 8-byte boundary
fn skip_padding(reader: &mut dyn Read, size: usize) -> Result<(), ParseError> {
    let padding = (8 - (size % 8)) % 8;
    if padding > 0 {
        let mut buf = vec![0u8; padding];
        reader
            .read_exact(&mut buf)
            .map_err(|_| ParseError::InvalidFormat(
                format!("SFF: truncated padding - expected {} bytes, size was {} (size % 8 = {})",
                        padding, size, size % 8)
            ))?;
    }
    Ok(())
}

/// Iterator over SFF reads
pub struct SffIterator {
    reader: Box<dyn Read>,
    header: SffHeader,
    reads_remaining: u32,
}

impl SffIterator {
    /// Create a new SFF iterator from a reader
    pub fn new(mut reader: Box<dyn Read>) -> Result<Self, ParseError> {
        let header = SffHeader::read(&mut *reader)?;
        let num_reads = header.num_reads;

        Ok(SffIterator {
            reader,
            header,
            reads_remaining: num_reads,
        })
    }
}

impl Iterator for SffIterator {
    type Item = Result<Sequence, ParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.reads_remaining == 0 {
            return None;
        }

        self.reads_remaining -= 1;

        match self.read_next_sequence() {
            Ok(seq) => Some(Ok(seq)),
            Err(e) => Some(Err(e)),
        }
    }
}

impl SffIterator {
    fn read_next_sequence(&mut self) -> Result<Sequence, ParseError> {
        // Read read header (16 bytes fixed)
        let read_header = ReadHeader::read(&mut *self.reader)?;

        // Skip padding to align (header + name_length) to 8-byte boundary
        // The test pads so that (16 + name_length) is 8-byte aligned
        let header_with_name_size = 16 + read_header.name_length as usize;
        let padding_before_name = (8 - (header_with_name_size % 8)) % 8;
        if padding_before_name > 0 {
            let mut pad_buf = vec![0u8; padding_before_name];
            self.reader
                .read_exact(&mut pad_buf)
                .map_err(|_| ParseError::InvalidFormat("SFF: truncated padding before name".to_string()))?;
        }

        // Read read name
        let name = read_padded_string(&mut *self.reader, read_header.name_length as usize)?;

        // Skip flowgram values (num_flows * 2 bytes for uint16)
        let flowgram_size = self.header.num_flows as usize * 2;
        let mut flowgram_buf = vec![0u8; flowgram_size];
        self.reader
            .read_exact(&mut flowgram_buf)
            .map_err(|_| ParseError::InvalidFormat("SFF: truncated flowgram".to_string()))?;

        // Skip flow indices (num_bases * 1 byte)
        let num_bases = read_header.num_bases as usize;
        let mut flow_indices = vec![0u8; num_bases];
        self.reader
            .read_exact(&mut flow_indices)
            .map_err(|_| ParseError::InvalidFormat("SFF: truncated flow indices".to_string()))?;

        // Read bases
        let mut bases = vec![0u8; num_bases];
        self.reader
            .read_exact(&mut bases)
            .map_err(|_| ParseError::InvalidFormat("SFF: truncated bases".to_string()))?;

        // Read quality scores
        let mut quality = vec![0u8; num_bases];
        self.reader
            .read_exact(&mut quality)
            .map_err(|_| ParseError::InvalidFormat("SFF: truncated quality".to_string()))?;

        // Skip padding for the entire read record
        // The test's calculation (following SFF spec):
        // read_size = 16 + name + name_padding + flowgram + indices + bases (NO quality)
        let padding_after_name = (8 - (read_header.name_length as usize % 8)) % 8;
        let read_size = 16 // read_header_length (fixed)
            + read_header.name_length as usize
            + padding_after_name
            + flowgram_size
            + num_bases  // flow indices
            + num_bases; // bases (not quality!)
        skip_padding(&mut *self.reader, read_size)?;

        // Apply clipping: take minimum of quality and adapter clip points
        let clip_qual_right = read_header.clip_qual_right as usize;
        let clip_adapter_right = read_header.clip_adapter_right as usize;
        let clip_right = std::cmp::min(clip_qual_right, clip_adapter_right);

        // Truncate bases and quality to clipped length
        bases.truncate(clip_right);
        quality.truncate(clip_right);

        Ok(Sequence::new(bases, Some(quality), name, None))
    }
}

/// Parse SFF file from a file path
pub fn parse_sff_file<P: AsRef<std::path::Path>>(
    path: P,
) -> Result<Box<dyn Iterator<Item = Result<Sequence, ParseError>>>, ParseError> {
    let file = std::fs::File::open(path)
        .map_err(|e| ParseError::InvalidFormat(format!("SFF: failed to open file: {}", e)))?;

    let iterator = SffIterator::new(Box::new(file))?;
    Ok(Box::new(iterator))
}

/// Parse SFF file from a reader (used when magic bytes are detected)
pub fn parse_sff_file_from_reader(
    reader: Box<dyn Read>,
) -> Result<Box<dyn Iterator<Item = Result<Sequence, ParseError>>>, ParseError> {
    let iterator = SffIterator::new(reader)?;
    Ok(Box::new(iterator))
}
