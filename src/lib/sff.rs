//! A reader of 454 SFF (Standard Flowgram Format) files, version 1, with
//! vsearch --sff_convert's checks.
//!
//! An SFF file holds a header, then for every read a header (name and clip
//! points) and the read's data (flowgram, bases and Phred scores), each
//! padded to 8 bytes, and optionally an index, between two reads or after
//! the last. Numbers are big-endian. Flowgrams and the index are skipped.

use std::fmt;
use std::io::BufRead;

/// The first four bytes of an SFF file.
pub const MAGIC: &[u8; 4] = b".sff";

/// The fixed part of the header, without flow characters and key.
const HEADER_LENGTH: usize = 31;
/// The fixed part of a read header, without the name.
const READ_HEADER_LENGTH: usize = 16;
/// Length of the index kind (".mft1.00" and the like) opening the index.
const INDEX_KIND_LENGTH: u64 = 8;
const KEY_LENGTH: u16 = 4;

/// An invalid or truncated SFF file.
#[derive(Debug)]
pub struct SffError(String);

impl fmt::Display for SffError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Invalid SFF file. {}", self.0)
    }
}

impl std::error::Error for SffError {}

fn invalid(message: &str) -> SffError {
    SffError(message.to_string())
}

fn truncated(section: &str) -> SffError {
    SffError(format!("Unable to read {section}. File may be truncated."))
}

/// rounds up to a multiple of 8, the alignment of SFF sections.
fn round_up_to_8(n: u64) -> u64 {
    n.div_ceil(8) * 8
}

/// A read: its name, bases and Phred scores, and the region the clip
/// points keep, `clip_start..clip_end` (0-based). The clip points may be
/// inverted, leaving `clip_start > clip_end`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SffRecord {
    pub name: String,
    pub bases: Vec<u8>,
    pub quals: Vec<u8>,
    pub clip_start: usize,
    pub clip_end: usize,
}

/// Reads SFF records one by one.
///
/// # Example
/// ```rust
/// # use rsearch::sff::Reader;
/// // the header of a file without reads: magic, version 1, no index,
/// // 0 reads, header length 40, key length 4, 4 flows, format 1, flow
/// // characters, key and padding
/// let mut file = b".sff\0\0\0\x01".to_vec();
/// file.extend([0; 16]);
/// file.extend([0, 40, 0, 4, 0, 4, 1]);
/// file.extend(b"TACGTCAG\0");
///
/// let mut reader = Reader::new(&file[..]).unwrap();
/// assert_eq!(reader.reads(), 0);
/// assert!(reader.read().unwrap().is_none());
/// ```
pub struct Reader<R> {
    inner: R,
    /// Bytes read so far, to find the index.
    position: u64,
    reads: u32,
    read_count: u32,
    flows: u16,
    index_offset: u64,
    index_length: u64,
    index_done: bool,
}

impl<R: BufRead> Reader<R> {
    /// reads and checks the header.
    pub fn new(inner: R) -> Result<Reader<R>, SffError> {
        let mut reader = Reader {
            inner,
            position: 0,
            reads: 0,
            read_count: 0,
            flows: 0,
            index_offset: 0,
            index_length: 0,
            index_done: false,
        };
        let header = reader
            .bytes(HEADER_LENGTH)
            .map_err(|_| invalid("Unable to read header. File may be truncated."))?;
        let u16_at = |i: usize| u16::from_be_bytes([header[i], header[i + 1]]);
        let u32_at = |i: usize| u32::from_be_bytes(header[i..i + 4].try_into().unwrap());
        if header[..4] != MAGIC[..] {
            return Err(invalid(
                "Incorrect magic number. Must be 0x2e736666 (.sff).",
            ));
        }
        if u32_at(4) != 1 {
            return Err(invalid("Incorrect version. Must be 1."));
        }
        let index_offset = u64::from_be_bytes(header[8..16].try_into().unwrap());
        let index_length = u32_at(16);
        let (header_length, key_length, flows) = (u16_at(24), u16_at(26), u16_at(28));
        if header[30] != 1 {
            return Err(invalid("Incorrect flowgram format code. Must be 1."));
        }
        let fixed = HEADER_LENGTH as u64 + u64::from(flows) + u64::from(key_length);
        if u64::from(header_length) != round_up_to_8(fixed) {
            return Err(invalid("Incorrect header length."));
        }
        if key_length != KEY_LENGTH {
            return Err(invalid("Incorrect key length. Must be 4."));
        }
        if index_length != 0 && u64::from(index_length) < INDEX_KIND_LENGTH {
            return Err(invalid("Incorrect index size. Must be at least 8."));
        }
        reader.skip(u64::from(flows), "flow characters")?;
        reader.skip(u64::from(key_length), "key sequence")?;
        reader.skip(u64::from(header_length) - fixed, "read padding")?;
        reader.reads = u32_at(20);
        reader.flows = flows;
        reader.index_offset = index_offset;
        reader.index_length = u64::from(index_length);
        reader.index_done = index_offset == 0 || index_length == 0;
        Ok(reader)
    }

    /// the number of reads the header announces.
    pub fn reads(&self) -> u32 {
        self.reads
    }

    fn bytes(&mut self, n: usize) -> Result<Vec<u8>, SffError> {
        let mut buffer = Vec::new();
        self.read_into(&mut buffer, n, "File may be truncated")?;
        Ok(buffer)
    }

    /// reads `n` bytes into `buffer`, replacing its contents.
    fn read_into(&mut self, buffer: &mut Vec<u8>, n: usize, section: &str) -> Result<(), SffError> {
        buffer.clear();
        buffer.resize(n, 0);
        self.inner
            .read_exact(buffer)
            .map_err(|_| truncated(section))?;
        self.position += n as u64;
        Ok(())
    }

    /// skips `n` bytes, returning how many there were.
    fn skip_some(&mut self, n: u64) -> u64 {
        let mut skipped = 0;
        while skipped < n {
            let available = match self.inner.fill_buf() {
                Ok([]) | Err(_) => break,
                Ok(buffer) => buffer.len() as u64,
            };
            let step = available.min(n - skipped);
            self.inner.consume(step as usize);
            skipped += step;
        }
        self.position += skipped;
        skipped
    }

    fn skip(&mut self, n: u64, section: &str) -> Result<(), SffError> {
        if self.skip_some(n) < n {
            return Err(truncated(section));
        }
        Ok(())
    }

    /// skips the index if it starts here; `last` tells that no read
    /// follows, where padding after the index may be missing.
    fn skip_index(&mut self, last: bool) -> Result<(), SffError> {
        if self.index_done || self.position != self.index_offset {
            return Ok(());
        }
        let kind = self.bytes(INDEX_KIND_LENGTH as usize);
        if kind.is_err() {
            return Err(SffError(
                "Unable to read index header. File may be truncated.".into(),
            ));
        }
        let padding = round_up_to_8(self.index_length) - self.index_length;
        let rest = self.index_length - INDEX_KIND_LENGTH;
        let size = if last { rest } else { rest + padding };
        if self.skip_some(size) != size {
            return Err(SffError(
                "Unable to read entire index. File may be truncated.".into(),
            ));
        }
        if last {
            self.skip_some(padding);
        }
        self.index_done = true;
        Ok(())
    }

    /// reads the next record, or `None` after the last.
    pub fn read(&mut self) -> Result<Option<SffRecord>, SffError> {
        let mut record = SffRecord::default();
        Ok(self.read_record(&mut record)?.then_some(record))
    }

    /// reads the next record into `record`, reusing its buffers; returns
    /// false after the last.
    pub fn read_record(&mut self, record: &mut SffRecord) -> Result<bool, SffError> {
        if self.read_count == self.reads {
            self.skip_index(true)?;
            return Ok(false);
        }
        self.skip_index(false)?;
        self.read_count += 1;

        let header = self
            .bytes(READ_HEADER_LENGTH)
            .map_err(|_| invalid("Unable to read read header. File may be truncated."))?;
        let u16_at = |i: usize| u16::from_be_bytes([header[i], header[i + 1]]);
        let (header_length, name_length) = (u16_at(0), u16_at(2));
        let bases = u32::from_be_bytes(header[4..8].try_into().unwrap());
        let (qual_left, qual_right) = (u16_at(8), u16_at(10));
        let (adapter_left, adapter_right) = (u16_at(12), u16_at(14));

        let fixed = (READ_HEADER_LENGTH + usize::from(name_length)) as u64;
        if u64::from(header_length) != round_up_to_8(fixed) {
            return Err(invalid("Incorrect read header length."));
        }
        for (clip, name) in [
            (qual_left, "clip_qual_left"),
            (adapter_left, "clip_adapter_left"),
            (qual_right, "clip_qual_right"),
            (adapter_right, "clip_adapter_right"),
        ] {
            if u32::from(clip) > bases {
                return Err(SffError(format!("Incorrect {name} value.")));
            }
        }

        let mut name = std::mem::take(&mut record.bases);
        self.read_into(&mut name, usize::from(name_length), "read name")?;
        record.name.clear();
        record.name.push_str(&String::from_utf8_lossy(&name));
        record.bases = name;
        self.skip(u64::from(header_length) - fixed, "read header padding")?;

        let n = bases as usize;
        let flowgram = 2 * u64::from(self.flows);
        self.skip(flowgram, "flowgram values")?;
        self.skip(bases.into(), "flow indices")?;
        let mut bases_read = std::mem::take(&mut record.bases);
        self.read_into(&mut bases_read, n, "read length")?;
        record.bases = bases_read;
        let mut quals = std::mem::take(&mut record.quals);
        self.read_into(&mut quals, n, "quality scores")?;
        record.quals = quals;
        let data = flowgram + 3 * u64::from(bases);
        self.skip(round_up_to_8(data) - data, "read data padding")?;

        // clip points are 1-based, and 0 means none
        let clip_start = usize::from(qual_left.max(adapter_left).max(1)) - 1;
        let right = |clip: u16| if clip == 0 { n } else { usize::from(clip) };
        let clip_end = right(qual_right).min(right(adapter_right));

        record.clip_start = clip_start;
        record.clip_end = clip_end;
        Ok(true)
    }
}
