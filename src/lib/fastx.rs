//! A FASTA and FASTQ reader and record, compatible with `bio::io::fasta`
//! and `bio::io::fastq`.
//!
//! `Record` has the accessors of bio's records (`id`, `desc`, `seq`,
//! `qual`, `check`, `is_empty`) and `Reader` reads like bio's readers, with
//! `read` into a reused record or `records` as an iterator. The format is
//! detected from the first byte of the input. Unlike bio, the record keeps
//! bytes in buffers that are reused from one read to the next. A read loop
//! thus allocates nothing per record. FASTQ reads with an empty sequence,
//! which adapter trimming commonly produces and vsearch accepts, are read
//! rather than rejected.
//!
//! Records are checked as vsearch checks them. A header ends at a carriage
//! return, may hold no control character but tab, and is split at its first
//! space or tab into an id and a description. Sequences, and FASTQ
//! qualities, may span several lines, and carriage returns in them are
//! dropped. A FASTQ sequence may hold only IUPAC nucleotide codes, and its
//! qualities only printable symbols. In a FASTA sequence, tabs, vertical
//! tabs and form feeds are dropped silently. '-', '.' and other control
//! characters are errors. Anything else that is not an IUPAC code is
//! dropped and counted, for `Reader::stripped_warning`.

use std::io::{self, BufRead};

use crate::simd::Level;

/// The format of a sequence file.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Format {
    Fasta,
    Fastq,
}

/// A FASTA or FASTQ record. A FASTA record has no qualities.
///
/// # Example
/// ```rust
/// # use rsearch::fastx::Record;
/// let record = Record::with_attrs("a", Some("desc"), b"ACGT", b"IIII");
/// assert_eq!(record.id(), "a");
/// assert_eq!(record.desc(), Some("desc"));
/// assert_eq!(record.seq(), b"ACGT");
/// assert!(record.check().is_ok());
///
/// let record = Record::fasta_with_attrs("b", None, b"ACGT");
/// assert_eq!(record.qual(), b"");
/// assert!(record.check().is_ok());
/// ```
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Record {
    /// The header without '>' or '@'.
    header: String,
    /// Where the id ends in the header.
    id_end: usize,
    /// Where the description starts in the header, if there is one.
    desc_start: Option<usize>,
    seq: Vec<u8>,
    qual: Vec<u8>,
    format: Option<Format>,
}

/// Why a record fails `Record::check`, as in bio.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CheckError {
    EmptyId,
    NonAsciiSequence,
    InvalidSequence,
    NonAsciiQualities,
    UnequalLength,
}

impl Record {
    /// returns an empty record.
    pub fn new() -> Record {
        Record::default()
    }

    fn build(id: &str, desc: Option<&str>, seq: &[u8], qual: &[u8], format: Format) -> Record {
        let mut header = id.to_owned();
        if let Some(desc) = desc {
            header.push(' ');
            header.push_str(desc);
        }
        Record {
            header,
            id_end: id.len(),
            desc_start: desc.map(|_| id.len() + 1),
            seq: seq.to_vec(),
            qual: qual.to_vec(),
            format: Some(format),
        }
    }

    /// returns a FASTQ record with the given attributes.
    pub fn with_attrs(id: &str, desc: Option<&str>, seq: &[u8], qual: &[u8]) -> Record {
        Record::build(id, desc, seq, qual, Format::Fastq)
    }

    /// returns a FASTA record with the given attributes.
    pub fn fasta_with_attrs(id: &str, desc: Option<&str>, seq: &[u8]) -> Record {
        Record::build(id, desc, seq, b"", Format::Fasta)
    }

    /// whether the record is empty, as after reading past the end.
    pub fn is_empty(&self) -> bool {
        self.format.is_none()
    }

    /// whether the record was read from FASTQ, and so has qualities.
    pub fn is_fastq(&self) -> bool {
        self.format == Some(Format::Fastq)
    }

    /// the whole header: the id and the description.
    pub fn header(&self) -> &str {
        &self.header
    }

    /// the identifier: the header up to the first space (FASTQ) or
    /// whitespace (FASTA).
    pub fn id(&self) -> &str {
        &self.header[..self.id_end]
    }

    /// the description: the header after the id and its separator.
    pub fn desc(&self) -> Option<&str> {
        self.desc_start.map(|start| &self.header[start..])
    }

    pub fn seq(&self) -> &[u8] {
        &self.seq
    }

    /// the quality string, empty for a FASTA record.
    pub fn qual(&self) -> &[u8] {
        &self.qual
    }

    /// checks the record as bio does. Both formats need a non-empty id and
    /// an ASCII sequence; FASTQ also needs an alphabetic (or '-', '.', '*')
    /// sequence, ASCII qualities, and as many qualities as bases.
    pub fn check(&self) -> Result<(), CheckError> {
        if self.id().is_empty() {
            return Err(CheckError::EmptyId);
        }
        if !self.seq.is_ascii() {
            return Err(CheckError::NonAsciiSequence);
        }
        if !self.is_fastq() {
            return Ok(());
        }
        if !self
            .seq
            .iter()
            .all(|b| b.is_ascii_alphabetic() || matches!(b, b'-' | b'.' | b'*'))
        {
            return Err(CheckError::InvalidSequence);
        }
        if !self.qual.is_ascii() {
            return Err(CheckError::NonAsciiQualities);
        }
        if self.seq.len() != self.qual.len() {
            return Err(CheckError::UnequalLength);
        }
        Ok(())
    }

    fn clear(&mut self) {
        self.header.clear();
        self.id_end = 0;
        self.desc_start = None;
        self.seq.clear();
        self.qual.clear();
        self.format = None;
    }

    /// sets the header, splitting it at the first space or tab.
    fn set_header(&mut self, header: &[u8]) -> io::Result<()> {
        let header =
            std::str::from_utf8(header).map_err(|_| invalid("header is not valid UTF-8"))?;
        self.header.push_str(header);
        match header.find([' ', '\t']) {
            Some(end) => {
                self.id_end = end;
                self.desc_start = Some(end + 1);
            }
            None => self.id_end = header.len(),
        }
        Ok(())
    }
}

/// What a reader does with a byte of a sequence or quality line.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Action {
    Accept,
    /// Dropped silently.
    Skip,
    /// Dropped and counted.
    Strip,
    Reject,
}

const fn is_iupac(byte: u8) -> bool {
    matches!(
        byte.to_ascii_uppercase(),
        b'A' | b'C'
            | b'G'
            | b'T'
            | b'U'
            | b'R'
            | b'Y'
            | b'S'
            | b'W'
            | b'K'
            | b'M'
            | b'D'
            | b'B'
            | b'H'
            | b'V'
            | b'N'
    )
}

const fn fasta_action(byte: u8) -> Action {
    match byte {
        _ if is_iupac(byte) => Action::Accept,
        b'\r' | b'\t' | 0x0b | 0x0c => Action::Skip,
        b'-' | b'.' => Action::Reject,
        _ if byte < b' ' => Action::Reject,
        _ => Action::Strip,
    }
}

const fn fastq_sequence_action(byte: u8) -> Action {
    match byte {
        _ if is_iupac(byte) => Action::Accept,
        b'\r' => Action::Skip,
        _ => Action::Reject,
    }
}

const fn fastq_quality_action(byte: u8) -> Action {
    match byte {
        b'!'..=b'~' => Action::Accept,
        b'\r' => Action::Skip,
        _ => Action::Reject,
    }
}

/// Which bytes a line may hold.
#[derive(Clone, Copy)]
enum Policy {
    Fasta,
    FastqSequence,
    FastqQuality,
}

/// the action of every byte under `policy`.
const fn table(policy: Policy) -> [Action; 256] {
    let mut table = [Action::Reject; 256];
    let mut byte = 0;
    while byte < 256 {
        table[byte] = match policy {
            Policy::Fasta => fasta_action(byte as u8),
            Policy::FastqSequence => fastq_sequence_action(byte as u8),
            Policy::FastqQuality => fastq_quality_action(byte as u8),
        };
        byte += 1;
    }
    table
}

static FASTA: [Action; 256] = table(Policy::Fasta);
static FASTQ_SEQUENCE: [Action; 256] = table(Policy::FastqSequence);
static FASTQ_QUALITY: [Action; 256] = table(Policy::FastqQuality);

/// whether a line holds anything but IUPAC codes, at the SIMD level
/// `simd::level` chose.
fn has_non_iupac(line: &[u8]) -> bool {
    has_non_iupac_at(crate::simd::level(), line)
}

/// whether a line holds anything but IUPAC codes, at `level`, by NEON or
/// AVX2 table lookups or else by byte comparisons the compiler vectorises.
fn has_non_iupac_at(level: Level, line: &[u8]) -> bool {
    match level {
        #[cfg(target_arch = "aarch64")]
        Level::Neon => has_non_iupac_neon(line),
        // SAFETY: the level is AVX2 only on a processor that has it
        #[cfg(target_arch = "x86_64")]
        Level::Avx2 => unsafe { has_non_iupac_avx2(line) },
        _ => has_non_iupac_portable(line),
    }
}

/// whether a line holds anything but IUPAC codes. Byte comparisons only, so
/// that it vectorizes. Clearing bit 5 maps exactly the letters to
/// 'A'..='Z', where the codes are six ranges.
fn has_non_iupac_portable(line: &[u8]) -> bool {
    let within = |b: u8, low: u8, high: u8| b.wrapping_sub(low) <= high - low;
    line.iter().fold(0u8, |other, &b| {
        let upper = b & !0x20;
        let iupac = within(upper, b'A', b'D')
            | within(upper, b'G', b'H')
            | (upper == b'K')
            | within(upper, b'M', b'N')
            | within(upper, b'R', b'W')
            | (upper == b'Y');
        other | u8::from(!iupac)
    }) != 0
}

/// The IUPAC codes by the low nibble of their byte: bit 0 for the rows
/// 0x4_ and 0x6_ ('A' to 'O'), bit 1 for the rows 0x5_ and 0x7_ ('P' to
/// 'Z'). A byte is a code when this and `IUPAC_HIGH` share a bit.
const IUPAC_LOW: [u8; 16] = {
    let mut table = [0; 16];
    let mut letter = b'A';
    while letter <= b'Z' {
        if is_iupac(letter) {
            table[(letter & 15) as usize] |= if letter < b'P' { 1 } else { 2 };
        }
        letter += 1;
    }
    table
};

/// The rows of letters by the high nibble of a byte; see `IUPAC_LOW`.
const IUPAC_HIGH: [u8; 16] = {
    let mut table = [0; 16];
    table[4] = 1;
    table[6] = 1;
    table[5] = 2;
    table[7] = 2;
    table
};

/// whether the bytes past the vector chunks hold anything but IUPAC codes,
/// by the same two tables.
fn has_non_iupac_rest(rest: &[u8]) -> bool {
    rest.iter()
        .any(|&b| IUPAC_LOW[usize::from(b & 15)] & IUPAC_HIGH[usize::from(b >> 4)] == 0)
}

/// whether a line holds anything but IUPAC codes: two 16-entry table
/// lookups per 16 bytes, with NEON's table instruction.
#[cfg(target_arch = "aarch64")]
fn has_non_iupac_neon(line: &[u8]) -> bool {
    use std::arch::aarch64::*;
    let (chunks, rest) = line.as_chunks::<16>();
    // SAFETY: NEON is part of every aarch64 target, and each load reads
    // 16 bytes from a chunk of exactly 16.
    let bad = unsafe {
        let low = vld1q_u8(IUPAC_LOW.as_ptr());
        let high = vld1q_u8(IUPAC_HIGH.as_ptr());
        let nibble = vdupq_n_u8(15);
        let mut bad = vdupq_n_u8(0);
        for chunk in chunks {
            let bytes = vld1q_u8(chunk.as_ptr());
            let by_low = vqtbl1q_u8(low, vandq_u8(bytes, nibble));
            let by_high = vqtbl1q_u8(high, vshrq_n_u8(bytes, 4));
            bad = vorrq_u8(bad, vceqzq_u8(vandq_u8(by_low, by_high)));
        }
        vmaxvq_u8(bad) != 0
    };
    bad || has_non_iupac_rest(rest)
}

/// whether a line holds anything but IUPAC codes, by two 16-entry table
/// lookups per 32 bytes with AVX2's byte shuffle. It looks up each 16-byte
/// half in its own copy of the table.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
fn has_non_iupac_avx2(line: &[u8]) -> bool {
    use std::arch::x86_64::*;
    let (chunks, rest) = line.as_chunks::<32>();
    // SAFETY: the tables are 16 bytes, and each load reads 32 bytes from a
    // chunk of exactly 32.
    let bad = unsafe {
        let low = _mm256_broadcastsi128_si256(_mm_loadu_si128(IUPAC_LOW.as_ptr().cast()));
        let high = _mm256_broadcastsi128_si256(_mm_loadu_si128(IUPAC_HIGH.as_ptr().cast()));
        let nibble = _mm256_set1_epi8(15);
        let zero = _mm256_setzero_si256();
        let mut bad = zero;
        for chunk in chunks {
            let bytes = _mm256_loadu_si256(chunk.as_ptr().cast());
            let by_low = _mm256_shuffle_epi8(low, _mm256_and_si256(bytes, nibble));
            // AVX2 has no byte shift. Shift the 16-bit words, then keep the
            // nibble.
            let high_nibbles = _mm256_and_si256(_mm256_srli_epi16::<4>(bytes), nibble);
            let by_high = _mm256_shuffle_epi8(high, high_nibbles);
            let codes = _mm256_and_si256(by_low, by_high);
            bad = _mm256_or_si256(bad, _mm256_cmpeq_epi8(codes, zero));
        }
        _mm256_movemask_epi8(bad) != 0
    };
    bad || has_non_iupac_rest(rest)
}

/// whether a line holds anything but printable quality symbols.
fn has_non_quality(line: &[u8]) -> bool {
    line.iter().fold(0u8, |other, &b| {
        other | u8::from(b.wrapping_sub(b'!') > b'~' - b'!')
    }) != 0
}

/// What a line holds, for error messages.
#[derive(Clone, Copy)]
enum Field {
    Header,
    Sequence,
    Quality,
}

/// vsearch's message about an illegal character.
fn illegal(field: Field, file: &str, byte: u8, line: u64) -> io::Error {
    let field = match field {
        Field::Header => "header",
        Field::Sequence => "sequence",
        Field::Quality => "quality",
    };
    let symbol = if (b' '..=b'~').contains(&byte) {
        format!("'{}'", char::from(byte))
    } else {
        format!("(unprintable, no {byte})")
    };
    invalid(&format!(
        "Illegal {field} character {symbol} on line {line} of {file} file"
    ))
}

/// vsearch's message about an invalid FASTQ line.
fn invalid_fastq(line: u64, message: &str) -> io::Error {
    invalid(&format!("Invalid line {line} in FASTQ file: {message}"))
}

/// Reads FASTA or FASTQ records from a buffered reader.
///
/// # Example
/// ```rust
/// # use rsearch::fastx::{Format, Reader, Record};
/// let input = b"@a desc\nACGT\n+\nIIII\n@empty\n\n+\n\n";
/// let mut reader = Reader::new(&input[..]);
/// assert_eq!(reader.format().unwrap(), Some(Format::Fastq));
/// let mut record = Record::new();
///
/// reader.read(&mut record).unwrap();
/// assert_eq!(record.id(), "a");
/// assert_eq!(record.desc(), Some("desc"));
///
/// reader.read(&mut record).unwrap();
/// assert_eq!(record.seq(), b"");
///
/// // at the end of the input, the record is empty
/// reader.read(&mut record).unwrap();
/// assert!(record.is_empty());
///
/// let mut reader = Reader::new(&b">b\tdesc\nAC\nGT\n"[..]);
/// reader.read(&mut record).unwrap();
/// assert_eq!((record.id(), record.seq()), ("b", &b"ACGT"[..]));
/// ```
pub struct Reader<R> {
    reader: R,
    /// The last line read, without its line feed.
    line: Vec<u8>,
    /// Whether the last line read ended with a line feed.
    line_complete: bool,
    /// Lines read so far.
    line_number: u64,
    /// The last header line read, as it was with its '>' or '@', and its
    /// line number.
    header: Vec<u8>,
    header_line: u64,
    /// How many times each byte was dropped from FASTA sequences.
    stripped: Box<[u64; 256]>,
}

fn invalid(message: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.to_string())
}

impl<R: BufRead> Reader<R> {
    pub fn new(reader: R) -> Reader<R> {
        Reader {
            reader,
            line: Vec::new(),
            line_complete: true,
            line_number: 0,
            header: Vec::new(),
            header_line: 0,
            stripped: Box::new([0; 256]),
        }
    }

    /// the format of the input, from its next byte: `None` at the end of
    /// the input.
    pub fn format(&mut self) -> io::Result<Option<Format>> {
        match self.peek()? {
            None => Ok(None),
            Some(b'>') => Ok(Some(Format::Fasta)),
            Some(b'@') => Ok(Some(Format::Fastq)),
            Some(_) => Err(invalid("File type not recognized.")),
        }
    }

    /// the next byte, if any, without reading it.
    fn peek(&mut self) -> io::Result<Option<u8>> {
        Ok(self.reader.fill_buf()?.first().copied())
    }

    /// the warning vsearch gives about the bytes dropped from FASTA
    /// sequences so far, if any were. The bytes are shown as they are. The
    /// warning may therefore not be UTF-8.
    ///
    /// # Example
    /// ```rust
    /// # use rsearch::fastx::{Reader, Record};
    /// let mut reader = Reader::new(&b">a\nAC1 G*T\n"[..]);
    /// let mut record = Record::new();
    /// reader.read(&mut record).unwrap();
    /// assert_eq!(record.seq(), b"ACGT");
    /// assert_eq!(
    ///     reader.stripped_warning().unwrap(),
    ///     b"3 invalid characters stripped from FASTA file:  (1) *(1) 1(1)"
    /// );
    /// ```
    pub fn stripped_warning(&self) -> Option<Vec<u8>> {
        let total: u64 = self.stripped.iter().sum();
        if total == 0 {
            return None;
        }
        let mut message =
            format!("{total} invalid characters stripped from FASTA file:").into_bytes();
        for (byte, &count) in self.stripped.iter().enumerate() {
            if count > 0 {
                message.push(b' ');
                message.push(byte as u8);
                message.extend(format!("({count})").into_bytes());
            }
        }
        Some(message)
    }

    /// reads the next line without its line feed; false at end of input.
    fn next_line(&mut self) -> io::Result<bool> {
        self.line.clear();
        if self.reader.read_until(b'\n', &mut self.line)? == 0 {
            return Ok(false);
        }
        self.line_complete = self.line.last() == Some(&b'\n');
        if self.line_complete {
            self.line.pop();
        }
        self.line_number += 1;
        Ok(true)
    }

    /// appends the bytes of the last line that `actions` accepts to `out`,
    /// dropping the others; returns the first rejected byte, if any.
    fn filter_line(&mut self, out: &mut Vec<u8>, actions: &[Action; 256]) -> Option<u8> {
        // most lines hold only accepted bytes and are copied whole
        let clean = if std::ptr::eq(actions, &FASTQ_QUALITY) {
            !has_non_quality(&self.line)
        } else {
            !has_non_iupac(&self.line)
        };
        if clean {
            out.extend_from_slice(&self.line);
            return None;
        }
        let mut rejected = None;
        for &byte in &self.line {
            match actions[usize::from(byte)] {
                Action::Accept => out.push(byte),
                Action::Skip => {}
                Action::Strip => self.stripped[usize::from(byte)] += 1,
                Action::Reject => {
                    rejected.get_or_insert(byte);
                }
            }
        }
        rejected
    }

    /// filters the last line into `out`, failing on a rejected byte.
    fn filter_sequence(
        &mut self,
        out: &mut Vec<u8>,
        actions: &[Action; 256],
        file: &str,
    ) -> io::Result<()> {
        match self.filter_line(out, actions) {
            Some(byte) => Err(illegal(Field::Sequence, file, byte, self.line_number)),
            None => Ok(()),
        }
    }

    /// keeps the last line as the header line, swapping buffers rather
    /// than copying it.
    fn read_header(&mut self) {
        std::mem::swap(&mut self.header, &mut self.line);
        self.header_line = self.line_number;
    }

    /// sets the record's header: cut at a carriage return or NUL, and free
    /// of other control characters but tab. vsearch checks it once the
    /// record is read.
    fn set_header(&self, record: &mut Record) -> io::Result<()> {
        let line = &self.header[1..];
        let end = line
            .iter()
            .position(|&b| b == b'\r' || b == 0)
            .unwrap_or(line.len());
        let header = &line[..end];
        if let Some(&byte) = header
            .iter()
            .find(|&&b| b == 0x7f || (b < b' ' && b != b'\t'))
        {
            return Err(illegal(
                Field::Header,
                "FASTA/FASTQ",
                byte,
                self.header_line,
            ));
        }
        record.set_header(header)
    }

    /// reads the next record into `record`, reusing its buffers. At the end
    /// of the input, `record` is left empty.
    pub fn read(&mut self, record: &mut Record) -> io::Result<()> {
        record.clear();
        match self.format()? {
            None => Ok(()),
            Some(Format::Fasta) => self.read_fasta(record),
            Some(Format::Fastq) => self.read_fastq(record),
        }
    }

    fn read_fasta(&mut self, record: &mut Record) -> io::Result<()> {
        self.next_line()?;
        if !self.line_complete {
            return Err(invalid(
                "Invalid FASTA - header must be terminated with newline",
            ));
        }
        self.read_header();
        self.set_header(record)?;
        // sequence lines, up to the next record
        while !matches!(self.peek()?, None | Some(b'>')) {
            self.next_line()?;
            self.filter_sequence(&mut record.seq, &FASTA, "FASTA")?;
        }
        record.format = Some(Format::Fasta);
        Ok(())
    }

    fn read_fastq(&mut self, record: &mut Record) -> io::Result<()> {
        // the line an unexpected end of the input is reported on
        let end_line = |reader: &Self| reader.line_number + u64::from(reader.line_complete);
        let unexpected_end =
            |reader: &Self| invalid_fastq(end_line(reader), "Unexpected end of file");

        self.next_line()?;
        if !self.line_complete {
            return Err(unexpected_end(self));
        }
        self.read_header();

        // sequence lines, up to a line starting with '+'
        loop {
            match self.peek()? {
                None => return Err(unexpected_end(self)),
                Some(b'+') if self.line_complete => break,
                Some(_) => {}
            }
            self.next_line()?;
            self.filter_sequence(&mut record.seq, &FASTQ_SEQUENCE, "FASTQ")?;
        }

        // the '+' line: empty, or the header again
        self.next_line()?;
        if !self.line_complete {
            return Err(unexpected_end(self));
        }
        let plus = &self.line[1..];
        if !(plus.is_empty() || plus == b"\r" || plus == &self.header[1..]) {
            return Err(invalid_fastq(
                self.line_number,
                "'+' line must be empty or identical to header",
            ));
        }

        // quality lines, until as long as the sequence and followed by the
        // next record or the end of the input
        let mut last_complete = false;
        let mut lines = 0;
        loop {
            match self.peek()? {
                None => break,
                Some(b'@') if last_complete && record.qual.len() == record.seq.len() => break,
                Some(_) => {}
            }
            self.next_line()?;
            lines += 1;
            let rejected = self.filter_line(&mut record.qual, &FASTQ_QUALITY);
            last_complete = self.line_complete;
            // a quality line too long is reported before an illegal symbol
            if record.qual.len() > record.seq.len() {
                break;
            }
            if let Some(byte) = rejected {
                return Err(illegal(Field::Quality, "FASTQ", byte, self.line_number));
            }
        }
        if record.qual.len() != record.seq.len() {
            let line = if lines == 0 {
                self.line_number + 1
            } else {
                self.line_number
            };
            return Err(invalid_fastq(
                line,
                "Sequence and quality lines must be equally long",
            ));
        }
        self.set_header(record)?;
        record.format = Some(Format::Fastq);
        Ok(())
    }

    /// returns an iterator over the records.
    pub fn records(self) -> Records<R> {
        Records { reader: self }
    }
}

/// An iterator over FASTA or FASTQ records.
pub struct Records<R> {
    reader: Reader<R>,
}

impl<R: BufRead> Iterator for Records<R> {
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut record = Record::new();
        match self.reader.read(&mut record) {
            Ok(()) if record.is_empty() => None,
            result => Some(result.map(|()| record)),
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn read(input: &[u8]) -> io::Result<Vec<Record>> {
        Reader::new(input).records().collect()
    }

    #[test]
    fn test_fastq() {
        let records = read(b"@a x y\nAC\nGT\n+a x y\nII\nII\n@b\nA\n+\n@\n").unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id(), "a");
        assert_eq!(records[0].desc(), Some("x y"));
        assert_eq!(records[0].header(), "a x y");
        assert_eq!(records[0].seq(), b"ACGT");
        assert_eq!(records[0].qual(), b"IIII");
        assert_eq!(records[1].desc(), None);
        // a quality line may start with '@'
        assert_eq!(records[1].qual(), b"@");
        assert!(records.iter().all(|r| r.is_fastq() && r.check().is_ok()));
    }

    #[test]
    fn test_fasta() {
        let records = read(b">a x\nAC\n\nGT \n>b\tdesc\n>c\nNN\n").unwrap();
        assert_eq!(records.len(), 3);
        assert_eq!(
            records[0],
            Record::fasta_with_attrs("a", Some("x"), b"ACGT")
        );
        // an id ends at a space or a tab
        assert_eq!((records[1].id(), records[1].desc()), ("b", Some("desc")));
        assert_eq!(records[1].seq(), b"");
        assert_eq!(records[2].qual(), b"");
        assert!(records.iter().all(|r| !r.is_fastq() && r.check().is_ok()));
    }

    #[test]
    fn test_format() {
        assert_eq!(Reader::new(&b""[..]).format().unwrap(), None);
        assert_eq!(
            Reader::new(&b">a\n"[..]).format().unwrap(),
            Some(Format::Fasta)
        );
        assert_eq!(
            Reader::new(&b"@a\n"[..]).format().unwrap(),
            Some(Format::Fastq)
        );
        assert!(Reader::new(&b"a\n"[..]).format().is_err());
        assert!(read(b"").unwrap().is_empty());
    }

    #[test]
    fn test_reuse() {
        let mut reader = Reader::new(&b"@long name\nACGTACGT\n+\nIIIIIIII\n@b\nA\n+\nI\n"[..]);
        let mut record = Record::new();
        reader.read(&mut record).unwrap();
        reader.read(&mut record).unwrap();
        // nothing is left over from the longer first record
        assert_eq!(record, Record::with_attrs("b", None, b"A", b"I"));
        reader.read(&mut record).unwrap();
        assert!(record.is_empty());
    }

    #[test]
    fn test_empty_reads() {
        let records = read(b"@a\n\n+\n\n@b\nA\n+\nI\n@c\n\n+\n").unwrap();
        assert_eq!(records.len(), 3);
        assert_eq!(records[0].seq(), b"");
        assert_eq!(records[1].seq(), b"A");
        assert_eq!(records[2].qual(), b"");
        assert!(records[0].check().is_ok());
    }

    #[test]
    fn test_crlf() {
        let records = read(b"@a d\r\nAC\r\n+\r\nII\r\n").unwrap();
        assert_eq!(records[0].desc(), Some("d"));
        assert_eq!(records[0].seq(), b"AC");
        let records = read(b">a d\r\nAC\r\nGT\r\n").unwrap();
        assert_eq!(records[0].seq(), b"ACGT");
    }

    #[test]
    fn test_fast_checks() {
        // the arithmetic checks agree with the tables on every byte
        for byte in 0..=255u8 {
            let line = [byte];
            let accepts = |table: &[Action; 256]| table[usize::from(byte)] == Action::Accept;
            assert_eq!(has_non_iupac(&line), !accepts(&FASTQ_SEQUENCE));
            assert_eq!(has_non_iupac_portable(&line), !accepts(&FASTQ_SEQUENCE));
            assert_eq!(has_non_iupac(&line), !accepts(&FASTA));
            assert_eq!(has_non_quality(&line), !accepts(&FASTQ_QUALITY));
        }
    }

    #[test]
    fn test_fast_checks_on_long_lines() {
        // every position of a 16-byte chunk, and the rest after the chunks
        for len in [15, 16, 17, 33, 48] {
            let line = vec![b'a'; len];
            assert!(!has_non_iupac(&line));
            for i in 0..len {
                let mut bad = line.clone();
                bad[i] = b'*';
                assert!(has_non_iupac(&bad), "{len} {i}");
            }
        }
    }

    #[test]
    fn test_finds_line_ends_at_every_level() {
        // each SIMD level finds the first newline anywhere in buffers that
        // end within, on, and after the 16- and 32-byte chunks, or none
        for level in crate::simd::available() {
            for len in [0, 1, 15, 16, 17, 31, 32, 33, 47, 48, 63, 64, 65, 100] {
                let none = vec![b'\r'; len];
                assert_eq!(find_newline_at(level, &none), None, "{level:?} {len}");
                for at in 0..len {
                    let mut bytes = none.clone();
                    bytes[at] = b'\n';
                    assert_eq!(find_newline_at(level, &bytes), Some(at), "{level:?} {len}");
                    if at + 1 < len {
                        bytes[len - 1] = b'\n';
                        assert_eq!(find_newline_at(level, &bytes), Some(at), "{level:?} {len}");
                    }
                }
            }
        }
    }

    #[test]
    fn test_reads_lines_at_every_level() {
        // lines across the reader's buffers, read as read_until reads them
        let text: Vec<u8> = (0..400u32)
            .flat_map(|i| {
                let mut line = vec![b'a' + (i % 26) as u8; (i * 7 % 70) as usize];
                line.push(b'\n');
                line
            })
            .chain(*b"last line without an end")
            .collect();
        for input in [&text[..], b"", b"\n", b"\n\n", b"x"] {
            let mut expected = Vec::new();
            let mut plain = input;
            loop {
                let mut line = Vec::new();
                let n = plain.read_until(b'\n', &mut line).unwrap();
                expected.push(line);
                if n == 0 {
                    break;
                }
            }
            for level in crate::simd::available() {
                for capacity in [1, 2, 3, 16, 33, 64, 8192] {
                    let mut buffered = io::BufReader::with_capacity(capacity, input);
                    for (i, want) in expected.iter().enumerate() {
                        let mut line = Vec::new();
                        let n = read_line_at(level, &mut buffered, &mut line).unwrap();
                        assert_eq!(n, want.len(), "{level:?} {capacity} line {i}");
                        assert_eq!(&line, want, "{level:?} {capacity} line {i}");
                    }
                }
            }
        }
    }

    #[test]
    fn test_fast_checks_at_every_level() {
        // each SIMD level this processor has agrees with the tables on every
        // byte, at the start, the middle, and the end of lines that end
        // within, on, and after the 16- and 32-byte chunks
        for level in crate::simd::available() {
            for len in [1, 15, 16, 17, 31, 32, 33, 47, 48, 63, 64, 65, 100] {
                let clean = vec![b'g'; len];
                assert!(!has_non_iupac_at(level, &clean), "{level:?} {len}");
                for byte in 0..=255u8 {
                    let rejected = FASTQ_SEQUENCE[usize::from(byte)] != Action::Accept;
                    for at in [0, len / 2, len - 1] {
                        let mut line = clean.clone();
                        line[at] = byte;
                        assert_eq!(
                            has_non_iupac_at(level, &line),
                            rejected,
                            "{level:?} {len} {at} {byte}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn test_check() {
        let check = |seq: &[u8], qual: &[u8]| Record::with_attrs("a", None, seq, qual).check();
        assert_eq!(check(b"AC-*.n", b"IIIIII"), Ok(()));
        assert_eq!(
            Record::with_attrs("", None, b"A", b"I").check(),
            Err(CheckError::EmptyId)
        );
        assert_eq!(check(b"A1", b"II"), Err(CheckError::InvalidSequence));
        assert_eq!(
            check(b"A\xc3\xa9", b"III"),
            Err(CheckError::NonAsciiSequence)
        );
        assert_eq!(check(b"A", b"\xc3"), Err(CheckError::NonAsciiQualities));
        assert_eq!(check(b"AC", b"I"), Err(CheckError::UnequalLength));
        // FASTA checks only the id and that the sequence is ASCII
        let fasta = |id, seq: &[u8]| Record::fasta_with_attrs(id, None, seq).check();
        assert_eq!(fasta("a", b"A1"), Ok(()));
        assert_eq!(fasta("", b"A"), Err(CheckError::EmptyId));
        assert_eq!(fasta("a", b"\xc3\xa9"), Err(CheckError::NonAsciiSequence));
    }

    #[test]
    fn test_errors() {
        assert!(read(b"a\nAC\n+\nII\n").is_err());
        assert!(read(b"@a\nAC\n").is_err());
        assert!(read(b"@a\nAC\n+\nI\n").is_err());
        assert!(read(b"@a\nAC\n+\nIII\n").is_err());
        assert!(read(b"@\xff\nA\n+\nI\n").is_err());
        assert!(read(b">\xff\nA\n").is_err());
    }
}
