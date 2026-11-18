//! A FASTQ reader and record, compatible with `bio::io::fastq`.
//!
//! `Record` has the accessors of bio's record (`id`, `desc`, `seq`,
//! `qual`, `check`, `is_empty`) and `Reader` reads like bio's reader, with
//! `read` into a reused record or `records` as an iterator. Unlike bio, the
//! record keeps bytes in buffers that are reused from one read to the next,
//! so a read loop allocates nothing per record; and reads with an empty
//! sequence, which adapter trimming commonly produces and vsearch accepts,
//! are read rather than rejected.
//!
//! Records are parsed as bio does: the header is split at the first space
//! into an id and a description, and sequences and qualities may span
//! several lines.

use std::io::{self, BufRead};

/// A FASTQ record.
///
/// # Example
/// ```rust
/// # use rsearch::fastq::Record;
/// let record = Record::with_attrs("a", Some("desc"), b"ACGT", b"IIII");
/// assert_eq!(record.id(), "a");
/// assert_eq!(record.desc(), Some("desc"));
/// assert_eq!(record.seq(), b"ACGT");
/// assert!(record.check().is_ok());
/// ```
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Record {
    /// The header without '@'.
    header: String,
    /// Where the id ends in the header: at the first space, or the end.
    id_end: usize,
    seq: Vec<u8>,
    qual: Vec<u8>,
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

    /// returns a record with the given attributes.
    pub fn with_attrs(id: &str, desc: Option<&str>, seq: &[u8], qual: &[u8]) -> Record {
        let mut header = id.to_owned();
        if let Some(desc) = desc {
            header.push(' ');
            header.push_str(desc);
        }
        Record {
            header,
            id_end: id.len(),
            seq: seq.to_vec(),
            qual: qual.to_vec(),
        }
    }

    /// whether the record is empty, as after reading past the end.
    pub fn is_empty(&self) -> bool {
        self.header.is_empty() && self.seq.is_empty() && self.qual.is_empty()
    }

    /// the whole header: the id and the description.
    pub fn header(&self) -> &str {
        &self.header
    }

    /// the identifier: the header up to the first space.
    pub fn id(&self) -> &str {
        &self.header[..self.id_end]
    }

    /// the description: the header after the first space.
    pub fn desc(&self) -> Option<&str> {
        self.header.get(self.id_end + 1..)
    }

    pub fn seq(&self) -> &[u8] {
        &self.seq
    }

    pub fn qual(&self) -> &[u8] {
        &self.qual
    }

    /// checks the record as bio does: a non-empty id, an alphabetic (or
    /// '-', '.', '*') ASCII sequence, ASCII qualities, and as many
    /// qualities as bases.
    pub fn check(&self) -> Result<(), CheckError> {
        if self.id().is_empty() {
            return Err(CheckError::EmptyId);
        }
        if !self.seq.is_ascii() {
            return Err(CheckError::NonAsciiSequence);
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
        self.seq.clear();
        self.qual.clear();
    }
}

/// Reads FASTQ records from a buffered reader.
///
/// # Example
/// ```rust
/// # use rsearch::fastq::{Reader, Record};
/// let input = b"@a desc\nACGT\n+\nIIII\n@empty\n\n+\n\n";
/// let mut reader = Reader::new(&input[..]);
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
/// ```
pub struct Reader<R> {
    reader: R,
    line: Vec<u8>,
}

fn invalid(message: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.to_string())
}

impl<R: BufRead> Reader<R> {
    pub fn new(reader: R) -> Reader<R> {
        Reader {
            reader,
            line: Vec::new(),
        }
    }

    /// reads the next line without its line ending; false at end of input.
    fn next_line(&mut self) -> io::Result<bool> {
        self.line.clear();
        if self.reader.read_until(b'\n', &mut self.line)? == 0 {
            return Ok(false);
        }
        while matches!(self.line.last(), Some(b'\n' | b'\r')) {
            self.line.pop();
        }
        Ok(true)
    }

    /// reads the next record into `record`, reusing its buffers. At the end
    /// of the input, `record` is left empty.
    pub fn read(&mut self, record: &mut Record) -> io::Result<()> {
        record.clear();
        if !self.next_line()? {
            return Ok(());
        }
        let Some(header) = self.line.strip_prefix(b"@") else {
            return Err(invalid("FASTQ record does not start with '@'"));
        };
        let header = std::str::from_utf8(header.trim_ascii_end())
            .map_err(|_| invalid("FASTQ header is not valid UTF-8"))?;
        record.header.push_str(header);
        record.id_end = header.find(' ').unwrap_or(header.len());

        // sequence lines, up to the '+' line
        loop {
            if !self.next_line()? {
                return Err(invalid("incomplete FASTQ record"));
            }
            if self.line.first() == Some(&b'+') {
                break;
            }
            record.seq.extend_from_slice(self.line.trim_ascii_end());
        }

        // quality lines: at least one, until as long as the sequence
        loop {
            if !self.next_line()? {
                if record.qual.is_empty() && record.seq.is_empty() {
                    break;
                }
                return Err(invalid("incomplete FASTQ record"));
            }
            record.qual.extend_from_slice(self.line.trim_ascii_end());
            if record.qual.len() >= record.seq.len() {
                break;
            }
        }
        if record.qual.len() != record.seq.len() {
            return Err(invalid("FASTQ sequence and quality lengths differ"));
        }
        Ok(())
    }

    /// returns an iterator over the records.
    pub fn records(self) -> Records<R> {
        Records { reader: self }
    }
}

/// An iterator over FASTQ records.
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
    fn test_records() {
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
        assert!(records.iter().all(|r| r.check().is_ok()));
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
    }

    #[test]
    fn test_errors() {
        assert!(read(b"a\nAC\n+\nII\n").is_err());
        assert!(read(b"@a\nAC\n").is_err());
        assert!(read(b"@a\nAC\n+\nI\n").is_err());
        assert!(read(b"@a\nAC\n+\nIII\n").is_err());
        assert!(read(b"@\xff\nA\n+\nI\n").is_err());
    }
}
