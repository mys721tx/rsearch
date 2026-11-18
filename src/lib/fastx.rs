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
//! Records are parsed as bio does. A FASTA header is split at the first
//! whitespace and a FASTQ header at the first space, into an id and a
//! description. Sequences, and FASTQ qualities, may span several lines.

use std::io::{self, BufRead};

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

    /// sets the header, splitting it at the first character matching
    /// `separator`.
    fn set_header(&mut self, header: &[u8], separator: impl Fn(char) -> bool) -> io::Result<()> {
        let header = std::str::from_utf8(header.trim_ascii_end())
            .map_err(|_| invalid("header is not valid UTF-8"))?;
        self.header.push_str(header);
        match header.char_indices().find(|&(_, c)| separator(c)) {
            Some((end, c)) => {
                self.id_end = end;
                self.desc_start = Some(end + c.len_utf8());
            }
            None => self.id_end = header.len(),
        }
        Ok(())
    }
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

    /// the format of the input, from its next byte: `None` at the end of
    /// the input.
    pub fn format(&mut self) -> io::Result<Option<Format>> {
        match self.reader.fill_buf()?.first() {
            None => Ok(None),
            Some(b'>') => Ok(Some(Format::Fasta)),
            Some(b'@') => Ok(Some(Format::Fastq)),
            Some(_) => Err(invalid("record does not start with '>' or '@'")),
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
        match self.format()? {
            None => Ok(()),
            Some(Format::Fasta) => self.read_fasta(record),
            Some(Format::Fastq) => self.read_fastq(record),
        }
    }

    fn read_fasta(&mut self, record: &mut Record) -> io::Result<()> {
        self.next_line()?;
        record.set_header(&self.line[1..], char::is_whitespace)?;
        // sequence lines, up to the next record
        while !matches!(self.reader.fill_buf()?.first(), None | Some(b'>')) {
            self.next_line()?;
            record.seq.extend_from_slice(self.line.trim_ascii_end());
        }
        record.format = Some(Format::Fasta);
        Ok(())
    }

    fn read_fastq(&mut self, record: &mut Record) -> io::Result<()> {
        self.next_line()?;
        record.set_header(&self.line[1..], |c| c == ' ')?;

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
        // a FASTA id ends at any whitespace
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
