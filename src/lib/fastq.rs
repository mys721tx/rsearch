//! A FASTQ reader yielding `bio::io::fastq::Record`s.
//!
//! bio's own reader rejects reads with an empty sequence, which adapter
//! trimming commonly produces and vsearch accepts. This reader accepts
//! them, and otherwise parses records as bio does: the header is split at
//! the first space into an id and a description, and sequences and
//! qualities may span several lines.

use std::io::{self, BufRead};

use bio::io::fastq::Record;

/// Reads FASTQ records from a buffered reader.
///
/// # Example
/// ```rust
/// # use rsearch::fastq::Reader;
/// let input = b"@a desc\nACGT\n+\nIIII\n@empty\n\n+\n\n";
/// let records: Vec<_> = Reader::new(&input[..]).records().map(Result::unwrap).collect();
///
/// assert_eq!(records[0].id(), "a");
/// assert_eq!(records[0].desc(), Some("desc"));
/// assert_eq!(records[0].qual(), b"IIII");
/// assert_eq!(records[1].seq(), b"");
/// ```
pub struct Reader<R> {
    reader: R,
    line: String,
}

fn invalid(message: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.to_string())
}

impl<R: BufRead> Reader<R> {
    pub fn new(reader: R) -> Reader<R> {
        Reader {
            reader,
            line: String::new(),
        }
    }

    /// reads the next line without its line ending; false at end of input.
    fn next_line(&mut self) -> io::Result<bool> {
        self.line.clear();
        if self.reader.read_line(&mut self.line)? == 0 {
            return Ok(false);
        }
        let len = self.line.trim_end_matches(['\n', '\r']).len();
        self.line.truncate(len);
        Ok(true)
    }

    /// reads the next record, or `None` at end of input.
    pub fn read(&mut self) -> io::Result<Option<Record>> {
        if !self.next_line()? {
            return Ok(None);
        }
        let Some(header) = self.line.strip_prefix('@') else {
            return Err(invalid("FASTQ record does not start with '@'"));
        };
        let mut fields = header.trim_end().splitn(2, ' ');
        let id = fields.next().unwrap_or_default().to_owned();
        let desc = fields.next().map(str::to_owned);

        // sequence lines, up to the '+' line
        let mut seq = String::new();
        loop {
            if !self.next_line()? {
                return Err(invalid("incomplete FASTQ record"));
            }
            if self.line.starts_with('+') {
                break;
            }
            seq.push_str(self.line.trim_end());
        }

        // quality lines: at least one, until as long as the sequence
        let mut qual = String::new();
        loop {
            if !self.next_line()? {
                if qual.is_empty() && seq.is_empty() {
                    break;
                }
                return Err(invalid("incomplete FASTQ record"));
            }
            qual.push_str(self.line.trim_end());
            if qual.len() >= seq.len() {
                break;
            }
        }
        if qual.len() != seq.len() {
            return Err(invalid("FASTQ sequence and quality lengths differ"));
        }

        Ok(Some(Record::with_attrs(
            &id,
            desc.as_deref(),
            seq.as_bytes(),
            qual.as_bytes(),
        )))
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
        self.reader.read().transpose()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn read(input: &str) -> io::Result<Vec<Record>> {
        Reader::new(input.as_bytes()).records().collect()
    }

    #[test]
    fn test_records() {
        let records = read("@a x y\nAC\nGT\n+a x y\nII\nII\n@b\nA\n+\n@\n").unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].desc(), Some("x y"));
        assert_eq!(records[0].seq(), b"ACGT");
        assert_eq!(records[0].qual(), b"IIII");
        // a quality line may start with '@'
        assert_eq!(records[1].qual(), b"@");
        assert!(records.iter().all(|r| r.check().is_ok()));
    }

    #[test]
    fn test_empty_reads() {
        let records = read("@a\n\n+\n\n@b\nA\n+\nI\n@c\n\n+\n").unwrap();
        assert_eq!(records.len(), 3);
        assert_eq!(records[0].seq(), b"");
        assert_eq!(records[1].seq(), b"A");
        assert_eq!(records[2].qual(), b"");
        assert!(records[0].check().is_ok());
    }

    #[test]
    fn test_crlf() {
        let records = read("@a d\r\nAC\r\n+\r\nII\r\n").unwrap();
        assert_eq!(records[0].desc(), Some("d"));
        assert_eq!(records[0].seq(), b"AC");
    }

    #[test]
    fn test_errors() {
        assert!(read("a\nAC\n+\nII\n").is_err());
        assert!(read("@a\nAC\n").is_err());
        assert!(read("@a\nAC\n+\nI\n").is_err());
        assert!(read("@a\nAC\n+\nIII\n").is_err());
    }
}
