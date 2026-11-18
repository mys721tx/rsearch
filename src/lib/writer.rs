//! Writing FASTA and FASTQ records.

use std::io::{self, Write};

use bio::io::{fasta, fastq};

/// Writes records as FASTA or FASTQ, with bio's writers.
///
/// # Example
/// ```rust
/// # use rsearch::writer::Writer;
/// let mut out = Vec::new();
/// let mut writer = Writer::new(&mut out, true, 2);
/// writer.write("a", None, b"ACGT", b"IIII").unwrap();
/// writer.flush().unwrap();
/// drop(writer);
/// assert_eq!(out, b">a\nAC\nGT\n");
/// ```
pub enum Writer<W: Write> {
    Fasta(fasta::Writer<W>),
    Fastq(fastq::Writer<W>),
}

impl<W: Write> Writer<W> {
    /// writes FASTA if `as_fasta`, FASTQ otherwise. `width` wraps FASTA
    /// sequence lines; 0 writes each on one line.
    pub fn new(out: W, as_fasta: bool, width: usize) -> Writer<W> {
        if as_fasta {
            let mut writer = fasta::Writer::new(out);
            writer.set_linewrap((width > 0).then_some(width));
            Writer::Fasta(writer)
        } else {
            Writer::Fastq(fastq::Writer::new(out))
        }
    }

    /// writes a record; `qual` is ignored when writing FASTA.
    pub fn write(
        &mut self,
        id: &str,
        desc: Option<&str>,
        seq: &[u8],
        qual: &[u8],
    ) -> io::Result<()> {
        match self {
            Writer::Fasta(w) => w.write(id, desc, seq),
            Writer::Fastq(w) => w.write(id, desc, seq, qual),
        }
    }

    pub fn flush(&mut self) -> io::Result<()> {
        match self {
            Writer::Fasta(w) => w.flush(),
            Writer::Fastq(w) => w.flush(),
        }
    }
}
