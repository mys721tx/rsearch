//! UDB files, as vsearch writes them with `--makeudb_usearch`.
//!
//! Every number is a little-endian 32-bit word, but the counts of residues
//! and of header bytes, which take two (low then high). In order:
//!
//! - a header of 50 words: the magic `UDBF`, 32 (bits), the word length,
//!   1 (the step between database words), 100 (the share of words indexed),
//!   the number of sequences, the alphabet `nt`, and the magic `UDBf`;
//! - for each of the 4^k words, the number of sequences holding it;
//! - the magic `UDB3`, then the numbers of the sequences holding each word,
//!   in word order, each list increasing;
//! - a second header of 8 words: the magic `UDB4`, 0x005e0db3, the number of
//!   sequences, the number of residues, the number of header bytes with
//!   their terminating zeros, and 0x005e0db4;
//! - the offset of each header, the headers each ended by a zero, the
//!   length of each sequence, and the sequences, unpadded.

use std::io::{self, Write};

use crate::words::Index;

const MAGIC_HEADER: u32 = 0x5544_4246;
const MAGIC_HEADER_END: u32 = 0x5544_4266;
const MAGIC_WORDS: u32 = 0x5544_4233;
const MAGIC_SEQUENCES: u32 = 0x5544_4234;
const MAGIC_SEQUENCES_2: u32 = 0x005e_0db3;
const MAGIC_SEQUENCES_3: u32 = 0x005e_0db4;
/// `nt`, the nucleotide alphabet.
const ALPHABET: u32 = 0x0000_746e;
/// The words in the first header.
const HEADER_WORDS: usize = 50;

fn put(out: &mut impl Write, values: &[u32]) -> io::Result<()> {
    let mut bytes = Vec::with_capacity(4 * values.len());
    for value in values {
        bytes.extend_from_slice(&value.to_le_bytes());
    }
    out.write_all(&bytes)
}

/// writes the UDB file of the sequences `records`, (header, sequence)
/// pairs, whose words are indexed in `index`.
///
/// # Example
/// ```rust
/// # use rsearch::{udb::write, words::Index};
/// let records: [(&[u8], &[u8]); 1] = [(b"a", b"ACGT")];
/// let index = Index::new(3, records.iter().map(|r| r.1));
/// let mut out = Vec::new();
/// write(&mut out, &index, &records).unwrap();
/// // the headers, the word counts, two lists, the second header, one
/// // offset, "a\0", one length and four residues
/// assert_eq!(out.len(), 4 * (50 + 64 + 1 + 2 + 8 + 1) + 2 + 4 + 4);
/// assert_eq!(&out[..4], b"FBDU");
/// ```
pub fn write<H: AsRef<[u8]>, S: AsRef<[u8]>>(
    out: &mut impl Write,
    index: &Index,
    records: &[(H, S)],
) -> io::Result<()> {
    let count = records.len() as u32;
    let mut header = [0u32; HEADER_WORDS];
    header[0] = MAGIC_HEADER;
    header[2] = 32;
    header[4] = index.word_len();
    header[5] = 1;
    header[6] = 100;
    header[13] = count;
    header[17] = ALPHABET;
    header[49] = MAGIC_HEADER_END;
    put(out, &header)?;
    put(out, &index.counts())?;
    put(out, &[MAGIC_WORDS])?;
    put(out, index.all_sequences())?;

    let residues: u64 = records.iter().map(|r| r.1.as_ref().len() as u64).sum();
    let header_bytes: u64 = records.iter().map(|r| r.0.as_ref().len() as u64 + 1).sum();
    put(
        out,
        &[
            MAGIC_SEQUENCES,
            MAGIC_SEQUENCES_2,
            count,
            residues as u32,
            (residues >> 32) as u32,
            header_bytes as u32,
            (header_bytes >> 32) as u32,
            MAGIC_SEQUENCES_3,
        ],
    )?;
    let mut offset = 0u32;
    let offsets: Vec<u32> = records
        .iter()
        .map(|r| {
            let this = offset;
            offset = offset.wrapping_add(r.0.as_ref().len() as u32 + 1);
            this
        })
        .collect();
    put(out, &offsets)?;
    for (h, _) in records {
        out.write_all(h.as_ref())?;
        out.write_all(&[0])?;
    }
    let lengths: Vec<u32> = records.iter().map(|r| r.1.as_ref().len() as u32).collect();
    put(out, &lengths)?;
    for (_, seq) in records {
        out.write_all(seq.as_ref())?;
    }
    Ok(())
}

/// Why a UDB file was refused.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum UdbError {
    /// The file is shorter than the first header.
    Short,
    /// A field of the first header is wrong, or the file is too short for
    /// its number of sequences.
    Invalid,
}

impl std::fmt::Display for UdbError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            UdbError::Short => write!(f, "invalid UDB file: shorter than its header"),
            UdbError::Invalid => write!(f, "invalid UDB file"),
        }
    }
}

/// The first header of a UDB file, checked as vsearch does.
#[derive(Clone, Debug)]
pub struct Header {
    words: [u32; HEADER_WORDS],
}

impl Header {
    /// reads the header at the start of `bytes`, of a file of `file_size`
    /// bytes. Like vsearch, it checks the magic numbers, the bits, the word
    /// length (3 to 15), the alphabet, that there are sequences, and that
    /// the file holds 4 bytes for each.
    ///
    /// # Example
    /// ```rust
    /// # use rsearch::udb::{Header, UdbError};
    /// let mut words = [0u32; 50];
    /// (words[0], words[2], words[4], words[13]) = (0x5544_4246, 32, 8, 1);
    /// (words[17], words[49]) = (0x746e, 0x5544_4266);
    /// let bytes: Vec<u8> = words.iter().flat_map(|w| w.to_le_bytes()).collect();
    /// let header = Header::parse(&bytes, 200).unwrap();
    /// assert_eq!((header.sequences(), header.word_len()), (1, 8));
    /// assert_eq!(Header::parse(&bytes[..199], 199).unwrap_err(), UdbError::Short);
    /// ```
    pub fn parse(bytes: &[u8], file_size: u64) -> Result<Header, UdbError> {
        if bytes.len() < 4 * HEADER_WORDS {
            return Err(UdbError::Short);
        }
        let mut words = [0u32; HEADER_WORDS];
        for (i, word) in words.iter_mut().enumerate() {
            *word = u32::from_le_bytes(bytes[4 * i..4 * i + 4].try_into().unwrap());
        }
        let header = Header { words };
        if words[0] != MAGIC_HEADER
            || words[2] != 32
            || !(3..=15).contains(&words[4])
            || words[13] == 0
            || words[17] != ALPHABET
            || words[49] != MAGIC_HEADER_END
            || u64::from(words[13]) > file_size / 4
        {
            return Err(UdbError::Invalid);
        }
        Ok(header)
    }

    /// the number of bits of a sequence number: always 32.
    pub fn bits(&self) -> u32 {
        self.words[2]
    }

    pub fn word_len(&self) -> u32 {
        self.words[4]
    }

    /// the step between the database positions indexed.
    pub fn db_step(&self) -> u32 {
        self.words[5]
    }

    /// the share of database words indexed, in percent.
    pub fn db_accel(&self) -> u32 {
        self.words[6]
    }

    /// the number of slots, 0 when the index holds every word.
    pub fn slots(&self) -> u32 {
        self.words[11]
    }

    pub fn sequences(&self) -> u32 {
        self.words[13]
    }

    /// vsearch's `--udbinfo` report of the header.
    ///
    /// # Example
    /// ```rust
    /// # use rsearch::udb::Header;
    /// # let mut words = [0u32; 50];
    /// # (words[0], words[2], words[4], words[13]) = (0x5544_4246, 32, 8, 1);
    /// # (words[17], words[49]) = (0x746e, 0x5544_4266);
    /// # let bytes: Vec<u8> = words.iter().flat_map(|w| w.to_le_bytes()).collect();
    /// let report = Header::parse(&bytes, 200).unwrap().report();
    /// assert!(report.contains("      Dict size  65536 (65.5k)\n"));
    /// ```
    pub fn report(&self) -> String {
        let dict = 1u64 << (2 * self.word_len());
        format!(
            "           Seqs  {}\n     SeqIx bits  {}\n          Alpha  nt (4)\n     \
             Word width  {}\n          Slots  {}\n      Dict size  {dict} ({:.1}k)\n         \
             DBstep  {}\n        DBAccel  {}%\n",
            self.sequences(),
            self.bits(),
            self.word_len(),
            self.slots(),
            dict as f64 / 1000.0,
            self.db_step(),
            self.db_accel(),
        )
    }
}
