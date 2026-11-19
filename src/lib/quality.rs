//! Converting quality scores between encodings, as vsearch --fastq_convert
//! does: between ASCII offsets, from the Solexa scale to Phred, and into an
//! output range.

use crate::filter::QualityError;

/// The lowest and highest printable ASCII symbols.
const PRINTABLE: (i16, i16) = (33, 126);

/// returns the Phred score of a Solexa score.
///
/// A Solexa score is -10 log10(p / (1 - p)) where a Phred score is -10
/// log10(p). The Phred score is therefore 10 log10(10^(Q/10) + 1), rounded.
/// The two agree from 10 upwards.
///
/// # Example
/// ```rust
/// # use rsearch::quality::solexa_to_phred;
/// assert_eq!(solexa_to_phred(-5), 1);
/// assert_eq!(solexa_to_phred(0), 3);
/// assert_eq!(solexa_to_phred(10), 10);
/// ```
pub fn solexa_to_phred(score: i16) -> i16 {
    (10.0 * (10f64.powf(f64::from(score) / 10.0) + 1.0).log10()).round() as i16
}

/// The input and output encodings of a conversion.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Conversion {
    /// ASCII offset of the input scores.
    pub offset: u8,
    /// Lowest and highest input score accepted.
    pub min: i16,
    pub max: i16,
    /// Input scores are on the Solexa scale.
    pub solexa: bool,
    /// ASCII offset of the output scores.
    pub offset_out: u8,
    /// Lowest and highest score written.
    pub min_out: i16,
    pub max_out: i16,
}

impl Default for Conversion {
    /// Phred+33 in and out, with every score the encoding can carry.
    fn default() -> Conversion {
        Conversion {
            offset: 33,
            min: 0,
            max: 93,
            solexa: false,
            offset_out: 33,
            min_out: 0,
            max_out: 93,
        }
    }
}

impl Conversion {
    /// returns a converter with the output symbol of every input symbol
    /// computed once.
    pub fn converter(&self) -> Converter {
        let symbols = std::array::from_fn(|symbol| {
            let score = symbol as i16 - i16::from(self.offset);
            if score < self.min || score > self.max {
                return REJECTED;
            }
            let score = if self.solexa {
                solexa_to_phred(score)
            } else {
                score
            };
            let symbol = score.max(self.min_out).min(self.max_out) + i16::from(self.offset_out);
            symbol.clamp(PRINTABLE.0, PRINTABLE.1) as u8
        });
        Converter {
            conversion: *self,
            symbols,
        }
    }
}

/// Converts quality strings, symbol by symbol.
///
/// # Example
/// ```rust
/// # use rsearch::quality::Conversion;
/// // Phred+64 to Phred+33, with scores limited to 40
/// let converter = Conversion {
///     offset: 64,
///     max: 62,
///     max_out: 40,
///     ..Conversion::default()
/// }
/// .converter();
/// let mut out = Vec::new();
/// converter.convert(b"@Jhi", &mut out).unwrap();
/// assert_eq!(out, b"!+II");
///
/// // a score below the minimum is refused
/// assert!(converter.convert(b"?", &mut out).is_err());
/// ```
pub struct Converter {
    conversion: Conversion,
    /// The output symbol of every input symbol, or `REJECTED`.
    symbols: [u8; 256],
}

/// Marks an input symbol out of range. No output symbol is 0.
const REJECTED: u8 = 0;

impl Converter {
    /// writes the converted scores into `out`, replacing its contents.
    /// Fails on the first score outside the accepted range.
    pub fn convert(&self, qual: &[u8], out: &mut Vec<u8>) -> Result<(), QualityError> {
        out.clear();
        out.extend(qual.iter().map(|&q| self.symbols[usize::from(q)]));
        // looking for a rejected symbol only afterwards keeps both loops
        // free of branches
        if !out
            .iter()
            .fold(false, |rejected, &q| rejected | (q == REJECTED))
        {
            return Ok(());
        }
        let i = out.iter().position(|&q| q == REJECTED).unwrap();
        Err(QualityError {
            score: i16::from(qual[i]) - i16::from(self.conversion.offset),
            min_score: self.conversion.min,
            max_score: self.conversion.max,
        })
    }
}

/// The quality encodings vsearch tells apart by the lowest and highest
/// score symbols of a file.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Encoding {
    /// Phred+33, scores up to 40 ('I').
    Sanger,
    /// Phred+33, scores above 40.
    Illumina18,
    /// Solexa+64, scores below 0.
    Solexa,
    /// Phred+64, scores from 0.
    Illumina13,
    /// Phred+64, scores from 2 ('B').
    Illumina15,
}

impl Encoding {
    /// guesses the encoding from the lowest and highest score symbols, as
    /// vsearch's `--fastq_chars` does.
    ///
    /// # Example
    /// ```rust
    /// # use rsearch::quality::Encoding;
    /// assert_eq!(Encoding::guess(b'!', b'I'), Encoding::Sanger);
    /// assert_eq!(Encoding::guess(b'I', b'J'), Encoding::Illumina18);
    /// assert_eq!(Encoding::guess(b';', b'h'), Encoding::Solexa);
    /// assert_eq!(Encoding::guess(b'@', b'h'), Encoding::Illumina13);
    /// assert_eq!(Encoding::guess(b'B', b'h'), Encoding::Illumina15);
    /// ```
    pub fn guess(lowest: u8, highest: u8) -> Encoding {
        // a low lowest or a low highest rules out the +64 offset
        if lowest < b';' || highest < b'K' {
            if highest > b'I' {
                Encoding::Illumina18
            } else {
                Encoding::Sanger
            }
        } else if lowest < b'@' {
            Encoding::Solexa
        } else if lowest < b'B' {
            Encoding::Illumina13
        } else {
            Encoding::Illumina15
        }
    }

    /// the ASCII offset of the encoding.
    pub fn offset(self) -> u8 {
        match self {
            Encoding::Sanger | Encoding::Illumina18 => 33,
            Encoding::Solexa | Encoding::Illumina13 | Encoding::Illumina15 => 64,
        }
    }
}
