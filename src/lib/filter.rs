//! Trimming and filtering of FASTQ reads, after vsearch's `--fastq_filter`.

use std::fmt;
use std::sync::LazyLock;

/// ASCII offset of Phred+33 quality scores.
const PHRED_OFFSET: u8 = 33;

/// The error probability `10^(-Q/10)` of every quality symbol, for a
/// given ASCII offset.
fn probabilities(offset: u8) -> [f64; 256] {
    std::array::from_fn(|symbol| {
        let q = symbol as f64 - f64::from(offset);
        10f64.powf(-q / 10.0)
    })
}

static ERROR_PROBABILITY: LazyLock<[f64; 256]> = LazyLock::new(|| probabilities(PHRED_OFFSET));
static ERROR_PROBABILITY_64: LazyLock<[f64; 256]> = LazyLock::new(|| probabilities(64));

/// returns the expected number of errors in a quality string, the sum of
/// the error probabilities of its symbols.
///
/// # Example
/// ```rust
/// # use rsearch::filter::expected_errors;
/// assert!((expected_errors(b"++5") - 0.21).abs() < 1e-12);
/// ```
pub fn expected_errors(qual: &[u8]) -> f64 {
    let probability = &*ERROR_PROBABILITY;
    qual.iter()
        .fold(0.0, |ee, &q| ee + probability[usize::from(q)])
}

/// returns the percentage of residues of `seq` that `counted` accepts,
/// computed as vsearch's `--fastx_mask` does: `NaN` for an empty sequence,
/// which no limit then discards.
fn percent(seq: &[u8], counted: impl Fn(u8) -> bool) -> f64 {
    let count = seq.iter().filter(|&&c| counted(c)).count();
    100.0 * count as f64 / seq.len() as f64
}

/// returns the percentage of upper-case residues, those a soft mask left
/// unmasked.
///
/// # Example
/// ```rust
/// # use rsearch::filter::upper_pct;
/// assert_eq!(upper_pct(b"ACGTacgt"), 50.0);
/// assert_eq!(upper_pct(b"NnACGTAC"), 87.5);
/// ```
pub fn upper_pct(seq: &[u8]) -> f64 {
    percent(seq, |c| c.is_ascii_uppercase())
}

/// returns the percentage of residues other than N (or n), those a hard
/// mask left unmasked.
///
/// # Example
/// ```rust
/// # use rsearch::filter::non_n_pct;
/// assert_eq!(non_n_pct(b"ACGTNNNN"), 50.0);
/// assert_eq!(non_n_pct(b"NnACGTAC"), 75.0);
/// ```
pub fn non_n_pct(seq: &[u8]) -> f64 {
    percent(seq, |c| c != b'N' && c != b'n')
}

/// Trimming and filtering settings. Trimming happens first. The ends are
/// stripped, then the read is truncated by length, quality and expected
/// errors. The filters then apply to what is left.
#[derive(Clone, Debug)]
pub struct Filter {
    /// Bases removed from the 5' end.
    pub strip_left: usize,
    /// Bases removed from the 3' end.
    pub strip_right: usize,
    /// Truncate to this length, keeping shorter reads.
    pub trunc_keep: Option<usize>,
    /// Truncate at the first base with a quality score at or below this.
    pub trunc_qual: Option<i16>,
    /// Truncate before the expected errors exceed this.
    pub trunc_ee: Option<f64>,
    /// Truncate before the expected errors per base exceed this.
    pub trunc_rate: Option<f64>,
    /// Discard reads shorter than this.
    pub min_len: usize,
    /// Discard reads longer than this.
    pub max_len: Option<usize>,
    /// Discard reads with more N bases than this.
    pub max_ns: Option<usize>,
    /// Discard reads with more expected errors than this.
    pub max_ee: Option<f64>,
    /// Discard reads with more expected errors per base than this.
    pub max_rate: Option<f64>,
    /// Discard reads with any quality score below this.
    pub min_qual: i16,
    /// Discard reads with an abundance below this.
    pub min_size: u64,
    /// Discard reads with an abundance above this.
    pub max_size: Option<u64>,
    /// ASCII offset of quality symbols: 33 (Phred+33) or 64 (Phred+64);
    /// any other value is read as 64.
    pub offset: u8,
    /// Lowest quality score accepted in the input.
    pub min_score: i16,
    /// Highest quality score accepted in the input.
    pub max_score: i16,
}

/// A quality score outside the accepted range.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct QualityError {
    pub score: i16,
    pub min_score: i16,
    pub max_score: i16,
}

impl fmt::Display for QualityError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.score < self.min_score {
            write!(
                f,
                "quality score {} below the minimum {}",
                self.score, self.min_score
            )
        } else {
            write!(
                f,
                "quality score {} above the maximum {}",
                self.score, self.max_score
            )
        }
    }
}

impl std::error::Error for QualityError {}

impl Default for Filter {
    /// No trimming; only empty reads are discarded.
    fn default() -> Filter {
        Filter {
            strip_left: 0,
            strip_right: 0,
            trunc_keep: None,
            trunc_qual: None,
            trunc_ee: None,
            trunc_rate: None,
            min_len: 1,
            max_len: None,
            max_ns: None,
            max_ee: None,
            max_rate: None,
            min_qual: 0,
            min_size: 1,
            max_size: None,
            offset: PHRED_OFFSET,
            min_score: 0,
            max_score: 93,
        }
    }
}

/// The result of trimming and filtering a read.
#[derive(Clone, Debug, PartialEq)]
pub struct Outcome {
    /// The kept bases are `start..start + len`.
    pub start: usize,
    pub len: usize,
    /// Expected errors of the kept bases, `None` for a FASTA read.
    pub ee: Option<f64>,
    /// Whether bases were removed.
    pub truncated: bool,
    /// Whether the read fails a filter.
    pub discarded: bool,
}

impl Outcome {
    /// The kept part of a sequence or quality string.
    pub fn trim<'a>(&self, s: &'a [u8]) -> &'a [u8] {
        &s[self.start..self.start + self.len]
    }
}

impl Filter {
    /// whether the filter depends on abundance, so that callers can skip
    /// finding it.
    pub fn uses_size(&self) -> bool {
        self.min_size > 1 || self.max_size.is_some()
    }

    /// the score of a quality symbol, or an error outside the accepted range.
    fn score(&self, symbol: u8) -> Result<i16, QualityError> {
        let score = i16::from(symbol) - i16::from(self.offset);
        if score < self.min_score || score > self.max_score {
            return Err(QualityError {
                score,
                min_score: self.min_score,
                max_score: self.max_score,
            });
        }
        Ok(score)
    }

    /// trims and filters a read of abundance `size`. `qual` is `None` for a
    /// FASTA read, which only the length, N and abundance filters apply to.
    /// As in vsearch, the scores of the bases left after trimming must be in
    /// the accepted range.
    ///
    /// # Example
    /// ```rust
    /// # use rsearch::filter::Filter;
    /// let filter = Filter {
    ///     trunc_qual: Some(10),
    ///     min_len: 3,
    ///     ..Filter::default()
    /// };
    ///
    /// // truncated before the first base with a score of 10 or less ('+')
    /// let outcome = filter.apply(b"ACGTA", Some(b"IIII+"), 1).unwrap();
    /// assert_eq!(outcome.trim(b"ACGTA"), b"ACGT");
    /// assert!(outcome.truncated && !outcome.discarded);
    ///
    /// let outcome = filter.apply(b"ACGTA", Some(b"II+II"), 1).unwrap();
    /// assert_eq!(outcome.len, 2);
    /// assert!(outcome.discarded);
    ///
    /// // '~' is score 93, above an accepted maximum of 41
    /// let strict = Filter { max_score: 41, ..Filter::default() };
    /// let error = strict.apply(b"A", Some(b"~"), 1).unwrap_err();
    /// assert_eq!((error.score, error.max_score), (93, 41));
    /// ```
    pub fn apply(
        &self,
        seq: &[u8],
        qual: Option<&[u8]>,
        size: u64,
    ) -> Result<Outcome, QualityError> {
        let seqlen = seq.len();
        let start = self.strip_left.min(seqlen);
        let mut len = (seqlen - start).saturating_sub(self.strip_right);
        if let Some(n) = self.trunc_keep {
            len = len.min(n);
        }

        let mut discarded = false;
        let mut ee = 0.0;
        let trunc_ee = self.trunc_ee.unwrap_or(f64::MAX);
        let trunc_rate = self.trunc_rate.unwrap_or(f64::MAX);
        let per_base = self.trunc_qual.is_some()
            || self.trunc_ee.is_some()
            || self.trunc_rate.is_some()
            || self.min_qual > self.min_score;
        let window = qual.map(|qual| &qual[start..start + len]);
        // indexed by symbol, which spares lookups their bounds checks
        let probability = match self.offset {
            PHRED_OFFSET => &*ERROR_PROBABILITY,
            _ => &*ERROR_PROBABILITY_64,
        };
        if let (true, Some(window)) = (per_base, window) {
            let trunc_qual = self.trunc_qual.unwrap_or(i16::MIN);
            for (i, &symbol) in window.iter().enumerate() {
                let score = self.score(symbol)?;
                let error = probability[usize::from(symbol)];
                ee += error;
                // as vsearch: add, then take back the base that crosses a limit
                if score <= trunc_qual || ee > trunc_ee || ee > trunc_rate * (i + 1) as f64 {
                    ee -= error;
                    len = i;
                    break;
                }
                if score < self.min_qual {
                    discarded = true;
                }
            }
        } else if let Some(window) = window {
            // nothing to check base by base: add up the errors, tracking the
            // lowest and highest symbols to check the range in the same pass
            let (mut low, mut high) = (u8::MAX, u8::MIN);
            for &q in window {
                ee += probability[usize::from(q)];
                low = low.min(q);
                high = high.max(q);
            }
            let lowest = i16::from(self.offset) + self.min_score;
            let highest = i16::from(self.offset) + self.max_score;
            if !window.is_empty() && (i16::from(low) < lowest || i16::from(high) > highest) {
                // report the first symbol out of range, as vsearch does
                let first = window
                    .iter()
                    .find(|&&q| !(lowest..=highest).contains(&i16::from(q)));
                if let Some(&symbol) = first {
                    self.score(symbol)?;
                }
            }
        }

        let kept = &seq[start..start + len];
        let ee = qual.map(|_| ee);
        discarded |= ee.is_some_and(|ee| self.max_ee.is_some_and(|max| ee > max))
            || ee.is_some_and(|ee| {
                len > 0 && self.max_rate.is_some_and(|max| ee / len as f64 > max)
            })
            || len < self.min_len
            || self.max_len.is_some_and(|max| len > max)
            || self
                .max_ns
                .is_some_and(|max| kept.iter().filter(|&&c| c == b'N' || c == b'n').count() > max)
            || size < self.min_size
            || self.max_size.is_some_and(|max| size > max);

        Ok(Outcome {
            start,
            len,
            ee,
            truncated: len < seqlen,
            discarded,
        })
    }
}
