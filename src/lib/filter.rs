//! Trimming and filtering of FASTQ reads, after vsearch's `--fastq_filter`.

use std::sync::LazyLock;

/// ASCII offset of Phred quality scores.
const PHRED_OFFSET: u8 = 33;

/// The error probability of every quality symbol, `10^(-Q/10)`.
static ERROR_PROBABILITY: LazyLock<[f64; 256]> = LazyLock::new(|| {
    std::array::from_fn(|symbol| {
        let q = symbol as f64 - f64::from(PHRED_OFFSET);
        10f64.powf(-q / 10.0)
    })
});

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

/// Trimming and filtering settings. Trimming happens first. The ends are
/// stripped, then the read is truncated by length, quality and expected
/// errors. The filters then apply to what is left.
#[derive(Clone, Debug)]
pub struct Filter {
    /// Bases removed from the 5' end.
    pub strip_left: usize,
    /// Bases removed from the 3' end.
    pub strip_right: usize,
    /// Truncate to this length, and discard shorter reads.
    pub trunc_len: Option<usize>,
    /// Truncate to this length, keeping shorter reads.
    pub trunc_keep: Option<usize>,
    /// Truncate at the first base with a quality score at or below this.
    pub trunc_qual: Option<u8>,
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
    pub min_qual: u8,
    /// Discard reads with an abundance below this.
    pub min_size: u64,
    /// Discard reads with an abundance above this.
    pub max_size: Option<u64>,
}

impl Default for Filter {
    /// No trimming; only empty reads are discarded.
    fn default() -> Filter {
        Filter {
            strip_left: 0,
            strip_right: 0,
            trunc_len: None,
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

    /// trims and filters a read of abundance `size`. `qual` is `None` for a
    /// FASTA read, which only the length, N and abundance filters apply to.
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
    /// let outcome = filter.apply(b"ACGTA", Some(b"IIII+"), 1);
    /// assert_eq!(outcome.trim(b"ACGTA"), b"ACGT");
    /// assert!(outcome.truncated && !outcome.discarded);
    ///
    /// let outcome = filter.apply(b"ACGTA", Some(b"II+II"), 1);
    /// assert_eq!(outcome.len, 2);
    /// assert!(outcome.discarded);
    /// ```
    pub fn apply(&self, seq: &[u8], qual: Option<&[u8]>, size: u64) -> Outcome {
        let seqlen = seq.len();
        let start = self.strip_left.min(seqlen);
        let mut len = (seqlen - start).saturating_sub(self.strip_right);
        if let Some(n) = self.trunc_len {
            len = len.min(n);
        }
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
            || self.min_qual > 0;
        let window = qual.map(|qual| &qual[start..start + len]);
        if let (true, Some(window)) = (per_base, window) {
            let probability = &*ERROR_PROBABILITY;
            // scores are at least 0, so -1 never truncates
            let trunc_qual = self.trunc_qual.map_or(-1, i16::from);
            for (i, &symbol) in window.iter().enumerate() {
                let score = symbol.saturating_sub(PHRED_OFFSET);
                let error = probability[usize::from(symbol)];
                ee += error;
                // as vsearch: add, then take back the base that crosses a limit
                if i16::from(score) <= trunc_qual
                    || ee > trunc_ee
                    || ee > trunc_rate * (i + 1) as f64
                {
                    ee -= error;
                    len = i;
                    break;
                }
                if score < self.min_qual {
                    discarded = true;
                }
            }
        } else if let Some(window) = window {
            // nothing to check base by base: just add up the errors
            ee = expected_errors(window);
        }

        let kept = &seq[start..start + len];
        let ee = qual.map(|_| ee);
        discarded |= ee.is_some_and(|ee| self.max_ee.is_some_and(|max| ee > max))
            || ee.is_some_and(|ee| {
                len > 0 && self.max_rate.is_some_and(|max| ee / len as f64 > max)
            })
            || self.trunc_len.is_some_and(|n| len < n)
            || len < self.min_len
            || self.max_len.is_some_and(|max| len > max)
            || self
                .max_ns
                .is_some_and(|max| kept.iter().filter(|&&c| c == b'N' || c == b'n').count() > max)
            || size < self.min_size
            || self.max_size.is_some_and(|max| size > max);

        Outcome {
            start,
            len,
            ee,
            truncated: len < seqlen,
            discarded,
        }
    }
}
