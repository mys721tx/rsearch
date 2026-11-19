//! Merging of overlapping read pairs, as vsearch's `--fastq_mergepairs`.
//!
//! The overlap is found by 5-mers the forward read shares with the reverse
//! complement of its mate. Each diagonal with enough of them is scored by
//! ungapped alignment, with match and mismatch scores in bits from the
//! error probabilities of both bases. A pair merges at its best diagonal
//! when it is the only one scoring high enough. The merged bases get the
//! posterior scores of Edgar & Flyvbjerg (2015).

use crate::cluster::complement;
use crate::filter::QualityError;

/// The length of the words that seed the overlap.
const K: usize = 5;
/// The largest drop in score from a high, in bits, an alignment may have.
const DROP_MAX: f64 = 16.0;
/// The highest score of a mismatch, in bits.
const MISMATCH_MAX: f64 = -4.0;
/// The longest merged read, vsearch's default `--fastq_maxmergelen`.
const MAX_MERGE_LENGTH: usize = 1_000_000;

/// How pairs are merged, and the scores read and written.
#[derive(Clone, Copy, Debug)]
pub struct Settings {
    /// The shortest overlap, at least 5.
    pub min_overlap: usize,
    /// The most differences in the overlap.
    pub max_diffs: usize,
    /// The highest percentage of differences in the overlap.
    pub max_diff_pct: f64,
    /// Whether the reverse read may run past the start of the forward read.
    pub stagger: bool,
    /// The ASCII offset of scores, read and written.
    pub offset: u8,
    /// The range of scores accepted in the input.
    pub min_score: i16,
    pub max_score: i16,
    /// The range of merged scores.
    pub min_score_out: i16,
    pub max_score_out: i16,
}

/// Why a pair was not merged, in the order vsearch reports them.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Reason {
    /// A read is empty.
    TooShort,
    /// No diagonal has enough shared 5-mers.
    NoKmers,
    /// More than one diagonal scores high enough.
    Repeat,
    MaxDiffs,
    MaxDiffPct,
    /// The best alignment scores too low, or drops too far.
    MinScore,
    MinOverlap,
    /// The merged read would be longer than a million bases.
    MaxMergeLength,
    Staggered,
}

impl Reason {
    /// Every reason, in the order of vsearch's report.
    pub const ALL: [Reason; 9] = [
        Reason::TooShort,
        Reason::NoKmers,
        Reason::Repeat,
        Reason::MaxDiffs,
        Reason::MaxDiffPct,
        Reason::MinScore,
        Reason::MinOverlap,
        Reason::MaxMergeLength,
        Reason::Staggered,
    ];

    /// vsearch's description of the reason in its report.
    pub fn description(self) -> &'static str {
        match self {
            Reason::TooShort => "reads too short (after truncation)",
            Reason::NoKmers => "too few kmers found on same diagonal",
            Reason::Repeat => "multiple potential alignments",
            Reason::MaxDiffs => "too many differences",
            Reason::MaxDiffPct => "too high percentage of differences",
            Reason::MinScore => "alignment score too low, or score drop too high",
            Reason::MinOverlap => "overlap too short",
            Reason::MaxMergeLength => "merged fragment too long",
            Reason::Staggered => "staggered read pairs",
        }
    }
}

/// A merged pair: the read and its errors.
#[derive(Clone, Debug, Default)]
pub struct Merged {
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
    /// The expected errors of the merged read, and of the overlap in each
    /// read.
    pub ee_merged: f64,
    pub ee_forward: f64,
    pub ee_reverse: f64,
    /// The bases of the overlap where each read differs from the merged
    /// read.
    pub forward_errors: usize,
    pub reverse_errors: usize,
}

/// Scores per pair of quality symbols, indexed by symbol.
struct Tables {
    /// The error probability of each symbol.
    error: [f64; 128],
    /// The merged symbol when two bases agree, and when they differ (the
    /// higher symbol first).
    same: Vec<u8>,
    diff: Vec<u8>,
    /// The alignment score of a match and of a mismatch, in bits.
    match_score: Vec<f64>,
    mismatch_score: Vec<f64>,
}

/// The error probability of a symbol: 10^(-Q/10), and 0.75 below Q2.
fn probability(symbol: u8, offset: u8) -> f64 {
    let score = i32::from(symbol) - i32::from(offset);
    if score < 2 {
        0.75
    } else {
        10f64.powf(-f64::from(score) / 10.0)
    }
}

/// A probability as a symbol, bounded as C++'s std::min and std::max.
fn symbol(p: f64, settings: &Settings) -> u8 {
    let mut q = (-10.0 * p.log10()).round();
    let (low, high) = (
        f64::from(settings.min_score_out),
        f64::from(settings.max_score_out),
    );
    if high < q {
        q = high;
    }
    if q < low {
        q = low;
    }
    (f64::from(settings.offset) + q) as u8
}

impl Tables {
    fn new(settings: &Settings) -> Tables {
        let mut tables = Tables {
            error: [0.0; 128],
            same: vec![0; 128 * 128],
            diff: vec![0; 128 * 128],
            match_score: vec![0.0; 128 * 128],
            mismatch_score: vec![0.0; 128 * 128],
        };
        for x in 33..=126u8 {
            let px = probability(x, settings.offset);
            tables.error[usize::from(x)] = px;
            for y in 33..=126u8 {
                let py = probability(y, settings.offset);
                let i = usize::from(x) * 128 + usize::from(y);
                let p = px * py / 3.0 / (1.0 - px - py + (4.0 * px * py / 3.0));
                tables.same[i] = symbol(p, settings);
                let p = px * (1.0 - (py / 3.0)) / (px + py - (4.0 * px * py / 3.0));
                tables.diff[i] = symbol(p, settings);
                // the chance of observing a match of two bases that agree
                let p = 1.0 - px - py + (px * py * 4.0 / 3.0);
                tables.match_score[i] = (p / 0.25).log2();
                tables.mismatch_score[i] = ((1.0 - p) / 0.75).log2().min(MISMATCH_MAX);
            }
        }
        tables
    }
}

/// The 2-bit code of a base: A 0, C 1, G 2, T and U 3, others 0.
fn two_bit(base: u8) -> u32 {
    match base {
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' | b'U' | b'u' => 3,
        _ => 0,
    }
}

/// Whether a base is ambiguous, so no 5-mer covers it.
fn ambiguous(base: u8) -> u32 {
    u32::from(!matches!(
        base,
        b'A' | b'C' | b'G' | b'T' | b'U' | b'a' | b'c' | b'g' | b't' | b'u'
    ))
}

/// Merges pairs, reusing its buffers from pair to pair.
pub struct Merger {
    settings: Settings,
    tables: Tables,
    /// The pair being merged, upper-cased, N scored 0.
    forward: (Vec<u8>, Vec<u8>),
    reverse: (Vec<u8>, Vec<u8>),
    /// The forward positions of each 5-mer, chained and biased by one.
    chain_head: Vec<u32>,
    chain_next: Vec<u32>,
    /// The shared 5-mers on each diagonal.
    diagonals: Vec<u32>,
    merged: Merged,
}

impl Merger {
    pub fn new(settings: Settings) -> Merger {
        Merger {
            tables: Tables::new(&settings),
            settings,
            forward: (Vec::new(), Vec::new()),
            reverse: (Vec::new(), Vec::new()),
            chain_head: vec![0; 1 << (2 * K)],
            chain_next: Vec::new(),
            diagonals: Vec::new(),
            merged: Merged::default(),
        }
    }

    /// The forward read of the last pair, as vsearch writes it unmerged:
    /// upper-cased, and with N scored 0 when the pair had no empty read.
    pub fn forward(&self) -> (&[u8], &[u8]) {
        (&self.forward.0, &self.forward.1)
    }

    /// The reverse read of the last pair, as [`Merger::forward`].
    pub fn reverse(&self) -> (&[u8], &[u8]) {
        (&self.reverse.0, &self.reverse.1)
    }

    /// The last pair merged.
    pub fn merged(&self) -> &Merged {
        &self.merged
    }

    /// merges a pair, or says why it cannot be; fails at the first score
    /// out of range, in the forward read first.
    pub fn merge(
        &mut self,
        forward: (&[u8], &[u8]),
        reverse: (&[u8], &[u8]),
    ) -> Result<Option<Reason>, QualityError> {
        for ((seq, qual), (from_seq, from_qual)) in
            [(&mut self.forward, forward), (&mut self.reverse, reverse)]
        {
            seq.clear();
            seq.extend(from_seq.iter().map(u8::to_ascii_uppercase));
            qual.clear();
            qual.extend_from_slice(from_qual);
        }
        if self.forward.0.is_empty() || self.reverse.0.is_empty() {
            return Ok(Some(Reason::TooShort));
        }
        self.check(true)?;
        self.check(false)?;
        let offset = self.settings.offset;
        for (seq, qual) in [&mut self.forward, &mut self.reverse] {
            for (base, symbol) in seq.iter().zip(qual.iter_mut()) {
                if *base == b'N' {
                    *symbol = offset;
                }
            }
        }
        match self.overlap() {
            Ok(overlap) => {
                self.join(overlap);
                Ok(None)
            }
            Err(reason) => Ok(Some(reason)),
        }
    }

    /// fails at the first score out of range in a read.
    fn check(&self, forward: bool) -> Result<(), QualityError> {
        let qual = if forward {
            &self.forward.1
        } else {
            &self.reverse.1
        };
        let offset = i16::from(self.settings.offset);
        let (min, max) = (self.settings.min_score, self.settings.max_score);
        let (low, high) = qual
            .iter()
            .fold((u8::MAX, 0), |(low, high), &q| (low.min(q), high.max(q)));
        if i16::from(low) - offset >= min && i16::from(high) - offset <= max {
            return Ok(());
        }
        let score = qual
            .iter()
            .map(|&q| i16::from(q) - offset)
            .find(|score| !(min..=max).contains(score))
            .unwrap();
        Err(QualityError {
            score,
            min_score: min,
            max_score: max,
        })
    }

    /// counts the 5-mers shared on each diagonal.
    fn count_diagonals(&mut self) {
        let (fwd, rev) = (&self.forward.0, &self.reverse.0);
        let mask = (1u32 << (2 * K)) - 1;
        self.chain_head.fill(0);
        if self.chain_next.len() <= fwd.len() {
            self.chain_next.resize(fwd.len() + 1, 0);
        }
        let (mut bad, mut kmer) = (mask, 0u32);
        for (pos, &base) in fwd.iter().enumerate() {
            bad = ((bad << 2) | ambiguous(base)) & mask;
            kmer = ((kmer << 2) | two_bit(base)) & mask;
            if bad == 0 {
                // the start of the 5-mer, biased by one
                let entry = (pos + 2 - K) as u32;
                self.chain_next[entry as usize] = self.chain_head[kmer as usize];
                self.chain_head[kmer as usize] = entry;
            }
        }

        self.diagonals.clear();
        self.diagonals.resize(fwd.len() + rev.len(), 0);
        let length = rev.len() as i64;
        let (mut bad, mut kmer) = (mask, 0u32);
        for (pos, &base) in rev.iter().rev().enumerate() {
            bad = ((bad << 2) | ambiguous(base)) & mask;
            kmer = ((kmer << 2) | two_bit(complement(base))) & mask;
            if bad == 0 {
                let base_diagonal = length - (pos as i64 + 1 - K as i64);
                let mut entry = self.chain_head[kmer as usize];
                while entry != 0 {
                    let diagonal = base_diagonal + i64::from(entry) - 1;
                    if diagonal >= 0 {
                        self.diagonals[diagonal as usize] += 1;
                    }
                    entry = self.chain_next[entry as usize];
                }
            }
        }
    }

    /// The best overlap, scoring each diagonal with enough shared 5-mers.
    fn overlap(&mut self) -> Result<usize, Reason> {
        let min_overlap = self.settings.min_overlap;
        let (min_count, min_score) = if min_overlap < 9 {
            ((min_overlap - 4) as u32, 1.6 * min_overlap as f64)
        } else {
            (4, 16.0)
        };
        self.count_diagonals();
        let (fwd, fwd_qual) = (&self.forward.0, &self.forward.1);
        let (rev, rev_qual) = (&self.reverse.0, &self.reverse.1);
        let (fwd_len, rev_len) = (fwd.len(), rev.len());

        let (mut best_score, mut best_i, mut best_diffs) = (0.0, 0, 0);
        let (mut hits, mut kmers) = (0, false);
        for i in 1..fwd_len + rev_len {
            if self.diagonals[fwd_len + rev_len - i] < min_count {
                continue;
            }
            kmers = true;
            let fwd_overhang = i.saturating_sub(rev_len);
            let rev_overhang = i.saturating_sub(fwd_len);
            let overlap = i - fwd_overhang - rev_overhang;
            let fwd_start = fwd_len - fwd_overhang - 1;
            let rev_start = rev_len - rev_overhang - overlap;
            let (mut score, mut high, mut drop) = (0.0, 0.0f64, 0.0);
            let mut diffs = 0;
            for j in 0..overlap {
                let (f, r) = (fwd_start - j, rev_start + j);
                let cell = usize::from(fwd_qual[f]) * 128 + usize::from(rev_qual[r]);
                if fwd[f] == complement(rev[r]) {
                    score += self.tables.match_score[cell];
                    high = high.max(score);
                } else {
                    score += self.tables.mismatch_score[cell];
                    diffs += 1;
                    if score < high - drop {
                        drop = high - score;
                    }
                }
            }
            if drop >= DROP_MAX {
                score = 0.0;
            }
            if score >= min_score {
                hits += 1;
            }
            if score > best_score {
                (best_score, best_i, best_diffs) = (score, i, diffs);
            }
        }

        if hits > 1 {
            Err(Reason::Repeat)
        } else if !self.settings.stagger && best_i > fwd_len {
            Err(Reason::Staggered)
        } else if best_diffs > self.settings.max_diffs {
            Err(Reason::MaxDiffs)
        } else if 100.0 * best_diffs as f64 / best_i as f64 > self.settings.max_diff_pct {
            Err(Reason::MaxDiffPct)
        } else if !kmers {
            Err(Reason::NoKmers)
        } else if best_score < min_score {
            Err(Reason::MinScore)
        } else if best_i < min_overlap {
            Err(Reason::MinOverlap)
        } else if fwd_len + rev_len - best_i > MAX_MERGE_LENGTH {
            Err(Reason::MaxMergeLength)
        } else {
            Ok(best_i)
        }
    }

    /// builds the merged read from the overlap.
    fn join(&mut self, overlap: usize) {
        let (fwd, fwd_qual) = (&self.forward.0, &self.forward.1);
        let (rev, rev_qual) = (&self.reverse.0, &self.reverse.1);
        let tables = &self.tables;
        let error = |symbol: u8| tables.error[usize::from(symbol)];
        let m = &mut self.merged;
        m.seq.clear();
        m.qual.clear();
        (m.ee_merged, m.ee_forward, m.ee_reverse) = (0.0, 0.0, 0.0);
        (m.forward_errors, m.reverse_errors) = (0, 0);

        // the forward read before the overlap
        let fwd_overhang = fwd.len().saturating_sub(overlap);
        for f in 0..fwd_overhang {
            m.seq.push(fwd[f]);
            m.qual.push(fwd_qual[f]);
            m.ee_merged += error(fwd_qual[f]);
            m.ee_forward += error(fwd_qual[f]);
        }

        // the overlap, after any reverse bases past the end of the forward
        // read
        let rev_overhang = overlap.saturating_sub(fwd.len());
        let mut r = rev.len() as i64 - 1 - rev_overhang as i64;
        let mut f = fwd_overhang;
        while f < fwd.len() && r >= 0 {
            let (fwd_base, fq) = (fwd[f], fwd_qual[f]);
            let (rev_base, rq) = (complement(rev[r as usize]), rev_qual[r as usize]);
            let (base, symbol) = if rev_base == b'N' {
                (fwd_base, fq)
            } else if fwd_base == b'N' {
                (rev_base, rq)
            } else if fwd_base == rev_base {
                (
                    fwd_base,
                    tables.same[usize::from(fq) * 128 + usize::from(rq)],
                )
            } else if fq > rq {
                (
                    fwd_base,
                    tables.diff[usize::from(fq) * 128 + usize::from(rq)],
                )
            } else {
                (
                    rev_base,
                    tables.diff[usize::from(rq) * 128 + usize::from(fq)],
                )
            };
            if base != fwd_base {
                m.forward_errors += 1;
            }
            if base != rev_base {
                m.reverse_errors += 1;
            }
            m.seq.push(base);
            m.qual.push(symbol);
            m.ee_merged += error(symbol);
            m.ee_forward += error(fq);
            m.ee_reverse += error(rq);
            f += 1;
            r -= 1;
        }

        // the reverse read after the overlap
        while r >= 0 {
            let symbol = rev_qual[r as usize];
            m.seq.push(complement(rev[r as usize]));
            m.qual.push(symbol);
            m.ee_merged += error(symbol);
            m.ee_reverse += error(symbol);
            r -= 1;
        }
    }
}
