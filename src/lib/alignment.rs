//! Alignments as CIGAR strings, as vsearch writes them in its hits: runs
//! of `M` (a column of both sequences), `D` (a query residue facing a gap
//! in the target) and `I` (a target residue facing a gap in the query), a
//! run of one without its length; and the statistics vsearch derives from
//! them (`alignstats` and `align_trim`).

use crate::scoring::{End, FOUR_BIT, Scoring, Substitutions};

/// parses a CIGAR into (operation, run) pairs.
///
/// # Example
/// ```rust
/// # use rsearch::alignment::parse_cigar;
/// assert_eq!(parse_cigar("10I80MD").unwrap(), vec![(b'I', 10), (b'M', 80), (b'D', 1)]);
/// assert!(parse_cigar("10X").is_err());
/// ```
pub fn parse_cigar(cigar: &str) -> Result<Vec<(u8, usize)>, String> {
    let mut ops = Vec::new();
    let mut run = 0usize;
    let mut digits = false;
    for &c in cigar.as_bytes() {
        if c.is_ascii_digit() {
            run = run
                .checked_mul(10)
                .and_then(|r| r.checked_add(usize::from(c - b'0')))
                .ok_or("CIGAR run too long")?;
            digits = true;
        } else if matches!(c, b'M' | b'D' | b'I') {
            ops.push((c, if digits { run.max(1) } else { 1 }));
            run = 0;
            digits = false;
        } else {
            return Err(format!("invalid CIGAR operation '{}'", c as char));
        }
    }
    if digits {
        return Err("CIGAR ends with a run length".to_string());
    }
    Ok(ops)
}

/// What vsearch knows of a hit's alignment.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Stats {
    pub score: i64,
    pub length: i64,
    pub matches: i64,
    pub mismatches: i64,
    /// Gap runs, the terminal ones included.
    pub gaps: i64,
    pub indels: i64,
    pub shortest: i64,
    pub longest: i64,
    /// Residues of each sequence facing terminal gaps, and the CIGAR
    /// characters of the terminal gap runs.
    pub trim_q_left: i64,
    pub trim_q_right: i64,
    pub trim_t_left: i64,
    pub trim_t_right: i64,
    pub trim_aln_left: usize,
    pub trim_aln_right: usize,
    /// The same without terminal gaps.
    pub internal_length: i64,
    pub internal_indels: i64,
    pub internal_gaps: i64,
    /// The five identity definitions.
    pub id: [f64; 5],
}

impl Stats {
    /// derives the statistics of `cigar` aligning `query` (on the strand it
    /// was aligned) with `target`, as vsearch's `alignstats` and
    /// `align_trim`; fails when the CIGAR does not span both sequences.
    ///
    /// # Example
    /// ```rust
    /// # use rsearch::{alignment::Stats, scoring::Scoring};
    /// let scoring = Scoring::default();
    /// let table = scoring.substitutions();
    /// let stats = Stats::new("2I4M", b"ACGT", b"TTACGA", &scoring, &table).unwrap();
    /// assert_eq!((stats.matches, stats.mismatches, stats.gaps), (3, 1, 1));
    /// assert_eq!((stats.trim_t_left, stats.internal_length), (2, 4));
    /// // 3 matches, a mismatch and a left gap of 2 in the query: its
    /// // opening penalty covers its first residue
    /// assert_eq!(stats.score, 3 * 2 - 4 - (2 + 1));
    /// assert_eq!(stats.id[2], 75.0);
    /// ```
    pub fn new(
        cigar: &str,
        query: &[u8],
        target: &[u8],
        scoring: &Scoring,
        table: &Substitutions,
    ) -> Result<Stats, String> {
        let ops = parse_cigar(cigar)?;
        let mut s = Stats::default();
        let (mut a, mut b) = (0usize, 0usize);
        let last = ops.len().saturating_sub(1);
        for (k, &(op, run)) in ops.iter().enumerate() {
            let end = if a == 0 && b == 0 {
                End::Left
            } else if k == last {
                End::Right
            } else {
                End::Interior
            };
            match op {
                b'M' => {
                    if a + run > query.len() || b + run > target.len() {
                        return Err("CIGAR longer than the sequences".to_string());
                    }
                    s.length += run as i64;
                    for (&x, &y) in query[a..a + run].iter().zip(&target[b..b + run]) {
                        let i = usize::from(FOUR_BIT[usize::from(x)]);
                        let j = usize::from(FOUR_BIT[usize::from(y)]);
                        s.score += table.score[i][j];
                        s.matches += i64::from(table.is_match[i][j]);
                    }
                    s.mismatches += run as i64;
                    a += run;
                    b += run;
                }
                b'I' => {
                    if b + run > target.len() {
                        return Err("CIGAR longer than the sequences".to_string());
                    }
                    let e = end as usize;
                    s.score -= scoring.gap_penalty(true, e, run);
                    s.gaps += 1;
                    s.length += run as i64;
                    b += run;
                }
                _ => {
                    if a + run > query.len() {
                        return Err("CIGAR longer than the sequences".to_string());
                    }
                    let e = end as usize;
                    s.score -= scoring.gap_penalty(false, e, run);
                    s.gaps += 1;
                    s.length += run as i64;
                    a += run;
                }
            }
        }
        if a != query.len() || b != target.len() {
            return Err("CIGAR shorter than the sequences".to_string());
        }
        s.mismatches -= s.matches;
        s.indels = s.length - s.matches - s.mismatches;
        s.shortest = query.len().min(target.len()) as i64;
        s.longest = query.len().max(target.len()) as i64;
        s.trim(cigar, &ops);
        Ok(s)
    }

    /// the statistics of an alignment from its aligner's counts, as vsearch
    /// keeps those of its 16-bit aligner (which counts every column of an
    /// empty query's alignment as a gap), with `align_trim`'s derived ones.
    #[allow(clippy::too_many_arguments)]
    pub fn from_counts(
        cigar: &str,
        qlen: usize,
        tlen: usize,
        length: i64,
        matches: i64,
        mismatches: i64,
        gaps: i64,
        score: i64,
    ) -> Result<Stats, String> {
        let ops = parse_cigar(cigar)?;
        let mut s = Stats {
            score,
            length,
            matches,
            mismatches,
            gaps,
            indels: length - matches - mismatches,
            shortest: qlen.min(tlen) as i64,
            longest: qlen.max(tlen) as i64,
            ..Stats::default()
        };
        s.trim(cigar, &ops);
        Ok(s)
    }

    /// vsearch's `align_trim`: the terminal gaps, the internal counts and
    /// the identities.
    fn trim(&mut self, cigar: &str, ops: &[(u8, usize)]) {
        if let Some(&(op, run)) = ops.first()
            && op != b'M'
        {
            let digits = cigar.bytes().take_while(u8::is_ascii_digit).count();
            self.trim_aln_left = 1 + digits;
            if op == b'D' {
                self.trim_q_left = run as i64;
            } else {
                self.trim_t_left = run as i64;
            }
        }
        if let Some(&(op, run)) = ops.last()
            && op != b'M'
        {
            let body = &cigar.as_bytes()[..cigar.len() - 1];
            let digits = body.iter().rev().take_while(|c| c.is_ascii_digit()).count();
            self.trim_aln_right = 1 + digits;
            if op == b'D' {
                self.trim_q_right = run as i64;
            } else {
                self.trim_t_right = run as i64;
            }
        }
        if self.trim_q_left >= self.length {
            self.trim_q_right = 0;
        }
        if self.trim_t_left >= self.length {
            self.trim_t_right = 0;
        }
        let trimmed = self.trim_q_left + self.trim_t_left + self.trim_q_right + self.trim_t_right;
        self.internal_length = self.length - trimmed;
        self.internal_indels = self.indels - trimmed;
        self.internal_gaps = self.gaps
            - i64::from(self.trim_q_left + self.trim_t_left > 0)
            - i64::from(self.trim_q_right + self.trim_t_right > 0);
        let m = self.matches as f64;
        let share = |of: i64| if of > 0 { 100.0 * m / of as f64 } else { 0.0 };
        self.id = [
            share(self.shortest),
            share(self.length),
            share(self.internal_length),
            (100.0 * (1.0 - (1.0 * (self.mismatches + self.gaps) as f64) / self.longest as f64))
                .max(0.0),
            share(self.length),
        ];
    }
}

/// the row of `seq` in the alignment: its residues, and gaps where the
/// operations it does not take part in are; `query` says which sequence.
pub fn row(seq: &[u8], ops: &[(u8, usize)], query: bool) -> Vec<u8> {
    let own = if query { b'D' } else { b'I' };
    let mut row = Vec::new();
    let mut at = 0;
    for &(op, run) in ops {
        if op == b'M' || op == own {
            row.extend_from_slice(&seq[at..at + run]);
            at += run;
        } else {
            row.extend(std::iter::repeat_n(b'-', run));
        }
    }
    row
}
