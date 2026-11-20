//! Alignment scores, as vsearch's: a match and a mismatch score for
//! nucleotides (ambiguity codes score 0, and N a mismatch with
//! `--n_mismatch`), and gap opening and extension penalties for gaps in
//! the query or the target, at the left end, inside, or at the right end.

/// The 4-bit code of each byte: a bit for each of A, C, G and T or U, so
/// that IUPAC codes hold the bases they stand for; 0 for other bytes.
pub const FOUR_BIT: [u8; 256] = {
    let mut code = [0u8; 256];
    let codes: [(u8, u8); 16] = [
        (b'A', 1),
        (b'B', 14),
        (b'C', 2),
        (b'D', 13),
        (b'G', 4),
        (b'H', 11),
        (b'K', 12),
        (b'M', 3),
        (b'N', 15),
        (b'R', 5),
        (b'S', 6),
        (b'T', 8),
        (b'U', 8),
        (b'V', 7),
        (b'W', 9),
        (b'Y', 10),
    ];
    let mut i = 0;
    while i < codes.len() {
        let (letter, value) = codes[i];
        code[letter as usize] = value;
        code[letter.to_ascii_lowercase() as usize] = value;
        i += 1;
    }
    code
};

/// The code of N.
pub const N: u8 = 15;

/// whether a 4-bit code stands for anything but a single base.
pub fn is_ambiguous(code: u8) -> bool {
    !matches!(code, 1 | 2 | 4 | 8)
}

/// Where a gap is.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum End {
    Left = 0,
    Interior = 1,
    Right = 2,
}

/// The penalty vsearch gives an infinite gap penalty, `*`.
pub const INFINITE: i64 = i32::MAX as i64;

/// The largest finite gap penalty vsearch takes.
pub const MAX_GAP_PENALTY: i64 = 32767 / 5;

/// Scores of an alignment.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Scoring {
    pub match_score: i64,
    pub mismatch: i64,
    pub n_mismatch: bool,
    /// Opening penalties of gaps in the query, by end, and in the target.
    pub open_query: [i64; 3],
    pub open_target: [i64; 3],
    /// Extension penalties, likewise.
    pub extend_query: [i64; 3],
    pub extend_target: [i64; 3],
}

impl Default for Scoring {
    /// vsearch's defaults: match 2, mismatch -4, `--gapopen 20I/2E`,
    /// `--gapext 2I/1E`.
    fn default() -> Scoring {
        Scoring {
            match_score: 2,
            mismatch: -4,
            n_mismatch: false,
            open_query: [2, 20, 2],
            open_target: [2, 20, 2],
            extend_query: [1, 2, 1],
            extend_target: [1, 2, 1],
        }
    }
}

impl Scoring {
    /// the score of a column of `a` against `b`.
    pub fn substitution(&self, a: u8, b: u8) -> i64 {
        let (a, b) = (FOUR_BIT[usize::from(a)], FOUR_BIT[usize::from(b)]);
        if self.n_mismatch && (a == N || b == N) {
            self.mismatch
        } else if is_ambiguous(a) || is_ambiguous(b) {
            0
        } else if a == b {
            self.match_score
        } else {
            self.mismatch
        }
    }

    /// the scores and match classes of every pair of 4-bit codes.
    pub fn substitutions(&self) -> Substitutions {
        let mut table = Substitutions {
            score: [[0; 16]; 16],
            is_match: [[false; 16]; 16],
        };
        for a in 0..16u8 {
            for b in 0..16u8 {
                let (i, j) = (usize::from(a), usize::from(b));
                table.score[i][j] = if self.n_mismatch && (a == N || b == N) {
                    self.mismatch
                } else if is_ambiguous(a) || is_ambiguous(b) {
                    0
                } else if a == b {
                    self.match_score
                } else {
                    self.mismatch
                };
                table.is_match[i][j] = !(self.n_mismatch && (a == N || b == N)) && a & b != 0;
            }
        }
        table
    }

    /// the penalty of a gap of `run` residues in the query (`query`) or the
    /// target, at end `end` (0 left, 1 interior, 2 right): as vsearch, whose
    /// opening penalty covers the first residue, the opening penalty plus
    /// the extension penalty of each further residue.
    pub fn gap_penalty(&self, query: bool, end: usize, run: usize) -> i64 {
        let (open, extend) = if query {
            (self.open_query[end], self.extend_query[end])
        } else {
            (self.open_target[end], self.extend_target[end])
        };
        open - extend + run as i64 * extend
    }

    /// applies a gap penalty string of vsearch's `--gapopen` (`open`) or
    /// `--gapext`: tokens separated by `/`, each a number from 0 to 6553 or
    /// `*`, followed by the ends (`L`, `I`, `R`, or `E` for both ends) and
    /// the sequences (`Q`, `T`) it applies to, all of them by default.
    ///
    /// # Example
    /// ```rust
    /// # use rsearch::scoring::Scoring;
    /// let mut scoring = Scoring::default();
    /// scoring.set_gaps("3Q/5T/*LQ", true).unwrap();
    /// assert_eq!(scoring.open_query, [i32::MAX as i64, 3, 3]);
    /// assert_eq!(scoring.open_target, [5, 5, 5]);
    /// assert!(scoring.set_gaps("2EL", true).is_err());
    /// ```
    pub fn set_gaps(&mut self, text: &str, open: bool) -> Result<(), String> {
        let bytes = text.as_bytes();
        let mut i = 0;
        while i < bytes.len() {
            let start = i;
            let penalty = if bytes[i] == b'*' {
                i += 1;
                INFINITE
            } else {
                let digits = bytes[i..]
                    .iter()
                    .enumerate()
                    .take_while(|&(j, &c)| c.is_ascii_digit() || (j == 0 && c == b'+'))
                    .count();
                let number = &text[i..i + digits];
                if digits == 0 || number == "+" {
                    return Err(format!("invalid gap penalty argument ({})", &text[i..]));
                }
                i += digits;
                match number.parse::<i64>() {
                    Ok(p) if p <= MAX_GAP_PENALTY => p,
                    _ => {
                        return Err(format!(
                            "a finite gap penalty must be in the range 0 to {MAX_GAP_PENALTY}"
                        ));
                    }
                }
            };
            let (mut e, mut l, mut m, mut r, mut q, mut t) =
                (false, false, false, false, false, false);
            while i < bytes.len() && bytes[i] != b'/' {
                match bytes[i] {
                    b'E' => e = true,
                    b'I' => m = true,
                    b'L' => l = true,
                    b'R' => r = true,
                    b'Q' => q = true,
                    b'T' => t = true,
                    c => return Err(format!("invalid character '{}' in gap penalty", c as char)),
                }
                i += 1;
            }
            if i < bytes.len() {
                i += 1;
            }
            if e && (l || r) {
                return Err(format!(
                    "invalid gap penalty (E with L or R) '{}'",
                    &text[start..]
                ));
            }
            if e {
                (l, r) = (true, true);
            }
            if !l && !m && !r {
                (l, m, r) = (true, true, true);
            }
            if !q && !t {
                (q, t) = (true, true);
            }
            let (query, target) = if open {
                (&mut self.open_query, &mut self.open_target)
            } else {
                (&mut self.extend_query, &mut self.extend_target)
            };
            for (applies, penalties) in [(q, query), (t, target)] {
                if applies {
                    for (end, set) in [l, m, r].into_iter().enumerate() {
                        if set {
                            penalties[end] = penalty;
                        }
                    }
                }
            }
        }
        Ok(())
    }
}

/// The score of each pair of 4-bit codes, and whether vsearch counts it as
/// a match: codes sharing a base, but N with `--n_mismatch`.
#[derive(Clone, Debug)]
pub struct Substitutions {
    pub score: [[i64; 16]; 16],
    pub is_match: [[bool; 16]; 16],
}
