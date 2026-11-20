//! vsearch's linear-memory global aligner (`linmemalign.cpp`), ported step
//! for step so that its CIGARs and ties are vsearch's. It is Hirschberg's
//! divide and conquer over affine gaps with separate terminal and interior
//! penalties. vsearch uses it when its 16-bit aligner cannot score a pair.

use crate::scoring::{FOUR_BIT, Scoring, Substitutions};

/// Which ends of the two sequences a sub-problem still includes.
#[derive(Clone, Copy)]
struct Ends {
    a_left: bool,
    a_right: bool,
    b_left: bool,
    b_right: bool,
}

/// The linear-memory aligner and its CIGAR being built.
pub struct LinearMemoryAligner {
    table: Substitutions,
    // penalties, opening less the extension as vsearch's: q gaps in the
    // query (I), t gaps in the target (D)
    go_q_l: i64,
    go_q_i: i64,
    go_q_r: i64,
    go_t_l: i64,
    go_t_i: i64,
    go_t_r: i64,
    ge_q_l: i64,
    ge_q_i: i64,
    ge_q_r: i64,
    ge_t_l: i64,
    ge_t_i: i64,
    ge_t_r: i64,
    hh: Vec<i64>,
    ee: Vec<i64>,
    xx: Vec<i64>,
    yy: Vec<i64>,
    cigar: String,
    op: u8,
    run: i64,
}

impl LinearMemoryAligner {
    pub fn new(scoring: &Scoring) -> LinearMemoryAligner {
        let go = |open: i64, ext: i64| open - ext;
        LinearMemoryAligner {
            table: scoring.substitutions(),
            go_q_l: go(scoring.open_query[0], scoring.extend_query[0]),
            go_q_i: go(scoring.open_query[1], scoring.extend_query[1]),
            go_q_r: go(scoring.open_query[2], scoring.extend_query[2]),
            go_t_l: go(scoring.open_target[0], scoring.extend_target[0]),
            go_t_i: go(scoring.open_target[1], scoring.extend_target[1]),
            go_t_r: go(scoring.open_target[2], scoring.extend_target[2]),
            ge_q_l: scoring.extend_query[0],
            ge_q_i: scoring.extend_query[1],
            ge_q_r: scoring.extend_query[2],
            ge_t_l: scoring.extend_target[0],
            ge_t_i: scoring.extend_target[1],
            ge_t_r: scoring.extend_target[2],
            hh: Vec::new(),
            ee: Vec::new(),
            xx: Vec::new(),
            yy: Vec::new(),
            cigar: String::new(),
            op: 0,
            run: 0,
        }
    }

    fn flush(&mut self) {
        if self.run <= 0 {
            return;
        }
        if self.run > 1 {
            self.cigar.push_str(&self.run.to_string());
        }
        self.cigar.push(self.op as char);
    }

    fn add(&mut self, op: u8, run: i64) {
        if self.op == op {
            self.run += run;
        } else {
            self.flush();
            self.op = op;
            self.run = run;
        }
    }

    fn subst(&self, a: u8, b: u8) -> i64 {
        self.table.score[usize::from(FOUR_BIT[usize::from(b)])]
            [usize::from(FOUR_BIT[usize::from(a)])]
    }

    /// aligns `a` (the query) with `b` (the target), and returns the CIGAR.
    ///
    /// # Example
    /// ```rust
    /// # use rsearch::{linmem::LinearMemoryAligner, scoring::Scoring};
    /// let mut aligner = LinearMemoryAligner::new(&Scoring::default());
    /// // terminal gaps cost less than an interior one: 8 - 5 - 6 beats 16 - 20
    /// assert_eq!(aligner.align(b"ACGTACGT", b"ACGTTACGT"), "4D4M5I");
    /// assert_eq!(aligner.align(b"ACGTAACGT", b"ACGTAACGT"), "9M");
    /// ```
    pub fn align(&mut self, a: &[u8], b: &[u8]) -> String {
        self.cigar.clear();
        self.op = 0;
        self.run = 0;
        if self.hh.len() < b.len() + 1 {
            for v in [&mut self.hh, &mut self.ee, &mut self.xx, &mut self.yy] {
                v.resize(b.len() + 1, 0);
            }
        }
        let whole = Ends {
            a_left: true,
            a_right: true,
            b_left: true,
            b_right: true,
        };
        self.diff(a, b, false, false, whole);
        self.flush();
        std::mem::take(&mut self.cigar)
    }

    /// vsearch's `diff`: `gap_b_left` and `gap_b_right` say whether a gap in
    /// `b` is already open at that end.
    fn diff(&mut self, a: &[u8], b: &[u8], gap_b_left: bool, gap_b_right: bool, ends: Ends) {
        const MIN: i64 = i64::MIN;
        let (a_len, b_len) = (a.len() as i64, b.len() as i64);
        if b_len == 0 {
            if a_len > 0 {
                self.add(b'D', a_len);
            }
            return;
        }
        if a_len == 0 {
            self.add(b'I', b_len);
            return;
        }
        if a_len == 1 {
            let mut score = 0i64;
            if !gap_b_left {
                score -= if ends.b_left {
                    self.go_t_l
                } else {
                    self.go_t_i
                };
            }
            score -= if ends.b_left {
                self.ge_t_l
            } else {
                self.ge_t_i
            };
            score -= if ends.a_right {
                self.go_q_r + b_len * self.ge_q_r
            } else {
                self.go_q_i + b_len * self.ge_q_i
            };
            let mut max_score = score;
            let mut best = -1i64;
            score -= if ends.a_left {
                self.go_q_l + b_len * self.ge_q_l
            } else {
                self.go_q_i + b_len * self.ge_q_i
            };
            if !gap_b_right {
                score -= if ends.b_right {
                    self.go_t_r
                } else {
                    self.go_t_i
                };
            }
            score -= if ends.b_right {
                self.ge_t_r
            } else {
                self.ge_t_i
            };
            if score > max_score {
                max_score = score;
                best = b_len;
            }
            for i in 0..b_len {
                let mut score = 0i64;
                if i > 0 {
                    score -= if ends.a_left {
                        self.go_q_l + i * self.ge_q_l
                    } else {
                        self.go_q_i + i * self.ge_q_i
                    };
                }
                score += self.subst(a[0], b[i as usize]);
                if i < b_len - 1 {
                    score -= if ends.a_right {
                        self.go_q_r + (b_len - 1 - i) * self.ge_q_r
                    } else {
                        self.go_q_i + (b_len - 1 - i) * self.ge_q_i
                    };
                }
                if score > max_score {
                    max_score = score;
                    best = i;
                }
            }
            if best == -1 {
                self.add(b'D', 1);
                self.add(b'I', b_len);
            } else if best == b_len {
                self.add(b'I', b_len);
                self.add(b'D', 1);
            } else {
                if best > 0 {
                    self.add(b'I', best);
                }
                self.add(b'M', 1);
                if best < b_len - 1 {
                    self.add(b'I', b_len - 1 - best);
                }
            }
            return;
        }

        let half = a_len / 2;
        let bl = b_len as usize;
        // forward over the first half of a
        self.hh[0] = 0;
        self.ee[0] = 0;
        for i in 1..=bl {
            let i64i = i as i64;
            self.hh[i] = -(if ends.a_left {
                self.go_q_l + i64i * self.ge_q_l
            } else {
                self.go_q_i + i64i * self.ge_q_i
            });
            self.ee[i] = MIN;
        }
        for i in 1..=half {
            let mut p = self.hh[0];
            let mut h = -(if ends.b_left {
                (if gap_b_left { 0 } else { self.go_t_l }) + i * self.ge_t_l
            } else {
                (if gap_b_left { 0 } else { self.go_t_i }) + i * self.ge_t_i
            });
            self.hh[0] = h;
            let mut f = MIN;
            let a_code = usize::from(FOUR_BIT[usize::from(a[(i - 1) as usize])]);
            for j in 1..=bl {
                f = f.max(h.wrapping_sub(self.go_q_i)).wrapping_sub(self.ge_q_i);
                if ends.b_right && j == bl {
                    self.ee[j] = self.ee[j]
                        .max(self.hh[j].wrapping_sub(self.go_t_r))
                        .wrapping_sub(self.ge_t_r);
                } else {
                    self.ee[j] = self.ee[j]
                        .max(self.hh[j].wrapping_sub(self.go_t_i))
                        .wrapping_sub(self.ge_t_i);
                }
                let b_code = usize::from(FOUR_BIT[usize::from(b[j - 1])]);
                h = p.wrapping_add(self.table.score[b_code][a_code]);
                h = f.max(h);
                h = self.ee[j].max(h);
                p = self.hh[j];
                self.hh[j] = h;
            }
        }
        self.ee[0] = self.hh[0];
        // backward over the second half of a
        self.xx[0] = 0;
        self.yy[0] = 0;
        for i in 1..=bl {
            let i64i = i as i64;
            self.xx[i] = -(if ends.a_right {
                self.go_q_r + i64i * self.ge_q_r
            } else {
                self.go_q_i + i64i * self.ge_q_i
            });
            self.yy[i] = MIN;
        }
        for i in 1..=(a_len - half) {
            let mut p = self.xx[0];
            let mut h = -(if ends.b_right {
                (if gap_b_right { 0 } else { self.go_t_r }) + i * self.ge_t_r
            } else {
                (if gap_b_right { 0 } else { self.go_t_i }) + i * self.ge_t_i
            });
            self.xx[0] = h;
            let mut f = MIN;
            let a_code = usize::from(FOUR_BIT[usize::from(a[(a_len - i) as usize])]);
            for j in 1..=bl {
                f = f.max(h.wrapping_sub(self.go_q_i)).wrapping_sub(self.ge_q_i);
                if ends.b_left && j == bl {
                    self.yy[j] = self.yy[j]
                        .max(self.xx[j].wrapping_sub(self.go_t_l))
                        .wrapping_sub(self.ge_t_l);
                } else {
                    self.yy[j] = self.yy[j]
                        .max(self.xx[j].wrapping_sub(self.go_t_i))
                        .wrapping_sub(self.ge_t_i);
                }
                let b_code = usize::from(FOUR_BIT[usize::from(b[bl - j])]);
                h = p.wrapping_add(self.table.score[b_code][a_code]);
                h = f.max(h);
                h = self.yy[j].max(h);
                p = self.xx[j];
                self.xx[j] = h;
            }
        }
        self.yy[0] = self.xx[0];
        let (mut max0, mut best0) = (MIN, -1i64);
        for i in 0..=bl {
            let score = self.hh[i].wrapping_add(self.xx[bl - i]);
            if score > max0 {
                max0 = score;
                best0 = i as i64;
            }
        }
        let (mut max1, mut best1) = (MIN, -1i64);
        for i in 0..=bl {
            let g = if ends.b_left && i == 0 {
                self.go_t_l
            } else if ends.b_right && i == bl {
                self.go_t_r
            } else {
                self.go_t_i
            };
            let score = self.ee[i].wrapping_add(self.yy[bl - i]).wrapping_add(g);
            if score > max1 {
                max1 = score;
                best1 = i as i64;
            }
        }
        let (parted, best) = if max0 > max1 {
            (false, best0)
        } else if max1 > max0 {
            (true, best1)
        } else if best0 <= best1 {
            (false, best0)
        } else {
            (true, best1)
        };
        let (half, best) = (half as usize, best as usize);
        let left_ends = Ends {
            a_left: ends.a_left,
            a_right: false,
            b_left: ends.b_left,
            b_right: ends.b_right && best == bl,
        };
        let right_ends = Ends {
            a_left: false,
            a_right: ends.a_right,
            b_left: ends.b_left && best == 0,
            b_right: ends.b_right,
        };
        if !parted {
            self.diff(&a[..half], &b[..best], gap_b_left, false, left_ends);
            self.diff(&a[half..], &b[best..], false, gap_b_right, right_ends);
        } else {
            self.diff(&a[..half - 1], &b[..best], gap_b_left, true, left_ends);
            self.add(b'D', 2);
            self.diff(&a[half + 1..], &b[best..], true, gap_b_right, right_ends);
        }
    }
}
