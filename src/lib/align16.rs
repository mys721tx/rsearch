//! vsearch's 16-bit global aligner (`align_simd.cpp`), ported lane for
//! lane. Eight targets are aligned at once against one query, four target
//! residues a step, in saturating 16-bit cells. The profile, direction
//! bits, channel scheduling and backtrack are the same, so that scores,
//! CIGARs and ties are vsearch's. A pair it cannot score in 16 bits gets
//! the score `i16::MAX` and no alignment, for the linear-memory aligner.

use crate::scoring::{FOUR_BIT, N, Scoring, is_ambiguous};

const CHANNELS: usize = 8;
const CDEPTH: usize = 4;
const MAX_LEN_PRODUCT: i64 = 25_000_000;
const MAX_LEN_SUM: i64 = 65_535;

type V = [i16; CHANNELS];

#[inline]
fn add(a: V, b: V) -> V {
    std::array::from_fn(|i| a[i].saturating_add(b[i]))
}
#[inline]
fn sub(a: V, b: V) -> V {
    std::array::from_fn(|i| a[i].saturating_sub(b[i]))
}
#[inline]
fn sub_unsigned(a: V, b: V) -> V {
    std::array::from_fn(|i| (a[i] as u16).saturating_sub(b[i] as u16) as i16)
}
#[inline]
fn max(a: V, b: V) -> V {
    std::array::from_fn(|i| a[i].max(b[i]))
}
#[inline]
fn min(a: V, b: V) -> V {
    std::array::from_fn(|i| a[i].min(b[i]))
}
#[inline]
fn and(a: V, b: V) -> V {
    std::array::from_fn(|i| a[i] & b[i])
}
#[inline]
fn xor(a: V, b: V) -> V {
    std::array::from_fn(|i| a[i] ^ b[i])
}
#[inline]
fn dup(x: i16) -> V {
    [x; CHANNELS]
}
#[inline]
fn shift_left(a: V) -> V {
    std::array::from_fn(|i| if i == 0 { 0 } else { a[i - 1] })
}
/// two bits a channel where `a` is above `b`, as `_mm_movemask_epi8` of
/// `_mm_cmpgt_epi16`.
#[inline]
fn mask_gt(a: V, b: V) -> u16 {
    let mut mask = 0u16;
    for i in 0..CHANNELS {
        if a[i] > b[i] {
            mask |= 3 << (2 * i);
        }
    }
    mask
}

/// whether vsearch's 16-bit aligner takes sequences of these lengths.
pub fn fits(qlen: usize, dlen: usize) -> bool {
    (qlen + dlen) as i64 <= MAX_LEN_SUM && (qlen as i64) * (dlen as i64) <= MAX_LEN_PRODUCT
}

/// The alignment of a query with a target: its score, or `i16::MAX` when
/// the 16-bit aligner could not score it; the columns, matches,
/// mismatches and gap runs; and the CIGAR.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Aligned {
    pub score: i16,
    pub length: u16,
    pub matches: u16,
    pub mismatches: u16,
    pub gaps: u16,
    pub cigar: String,
}

impl Aligned {
    fn overflow() -> Aligned {
        Aligned {
            score: i16::MAX,
            ..Aligned::default()
        }
    }
}

/// vsearch's `s16info_s`: the scores and penalties in 16-bit cells.
pub struct Search16 {
    /// The score of each pair of 4-bit codes.
    matrix: [[i16; 16]; 16],
    /// Penalties (opening already less the extension, as vsearch's):
    /// [query, target] x [left, interior, right].
    open: [[i16; 3]; 2],
    extend: [[i16; 3]; 2],
    n_mismatch: bool,
    /// A score or penalty beyond the cells: every pair to the fallback.
    fallback: bool,
    hearray: Vec<V>,
    dir: Vec<u16>,
    qcodes: Vec<u8>,
    qseq: Vec<u8>,
}

/// the value of `value` in a cell, as vsearch's `clamp_to_cell`.
fn clamp(value: i64, limit: i64, fallback: &mut bool) -> i16 {
    if value > limit {
        *fallback = true;
        limit as i16
    } else if value < -limit {
        *fallback = true;
        (-limit) as i16
    } else {
        value as i16
    }
}

impl Search16 {
    pub fn new(scoring: &Scoring) -> Search16 {
        let mut fallback = false;
        let score_limit = i64::from(i16::MAX);
        let penalty_limit = i64::from(i16::MAX) / (1 + CDEPTH as i64);
        let m = clamp(scoring.match_score, score_limit, &mut fallback);
        let mm = clamp(scoring.mismatch, score_limit, &mut fallback);
        let mut matrix = [[0i16; 16]; 16];
        for (a, row) in matrix.iter_mut().enumerate() {
            for (b, cell) in row.iter_mut().enumerate() {
                let (a, b) = (a as u8, b as u8);
                *cell = if scoring.n_mismatch && (a == N || b == N) {
                    mm
                } else if is_ambiguous(a) || is_ambiguous(b) {
                    0
                } else if a == b {
                    m
                } else {
                    mm
                };
            }
        }
        let mut open = [[0i16; 3]; 2];
        let mut extend = [[0i16; 3]; 2];
        for e in 0..3 {
            for (s, (opens, exts)) in [
                (&scoring.open_query, &scoring.extend_query),
                (&scoring.open_target, &scoring.extend_target),
            ]
            .into_iter()
            .enumerate()
            {
                open[s][e] = clamp(opens[e] - exts[e], penalty_limit, &mut fallback);
                extend[s][e] = clamp(exts[e], penalty_limit, &mut fallback);
            }
        }
        Search16 {
            matrix,
            open,
            extend,
            n_mismatch: scoring.n_mismatch,
            fallback,
            hearray: Vec::new(),
            dir: Vec::new(),
            qcodes: Vec::new(),
            qseq: Vec::new(),
        }
    }

    /// vsearch's `search16_qprep`.
    pub fn set_query(&mut self, query: &[u8]) {
        self.qseq = query.to_vec();
        self.qcodes = query.iter().map(|&c| FOUR_BIT[usize::from(c)]).collect();
        self.hearray.clear();
        self.hearray.resize(2 * query.len(), [0; CHANNELS]);
    }

    fn score_min(&self) -> i16 {
        let mut most = 0i32;
        for s in 0..2 {
            for e in 0..3 {
                most = most.max(i32::from(self.open[s][e]) + i32::from(self.extend[s][e]));
            }
        }
        (i32::from(i16::MIN) + most) as i16
    }

    /// aligns the query with each of `targets`, as vsearch's `search16`.
    pub fn align(&mut self, targets: &[&[u8]]) -> Vec<Aligned> {
        let n = targets.len();
        let mut out = vec![Aligned::default(); n];
        let qlen = self.qseq.len();
        if self.fallback {
            return vec![Aligned::overflow(); n];
        }
        let (q, t) = (0usize, 1usize);
        let (left, interior, right) = (0usize, 1usize, 2usize);
        if qlen == 0 {
            for (slot, target) in out.iter_mut().zip(targets) {
                let length = target.len() as i64;
                if !fits(0, target.len()) {
                    *slot = Aligned::overflow();
                    continue;
                }
                slot.length = length as u16;
                slot.gaps = length as u16;
                slot.score = if length == 0 {
                    0
                } else {
                    let l =
                        -i64::from(self.open[t][left]) - length * i64::from(self.extend[t][left]);
                    let r =
                        -i64::from(self.open[t][right]) - length * i64::from(self.extend[t][right]);
                    l.max(r) as i16
                };
                if length > 0 {
                    slot.cigar = format!("{length}I");
                }
            }
            return out;
        }
        let mut maxdlen = 0usize;
        for target in targets {
            if fits(qlen, target.len()) {
                maxdlen = maxdlen.max(target.len());
            }
        }
        maxdlen = 4 * maxdlen.div_ceil(4);
        let dirsize = qlen * maxdlen * 4;
        if self.dir.len() < dirsize {
            self.dir.resize(dirsize, 0);
        }

        let pen = |x: i16| dup(x);
        let qr = |s: usize, e: usize| dup(self.open[s][e].wrapping_add(self.extend[s][e]));
        let t0: V = [-1, 0, 0, 0, 0, 0, 0, 0];
        let r_query_left = pen(self.extend[q][left]);
        let qr_query_interior = qr(q, interior);
        let r_query_interior = pen(self.extend[q][interior]);
        let qr_query_right = qr(q, right);
        let r_query_right = pen(self.extend[q][right]);
        let qr_target_left = qr(t, left);
        let r_target_left = pen(self.extend[t][left]);
        let qr_target_interior = qr(t, interior);
        let r_target_interior = pen(self.extend[t][interior]);
        let qr_target_right = qr(t, right);
        let r_target_right = pen(self.extend[t][right]);

        // per channel: the target, where its next residue is, its offset in
        // the direction buffer, and whether it overflowed
        let mut begin = [0usize; CHANNELS];
        let mut end = [0usize; CHANNELS];
        let mut offset = [0usize; CHANNELS];
        let mut seq_id = [-1i64; CHANNELS];
        let mut overflow = [false; CHANNELS];
        let mut dseq = [[0u8; CHANNELS]; CDEPTH];
        let mut s_vec: [V; 4] = [[0; CHANNELS]; 4];
        let score_min = self.score_min();
        let score_max = i16::MAX;
        let (mut h0, mut h1, mut h2, mut h3) =
            ([0; CHANNELS], [0; CHANNELS], [0; CHANNELS], [0; CHANNELS]);
        let (mut f0, mut f1, mut f2, mut f3) =
            ([0; CHANNELS], [0; CHANNELS], [0; CHANNELS], [0; CHANNELS]);
        let mut easy = false;
        let mut dir = 0usize;
        let mut next_id = 0usize;
        let mut done = 0usize;
        let go_q_l = i64::from(self.open[q][left]);
        let ge_q_l = i64::from(self.extend[q][left]);

        loop {
            let mut masked: Option<V> = None;
            if easy {
                for c in 0..CHANNELS {
                    for depth in dseq.iter_mut() {
                        depth[c] = if begin[c] < end[c] {
                            let code = FOUR_BIT[usize::from(targets[seq_id[c] as usize][begin[c]])];
                            begin[c] += 1;
                            code
                        } else {
                            0
                        };
                    }
                    if begin[c] == end[c] {
                        easy = false;
                    }
                }
            } else {
                easy = true;
                let mut m: V = [0; CHANNELS];
                let mut tt = t0;
                for c in 0..CHANNELS {
                    if begin[c] < end[c] {
                        let target = targets[seq_id[c] as usize];
                        for depth in dseq.iter_mut() {
                            depth[c] = if begin[c] < end[c] {
                                let code = FOUR_BIT[usize::from(target[begin[c]])];
                                begin[c] += 1;
                                code
                            } else {
                                0
                            };
                        }
                        if begin[c] == end[c] {
                            easy = false;
                        }
                    } else {
                        m = xor(m, tt);
                        let cand = seq_id[c];
                        if cand >= 0 {
                            let slot = cand as usize;
                            let dlen = targets[slot].len();
                            let z = (dlen + 3) % 4;
                            if overflow[c] {
                                out[slot] = Aligned::overflow();
                            } else {
                                let score = s_vec[z][c];
                                out[slot] = self.backtrack(targets[slot], offset[c], c, maxdlen);
                                out[slot].score = score;
                            }
                            done += 1;
                        }
                        let mut length = 0usize;
                        let mut cand_id = 0usize;
                        while length == 0 && next_id < n {
                            cand_id = next_id;
                            next_id += 1;
                            length = targets[cand_id].len();
                            if length == 0 || !fits(qlen, length) {
                                out[cand_id] = Aligned::overflow();
                                length = 0;
                                done += 1;
                            }
                        }
                        if length > 0 {
                            seq_id[c] = cand_id as i64;
                            begin[c] = 0;
                            end[c] = length;
                            offset[c] = dir;
                            overflow[c] = false;
                            h0[c] = 0;
                            let edge = |k: i64| (-go_q_l - k * ge_q_l) as i16;
                            h1[c] = edge(1);
                            h2[c] = edge(2);
                            h3[c] = edge(3);
                            f0[c] = edge(1);
                            f1[c] = edge(2);
                            f2[c] = edge(3);
                            f3[c] = edge(4);
                            let target = targets[cand_id];
                            for depth in dseq.iter_mut() {
                                depth[c] = if begin[c] < end[c] {
                                    let code = FOUR_BIT[usize::from(target[begin[c]])];
                                    begin[c] += 1;
                                    code
                                } else {
                                    0
                                };
                            }
                            if begin[c] == end[c] {
                                easy = false;
                            }
                        } else {
                            seq_id[c] = -1;
                            begin[c] = 0;
                            end[c] = 0;
                            offset[c] = 0;
                            for depth in dseq.iter_mut() {
                                depth[c] = 0;
                            }
                        }
                    }
                    tt = shift_left(tt);
                }
                if done == n {
                    break;
                }
                masked = Some(m);
            }

            // the profile: the score of each query code against each target
            // residue of the step
            let mut profile = [[[0i16; CHANNELS]; CDEPTH]; 16];
            for (r, rows) in profile.iter_mut().enumerate() {
                for (j, lanes) in rows.iter_mut().enumerate() {
                    for c in 0..CHANNELS {
                        lanes[c] = self.matrix[usize::from(dseq[j][c])][r];
                    }
                }
            }
            let mut qr_target = [qr_target_interior; 4];
            let mut r_target = [r_target_interior; 4];
            if !easy {
                let qr_diff = sub(qr_target_right, qr_target_interior);
                let r_diff = sub(r_target_right, r_target_interior);
                for j in 0..CDEPTH {
                    let mut mm: V = [0; CHANNELS];
                    let mut tt = t0;
                    for c in 0..CHANNELS {
                        let dlen = if seq_id[c] >= 0 { end[c] } else { 0 };
                        if begin[c] == end[c] && j >= (dlen + 3) % 4 {
                            mm = xor(mm, tt);
                        }
                        tt = shift_left(tt);
                    }
                    qr_target[j] = add(qr_target_interior, and(qr_diff, mm));
                    r_target[j] = add(r_target_interior, and(r_diff, mm));
                }
            }
            let mask = masked.map(|m| Masked {
                mm: m,
                qr_t_left: and(m, qr_target_left),
                r_t_left: and(m, r_target_left),
                qr_q_interior: and(m, qr_query_interior),
                qr_q_right: and(m, qr_query_right),
            });
            let (h_min, h_max) = self.align_columns(
                &mut s_vec,
                &profile,
                [
                    qr_query_interior,
                    r_query_interior,
                    qr_query_right,
                    r_query_right,
                ],
                qr_target,
                r_target,
                [h0, h1, h2, h3],
                [f0, f1, f2, f3],
                mask,
                dir,
            );
            for c in 0..CHANNELS {
                if !overflow[c] && (h_min[c] <= score_min || h_max[c] >= score_max) {
                    overflow[c] = true;
                }
            }
            h0 = sub(h3, r_query_left);
            h1 = sub(h0, r_query_left);
            h2 = sub(h1, r_query_left);
            h3 = sub(h2, r_query_left);
            f0 = sub(f3, r_query_left);
            f1 = sub(f0, r_query_left);
            f2 = sub(f1, r_query_left);
            f3 = sub(f2, r_query_left);
            dir += 16 * qlen;
            if dir >= dirsize {
                dir -= dirsize;
            }
        }
        out
    }

    /// vsearch's `aligncolumns`: four target columns down the query.
    #[allow(clippy::too_many_arguments)]
    fn align_columns(
        &mut self,
        s_vec: &mut [V; 4],
        profile: &[[V; CDEPTH]; 16],
        query: [V; 4],
        qr_t: [V; 4],
        r_t: [V; 4],
        h: [V; 4],
        f: [V; 4],
        mut mask: Option<Masked>,
        dir: usize,
    ) -> (V, V) {
        let [qr_q_i, r_q_i, qr_q_r, r_q_r] = query;
        let ql = self.qcodes.len();
        let (mut h0, mut h1, mut h2, mut h3) = (h[0], h[1], h[2], h[3]);
        let mut f = [
            sub(f[0], qr_t[0]),
            sub(f[1], qr_t[1]),
            sub(f[2], qr_t[2]),
            sub(f[3], qr_t[3]),
        ];
        let mut h_min: V = [0; CHANNELS];
        let mut h_max: V = [0; CHANNELS];
        let (mut h5, mut h6, mut h7, mut h8) =
            ([0; CHANNELS], [0; CHANNELS], [0; CHANNELS], [0; CHANNELS]);
        for i in 0..ql {
            let last = i == ql - 1;
            let vp = &profile[usize::from(self.qcodes[i])];
            let mut h4 = self.hearray[2 * i];
            let mut e = self.hearray[2 * i + 1];
            let (qr_q, r_q) = if last {
                (qr_q_r, r_q_r)
            } else {
                (qr_q_i, r_q_i)
            };
            if let Some(m) = mask.as_mut() {
                if last {
                    e = sub_unsigned(e, m.mm);
                    e = sub(e, m.qr_t_left);
                    e = sub(e, m.qr_q_right);
                } else {
                    h4 = sub_unsigned(h4, m.mm);
                    h4 = sub(h4, m.qr_t_left);
                    e = sub_unsigned(e, m.mm);
                    e = sub(e, m.qr_t_left);
                    e = sub(e, m.qr_q_interior);
                    m.qr_t_left = add(m.qr_t_left, m.r_t_left);
                }
            }
            let at = dir + 16 * i;
            let (dir_buf, hmin, hmax) = (&mut self.dir, &mut h_min, &mut h_max);
            let mut step = |hh: &mut V, nn: &mut V, ff: &mut V, v: V, k: usize, e: &mut V| {
                let p = &mut dir_buf[at + 4 * k..at + 4 * k + 4];
                *hh = add(*hh, v);
                p[0] = mask_gt(*ff, *hh);
                *hh = max(*hh, *ff);
                p[1] = mask_gt(*e, *hh);
                *hh = max(*hh, *e);
                *hmin = min(*hmin, *hh);
                *hmax = max(*hmax, *hh);
                *nn = *hh;
                let hf = sub(*hh, qr_t[k]);
                *ff = sub(*ff, r_t[k]);
                p[2] = mask_gt(*ff, hf);
                *ff = max(*ff, hf);
                let he = sub(*hh, qr_q);
                *e = sub(*e, r_q);
                p[3] = mask_gt(*e, he);
                *e = max(*e, he);
            };
            let [f0, f1, f2, f3] = &mut f;
            step(&mut h0, &mut h5, f0, vp[0], 0, &mut e);
            step(&mut h1, &mut h6, f1, vp[1], 1, &mut e);
            step(&mut h2, &mut h7, f2, vp[2], 2, &mut e);
            step(&mut h3, &mut h8, f3, vp[3], 3, &mut e);
            self.hearray[2 * i] = h8;
            self.hearray[2 * i + 1] = e;
            h0 = h4;
            h1 = h5;
            h2 = h6;
            h3 = h7;
        }
        *s_vec = [h5, h6, h7, h8];
        (h_min, h_max)
    }

    /// vsearch's `backtrack16`: the alignment of a target from the
    /// direction bits of its channel.
    fn backtrack(&self, dseq: &[u8], offset: usize, channel: usize, maxdlen: usize) -> Aligned {
        let qlen = self.qseq.len();
        let dlen = dseq.len();
        let dirsize = qlen * maxdlen * 4;
        let shift = 2 * channel;
        let (mut aligned, mut matches, mut mismatches, mut gaps) = (0u16, 0u16, 0u16, 0u16);
        let (mut i, mut j) = (qlen as i64 - 1, dlen as i64 - 1);
        let mut ops: Vec<(u8, usize)> = Vec::new();
        let mut op = 0u8;
        let push = |ops: &mut Vec<(u8, usize)>, op: &mut u8, new: u8| {
            if let Some(last) = ops.last_mut()
                && last.0 == new
            {
                last.1 += 1;
            } else {
                ops.push((new, 1));
            }
            *op = new;
        };
        while i >= 0 && j >= 0 {
            aligned = aligned.wrapping_add(1);
            let at =
                (offset + 16 * qlen * (j as usize / 4) + 16 * i as usize + 4 * (j as usize & 3))
                    % dirsize;
            let d = &self.dir[at..at + 4];
            let bit = |w: u16| (w >> shift) & 3 != 0;
            if op == b'I' && bit(d[3]) {
                j -= 1;
                push(&mut ops, &mut op, b'I');
            } else if op == b'D' && bit(d[2]) {
                i -= 1;
                push(&mut ops, &mut op, b'D');
            } else if bit(d[1]) {
                if op != b'I' {
                    gaps = gaps.wrapping_add(1);
                }
                j -= 1;
                push(&mut ops, &mut op, b'I');
            } else if bit(d[0]) {
                if op != b'D' {
                    gaps = gaps.wrapping_add(1);
                }
                i -= 1;
                push(&mut ops, &mut op, b'D');
            } else {
                let (a, b) = (
                    FOUR_BIT[usize::from(self.qseq[i as usize])],
                    FOUR_BIT[usize::from(dseq[j as usize])],
                );
                if a & b != 0 {
                    if self.n_mismatch && (a == N || b == N) {
                        mismatches = mismatches.wrapping_add(1);
                    } else {
                        matches = matches.wrapping_add(1);
                    }
                } else {
                    mismatches = mismatches.wrapping_add(1);
                }
                i -= 1;
                j -= 1;
                push(&mut ops, &mut op, b'M');
            }
        }
        while i >= 0 {
            aligned = aligned.wrapping_add(1);
            if op != b'D' {
                gaps = gaps.wrapping_add(1);
            }
            i -= 1;
            push(&mut ops, &mut op, b'D');
        }
        while j >= 0 {
            aligned = aligned.wrapping_add(1);
            if op != b'I' {
                gaps = gaps.wrapping_add(1);
            }
            j -= 1;
            push(&mut ops, &mut op, b'I');
        }
        let mut cigar = String::new();
        for &(op, run) in ops.iter().rev() {
            if run > 1 {
                cigar.push_str(&run.to_string());
            }
            cigar.push(op as char);
        }
        Aligned {
            score: 0,
            length: aligned,
            matches,
            mismatches,
            gaps,
            cigar,
        }
    }
}

/// The masks of vsearch's `Masked` columns: the channels starting a target.
struct Masked {
    mm: V,
    qr_t_left: V,
    r_t_left: V,
    qr_q_interior: V,
    qr_q_right: V,
}
