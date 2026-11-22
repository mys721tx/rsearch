//! vsearch's 16-bit global aligner (`align_simd.cpp`), ported lane for
//! lane. Eight targets are aligned at once against one query, sixteen with
//! AVX2, four target residues a step, in saturating 16-bit cells. The
//! profile, direction bits, channel scheduling, and backtrack are the same,
//! so that scores, CIGARs, and ties are vsearch's. Each target has a
//! channel of its own, so the number of lanes changes none of them. A pair
//! it cannot score in 16 bits gets the score `i16::MAX` and no alignment,
//! for the linear-memory aligner.

use std::fmt::Write as _;

use crate::lanes::{Lanes, MAX_LANES, Scalar};
use crate::scoring::{FOUR_BIT, N, Scoring, is_ambiguous};
use crate::simd::{Level, level};

const CDEPTH: usize = 4;
const MAX_LEN_PRODUCT: i64 = 25_000_000;
const MAX_LEN_SUM: i64 = 65_535;

#[cfg(target_arch = "aarch64")]
use crate::lanes::Neon;
#[cfg(target_arch = "x86_64")]
use crate::lanes::{Avx2, Sse2};

// The lanes' operations as functions, as vsearch writes its intrinsics.
#[inline(always)]
fn add<L: Lanes>(a: L, b: L) -> L {
    a.add(b)
}
#[inline(always)]
fn sub<L: Lanes>(a: L, b: L) -> L {
    a.sub(b)
}
#[inline(always)]
fn sub_unsigned<L: Lanes>(a: L, b: L) -> L {
    a.sub_unsigned(b)
}
#[inline(always)]
fn max<L: Lanes>(a: L, b: L) -> L {
    a.max(b)
}
#[inline(always)]
fn min<L: Lanes>(a: L, b: L) -> L {
    a.min(b)
}
#[inline(always)]
fn and<L: Lanes>(a: L, b: L) -> L {
    a.and(b)
}
#[inline(always)]
fn xor<L: Lanes>(a: L, b: L) -> L {
    a.xor(b)
}
#[inline(always)]
fn mask_gt<L: Lanes>(a: L, b: L) -> L::Bits {
    a.mask_gt(b)
}
#[inline(always)]
fn dup<L: Lanes>(x: i16) -> L {
    L::dup(x)
}
#[inline(always)]
fn zero<L: Lanes>() -> L {
    L::zero()
}
#[inline(always)]
fn set<L: Lanes>(v: &mut L, c: usize, x: i16) {
    *v = v.with(c, x);
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

/// vsearch's 16-bit aligner, on the lanes [`crate::simd::level`] chose.
pub struct Search16 {
    kernel: Kernel,
}

/// The aligner on each set of lanes.
enum Kernel {
    Scalar(Engine<Scalar>),
    #[cfg(target_arch = "x86_64")]
    Sse2(Engine<Sse2>),
    /// sixteen lanes, and eight for a call of eight targets or fewer,
    /// which would leave half of sixteen idle
    #[cfg(target_arch = "x86_64")]
    Avx2 {
        wide: Engine<Avx2>,
        narrow: Engine<Sse2>,
    },
    #[cfg(target_arch = "aarch64")]
    Neon(Engine<Neon>),
}

impl Search16 {
    pub fn new(scoring: &Scoring) -> Search16 {
        let kernel = match level() {
            #[cfg(target_arch = "x86_64")]
            Level::Avx2 => Kernel::Avx2 {
                wide: Engine::new(scoring),
                narrow: Engine::new(scoring),
            },
            #[cfg(target_arch = "x86_64")]
            Level::Sse2 => Kernel::Sse2(Engine::new(scoring)),
            #[cfg(target_arch = "aarch64")]
            Level::Neon => Kernel::Neon(Engine::new(scoring)),
            _ => Kernel::Scalar(Engine::new(scoring)),
        };
        Search16 { kernel }
    }

    /// vsearch's `search16_qprep`.
    pub fn set_query(&mut self, query: &[u8]) {
        match &mut self.kernel {
            Kernel::Scalar(e) => e.set_query(query),
            #[cfg(target_arch = "x86_64")]
            Kernel::Sse2(e) => e.set_query(query),
            #[cfg(target_arch = "x86_64")]
            Kernel::Avx2 { wide, narrow } => {
                narrow.set_query(query);
                // SAFETY: the kernel is AVX2 only when `simd::level` found
                // the processor has it
                unsafe { set_query_avx2(wide, query) };
            }
            #[cfg(target_arch = "aarch64")]
            Kernel::Neon(e) => e.set_query(query),
        }
    }

    /// the targets aligned at once in the widest calls.
    pub fn lanes(&self) -> usize {
        match &self.kernel {
            Kernel::Scalar(_) => Scalar::N,
            #[cfg(target_arch = "x86_64")]
            Kernel::Sse2(_) => Sse2::N,
            #[cfg(target_arch = "x86_64")]
            Kernel::Avx2 { .. } => Avx2::N,
            #[cfg(target_arch = "aarch64")]
            Kernel::Neon(_) => Neon::N,
        }
    }

    /// aligns the query with each of `targets`, as vsearch's `search16`.
    pub fn align(&mut self, targets: &[&[u8]]) -> Vec<Aligned> {
        match &mut self.kernel {
            Kernel::Scalar(e) => e.align(targets),
            #[cfg(target_arch = "x86_64")]
            Kernel::Sse2(e) => e.align(targets),
            #[cfg(target_arch = "x86_64")]
            Kernel::Avx2 { narrow, .. } if targets.len() <= Sse2::N => narrow.align(targets),
            // SAFETY: the kernel is AVX2 only when `simd::level` found the
            // processor has it
            #[cfg(target_arch = "x86_64")]
            Kernel::Avx2 { wide, .. } => unsafe { align_avx2(wide, targets) },
            #[cfg(target_arch = "aarch64")]
            Kernel::Neon(e) => e.align(targets),
        }
    }
}

/// [`Engine::set_query`] compiled with AVX2, so that the lanes'
/// instructions inline into it.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
fn set_query_avx2(engine: &mut Engine<Avx2>, query: &[u8]) {
    engine.set_query(query)
}

/// [`Engine::align`] compiled with AVX2, so that the lanes' instructions
/// inline into it.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
fn align_avx2(engine: &mut Engine<Avx2>, targets: &[&[u8]]) -> Vec<Aligned> {
    engine.align(targets)
}

/// vsearch's `s16info_s`: the scores and penalties in 16-bit cells, and
/// the aligner's buffers, on lanes `L`.
struct Engine<L: Lanes> {
    /// The score of each pair of 4-bit codes.
    matrix: [[i16; 16]; 16],
    /// Penalties (opening already less the extension, as vsearch's):
    /// [query, target] x [left, interior, right].
    open: [[i16; 3]; 2],
    extend: [[i16; 3]; 2],
    n_mismatch: bool,
    /// A score or penalty beyond the cells: every pair to the fallback.
    fallback: bool,
    hearray: Vec<L>,
    dir: Vec<L::Bits>,
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

impl<L: Lanes> Engine<L> {
    fn new(scoring: &Scoring) -> Engine<L> {
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
        Engine {
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
    #[inline(always)]
    fn set_query(&mut self, query: &[u8]) {
        self.qseq = query.to_vec();
        self.qcodes = query.iter().map(|&c| FOUR_BIT[usize::from(c)]).collect();
        self.hearray.clear();
        self.hearray.resize(2 * query.len(), zero());
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
    #[inline(always)]
    fn align(&mut self, targets: &[&[u8]]) -> Vec<Aligned> {
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
            self.dir.resize(dirsize, L::Bits::default());
        }

        let pen = |x: i16| dup::<L>(x);
        let qr = |s: usize, e: usize| dup::<L>(self.open[s][e].wrapping_add(self.extend[s][e]));
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
        let mut begin = [0usize; MAX_LANES];
        let mut end = [0usize; MAX_LANES];
        let mut offset = [0usize; MAX_LANES];
        let mut seq_id = [-1i64; MAX_LANES];
        let mut overflow = [false; MAX_LANES];
        let mut dseq = [[0u8; MAX_LANES]; CDEPTH];
        let mut s_vec: [L; 4] = [zero(); 4];
        let score_min = self.score_min();
        let score_max = i16::MAX;
        let (mut h0, mut h1, mut h2, mut h3): (L, L, L, L) = (zero(), zero(), zero(), zero());
        let (mut f0, mut f1, mut f2, mut f3): (L, L, L, L) = (zero(), zero(), zero(), zero());
        let mut easy = false;
        let mut dir = 0usize;
        let mut next_id = 0usize;
        let mut done = 0usize;
        let go_q_l = i64::from(self.open[q][left]);
        let ge_q_l = i64::from(self.extend[q][left]);

        loop {
            let mut masked: Option<L> = None;
            if easy {
                for c in 0..L::N {
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
                let mut m: L = zero();
                for c in 0..L::N {
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
                        m = xor(m, L::lane(c));
                        let cand = seq_id[c];
                        if cand >= 0 {
                            let slot = cand as usize;
                            let dlen = targets[slot].len();
                            let z = (dlen + 3) % 4;
                            if overflow[c] {
                                out[slot] = Aligned::overflow();
                            } else {
                                let score = s_vec[z].store()[c];
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
                            set(&mut h0, c, 0);
                            let edge = |k: i64| (-go_q_l - k * ge_q_l) as i16;
                            set(&mut h1, c, edge(1));
                            set(&mut h2, c, edge(2));
                            set(&mut h3, c, edge(3));
                            set(&mut f0, c, edge(1));
                            set(&mut f1, c, edge(2));
                            set(&mut f2, c, edge(3));
                            set(&mut f3, c, edge(4));
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
                }
                if done == n {
                    break;
                }
                masked = Some(m);
            }

            // the profile: the score of each query code against each target
            // residue of the step
            let mut profile: [[L; CDEPTH]; 16] = [[zero(); CDEPTH]; 16];
            for (r, rows) in profile.iter_mut().enumerate() {
                for (j, lanes) in rows.iter_mut().enumerate() {
                    *lanes = L::load(&std::array::from_fn(|c| {
                        self.matrix[usize::from(dseq[j][c])][r]
                    }));
                }
            }
            let mut qr_target = [qr_target_interior; 4];
            let mut r_target = [r_target_interior; 4];
            if !easy {
                let qr_diff = sub(qr_target_right, qr_target_interior);
                let r_diff = sub(r_target_right, r_target_interior);
                for j in 0..CDEPTH {
                    let mut mm: L = zero();
                    for c in 0..L::N {
                        let dlen = if seq_id[c] >= 0 { end[c] } else { 0 };
                        if begin[c] == end[c] && j >= (dlen + 3) % 4 {
                            mm = xor(mm, L::lane(c));
                        }
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
            let (h_min, h_max) = (h_min.store(), h_max.store());
            for c in 0..L::N {
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
    /// Inlined, as everything `align` calls in its loop is. With AVX2 the
    /// lanes' instructions must land in `align_avx2`, which has the
    /// feature. A function without it would call each one.
    #[allow(clippy::too_many_arguments)]
    #[inline(always)]
    fn align_columns(
        &mut self,
        s_vec: &mut [L; 4],
        profile: &[[L; CDEPTH]; 16],
        query: [L; 4],
        qr_t: [L; 4],
        r_t: [L; 4],
        h: [L; 4],
        f: [L; 4],
        mask: Option<Masked<L>>,
        dir: usize,
    ) -> (L, L) {
        match mask {
            Some(mask) => self.columns::<true>(s_vec, profile, query, qr_t, r_t, h, f, mask, dir),
            None => {
                let none = Masked {
                    mm: zero(),
                    qr_t_left: zero(),
                    r_t_left: zero(),
                    qr_q_interior: zero(),
                    qr_q_right: zero(),
                };
                self.columns::<false>(s_vec, profile, query, qr_t, r_t, h, f, none, dir)
            }
        }
    }

    /// the columns, with the channels starting a target masked or not, as
    /// vsearch's `Masked` and `Unmasked`.
    #[allow(clippy::too_many_arguments)]
    #[inline(always)]
    fn columns<const MASKED: bool>(
        &mut self,
        s_vec: &mut [L; 4],
        profile: &[[L; CDEPTH]; 16],
        query: [L; 4],
        qr_t: [L; 4],
        r_t: [L; 4],
        h: [L; 4],
        f: [L; 4],
        mut m: Masked<L>,
        dir: usize,
    ) -> (L, L) {
        let [qr_q_i, r_q_i, qr_q_r, r_q_r] = query;
        let ql = self.qcodes.len();
        let [mut h0, mut h1, mut h2, mut h3] = h;
        let mut f0 = sub(f[0], qr_t[0]);
        let mut f1 = sub(f[1], qr_t[1]);
        let mut f2 = sub(f[2], qr_t[2]);
        let mut f3 = sub(f[3], qr_t[3]);
        let mut h_min = zero();
        let mut h_max = zero();
        let (mut h5, mut h6, mut h7, mut h8) = (zero(), zero(), zero(), zero());
        let rows = self.dir[dir..dir + 16 * ql]
            .as_chunks_mut::<16>()
            .0
            .iter_mut();
        let cells = self.hearray.as_chunks_mut::<2>().0.iter_mut();
        for (i, ((paths, he), &code)) in rows.zip(cells).zip(&self.qcodes).enumerate() {
            let last = i + 1 == ql;
            let vp = &profile[usize::from(code)];
            if !last {
                let mut h4 = he[0];
                let mut e = he[1];
                if MASKED {
                    h4 = sub_unsigned(h4, m.mm);
                    h4 = sub(h4, m.qr_t_left);
                    e = sub_unsigned(e, m.mm);
                    e = sub(e, m.qr_t_left);
                    e = sub(e, m.qr_q_interior);
                    m.qr_t_left = add(m.qr_t_left, m.r_t_left);
                }
                let (qr_q, r_q) = (qr_q_i, r_q_i);
                step(
                    &mut h0,
                    &mut h5,
                    &mut f0,
                    vp[0],
                    &mut paths[0..4],
                    &mut e,
                    qr_q,
                    r_q,
                    qr_t[0],
                    r_t[0],
                    &mut h_min,
                    &mut h_max,
                );
                step(
                    &mut h1,
                    &mut h6,
                    &mut f1,
                    vp[1],
                    &mut paths[4..8],
                    &mut e,
                    qr_q,
                    r_q,
                    qr_t[1],
                    r_t[1],
                    &mut h_min,
                    &mut h_max,
                );
                step(
                    &mut h2,
                    &mut h7,
                    &mut f2,
                    vp[2],
                    &mut paths[8..12],
                    &mut e,
                    qr_q,
                    r_q,
                    qr_t[2],
                    r_t[2],
                    &mut h_min,
                    &mut h_max,
                );
                step(
                    &mut h3,
                    &mut h8,
                    &mut f3,
                    vp[3],
                    &mut paths[12..16],
                    &mut e,
                    qr_q,
                    r_q,
                    qr_t[3],
                    r_t[3],
                    &mut h_min,
                    &mut h_max,
                );
                he[0] = h8;
                he[1] = e;
                h0 = h4;
                h1 = h5;
                h2 = h6;
                h3 = h7;
            } else {
                let mut e = he[1];
                if MASKED {
                    e = sub_unsigned(e, m.mm);
                    e = sub(e, m.qr_t_left);
                    e = sub(e, m.qr_q_right);
                }
                let (qr_q, r_q) = (qr_q_r, r_q_r);
                step(
                    &mut h0,
                    &mut h5,
                    &mut f0,
                    vp[0],
                    &mut paths[0..4],
                    &mut e,
                    qr_q,
                    r_q,
                    qr_t[0],
                    r_t[0],
                    &mut h_min,
                    &mut h_max,
                );
                step(
                    &mut h1,
                    &mut h6,
                    &mut f1,
                    vp[1],
                    &mut paths[4..8],
                    &mut e,
                    qr_q,
                    r_q,
                    qr_t[1],
                    r_t[1],
                    &mut h_min,
                    &mut h_max,
                );
                step(
                    &mut h2,
                    &mut h7,
                    &mut f2,
                    vp[2],
                    &mut paths[8..12],
                    &mut e,
                    qr_q,
                    r_q,
                    qr_t[2],
                    r_t[2],
                    &mut h_min,
                    &mut h_max,
                );
                step(
                    &mut h3,
                    &mut h8,
                    &mut f3,
                    vp[3],
                    &mut paths[12..16],
                    &mut e,
                    qr_q,
                    r_q,
                    qr_t[3],
                    r_t[3],
                    &mut h_min,
                    &mut h_max,
                );
                he[0] = h8;
                he[1] = e;
            }
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
        let mut ops: Vec<(u8, usize)> = Vec::with_capacity(64);
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
            let bit = |w: L::Bits| (L::bits(w) >> shift) & 3 != 0;
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
        let mut cigar = String::with_capacity(6 * ops.len());
        for &(op, run) in ops.iter().rev() {
            if run > 1 {
                let _ = write!(cigar, "{run}");
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

/// vsearch's `onestep`: one cell of four target columns, its direction
/// bits in `path`.
#[allow(clippy::too_many_arguments)]
#[inline(always)]
fn step<L: Lanes>(
    h: &mut L,
    n: &mut L,
    f: &mut L,
    v: L,
    path: &mut [L::Bits],
    e: &mut L,
    qr_q: L,
    r_q: L,
    qr_t: L,
    r_t: L,
    h_min: &mut L,
    h_max: &mut L,
) {
    *h = add(*h, v);
    path[0] = mask_gt(*f, *h);
    *h = max(*h, *f);
    path[1] = mask_gt(*e, *h);
    *h = max(*h, *e);
    *h_min = min(*h_min, *h);
    *h_max = max(*h_max, *h);
    *n = *h;
    let hf = sub(*h, qr_t);
    *f = sub(*f, r_t);
    path[2] = mask_gt(*f, hf);
    *f = max(*f, hf);
    let he = sub(*h, qr_q);
    *e = sub(*e, r_q);
    path[3] = mask_gt(*e, he);
    *e = max(*e, he);
}

/// The masks of vsearch's `Masked` columns: the channels starting a target.
struct Masked<L: Lanes> {
    mm: L,
    qr_t_left: L,
    r_t_left: L,
    qr_q_interior: L,
    qr_q_right: L,
}
