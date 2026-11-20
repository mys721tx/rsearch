//! vsearch's UCHIME chimera detection (`chimera.cpp`), ported step for step
//! so that scores, verdicts, tables and alignments are vsearch's. The four
//! parts of a query are searched against the parents, and the hits are
//! aligned with the whole query. The two candidates winning the most
//! windows are its parents. The votes of the columns where the query
//! follows one parent and not the other give its score.

use std::fmt::Write;

use crate::accept::Criteria;
use crate::align16::Search16;
use crate::alignment::parse_cigar;
use crate::linmem::LinearMemoryAligner;
use crate::scoring::{FOUR_BIT, Scoring, is_ambiguous};
use crate::search::{Database, Limits, MIN_WORD_MATCHES, Searcher, order};

/// The parts a query is searched in.
const PARTS: usize = 4;
/// The hits kept of all parts.
const MAX_CANDIDATES: usize = 400;
/// The window over which matches are counted.
const WINDOW: usize = 32;

/// UCHIME's options.
#[derive(Clone, Copy)]
pub struct Options {
    /// UCHIME2's verdict: a chimera when its model matches it perfectly.
    pub uchime2: bool,
    /// de novo: the parents are at least this much more abundant, and
    /// neither the query nor identical to it.
    pub abskew: Option<f64>,
    pub xn: f64,
    pub dn: f64,
    pub minh: f64,
    pub mindiv: f64,
    pub mindiffs: i32,
    /// the width of the alignments' rows, 0 for the whole.
    pub width: usize,
}

/// A query's verdict, in vsearch's order.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Status {
    NoParents,
    NoAlignment,
    LowScore,
    Suspicious,
    Chimeric,
}

/// A query's verdict, its line of the table, and its alignment if it is a
/// chimera.
pub struct Verdict {
    pub status: Status,
    pub table: String,
    pub alignment: Option<String>,
}

/// The four-bit code of a residue, 0 for a gap.
fn code(c: u8) -> u8 {
    FOUR_BIT[usize::from(c)]
}

/// The chimera detector.
pub struct Uchime {
    options: Options,
    searcher: Searcher,
    simd: Search16,
    fallback: LinearMemoryAligner,
}

impl Uchime {
    /// a detector searching `db`, which de novo grows with the non-chimeras.
    pub fn new(db: &Database, options: Options) -> Uchime {
        let scoring = Scoring::default();
        let mut criteria = Criteria {
            id: 0.55,
            weak_id: 0.55,
            ..Criteria::default()
        };
        if let Some(abskew) = options.abskew {
            criteria.self_label = true;
            criteria.self_id = true;
            criteria.sizes.max_size_ratio = 1.0 / abskew;
        }
        let limits = Limits {
            max_accepts: 4,
            max_rejects: 16,
            // vsearch takes the default of its word length, 8
            min_word_matches: MIN_WORD_MATCHES[8],
            top_hits: Some(20),
        };
        Uchime {
            options,
            searcher: Searcher::new(db, criteria, scoring, limits),
            simd: Search16::new(&scoring),
            fallback: LinearMemoryAligner::new(&scoring),
        }
    }

    /// the candidate parents of a query: the accepted hits of its parts,
    /// best first, each once.
    fn candidates(&mut self, db: &Database, label: &str, size: i64, query: &[u8]) -> Vec<usize> {
        let mut candidates: Vec<usize> = Vec::new();
        if query.len() < PARTS {
            return candidates;
        }
        let mut hits_count = 0;
        let (mut rest, mut start) = (query.len(), 0);
        for i in 0..PARTS {
            let len = rest.div_ceil(PARTS - i);
            let part = &query[start..start + len];
            let mut hits = self.searcher.search_part(db, label, size, part, query);
            order(&mut hits, 2);
            for h in hits.iter().filter(|h| !h.weak) {
                if hits_count < MAX_CANDIDATES {
                    hits_count += 1;
                    if !candidates.contains(&h.target) {
                        candidates.push(h.target);
                    }
                }
            }
            rest -= len;
            start += len;
        }
        candidates
    }

    /// judges `query` against the parents of `db`.
    pub fn evaluate(&mut self, db: &Database, label: &str, size: i64, query: &[u8]) -> Verdict {
        let candidates = self.candidates(db, label, size, query);
        let mut cigars = Vec::with_capacity(candidates.len());
        if !candidates.is_empty() {
            self.simd.set_query(query);
            let targets: Vec<&[u8]> = candidates.iter().map(|&t| db.seqs[t].as_slice()).collect();
            for (&t, aligned) in candidates.iter().zip(self.simd.align(&targets)) {
                cigars.push(if aligned.score == i16::MAX {
                    self.fallback.align(query, &db.seqs[t])
                } else {
                    aligned.cigar
                });
            }
        }
        let ops: Vec<Vec<(u8, usize)>> = cigars
            .iter()
            .map(|c| parse_cigar(c).expect("the aligners' CIGARs are valid"))
            .collect();
        match best_parents(db, query, &candidates, &ops) {
            Some(parents) => self.eval_parents(db, label, query, &candidates, &ops, parents),
            None => Verdict {
                status: Status::NoParents,
                table: star_row(label),
                alignment: None,
            },
        }
    }

    /// vsearch's `eval_parents`: the votes and the score of the query
    /// against its two parents.
    fn eval_parents(
        &self,
        db: &Database,
        label: &str,
        query: &[u8],
        candidates: &[usize],
        ops: &[Vec<(u8, usize)>],
        parents: [usize; 2],
    ) -> Verdict {
        let o = &self.options;
        let qlen = query.len();
        // the longest insertion of a parent before each query position
        let mut maxi = vec![0usize; qlen + 1];
        for &p in &parents {
            let mut pos = 0;
            for &(op, run) in &ops[p] {
                match op {
                    b'I' => maxi[pos] = maxi[pos].max(run),
                    _ => pos += run,
                }
            }
        }
        let alnlen = qlen + maxi.iter().sum::<usize>();
        let mut paln: [Vec<u8>; 2] = [Vec::with_capacity(alnlen), Vec::with_capacity(alnlen)];
        for (row, &p) in paln.iter_mut().zip(&parents) {
            let target = &db.seqs[candidates[p]];
            let (mut qpos, mut tpos) = (0usize, 0usize);
            let mut inserted = false;
            for &(op, run) in &ops[p] {
                if op == b'I' {
                    for j in 0..maxi[qpos] {
                        if j < run {
                            row.push(target[tpos].to_ascii_uppercase());
                            tpos += 1;
                        } else {
                            row.push(b'-');
                        }
                    }
                    inserted = true;
                } else {
                    for _ in 0..run {
                        if !inserted {
                            row.extend(std::iter::repeat_n(b'-', maxi[qpos]));
                        }
                        if op == b'M' {
                            row.push(target[tpos].to_ascii_uppercase());
                            tpos += 1;
                        } else {
                            row.push(b'-');
                        }
                        qpos += 1;
                        inserted = false;
                    }
                }
            }
            if !inserted {
                row.extend(std::iter::repeat_n(b'-', maxi[qpos]));
            }
        }
        let mut qaln = Vec::with_capacity(alnlen);
        for (i, &c) in query.iter().enumerate() {
            qaln.extend(std::iter::repeat_n(b'-', maxi[i]));
            qaln.push(c.to_ascii_uppercase());
        }
        qaln.extend(std::iter::repeat_n(b'-', maxi[qlen]));

        let mut ignore = vec![false; alnlen];
        let mut diffs = vec![b' '; alnlen];
        for i in 0..alnlen {
            let (qsym, p1sym, p2sym) = (code(qaln[i]), code(paln[0][i]), code(paln[1][i]));
            if qsym == 0 || p1sym == 0 || p2sym == 0 {
                ignore[i] = true;
                if i > 0 {
                    ignore[i - 1] = true;
                }
                if i + 1 < alnlen {
                    ignore[i + 1] = true;
                }
            }
            if is_ambiguous(qsym) || is_ambiguous(p1sym) || is_ambiguous(p2sym) {
                ignore[i] = true;
            }
            if p1sym != 0 && p1sym != qsym {
                paln[0][i] = paln[0][i].to_ascii_lowercase();
            }
            if p2sym != 0 && p2sym != qsym {
                paln[1][i] = paln[1][i].to_ascii_lowercase();
            }
            diffs[i] = if qsym != 0 && p1sym != 0 && p2sym != 0 {
                if p1sym == p2sym {
                    if qsym == p1sym { b' ' } else { b'N' }
                } else if qsym == p1sym {
                    b'A'
                } else if qsym == p2sym {
                    b'B'
                } else {
                    b'?'
                }
            } else {
                b' '
            };
        }

        let (mut sum_a, mut sum_b, mut sum_n) = (0i32, 0i32, 0i32);
        for i in 0..alnlen {
            if ignore[i] {
                continue;
            }
            match diffs[i] {
                b'A' => sum_a += 1,
                b'B' => sum_b += 1,
                b' ' => {}
                _ => sum_n += 1,
            }
        }
        let (mut left_n, mut left_a, mut left_y) = (0i32, 0i32, 0i32);
        let (mut right_n, mut right_a, mut right_y) = (sum_a, sum_n, sum_b);
        let mut best_h = -1.0f64;
        let mut best_i: i64 = -1;
        let mut reverse = false;
        let mut best = [0i32; 6]; // left n, a, y; right n, a, y
        for i in 0..alnlen {
            if ignore[i] || diffs[i] == b' ' {
                continue;
            }
            match diffs[i] {
                b'A' => {
                    left_y += 1;
                    right_n -= 1;
                }
                b'B' => {
                    left_n += 1;
                    right_y -= 1;
                }
                _ => {
                    left_a += 1;
                    right_a -= 1;
                }
            }
            let score = |yes: i32, no: i32, abstain: i32| {
                f64::from(yes) / ((o.xn * (f64::from(no) + o.dn)) + f64::from(abstain))
            };
            if left_y > left_n && right_y > right_n {
                let h = score(left_y, left_n, left_a) * score(right_y, right_n, right_a);
                if h > best_h {
                    reverse = false;
                    best_h = h;
                    best_i = i as i64;
                    best = [left_n, left_a, left_y, right_n, right_a, right_y];
                }
            } else if left_n > left_y && right_n > right_y {
                let h = score(left_n, left_y, left_a) * score(right_n, right_y, right_a);
                if h > best_h {
                    reverse = true;
                    best_h = h;
                    best_i = i as i64;
                    best = [left_y, left_a, left_n, right_y, right_a, right_n];
                }
            }
        }
        if best_h < 0.0 {
            return Verdict {
                status: Status::NoAlignment,
                table: star_row(label),
                alignment: None,
            };
        }

        let mut status = Status::LowScore;
        if reverse {
            for d in diffs.iter_mut() {
                *d = match *d {
                    b'A' => b'B',
                    b'B' => b'A',
                    d => d,
                };
            }
        }
        let mut model = vec![b' '; alnlen];
        let mut votes = vec![b' '; alnlen];
        for i in 0..alnlen {
            let m = if (i as i64) <= best_i { b'A' } else { b'B' };
            model[i] = m;
            let mut v = b' ';
            if !ignore[i] {
                let d = diffs[i];
                if d == b'A' || d == b'B' {
                    v = if d == m { b'+' } else { b'!' };
                } else if d == b'N' || d == b'?' {
                    v = b'0';
                }
            }
            votes[i] = v;
            if v == b'!' {
                diffs[i] = diffs[i].to_ascii_lowercase();
            }
        }
        for i in (best_i + 1) as usize..alnlen {
            if diffs[i] == b' ' || diffs[i] == b'A' {
                model[i] = b'x';
            } else {
                break;
            }
        }

        let (index_a, index_b) = if reverse { (1, 0) } else { (0, 1) };
        let (mut match_qa, mut match_qb, mut match_ab, mut match_qm, mut cols) = (0, 0, 0, 0, 0);
        for i in 0..alnlen {
            if ignore[i] {
                continue;
            }
            cols += 1;
            let qsym = code(qaln[i]);
            let asym = code(paln[index_a][i]);
            let bsym = code(paln[index_b][i]);
            let msym = if (i as i64) <= best_i { asym } else { bsym };
            match_qa += i32::from(qsym == asym);
            match_qb += i32::from(qsym == bsym);
            match_ab += i32::from(asym == bsym);
            match_qm += i32::from(qsym == msym);
        }
        let seqno_a = candidates[parents[index_a]];
        let seqno_b = candidates[parents[index_b]];
        let pct = |n: i32| 100.0 * f64::from(n) / f64::from(cols);
        let (qa, qb, ab, qm) = (pct(match_qa), pct(match_qb), pct(match_ab), pct(match_qm));
        let qt = qa.max(qb);
        let divdiff = qm - qt;
        let divfrac = 100.0 * divdiff / qt;
        let [left_n, left_a, left_y, right_n, right_a, right_y] = best;
        let sum_l = left_n + left_a + left_y;
        let sum_r = right_n + right_a + right_y;
        if o.uchime2 {
            if match_qm == cols && qt < 100.0 {
                status = Status::Chimeric;
            }
        } else if best_h >= o.minh {
            status = Status::Suspicious;
            if divdiff >= o.mindiv && sum_l >= o.mindiffs && sum_r >= o.mindiffs {
                status = Status::Chimeric;
            }
        }

        let (label_a, label_b) = (&db.labels[seqno_a], &db.labels[seqno_b]);
        let mut alignment = None;
        if status == Status::Chimeric {
            let mut text = String::new();
            text.push('\n');
            text.push_str(&"-".repeat(72));
            text.push('\n');
            let _ = writeln!(text, "Query   ({qlen:5} nt) {label}");
            let _ = writeln!(text, "ParentA ({:5} nt) {label_a}", db.seqs[seqno_a].len());
            let _ = writeln!(
                text,
                "ParentB ({:5} nt) {label_b}\n",
                db.seqs[seqno_b].len()
            );
            let width = if o.width > 0 { o.width } else { alnlen };
            let (mut qpos, mut p1pos, mut p2pos) = (0usize, 0usize, 0usize);
            let mut i = 0;
            while i < alnlen {
                let w = width.min(alnlen - i);
                let residues = |row: &[u8]| row[i..i + w].iter().filter(|&&c| c != b'-').count();
                let (qnt, p1nt, p2nt) = (residues(&qaln), residues(&paln[0]), residues(&paln[1]));
                let row = |text: &mut String, name: char, start: usize, row: &[u8], n: usize| {
                    let _ = writeln!(
                        text,
                        "{name} {:5} {} {}",
                        start + 1,
                        String::from_utf8_lossy(&row[i..i + w]),
                        start + n
                    );
                };
                if !reverse {
                    row(&mut text, 'A', p1pos, &paln[0], p1nt);
                    row(&mut text, 'Q', qpos, &qaln, qnt);
                    row(&mut text, 'B', p2pos, &paln[1], p2nt);
                } else {
                    row(&mut text, 'A', p2pos, &paln[1], p2nt);
                    row(&mut text, 'Q', qpos, &qaln, qnt);
                    row(&mut text, 'B', p1pos, &paln[0], p1nt);
                }
                for (name, annotation) in [
                    ("Diffs   ", &diffs),
                    ("Votes   ", &votes),
                    ("Model   ", &model),
                ] {
                    text.push_str(name);
                    text.push_str(&String::from_utf8_lossy(&annotation[i..i + w]));
                    text.push('\n');
                }
                text.push('\n');
                qpos += qnt;
                p1pos += p1nt;
                p2pos += p2nt;
                i += width;
            }
            let _ = writeln!(
                text,
                "Ids.  QA {qa:.1}%, QB {qb:.1}%, AB {ab:.1}%, QModel {qm:.1}%, Div. {divfrac:+.1}%"
            );
            let _ = writeln!(
                text,
                "Diffs Left {sum_l}: N {left_n}, A {left_a}, Y {left_y} ({:.1}%); Right {sum_r}: N \
                 {right_n}, A {right_a}, Y {right_y} ({:.1}%), Score {best_h:.4}",
                100.0 * f64::from(left_y) / f64::from(sum_l),
                100.0 * f64::from(right_y) / f64::from(sum_r)
            );
            alignment = Some(text);
        }
        let closest = if qa >= qb { label_a } else { label_b };
        let flag = match status {
            Status::Chimeric => 'Y',
            Status::LowScore => 'N',
            _ => '?',
        };
        let table = format!(
            "{best_h:.4}\t{label}\t{label_a}\t{label_b}\t{closest}\t{qm:.1}\t{qa:.1}\t{qb:.1}\t\
             {ab:.1}\t{qt:.1}\t{left_y}\t{left_n}\t{left_a}\t{right_y}\t{right_n}\t{right_a}\t\
             {divdiff:.1}\t{flag}\n"
        );
        Verdict {
            status,
            table,
            alignment,
        }
    }
}

/// the table line of a query without parents.
fn star_row(label: &str) -> String {
    format!("0.0000\t{label}\t*\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0\t0\t0\t0\t*\tN\n")
}

/// vsearch's `find_best_parents`: the two candidates whose matches with
/// the query, counted over windows of 32, are the most often the best, the
/// second once the windows the first wins are wiped.
// indexed as vsearch indexes its arrays
#[allow(clippy::needless_range_loop)]
fn best_parents(
    db: &Database,
    query: &[u8],
    candidates: &[usize],
    ops: &[Vec<(u8, usize)>],
) -> Option<[usize; 2]> {
    let qlen = query.len();
    let n = candidates.len();
    let mut matches = vec![0i32; n * qlen];
    for (i, cand_ops) in ops.iter().enumerate() {
        let target = &db.seqs[candidates[i]];
        let (mut qpos, mut tpos) = (0usize, 0usize);
        for &(op, run) in cand_ops {
            match op {
                b'M' => {
                    for _ in 0..run {
                        if code(query[qpos]) & code(target[tpos]) != 0 {
                            matches[i * qlen + qpos] = 1;
                        }
                        qpos += 1;
                        tpos += 1;
                    }
                }
                b'I' => tpos += run,
                _ => qpos += run,
            }
        }
    }
    let mut smooth = vec![0i32; n * qlen];
    let mut maxsmooth = vec![0i32; qlen];
    let mut selected = vec![false; n];
    let mut best: [Option<usize>; 2] = [None, None];
    for f in 0..2 {
        if f > 0 {
            let prev = best[f - 1].unwrap();
            let mut wiped_upto = 0;
            for qpos in WINDOW.saturating_sub(1)..qlen {
                if smooth[prev * qlen + qpos] == maxsmooth[qpos] {
                    let first = (qpos + 1).saturating_sub(WINDOW).max(wiped_upto);
                    for j in 0..n {
                        matches[j * qlen + first..j * qlen + qpos + 1].fill(0);
                    }
                    wiped_upto = qpos + 1;
                }
            }
        }
        maxsmooth.fill(0);
        for (i, &chosen) in selected.iter().enumerate() {
            if chosen {
                continue;
            }
            let mut sum = 0;
            for qpos in 0..qlen {
                let z = i * qlen + qpos;
                sum += matches[z];
                if qpos >= WINDOW {
                    sum -= matches[z - WINDOW];
                }
                if qpos + 1 >= WINDOW {
                    smooth[z] = sum;
                    maxsmooth[qpos] = maxsmooth[qpos].max(sum);
                }
            }
        }
        let mut wins = vec![0i32; n];
        for qpos in WINDOW.saturating_sub(1)..qlen {
            if maxsmooth[qpos] != 0 {
                for i in 0..n {
                    if !selected[i] && smooth[i * qlen + qpos] == maxsmooth[qpos] {
                        wins[i] += 1;
                    }
                }
            }
        }
        let mut maxwins = 0;
        for (i, &w) in wins.iter().enumerate() {
            if w > maxwins {
                maxwins = w;
                best[f] = Some(i);
            }
        }
        match best[f] {
            Some(b) => selected[b] = true,
            None => break,
        }
    }
    match best {
        [Some(a), Some(b)] => Some([a, b]),
        _ => None,
    }
}
