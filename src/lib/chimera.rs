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

/// The search for a query's candidate parents, and their alignments with
/// it.
struct Finder {
    searcher: Searcher,
    simd: Search16,
    fallback: LinearMemoryAligner,
}

impl Finder {
    /// a finder searching `db`, which de novo grows with the non-chimeras.
    /// De novo, the parents are at least `abskew` times as abundant as the
    /// query, and neither the query nor identical to it.
    fn new(db: &Database, abskew: Option<f64>) -> Finder {
        let scoring = Scoring::default();
        let mut criteria = Criteria {
            id: 0.55,
            weak_id: 0.55,
            ..Criteria::default()
        };
        if let Some(abskew) = abskew {
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
        Finder {
            searcher: Searcher::new(db, criteria, scoring, limits),
            simd: Search16::new(&scoring),
            fallback: LinearMemoryAligner::new(&scoring),
        }
    }

    /// the candidate parents of a query searched in `parts`: the accepted
    /// hits of its parts, best first, each once.
    fn candidates(
        &mut self,
        db: &Database,
        label: &str,
        size: i64,
        query: &[u8],
        parts: usize,
    ) -> Vec<usize> {
        let mut candidates: Vec<usize> = Vec::new();
        if query.len() < parts {
            return candidates;
        }
        let mut hits_count = 0;
        let (mut rest, mut start) = (query.len(), 0);
        for i in 0..parts {
            let len = rest.div_ceil(parts - i);
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

    /// the alignments of `query` with the candidates, as operations.
    fn align(
        &mut self,
        db: &Database,
        query: &[u8],
        candidates: &[usize],
    ) -> Vec<Vec<(u8, usize)>> {
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
        cigars
            .iter()
            .map(|c| parse_cigar(c).expect("the aligners' CIGARs are valid"))
            .collect()
    }
}

/// The longest insertion of the parents before each query position.
fn insertions(qlen: usize, ops: &[Vec<(u8, usize)>], parents: &[usize]) -> Vec<usize> {
    let mut maxi = vec![0usize; qlen + 1];
    for &p in parents {
        let mut pos = 0;
        for &(op, run) in &ops[p] {
            match op {
                b'I' => maxi[pos] = maxi[pos].max(run),
                _ => pos += run,
            }
        }
    }
    maxi
}

/// the row of a parent in the alignment, upper case, gapped by `maxi`.
fn parent_row(target: &[u8], ops: &[(u8, usize)], maxi: &[usize]) -> Vec<u8> {
    let mut row = Vec::new();
    let (mut qpos, mut tpos) = (0usize, 0usize);
    let mut inserted = false;
    for &(op, run) in ops {
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
    row
}

/// the row of the query in the alignment, upper case, gapped by `maxi`.
fn query_row(query: &[u8], maxi: &[usize]) -> Vec<u8> {
    let mut row = Vec::new();
    for (i, &c) in query.iter().enumerate() {
        row.extend(std::iter::repeat_n(b'-', maxi[i]));
        row.push(c.to_ascii_uppercase());
    }
    row.extend(std::iter::repeat_n(b'-', maxi[query.len()]));
    row
}

/// The chimera detector.
pub struct Uchime {
    options: Options,
    finder: Finder,
}

impl Uchime {
    /// a detector searching `db`, which de novo grows with the non-chimeras.
    pub fn new(db: &Database, options: Options) -> Uchime {
        Uchime {
            options,
            finder: Finder::new(db, options.abskew),
        }
    }

    /// judges `query` against the parents of `db`.
    pub fn evaluate(&mut self, db: &Database, label: &str, size: i64, query: &[u8]) -> Verdict {
        let candidates = self.finder.candidates(db, label, size, query, PARTS);
        let ops = self.finder.align(db, query, &candidates);
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
        let maxi = insertions(qlen, ops, &parents);
        let alnlen = qlen + maxi.iter().sum::<usize>();
        let mut paln: [Vec<u8>; 2] =
            parents.map(|p| parent_row(&db.seqs[candidates[p]], &ops[p], &maxi));
        let qaln = query_row(query, &maxi);

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

/// The most parents vsearch gives a multi-parent chimera.
pub const MAX_PARENTS: usize = 20;

/// The options of multi-parent chimera detection.
#[derive(Clone, Copy)]
pub struct LongOptions {
    /// the parts a query is searched in, by default one a hundred bases
    pub parts: Option<usize>,
    pub max_parents: usize,
    /// the shortest segment of a parent
    pub min_segment: usize,
    /// the share of differences a segment may hold, in percent
    pub diff_pct: f64,
    pub abskew: f64,
    /// the width of the alignments' rows, 0 for the whole.
    pub width: usize,
}

/// A parent of a multi-parent chimera: a candidate and the query segment
/// it covers.
#[derive(Clone, Copy)]
struct Segment {
    cand: usize,
    start: usize,
    len: usize,
}

/// vsearch's `scan_matches`: the longest stretch of `matches` scoring at
/// least 0, a match `percentage` and a difference `percentage - 100`.
fn scan_matches(matches: &[i32], percentage: f64) -> Option<(usize, usize)> {
    let len = matches.len();
    let (score_match, score_mismatch) = (percentage, percentage - 100.0);
    let mut p = vec![0.0f64; len + 1];
    for i in 0..len {
        p[i + 1] = p[i]
            + if matches[i] != 0 {
                score_match
            } else {
                score_mismatch
            };
    }
    let mut q = vec![0.0f64; len + 1];
    q[len] = p[len];
    for i in (0..len).rev() {
        q[i] = q[i + 1].max(p[i]);
    }
    let (mut best_i, mut best_d, mut best_c) = (0usize, -1i64, -1.0f64);
    let (mut i, mut j) = (1usize, 1usize);
    while j <= len {
        let c = q[j] - p[i - 1];
        if c >= 0.0 {
            let d = j as i64 - i as i64 + 1;
            if d > best_d {
                best_i = i;
                best_d = d;
                best_c = c;
            }
            j += 1;
        } else {
            i += 1;
        }
    }
    (best_c >= 0.0).then(|| (best_i - 1, best_d as usize))
}

/// The multi-parent chimera detector, vsearch's `--chimeras_denovo`.
pub struct Mchime {
    options: LongOptions,
    finder: Finder,
}

impl Mchime {
    /// a detector searching `db`, which grows with the non-chimeras.
    pub fn new(db: &Database, options: LongOptions) -> Mchime {
        Mchime {
            options,
            finder: Finder::new(db, Some(options.abskew)),
        }
    }

    /// judges `query` against the parents of `db`: a chimera when the
    /// segments of two or more of them cover it.
    pub fn evaluate(&mut self, db: &Database, label: &str, size: i64, query: &[u8]) -> Verdict {
        let o = self.options;
        let parts = o.parts.unwrap_or(query.len().div_ceil(100)).clamp(2, 100);
        let candidates = self.finder.candidates(db, label, size, query, parts);
        let ops = self.finder.align(db, query, &candidates);
        match self.best_parents(db, query, &candidates, &ops) {
            Some(segments) => self.eval_parents(db, label, query, &candidates, &ops, &segments),
            None => Verdict {
                status: Status::NoParents,
                table: String::new(),
                alignment: None,
            },
        }
    }

    /// vsearch's `find_best_parents_long`. It takes in turn the longest
    /// segment of a candidate over positions not yet covered, without an
    /// insertion of the candidate inside it. It returns the parents, in
    /// query order, if they cover the whole query.
    fn best_parents(
        &self,
        db: &Database,
        query: &[u8],
        candidates: &[usize],
        ops: &[Vec<(u8, usize)>],
    ) -> Option<Vec<Segment>> {
        let o = &self.options;
        let qlen = query.len();
        let n = candidates.len();
        let mut matches = vec![0i32; n * qlen];
        let mut insert = vec![0usize; n * qlen];
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
                    b'I' => {
                        // vsearch's insertion past the end lands on the next
                        // candidate's first position, which no segment reads
                        if qpos < qlen {
                            insert[i * qlen + qpos] = run;
                        }
                        tpos += run;
                    }
                    _ => qpos += run,
                }
            }
        }
        let mut used = vec![false; qlen];
        let mut remaining = qlen;
        let mut segments: Vec<Segment> = Vec::new();
        for _ in 0..o.max_parents.min(MAX_PARENTS) {
            let mut best = Segment {
                cand: 0,
                start: 0,
                len: 0,
            };
            for i in 0..n {
                let mut j = 0;
                while j < qlen {
                    let start = j;
                    let mut len = 0;
                    while j < qlen && !used[j] && (len == 0 || insert[i * qlen + j] == 0) {
                        len += 1;
                        j += 1;
                    }
                    if len > best.len
                        && let Some((s, l)) = scan_matches(
                            &matches[i * qlen + start..i * qlen + start + len],
                            o.diff_pct,
                        )
                        && l > best.len
                    {
                        best = Segment {
                            cand: i,
                            start: start + s,
                            len: l,
                        };
                    }
                    j += 1;
                }
            }
            if best.len >= o.min_segment {
                segments.push(best);
                used[best.start..best.start + best.len].fill(true);
                remaining -= best.len;
            } else {
                break;
            }
        }
        segments.sort_by_key(|s| s.start);
        (segments.len() > 1 && remaining == 0).then_some(segments)
    }

    /// vsearch's `eval_parents_long`: the alignment of the query with its
    /// parents and their identities.
    fn eval_parents(
        &self,
        db: &Database,
        label: &str,
        query: &[u8],
        candidates: &[usize],
        ops: &[Vec<(u8, usize)>],
        segments: &[Segment],
    ) -> Verdict {
        let qlen = query.len();
        let parents: Vec<usize> = segments.iter().map(|s| s.cand).collect();
        let maxi = insertions(qlen, ops, &parents);
        let alnlen = qlen + maxi.iter().sum::<usize>();
        let mut paln: Vec<Vec<u8>> = parents
            .iter()
            .map(|&p| parent_row(&db.seqs[candidates[p]], &ops[p], &maxi))
            .collect();
        let qaln = query_row(query, &maxi);
        // the parent of each column: the segment covering its query position
        let mut model = Vec::with_capacity(alnlen);
        let mut nth = 0;
        for (qpos, &gaps) in maxi[..qlen].iter().enumerate() {
            if nth + 1 < segments.len() && qpos >= segments[nth].start + segments[nth].len {
                nth += 1;
            }
            let letter = b'A' + nth as u8;
            model.extend(std::iter::repeat_n(letter, gaps + 1));
        }
        model.extend(std::iter::repeat_n(b'A' + nth as u8, maxi[qlen]));

        let mut diffs = vec![b' '; alnlen];
        for i in 0..alnlen {
            let qsym = code(qaln[i]);
            let psym: Vec<u8> = paln.iter().map(|row| code(row[i])).collect();
            for (row, &p) in paln.iter_mut().zip(&psym) {
                if p != 0 && p != qsym {
                    row[i] = row[i].to_ascii_lowercase();
                }
            }
            if qsym != 0 && psym.iter().all(|&p| p != 0) {
                let mut found = 0;
                for (f, &p) in psym.iter().enumerate() {
                    if p == qsym {
                        diffs[i] = b'A' + f as u8;
                        found += 1;
                    }
                }
                if found > 1 {
                    diffs[i] = b' ';
                }
            }
        }
        let identities: Vec<f64> = paln
            .iter()
            .map(|row| {
                let same = (0..alnlen)
                    .filter(|&i| code(qaln[i]) == code(row[i]))
                    .count();
                100.0 * same as f64 / alnlen as f64
            })
            .collect();
        let qt = identities.iter().copied().fold(0.0, f64::max);
        let (qa, qb) = (identities[0], identities[1]);
        let qc = identities.get(2).copied().unwrap_or(0.0);
        let qm = 100.0;
        let divfrac = 100.0 * (qm - qt) / qt;
        let labels: Vec<&String> = parents.iter().map(|&p| &db.labels[candidates[p]]).collect();

        let mut text = String::new();
        text.push('\n');
        text.push_str(&"-".repeat(72));
        let _ = write!(text, "\nQuery   ({qlen:5} nt) {label}");
        for (f, &p) in parents.iter().enumerate() {
            let _ = write!(
                text,
                "\nParent{} ({:5} nt) {}",
                (b'A' + f as u8) as char,
                db.seqs[candidates[p]].len(),
                labels[f]
            );
        }
        text.push_str("\n\n");
        let width = if self.options.width > 0 {
            self.options.width
        } else {
            alnlen
        };
        let mut qpos = 0;
        let mut ppos = vec![0usize; parents.len()];
        let mut i = 0;
        while i < alnlen {
            let w = width.min(alnlen - i);
            let residues = |row: &[u8]| row[i..i + w].iter().filter(|&&c| c != b'-').count();
            let row = |text: &mut String, name: char, start: usize, row: &[u8], n: usize| {
                let _ = writeln!(
                    text,
                    "{name} {:5} {} {}",
                    start + 1,
                    String::from_utf8_lossy(&row[i..i + w]),
                    start + n
                );
            };
            let qnt = residues(&qaln);
            row(&mut text, 'Q', qpos, &qaln, qnt);
            for (f, prow) in paln.iter().enumerate() {
                let pnt = residues(prow);
                row(&mut text, (b'A' + f as u8) as char, ppos[f], prow, pnt);
                ppos[f] += pnt;
            }
            for (name, annotation) in [("Diffs   ", &diffs), ("Model   ", &model)] {
                text.push_str(name);
                text.push_str(&String::from_utf8_lossy(&annotation[i..i + w]));
                text.push('\n');
            }
            text.push('\n');
            qpos += qnt;
            i += width;
        }
        let _ = writeln!(
            text,
            "Ids.  QA {qa:.2}%, QB {qb:.2}%, QC {qc:.2}%, QT {qt:.2}%, QModel {qm:.2}%, Div. \
             {divfrac:+.2}%"
        );
        let third = labels.get(2).map_or("*", |l| l.as_str());
        let table = format!(
            "99.9999\t{label}\t{}\t{}\t{third}\t{qm:.2}\t{qa:.2}\t{qb:.2}\t{qc:.2}\t{qt:.2}\t0\t0\t0\t\
             0\t0\t0\t0.00\tY\n",
            labels[0], labels[1]
        );
        Verdict {
            status: Status::Chimeric,
            table,
            alignment: Some(text),
        }
    }
}
