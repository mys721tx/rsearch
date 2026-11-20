//! vsearch's search core (`searchcore.cpp`), ported step for step so that
//! the hits, their order and their ties are vsearch's. The database targets
//! sharing the most words with a query are aligned, best first and in
//! batches of eight, until enough are accepted or rejected.

use std::cmp::Reverse;
use std::collections::BinaryHeap;

use crate::accept::{Criteria, Verdict};
use crate::align16::Search16;
use crate::alignment::Stats;
use crate::linmem::LinearMemoryAligner;
use crate::scoring::{Scoring, Substitutions};
use crate::words::{Index, Uniquer};

/// The targets vsearch aligns at once.
const MAX_DELAYED: usize = 8;

/// The highest word count of a target, as vsearch's 16-bit counters.
const MAX_COUNT: u16 = i16::MAX as u16;

/// vsearch's default words a query and a target share to be a candidate,
/// by word length.
pub const MIN_WORD_MATCHES: [u32; 16] = [0, 0, 0, 18, 17, 16, 15, 14, 12, 11, 10, 9, 8, 7, 5, 3];

/// A database to search: its records, their abundances, and the index of
/// their words.
pub struct Database {
    pub labels: Vec<String>,
    pub seqs: Vec<Vec<u8>>,
    pub sizes: Vec<i64>,
    index: Index,
    /// the distinct words of each sequence.
    words: Vec<u32>,
}

impl Database {
    pub fn new(labels: Vec<String>, seqs: Vec<Vec<u8>>, sizes: Vec<i64>, index: Index) -> Database {
        let mut words = vec![0u32; seqs.len()];
        for &s in index.all_sequences() {
            words[s as usize] += 1;
        }
        Database {
            labels,
            seqs,
            sizes,
            index,
            words,
        }
    }

    pub fn len(&self) -> usize {
        self.seqs.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seqs.is_empty()
    }

    pub fn word_len(&self) -> u32 {
        self.index.word_len()
    }
}

/// The limits of a search: at most this many hits accepted or rejected (0
/// for every target), and the words a candidate shares with the query.
#[derive(Clone, Copy)]
pub struct Limits {
    pub max_accepts: usize,
    pub max_rejects: usize,
    pub min_word_matches: u32,
}

/// A hit kept: accepted, or weak.
pub struct Hit {
    pub target: usize,
    pub minus: bool,
    pub weak: bool,
    pub stats: Stats,
    pub cigar: String,
}

/// A candidate taken from the heap, rejected before alignment or not.
struct Candidate {
    target: usize,
    rejected: bool,
}

/// Searches queries against a database, one strand at a time.
pub struct Searcher<'a> {
    db: &'a Database,
    criteria: Criteria,
    scoring: Scoring,
    table: Substitutions,
    max_accepts: usize,
    max_rejects: usize,
    min_word_matches: u32,
    top_hits: usize,
    uniquer: Uniquer,
    counts: Vec<u16>,
    simd: Search16,
    fallback: LinearMemoryAligner,
}

impl<'a> Searcher<'a> {
    pub fn new(
        db: &'a Database,
        criteria: Criteria,
        scoring: Scoring,
        limits: Limits,
    ) -> Searcher<'a> {
        let n = db.len();
        let limit = |x: usize| if x == 0 || x > n { n } else { x };
        let (max_accepts, max_rejects) = (limit(limits.max_accepts), limit(limits.max_rejects));
        Searcher {
            db,
            table: scoring.substitutions(),
            simd: Search16::new(&scoring),
            fallback: LinearMemoryAligner::new(&scoring),
            criteria,
            scoring,
            max_accepts,
            max_rejects,
            min_word_matches: limits.min_word_matches,
            top_hits: (max_accepts + max_rejects + MAX_DELAYED).min(n),
            uniquer: Uniquer::new(db.word_len()),
            counts: vec![0; n],
        }
    }

    /// the candidates for `query`, best first: the targets sharing the
    /// most words, then the shorter, then the first; at most `top_hits`.
    fn candidates(&mut self, query: &[u8]) -> Vec<usize> {
        let db = self.db;
        let words = self.uniquer.words(query);
        let min_matches = self.min_word_matches.min(words.len() as u32);
        self.counts.fill(0);
        for &word in words {
            for &s in db.index.sequences_of(word) {
                let c = &mut self.counts[s as usize];
                if *c < MAX_COUNT {
                    *c += 1;
                }
            }
        }
        // the worst kept on top
        let mut heap: BinaryHeap<(Reverse<u16>, usize, usize)> = BinaryHeap::new();
        if self.top_hits == 0 {
            return Vec::new();
        }
        for (s, &count) in self.counts.iter().enumerate() {
            let c = u32::from(count);
            let kmers = db.words[s];
            let eligible = c >= min_matches
                || (kmers != 0 && kmers < self.min_word_matches && c >= min_matches.min(kmers));
            if !eligible {
                continue;
            }
            let key = (Reverse(count), db.seqs[s].len(), s);
            if heap.len() < self.top_hits {
                heap.push(key);
            } else if key < *heap.peek().unwrap() {
                heap.pop();
                heap.push(key);
            }
        }
        heap.into_sorted_vec()
            .into_iter()
            .map(|(_, _, s)| s)
            .collect()
    }

    /// searches one strand of a query (reverse-complemented for the minus
    /// strand), and returns the hits accepted and the weak ones, in the
    /// order found.
    pub fn search(&mut self, label: &str, size: i64, query: &[u8], minus: bool) -> Vec<Hit> {
        let db = self.db;
        let candidates = self.candidates(query);
        self.simd.set_query(query);
        let (max_accepts, max_rejects) = (self.max_accepts, self.max_rejects);
        let mut pending: Vec<Candidate> = Vec::new();
        let mut hits = Vec::new();
        let (mut finalized, mut delayed) = (0usize, 0usize);
        let (mut accepts, mut rejects) = (0usize, 0usize);
        let mut next = candidates.into_iter();
        while finalized + delayed < (max_accepts + max_rejects).saturating_sub(1)
            && rejects < max_rejects
            && accepts < max_accepts
        {
            let Some(t) = next.next() else { break };
            let ok = self.criteria.unaligned(
                label,
                query,
                size,
                &db.labels[t],
                &db.seqs[t],
                db.sizes[t],
            );
            if ok {
                delayed += 1;
            }
            pending.push(Candidate {
                target: t,
                rejected: !ok,
            });
            if delayed == MAX_DELAYED {
                self.align_delayed(
                    query,
                    minus,
                    &pending[finalized..],
                    &mut accepts,
                    &mut rejects,
                    &mut hits,
                );
                finalized = pending.len();
                delayed = 0;
            }
        }
        if delayed > 0 {
            self.align_delayed(
                query,
                minus,
                &pending[finalized..],
                &mut accepts,
                &mut rejects,
                &mut hits,
            );
        }
        hits
    }

    /// vsearch's `align_delayed`: aligns the candidates not yet rejected,
    /// then counts the batch's accepted and rejected ones up to the limits.
    fn align_delayed(
        &mut self,
        query: &[u8],
        minus: bool,
        batch: &[Candidate],
        accepts: &mut usize,
        rejects: &mut usize,
        hits: &mut Vec<Hit>,
    ) {
        let db = self.db;
        let targets: Vec<&[u8]> = batch
            .iter()
            .filter(|c| !c.rejected)
            .map(|c| db.seqs[c.target].as_slice())
            .collect();
        let mut aligned = if targets.is_empty() {
            Vec::new()
        } else {
            self.simd.align(&targets)
        }
        .into_iter();
        for c in batch {
            if *rejects >= self.max_rejects || *accepts >= self.max_accepts {
                break;
            }
            if c.rejected {
                *rejects += 1;
                continue;
            }
            let a = aligned.next().expect("an alignment for each candidate");
            let target = &db.seqs[c.target];
            // the 16-bit aligner's counts, or the fallback's statistics
            let (cigar, stats) = if a.score == i16::MAX {
                let cigar = self.fallback.align(query, target);
                let stats = Stats::new(&cigar, query, target, &self.scoring, &self.table);
                (cigar, stats)
            } else {
                let stats = Stats::from_counts(
                    &a.cigar,
                    query.len(),
                    target.len(),
                    i64::from(a.length),
                    i64::from(a.matches),
                    i64::from(a.mismatches),
                    i64::from(a.gaps),
                    i64::from(a.score),
                );
                (a.cigar, stats)
            };
            let stats = stats.expect("the aligners' CIGARs are valid");
            let verdict =
                self.criteria
                    .aligned(&stats, &cigar, query.len(), target.len(), &self.scoring);
            if verdict == Verdict::Accepted {
                *accepts += 1;
            } else {
                *rejects += 1;
            }
            if verdict != Verdict::Rejected {
                hits.push(Hit {
                    target: c.target,
                    minus,
                    weak: verdict == Verdict::Weak,
                    stats,
                    cigar,
                });
            }
        }
    }
}

/// orders the hits of both strands as vsearch reports them: the accepted
/// before the weak, by decreasing identity, then target.
pub fn order(hits: &mut [Hit], iddef: usize) {
    hits.sort_by(|a, b| {
        a.weak
            .cmp(&b.weak)
            .then(b.stats.id[iddef].total_cmp(&a.stats.id[iddef]))
            .then(a.target.cmp(&b.target))
    });
}
