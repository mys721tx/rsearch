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
use crate::words::{Index, MIN_WORD_LEN, Uniquer};

/// The targets vsearch aligns at once.
const MAX_DELAYED: usize = 8;

/// The highest word count of a target, as vsearch's 16-bit counters.
const MAX_COUNT: u16 = i16::MAX as u16;

/// vsearch's default words a query and a target share to be a candidate,
/// by word length.
pub const MIN_WORD_MATCHES: [u32; 16] = [0, 0, 0, 18, 17, 16, 15, 14, 12, 11, 10, 9, 8, 7, 5, 3];

/// The index of a database's words: for each word, the targets holding
/// it, by their place in the index; and for each target, its sequence and
/// the number of its distinct words. The targets are the database's
/// sequences, or those added as the index grows, as clustering adds its
/// centroids.
struct Postings {
    /// for each word, where its targets start in `entries`; one more ends
    /// the last.
    starts: Vec<u64>,
    /// for each word, its targets so far; empty when every sequence is a
    /// target.
    filled: Vec<u32>,
    entries: Vec<u32>,
    /// the sequence of each target.
    seqnos: Vec<u32>,
    /// the distinct words of each target.
    words: Vec<u32>,
}

impl Postings {
    /// room for the words of `counts`, the number of sequences holding each.
    fn with_room(counts: &[u32]) -> Postings {
        let mut starts = Vec::with_capacity(counts.len() + 1);
        let mut sum = 0u64;
        starts.push(0);
        for &c in counts {
            sum += u64::from(c);
            starts.push(sum);
        }
        Postings {
            starts,
            filled: vec![0; counts.len()],
            entries: vec![0; sum as usize],
            seqnos: Vec::new(),
            words: Vec::new(),
        }
    }

    fn sequences_of(&self, word: u32) -> &[u32] {
        let w = word as usize;
        let start = self.starts[w] as usize;
        let end = if self.filled.is_empty() {
            self.starts[w + 1] as usize
        } else {
            start + self.filled[w] as usize
        };
        &self.entries[start..end]
    }

    /// adds sequence `seqno` of `words` as the next target.
    fn add(&mut self, seqno: usize, words: &[u32]) {
        let target = self.seqnos.len() as u32;
        for &w in words {
            let w = w as usize;
            self.entries[(self.starts[w] + u64::from(self.filled[w])) as usize] = target;
            self.filled[w] += 1;
        }
        self.seqnos.push(seqno as u32);
        self.words.push(words.len() as u32);
    }
}

/// A database to search: its records, their abundances, and the index of
/// their words.
pub struct Database {
    pub labels: Vec<String>,
    pub seqs: Vec<Vec<u8>>,
    pub sizes: Vec<i64>,
    word_len: u32,
    postings: Postings,
    uniquer: Uniquer,
}

impl Database {
    /// a database of every sequence of `index`.
    pub fn new(labels: Vec<String>, seqs: Vec<Vec<u8>>, sizes: Vec<i64>, index: Index) -> Database {
        let word_len = index.word_len();
        let (starts, entries) = index.into_parts();
        let mut words = vec![0; seqs.len()];
        for &s in &entries {
            words[s as usize] += 1;
        }
        let postings = Postings {
            starts,
            filled: Vec::new(),
            entries,
            seqnos: (0..seqs.len() as u32).collect(),
            words,
        };
        Database {
            labels,
            seqs,
            sizes,
            word_len,
            postings,
            // only a growing database needs one
            uniquer: Uniquer::new(MIN_WORD_LEN),
        }
    }

    /// a database of `seqs` without targets yet, words of `word_len`, to
    /// which `add` adds them.
    pub fn growing(
        labels: Vec<String>,
        seqs: Vec<Vec<u8>>,
        sizes: Vec<i64>,
        word_len: u32,
    ) -> Database {
        let (counts, _) = crate::words::word_counts(word_len, seqs.iter().map(Vec::as_slice));
        Database {
            postings: Postings::with_room(&counts),
            labels,
            seqs,
            sizes,
            word_len,
            uniquer: Uniquer::new(word_len),
        }
    }

    /// makes sequence `seqno` a target.
    pub fn add(&mut self, seqno: usize) {
        let words = self.uniquer.words(&self.seqs[seqno]);
        self.postings.add(seqno, words);
    }

    pub fn len(&self) -> usize {
        self.seqs.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seqs.is_empty()
    }

    pub fn word_len(&self) -> u32 {
        self.word_len
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

/// The hits of a strand so far, and its accepted and rejected targets.
#[derive(Default)]
struct Tally {
    accepts: usize,
    rejects: usize,
    hits: Vec<Hit>,
}

/// Searches queries against a database, one strand at a time.
pub struct Searcher {
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

impl Searcher {
    /// a searcher of `db`, which may grow but keeps its sequences.
    pub fn new(db: &Database, criteria: Criteria, scoring: Scoring, limits: Limits) -> Searcher {
        let n = db.len();
        let limit = |x: usize| if x == 0 || x > n { n } else { x };
        let (max_accepts, max_rejects) = (limit(limits.max_accepts), limit(limits.max_rejects));
        Searcher {
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
            counts: Vec::new(),
        }
    }

    /// the candidates for `query`, best first: the targets sharing the
    /// most words, then the shorter, then the first; at most `top_hits`.
    fn candidates(&mut self, db: &Database, query: &[u8]) -> Vec<usize> {
        let postings = &db.postings;
        let words = self.uniquer.words(query);
        let min_matches = self.min_word_matches.min(words.len() as u32);
        self.counts.clear();
        self.counts.resize(postings.seqnos.len(), 0);
        for &word in words {
            for &s in postings.sequences_of(word) {
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
        for (t, &count) in self.counts.iter().enumerate() {
            let c = u32::from(count);
            let kmers = postings.words[t];
            let eligible = c >= min_matches
                || (kmers != 0 && kmers < self.min_word_matches && c >= min_matches.min(kmers));
            if !eligible {
                continue;
            }
            let s = postings.seqnos[t] as usize;
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
    pub fn search(
        &mut self,
        db: &Database,
        label: &str,
        size: i64,
        query: &[u8],
        minus: bool,
    ) -> Vec<Hit> {
        let candidates = self.candidates(db, query);
        self.simd.set_query(query);
        let (max_accepts, max_rejects) = (self.max_accepts, self.max_rejects);
        let mut pending: Vec<Candidate> = Vec::new();
        let mut tally = Tally::default();
        let (mut finalized, mut delayed) = (0usize, 0usize);
        let mut next = candidates.into_iter();
        while finalized + delayed < (max_accepts + max_rejects).saturating_sub(1)
            && tally.rejects < max_rejects
            && tally.accepts < max_accepts
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
                self.align_delayed(db, query, minus, &pending[finalized..], &mut tally);
                finalized = pending.len();
                delayed = 0;
            }
        }
        if delayed > 0 {
            self.align_delayed(db, query, minus, &pending[finalized..], &mut tally);
        }
        tally.hits
    }

    /// vsearch's `align_delayed`: aligns the candidates not yet rejected,
    /// then counts the batch's accepted and rejected ones up to the limits.
    fn align_delayed(
        &mut self,
        db: &Database,
        query: &[u8],
        minus: bool,
        batch: &[Candidate],
        tally: &mut Tally,
    ) {
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
            if tally.rejects >= self.max_rejects || tally.accepts >= self.max_accepts {
                break;
            }
            if c.rejected {
                tally.rejects += 1;
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
                tally.accepts += 1;
            } else {
                tally.rejects += 1;
            }
            if verdict != Verdict::Rejected {
                tally.hits.push(Hit {
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
