//! vsearch's search core (`searchcore.cpp`), ported step for step so that
//! the hits, their order and their ties are vsearch's. The database targets
//! sharing the most words with a query are aligned, best first and in
//! batches of eight, until enough are accepted or rejected.

use std::cmp::Reverse;
use std::collections::BinaryHeap;

use crate::accept::{Criteria, Verdict};
use crate::align16::{Aligned, Search16};
use crate::alignment::Stats;
use crate::linmem::LinearMemoryAligner;
use crate::scoring::{FOUR_BIT, Scoring, Substitutions};
use std::fs::File;
use std::io::{BufRead, BufReader, Read};

use crate::fastx::{Reader, Record};
use crate::label::checked_abundance;
use crate::udb::{self, Header, Udb};
use crate::words::{Index, MIN_WORD_LEN, Uniquer};

/// The targets vsearch aligns at once.
const MAX_DELAYED: usize = 8;

/// The highest word count of a target, as vsearch's 16-bit counters.
const MAX_COUNT: u16 = i16::MAX as u16;

/// vsearch's default words a query and a target share to be a candidate,
/// by word length.
pub const MIN_WORD_MATCHES: [u32; 16] = [0, 0, 0, 18, 17, 16, 15, 14, 12, 11, 10, 9, 8, 7, 5, 3];

/// A word held by at least one in this many of a database's sequences is
/// kept as a bitmap of the targets too, as vsearch's database index does.
const BITMAP_SHARE: usize = 8;

/// The targets of a bitmap word, and of a block of the scan for
/// candidates.
const BLOCK: usize = 64;

/// For each four bits of a bitmap, the four 16-bit counters they raise by one,
/// as one word to add. Bit `i` becomes counter `i`.
const EXPANSION: [u64; 16] = {
    let mut table = [0u64; 16];
    let mut bits = 0;
    while bits < 16 {
        let mut i = 0;
        while i < 4 {
            if bits & (1 << i) != 0 {
                table[bits] |= 1 << (16 * i);
            }
            i += 1;
        }
        bits += 1;
    }
    table
};

/// The high bit of each 16-bit counter of a word.
const HIGH_BITS: u64 = 0x8000_8000_8000_8000;
/// One, in each 16-bit counter of a word.
const ONES: u64 = 0x0001_0001_0001_0001;

/// The index of a database's words: for each word, the targets holding
/// it, by their place in the index; and for each target, its sequence and
/// the number of its distinct words. The targets are the database's
/// sequences, or those added as the index grows, as clustering adds its
/// centroids. The words held by many sequences have a bitmap of their
/// targets as well.
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
    /// the words with a bitmap, in increasing order, and their bitmaps, a
    /// bit per target, the first in the lowest bit; `bitmap_len` words
    /// each.
    common: Vec<u32>,
    bitmaps: Vec<u64>,
    bitmap_len: usize,
    /// for each block of targets, the fewest distinct words of a target
    /// holding some (`u32::MAX` for none).
    fewest: Vec<u32>,
}

impl Postings {
    /// room for the words of `counts`, the number of sequences holding
    /// each, among `n` sequences.
    fn with_room(counts: &[u32], n: usize) -> Postings {
        let mut starts = Vec::with_capacity(counts.len() + 1);
        let mut sum = 0u64;
        starts.push(0);
        for &c in counts {
            sum += u64::from(c);
            starts.push(sum);
        }
        let mut postings = Postings {
            starts,
            filled: vec![0; counts.len()],
            entries: vec![0; sum as usize],
            seqnos: Vec::new(),
            words: Vec::new(),
            common: Vec::new(),
            bitmaps: Vec::new(),
            bitmap_len: n.div_ceil(BLOCK),
            fewest: Vec::new(),
        };
        postings.make_bitmaps(Some(counts), n);
        postings
    }

    /// the bitmaps of the words held by enough of `n` sequences. A database
    /// that grows finds them by `counts`, the sequences holding each word, and
    /// fills them as targets are added. Any other finds them by its lists of
    /// targets, which fill them.
    fn make_bitmaps(&mut self, counts: Option<&[u32]>, n: usize) {
        let threshold = (n / BITMAP_SHARE).max(1);
        self.bitmap_len = n.div_ceil(BLOCK);
        for word in 0..self.starts.len() - 1 {
            let count = match counts {
                Some(counts) => u64::from(counts[word]),
                None => self.starts[word + 1] - self.starts[word],
            };
            if count >= threshold as u64 {
                let mut bitmap = vec![0u64; self.bitmap_len];
                if counts.is_none() {
                    for &t in self.sequences_of(word as u32) {
                        bitmap[t as usize / BLOCK] |= 1 << (t as usize % BLOCK);
                    }
                }
                self.bitmaps.extend_from_slice(&bitmap);
                self.common.push(word as u32);
            }
        }
    }

    /// the bitmap of `word`, if it has one.
    fn bitmap(&self, word: u32) -> Option<&[u64]> {
        let slot = self.common.binary_search(&word).ok()?;
        let start = slot * self.bitmap_len;
        Some(&self.bitmaps[start..start + self.bitmap_len])
    }

    /// notes the distinct words of the next target.
    fn note_words(&mut self, words: u32) {
        let block = self.words.len() / BLOCK;
        if block == self.fewest.len() {
            self.fewest.push(u32::MAX);
        }
        if words != 0 {
            self.fewest[block] = self.fewest[block].min(words);
        }
        self.words.push(words);
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
            let at = w as usize;
            self.entries[(self.starts[at] + u64::from(self.filled[at])) as usize] = target;
            self.filled[at] += 1;
            if let Ok(slot) = self.common.binary_search(&w) {
                let t = target as usize;
                self.bitmaps[slot * self.bitmap_len + t / BLOCK] |= 1 << (t % BLOCK);
            }
        }
        self.seqnos.push(seqno as u32);
        self.note_words(words.len() as u32);
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
        let mut postings = Postings {
            starts,
            filled: Vec::new(),
            entries,
            seqnos: (0..seqs.len() as u32).collect(),
            words: Vec::new(),
            common: Vec::new(),
            bitmaps: Vec::new(),
            bitmap_len: 0,
            fewest: Vec::new(),
        };
        for w in words {
            postings.note_words(w);
        }
        postings.make_bitmaps(None, seqs.len());
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
            postings: Postings::with_room(&counts, seqs.len()),
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

/// reads a database to search: a UDB file, whose word length wins over
/// `word_len`, or a FASTA or FASTQ file indexed at `word_len`; `warn` gets
/// the warnings.
pub fn read_database(
    path: &str,
    word_len: u32,
    warn: &mut dyn FnMut(&[u8]),
) -> Result<Database, String> {
    let file = File::open(path).map_err(|e| format!("cannot open {path}: {e}"))?;
    let file_size = file
        .metadata()
        .ok()
        .filter(|m| m.is_file())
        .map(|m| m.len());
    let mut fin = BufReader::new(file);
    let start = fin
        .fill_buf()
        .map_err(|e| format!("cannot read {path}: {e}"))?;
    if start.starts_with(udb::MAGIC) {
        let mut bytes = Vec::new();
        let file_size = match file_size {
            Some(s) => s,
            None => {
                fin.read_to_end(&mut bytes)
                    .map_err(|e| format!("cannot read {path}: {e}"))?;
                bytes.len() as u64
            }
        };
        // vsearch warns as it reads the header
        let first = if bytes.is_empty() {
            fin.fill_buf().unwrap_or_default()
        } else {
            &bytes
        };
        if let Ok(header) = Header::parse(first, file_size)
            && header.word_len() != word_len
        {
            warn(
                format!(
                    "the word length of the UDB file, {}, is used",
                    header.word_len()
                )
                .as_bytes(),
            );
        }
        let mut list = Vec::new();
        let keep = |_: &[u32], _: u64, block: &[u32]| list.extend_from_slice(block);
        let read = if bytes.is_empty() {
            Udb::read(&mut fin, file_size, keep)
        } else {
            Udb::read(&mut &bytes[..], file_size, keep)
        };
        let udb = read.map_err(|e| e.to_string())?;
        let labels: Vec<String> = udb
            .header_starts
            .windows(2)
            .map(|pair| {
                String::from_utf8_lossy(&udb.headers[pair[0] as usize..pair[1] as usize - 1])
                    .into_owned()
            })
            .collect();
        let sizes = labels
            .iter()
            .map(|l| {
                checked_abundance(l)
                    .map(|s| s as i64)
                    .map_err(|e| format!("{e} in the database"))
            })
            .collect::<Result<_, _>>()?;
        let mut seqs = Vec::with_capacity(labels.len());
        let mut at = 0;
        for &length in &udb.lengths {
            seqs.push(udb.sequences[at..at + length as usize].to_vec());
            at += length as usize;
        }
        let index = Index::from_parts(udb.header.word_len(), &udb.counts, list);
        return Ok(Database::new(labels, seqs, sizes, index));
    }

    let mut reader = Reader::new(fin);
    let (mut labels, mut seqs, mut sizes) = (Vec::new(), Vec::new(), Vec::new());
    let mut record = Record::new();
    loop {
        reader.read(&mut record).map_err(|e| e.to_string())?;
        if record.is_empty() {
            break;
        }
        let header = record.header();
        let abundance = checked_abundance(header)
            .map_err(|e| format!("{e} in read {header} of the database"))?;
        labels.push(header.to_string());
        seqs.push(record.seq().to_vec());
        sizes.push(abundance as i64);
    }
    if let Some(warning) = reader.stripped_warning() {
        warn(&warning);
    }
    let index = Index::new(word_len, seqs.iter().map(Vec::as_slice));
    if index.without_words() > 0 {
        warn(
            format!(
                "{} of {} sequences give no word to index: they are shorter than the word \
                 length ({word_len}), or each of their words holds a lower-case residue or a \
                 symbol other than A, C, G, T or U",
                index.without_words(),
                seqs.len()
            )
            .as_bytes(),
        );
    }
    Ok(Database::new(labels, seqs, sizes, index))
}

/// The limits of a search: at most this many hits accepted or rejected (0
/// for every target), and the words a candidate shares with the query.
#[derive(Clone, Copy)]
pub struct Limits {
    pub max_accepts: usize,
    pub max_rejects: usize,
    pub min_word_matches: u32,
    /// the candidates kept, by default the limits and a batch of eight
    pub top_hits: Option<usize>,
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

/// A batch of candidates taken by [`Searcher::pull`]: where it ends in the
/// candidates taken, how many are to be aligned and whether that is
/// `MAX_DELAYED`, the last taken, and whether the candidates ran out.
#[derive(Default)]
struct Batch {
    end: usize,
    delayed: usize,
    full: bool,
    last: Option<Key>,
    exhausted: bool,
}

/// The hits of a strand so far, and its accepted and rejected targets.
#[derive(Default)]
struct Tally {
    accepts: usize,
    rejects: usize,
    hits: Vec<Hit>,
}

/// A candidate's place in the order targets are aligned in: the most
/// words shared first, then the shorter, then the first.
type Key = (Reverse<u16>, usize, usize);

/// How far a search went down its candidates. A target added to the
/// database after the search can change it only if it would have come
/// before the last candidate the search took, or if the search ran out of
/// candidates before it ran out of room for them ([`Searcher::would_change`]).
#[derive(Clone, Debug, Default)]
pub struct Horizon {
    /// the distinct words of the query, sorted.
    words: Vec<u32>,
    /// the words a target had to share with the query.
    min_matches: u32,
    /// the key of the last candidate taken.
    last: Option<Key>,
    /// whether the search took every candidate, and there were fewer than
    /// it keeps.
    open: bool,
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
    /// a 16-bit counter per target, four to a word, target `t` in counter
    /// `t % 4` (from the lowest) of word `t / 4`.
    counts: Vec<u64>,
    simd: Search16,
    fallback: LinearMemoryAligner,
    /// UNOISE's alpha, when its abundance skew rule decides acceptance.
    unoise_alpha: Option<f64>,
    /// the abundance of the query searched.
    query_size: i64,
    /// the whole query, when a part of it is searched with `--self-id`.
    whole: Option<Vec<u8>>,
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
            top_hits: limits
                .top_hits
                .unwrap_or(max_accepts + max_rejects + MAX_DELAYED)
                .min(n),
            uniquer: Uniquer::new(db.word_len()),
            counts: Vec::new(),
            unoise_alpha: None,
            query_size: 1,
            whole: None,
        }
    }

    /// searches `part` of the query `whole`, as vsearch searches the parts
    /// of a query for chimeras. `--self-id` compares the whole query, not
    /// the part, with the targets.
    pub fn search_part(
        &mut self,
        db: &Database,
        label: &str,
        size: i64,
        part: &[u8],
        whole: &[u8],
    ) -> Vec<Hit> {
        // the part itself is not compared
        let self_id = std::mem::replace(&mut self.criteria.self_id, false);
        if self_id {
            self.whole = Some(whole.to_vec());
        }
        let hits = self.search(db, label, size, part, false);
        self.whole = None;
        self.criteria.self_id = self_id;
        hits
    }

    /// accepts the hits that pass the criteria only by UNOISE's rule, as
    /// vsearch --cluster_unoise. A target must be at most one in 2^(alpha d
    /// + 1) as abundant as the query, for d mismatches, or have no
    /// mismatches. The others are weak.
    pub fn set_unoise(&mut self, alpha: f64) {
        self.unoise_alpha = Some(alpha);
    }

    /// the candidates for `query`, best first: the targets sharing the
    /// most words, then the shorter, then the first; at most `top_hits`.
    /// With `words`, the distinct words of the query, sorted, and the words
    /// a candidate has to share.
    fn candidates(
        &mut self,
        db: &Database,
        query: &[u8],
        words_out: Option<&mut (Vec<u32>, u32)>,
    ) -> Vec<Key> {
        let postings = &db.postings;
        let words = self.uniquer.words(query);
        let min_matches = self.min_word_matches.min(words.len() as u32);
        if let Some((sorted, min)) = words_out {
            sorted.clear();
            sorted.extend_from_slice(words);
            sorted.sort_unstable();
            *min = min_matches;
        }
        let targets = postings.seqnos.len();
        let blocks = targets.div_ceil(BLOCK);
        self.counts.clear();
        self.counts.resize(blocks * BLOCK / 4, 0);
        if words.len() <= usize::from(MAX_COUNT) {
            // no count can reach the counters' limit
            for &word in words {
                if let Some(bitmap) = postings.bitmap(word) {
                    // four counters at once, four bits of the bitmap each
                    for (bits, counts) in bitmap[..blocks]
                        .iter()
                        .zip(self.counts.as_chunks_mut::<{ BLOCK / 4 }>().0)
                    {
                        if *bits == 0 {
                            continue;
                        }
                        for (i, count) in counts.iter_mut().enumerate() {
                            *count += EXPANSION[((bits >> (4 * i)) & 0xf) as usize];
                        }
                    }
                } else {
                    for &t in postings.sequences_of(word) {
                        self.counts[t as usize / 4] += 1 << (16 * (t % 4));
                    }
                }
            }
        } else {
            for &word in words {
                for &t in postings.sequences_of(word) {
                    let (at, shift) = (t as usize / 4, 16 * (t % 4));
                    if (self.counts[at] >> shift) & 0xffff < u64::from(MAX_COUNT) {
                        self.counts[at] += 1 << shift;
                    }
                }
            }
        }
        // the worst kept on top
        let mut heap: BinaryHeap<Key> = BinaryHeap::new();
        if self.top_hits == 0 {
            return Vec::new();
        }
        for (block, counts) in self
            .counts
            .as_chunks::<{ BLOCK / 4 }>()
            .0
            .iter()
            .enumerate()
        {
            // a candidate shares at least this many words. That is the fewest
            // any target of the block may, and once the heap is full, as many
            // as its worst
            let mut least = min_matches.min(postings.fewest[block]);
            if heap.len() == self.top_hits {
                least = least.max(u32::from(heap.peek().unwrap().0.0));
            }
            let least = u64::from(least.min(u32::from(MAX_COUNT)));
            for (at, &four) in counts.iter().enumerate() {
                // the counters reaching `least`. Adding 2^15 less it sets a
                // counter's high bit exactly when it does, and carries into no
                // other counter
                let mut reaching = four.wrapping_add(ONES * (0x8000 - least)) & HIGH_BITS;
                while reaching != 0 {
                    let lane = reaching.trailing_zeros() as usize / 16;
                    reaching &= reaching - 1;
                    let t = block * BLOCK + at * 4 + lane;
                    if t >= targets {
                        break;
                    }
                    let count = ((four >> (16 * lane)) & 0xffff) as u16;
                    let kmers = postings.words[t];
                    if !self.eligible(u32::from(count), kmers, min_matches) {
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
            }
        }
        heap.into_sorted_vec()
    }

    /// whether a target of `kmers` distinct words, sharing `shared` with a
    /// query, is a candidate for it.
    fn eligible(&self, shared: u32, kmers: u32, min_matches: u32) -> bool {
        shared >= min_matches
            || (kmers != 0 && kmers < self.min_word_matches && shared >= min_matches.min(kmers))
    }

    /// the distinct words of `seq`, sorted, as [`Searcher::would_change`]
    /// takes them.
    pub fn distinct_words(&mut self, seq: &[u8]) -> Vec<u32> {
        let mut words = self.uniquer.words(seq).to_vec();
        words.sort_unstable();
        words
    }

    /// whether making `seqno` of `db` a target, with the sorted distinct
    /// `words`, could change the search that went as far as `horizon`.
    /// When it cannot, the search finds the same with the target as without it.
    /// Then the target is no candidate, or it comes after every candidate the
    /// search took and the search stopped before its last.
    pub fn would_change(
        &self,
        db: &Database,
        horizon: &Horizon,
        seqno: usize,
        words: &[u32],
    ) -> bool {
        if self.top_hits == 0 {
            return false;
        }
        // the words the query and the target share, counted as candidates
        // counts them
        let (mut i, mut j, mut shared) = (0, 0, 0u32);
        while i < horizon.words.len() && j < words.len() {
            match horizon.words[i].cmp(&words[j]) {
                std::cmp::Ordering::Less => i += 1,
                std::cmp::Ordering::Greater => j += 1,
                std::cmp::Ordering::Equal => {
                    shared += 1;
                    i += 1;
                    j += 1;
                }
            }
        }
        let count = shared.min(u32::from(MAX_COUNT)) as u16;
        if !self.eligible(u32::from(count), words.len() as u32, horizon.min_matches) {
            return false;
        }
        let key = (Reverse(count), db.seqs[seqno].len(), seqno);
        horizon.open || horizon.last.is_none_or(|last| key < last)
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
        self.search_inner(db, label, size, query, minus, None)
    }

    /// searches as [`Searcher::search`] does, and returns how far it went
    /// down its candidates too.
    pub fn search_seen(
        &mut self,
        db: &Database,
        label: &str,
        size: i64,
        query: &[u8],
        minus: bool,
    ) -> (Vec<Hit>, Horizon) {
        let mut horizon = Horizon::default();
        let hits = self.search_inner(db, label, size, query, minus, Some(&mut horizon));
        (hits, horizon)
    }

    fn search_inner(
        &mut self,
        db: &Database,
        label: &str,
        size: i64,
        query: &[u8],
        minus: bool,
        mut horizon: Option<&mut Horizon>,
    ) -> Vec<Hit> {
        self.query_size = size;
        let mut words = horizon.as_ref().map(|_| (Vec::new(), 0));
        let candidates = self.candidates(db, query, words.as_mut());
        let kept = candidates.len();
        if let (Some(h), Some((sorted, min))) = (horizon.as_deref_mut(), words) {
            h.words = sorted;
            h.min_matches = min;
        }
        self.simd.set_query(query);
        let limit = (self.max_accepts + self.max_rejects).saturating_sub(1);
        // with sixteen lanes a batch and the next are aligned in one call;
        // the next is counted only if vsearch would have taken it
        let pair = self.simd.lanes() >= 2 * MAX_DELAYED;
        let mut pending: Vec<Candidate> = Vec::new();
        let mut tally = Tally::default();
        let mut next = candidates.into_iter();
        // the batch to pull starts after the last full batch
        let mut base = 0;
        let reached = |tally: &Tally, this: &Searcher| {
            tally.rejects >= this.max_rejects || tally.accepts >= this.max_accepts
        };
        let note = |horizon: &mut Option<&mut Horizon>, batch: &Batch, top_hits: usize| {
            if let Some(h) = horizon.as_deref_mut() {
                if batch.last.is_some() {
                    h.last = batch.last;
                }
                if batch.exhausted {
                    h.open = kept < top_hits;
                }
            }
        };
        while !reached(&tally, self) {
            let first = self.pull(&mut next, &mut pending, base, limit, db, label, query, size);
            let second = (pair && first.full).then(|| {
                self.pull(
                    &mut next,
                    &mut pending,
                    first.end,
                    limit,
                    db,
                    label,
                    query,
                    size,
                )
            });
            let both = second.as_ref().map_or(first.end, |b| b.end);
            let mut aligned = self.align_batch(db, &pending[base..both]).into_iter();
            note(&mut horizon, &first, self.top_hits);
            if first.delayed > 0 {
                self.align_delayed(
                    db,
                    query,
                    minus,
                    &pending[base..first.end],
                    &mut aligned,
                    &mut tally,
                );
            }
            if !first.full {
                break;
            }
            base = first.end;
            let Some(second) = second else {
                continue;
            };
            if reached(&tally, self) {
                // vsearch stops after the first batch, and the second is
                // dropped
                break;
            }
            note(&mut horizon, &second, self.top_hits);
            if second.delayed > 0 {
                self.align_delayed(
                    db,
                    query,
                    minus,
                    &pending[first.end..second.end],
                    &mut aligned,
                    &mut tally,
                );
            }
            if !second.full {
                break;
            }
            base = second.end;
        }
        tally.hits
    }

    /// takes the candidates of the next batch from `next` into `pending`,
    /// as vsearch's search loop takes them. It stops once `MAX_DELAYED` are
    /// to be aligned, the candidates run out, or those taken since the last
    /// full batch, at `base`, reach `limit`.
    #[allow(clippy::too_many_arguments)]
    fn pull(
        &self,
        next: &mut impl Iterator<Item = Key>,
        pending: &mut Vec<Candidate>,
        base: usize,
        limit: usize,
        db: &Database,
        label: &str,
        query: &[u8],
        size: i64,
    ) -> Batch {
        let mut batch = Batch::default();
        while base + batch.delayed < limit {
            let Some(key) = next.next() else {
                batch.exhausted = true;
                break;
            };
            batch.last = Some(key);
            let t = key.2;
            let mut ok = self.criteria.unaligned(
                label,
                query,
                size,
                &db.labels[t],
                &db.seqs[t],
                db.sizes[t],
            );
            // a part of a query is identical to the target when the whole is
            if let Some(whole) = &self.whole {
                let target = &db.seqs[t];
                ok = ok
                    && !(whole.len() == target.len()
                        && whole
                            .iter()
                            .zip(target)
                            .all(|(&a, &b)| FOUR_BIT[usize::from(a)] == FOUR_BIT[usize::from(b)]));
            }
            if ok {
                batch.delayed += 1;
            }
            pending.push(Candidate {
                target: t,
                rejected: !ok,
            });
            if batch.delayed == MAX_DELAYED {
                batch.full = true;
                break;
            }
        }
        batch.end = pending.len();
        batch
    }

    /// aligns the candidates of `batch` not yet rejected, in one call.
    fn align_batch(&mut self, db: &Database, batch: &[Candidate]) -> Vec<Aligned> {
        let targets: Vec<&[u8]> = batch
            .iter()
            .filter(|c| !c.rejected)
            .map(|c| db.seqs[c.target].as_slice())
            .collect();
        if targets.is_empty() {
            Vec::new()
        } else {
            self.simd.align(&targets)
        }
    }

    /// vsearch's `align_delayed`: counts the batch's accepted and rejected
    /// candidates up to the limits, with `aligned`, the alignments of those
    /// not yet rejected, in order.
    fn align_delayed(
        &mut self,
        db: &Database,
        query: &[u8],
        minus: bool,
        batch: &[Candidate],
        aligned: &mut impl Iterator<Item = Aligned>,
        tally: &mut Tally,
    ) {
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
            let mut verdict =
                self.criteria
                    .aligned(&stats, &cigar, query.len(), target.len(), &self.scoring);
            if let Some(alpha) = self.unoise_alpha
                && verdict == Verdict::Accepted
            {
                let mismatches = stats.mismatches;
                let skew = self.query_size as f64 / db.sizes[c.target] as f64;
                let beta = 1.0 / 2f64.powf(alpha * mismatches as f64 + 1.0);
                if !(skew <= beta || mismatches == 0) {
                    verdict = Verdict::Weak;
                }
            }
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
