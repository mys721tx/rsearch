//! SINTAX classification, as vsearch's `--sintax`. A hundred times, 32
//! words drawn at random from the query's distinct words pick the database
//! sequence sharing the most of them. With at least fifty such picks, each
//! taxonomic level takes the name most of the picks agree on, among the
//! picks that agreed at the levels above.

use crate::cluster::reverse_complement_keeping_case;
use crate::random::{SplitMix64, bounded};
use crate::taxonomy::{LEVELS, split};
use crate::words::{Index, Uniquer};

/// The words drawn for one pick, and the picks.
const SUBSET: usize = 32;
const BOOTSTRAPS: usize = 100;

/// A word held by at least one in this many database sequences is kept as
/// a bitmap of the sequences too, as vsearch's database index does.
const BITMAP_SHARE: usize = 8;

/// A database: its headers, sequence lengths and word index.
pub struct Database {
    pub headers: Vec<Vec<u8>>,
    pub lengths: Vec<u32>,
    pub index: Index,
    /// the words with a bitmap, in increasing order.
    common: Vec<u32>,
    /// the bitmaps of the common words, a bit per sequence, the first in
    /// the lowest bit; `bitmap_len` words each.
    bitmaps: Vec<u64>,
    bitmap_len: usize,
}

impl Database {
    /// a database of `headers` and `lengths`, indexed by `index`.
    pub fn new(headers: Vec<Vec<u8>>, lengths: Vec<u32>, index: Index) -> Database {
        let n = lengths.len();
        let bitmap_len = n.div_ceil(64);
        let threshold = (n / BITMAP_SHARE).max(1);
        let words = 1usize << (2 * index.word_len());
        let mut common = Vec::new();
        let mut bitmaps = Vec::new();
        for word in 0..words as u32 {
            let sequences = index.sequences_of(word);
            if sequences.len() >= threshold {
                let start = bitmaps.len();
                bitmaps.resize(start + bitmap_len, 0);
                for &s in sequences {
                    bitmaps[start + s as usize / 64] |= 1 << (s % 64);
                }
                common.push(word);
            }
        }
        Database {
            headers,
            lengths,
            index,
            common,
            bitmaps,
            bitmap_len,
        }
    }

    /// the bitmap of `word`, if it has one.
    fn bitmap(&self, word: u32) -> Option<&[u64]> {
        let slot = self.common.binary_search(&word).ok()?;
        let start = slot * self.bitmap_len;
        Some(&self.bitmaps[start..start + self.bitmap_len])
    }
}

/// For each byte of a bitmap, the eight byte counters it raises by one, as one
/// word to add. Bit `i` of the byte becomes byte `i` of the word.
const EXPANSION: [u64; 256] = {
    let mut table = [0u64; 256];
    let mut byte = 0;
    while byte < 256 {
        let mut bit = 0;
        while bit < 8 {
            if byte & (1 << bit) != 0 {
                table[byte] |= 1 << (8 * bit);
            }
            bit += 1;
        }
        byte += 1;
    }
    table
};

/// A byte with its high bit alone set, in each byte of a word.
const HIGH_BITS: u64 = 0x8080_8080_8080_8080;
/// One, in each byte of a word.
const ONES: u64 = 0x0101_0101_0101_0101;

/// The classification of a query.
pub struct Classification {
    /// The strand that won: `false` for plus, `true` for minus.
    pub minus: bool,
    /// The database sequences picked, on the winning strand.
    pub picks: Vec<u32>,
}

/// Classifies queries against a database.
pub struct Classifier<'a> {
    db: &'a Database,
    both_strands: bool,
    random_ties: bool,
    uniquer: Uniquer,
    /// a byte counter per database sequence, eight to a word, sequence `s` in
    /// byte `s % 8` (from the lowest) of word `s / 8`, padded to a whole bitmap
    /// word. A pick draws at most 32 words. No count thus carries into the next
    /// byte.
    counters: Vec<u64>,
    drawn: Vec<bool>,
    rc: Vec<u8>,
}

impl<'a> Classifier<'a> {
    pub fn new(db: &'a Database, both_strands: bool, random_ties: bool) -> Classifier<'a> {
        Classifier {
            db,
            both_strands,
            random_ties,
            uniquer: Uniquer::new(db.index.word_len()),
            counters: vec![0; db.bitmap_len * 8],
            drawn: Vec::new(),
            rc: Vec::new(),
        }
    }

    /// picks database sequences for `seq` on each strand searched, drawing
    /// from `rng`, and keeps the strand with the better best count, or the
    /// more picks.
    pub fn classify(&mut self, seq: &[u8], rng: &mut SplitMix64) -> Classification {
        let mut strands = [(Vec::new(), 0u32), (Vec::new(), 0u32)];
        let n = if self.both_strands { 2 } else { 1 };
        for (s, strand) in strands.iter_mut().enumerate().take(n) {
            if s == 1 {
                reverse_complement_keeping_case(&mut self.rc, seq);
            }
            let query: &[u8] = if s == 0 { seq } else { &self.rc };
            let words = self.uniquer.words(query).to_vec();
            if words.len() < SUBSET {
                continue;
            }
            let mut subset = Vec::with_capacity(SUBSET);
            for _ in 0..BOOTSTRAPS {
                // 32 draws, each word kept once
                subset.clear();
                self.drawn.clear();
                self.drawn.resize(words.len(), false);
                for _ in 0..SUBSET {
                    let x = bounded(rng, words.len() as u64) as usize;
                    if !self.drawn[x] {
                        self.drawn[x] = true;
                        subset.push(words[x]);
                    }
                }
                if let Some((seqno, count)) = self.best(&subset, rng) {
                    strand.0.push(seqno);
                    strand.1 = strand.1.max(count);
                }
            }
        }
        let minus = self.both_strands
            && (strands[1].1 > strands[0].1
                || (strands[1].1 == strands[0].1 && strands[1].0.len() > strands[0].0.len()));
        let [plus, minus_strand] = strands;
        Classification {
            minus,
            picks: if minus { minus_strand.0 } else { plus.0 },
        }
    }

    /// the database sequence sharing the most of `words`, and how many; the
    /// shortest, then the first, among ties, or one at random with
    /// `random_ties`; none unless it shares two or more.
    fn best(&mut self, words: &[u32], rng: &mut SplitMix64) -> Option<(u32, u32)> {
        let db = self.db;
        self.counters.fill(0);
        for &word in words {
            if let Some(bitmap) = db.bitmap(word) {
                // eight counters at once, a byte of the bitmap each
                for (bits, counters) in bitmap.iter().zip(self.counters.chunks_exact_mut(8)) {
                    if *bits == 0 {
                        continue;
                    }
                    for (i, counter) in counters.iter_mut().enumerate() {
                        *counter += EXPANSION[((bits >> (8 * i)) & 0xff) as usize];
                    }
                }
            } else {
                for &seqno in db.index.sequences_of(word) {
                    self.counters[seqno as usize / 8] += 1 << (8 * (seqno % 8));
                }
            }
        }
        let n = db.lengths.len();
        let (mut best_count, mut best_seqno, mut best_length) = (0u32, 0u32, 0u32);
        let mut ties = 0u64;
        for (at, &counters) in self.counters.iter().enumerate() {
            // the counters of the word that reach the best so far. No count
            // exceeds 32. Adding 128 less the best thus sets a byte's high bit
            // exactly when it does, and carries into no other byte
            let mut reaching =
                counters.wrapping_add(ONES * (128 - u64::from(best_count))) & HIGH_BITS;
            while reaching != 0 {
                let byte = reaching.trailing_zeros() as usize / 8;
                reaching &= reaching - 1;
                let seqno = at * 8 + byte;
                if seqno >= n {
                    break;
                }
                let count = ((self.counters[at] >> (8 * byte)) & 0xff) as u32;
                if count < best_count {
                    continue;
                }
                let length = db.lengths[seqno];
                let seqno = seqno as u32;
                if count > best_count {
                    (best_count, best_seqno, best_length) = (count, seqno, length);
                    ties = 1;
                } else if self.random_ties {
                    ties += 1;
                    if bounded(rng, ties) == 0 {
                        (best_seqno, best_length) = (seqno, length);
                    }
                } else if length < best_length {
                    (best_seqno, best_length) = (seqno, length);
                } else if length == best_length {
                    best_seqno = best_seqno.min(seqno);
                }
            }
        }
        (best_count > 1).then_some((best_seqno, best_count))
    }
}

/// vsearch's table line for a query, after its header and a tab: the
/// names at each level with the share of picks agreeing, the strand, and
/// with a `cutoff` above 0 the names whose share reaches it; or empty
/// fields with fewer than fifty picks.
pub fn table_fields(db: &Database, classification: &Classification, cutoff: f64) -> Vec<u8> {
    let picks = &classification.picks;
    let mut out = Vec::new();
    if picks.len() < BOOTSTRAPS.div_ceil(2) {
        out.extend_from_slice(if cutoff > 0.0 { b"\t\t" } else { b"\t" });
        return out;
    }
    let names: Vec<[&[u8]; 9]> = picks
        .iter()
        .map(|&seqno| {
            let header = &db.headers[seqno as usize];
            split(header).map(|range| &header[range])
        })
        .collect();
    let mut included = vec![true; picks.len()];
    let mut best = [0usize; 9];
    let mut agreeing = [0usize; 9];
    for level in 0..LEVELS.len() {
        // each pick joins the first included pick naming the same
        let mut group = vec![usize::MAX; picks.len()];
        let mut members = vec![0usize; picks.len()];
        for i in 0..picks.len() {
            if included[i] {
                let first = (0..=i)
                    .find(|&j| included[j] && names[i][level] == names[j][level])
                    .unwrap();
                group[i] = first;
                members[first] += 1;
            }
        }
        // the largest group, the first among equals
        let (mut winner, mut size) = (0, 0);
        for (i, &m) in members.iter().enumerate() {
            if m > size {
                (winner, size) = (i, m);
            }
        }
        best[level] = winner;
        agreeing[level] = size;
        for i in 0..picks.len() {
            if group[i] != winner {
                included[i] = false;
            }
        }
    }
    let share = |level: usize| agreeing[level] as f64 / picks.len() as f64;
    let mut first = true;
    for level in 0..LEVELS.len() {
        let name = names[best[level]][level];
        if !name.is_empty() {
            if !first {
                out.push(b',');
            }
            first = false;
            out.push(LEVELS[level]);
            out.push(b':');
            out.extend_from_slice(name);
            out.extend_from_slice(format!("({:.2})", share(level)).as_bytes());
        }
    }
    out.push(b'\t');
    out.push(if classification.minus { b'-' } else { b'+' });
    if cutoff > 0.0 {
        out.push(b'\t');
        let mut first = true;
        for level in 0..LEVELS.len() {
            let name = names[best[level]][level];
            if !name.is_empty() && share(level) >= cutoff {
                if !first {
                    out.push(b',');
                }
                first = false;
                out.push(LEVELS[level]);
                out.push(b':');
                out.extend_from_slice(name);
            }
        }
    }
    out
}

/// a generator for query number `query`, as vsearch's per-query streams.
pub fn query_rng(seed: u64, query: u64) -> SplitMix64 {
    SplitMix64::new(crate::random::substream(seed, query))
}
