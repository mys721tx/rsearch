//! Orientation by database words, as vsearch's `--orient`. Each distinct
//! word of a query votes for the strand on which it is held by more than
//! eight times as many database sequences as on the other. A query is
//! oriented when one strand has at least one vote and four times the votes
//! of the other.

use crate::words::{Uniquer, reverse_complement};

/// The strand a query is found on.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
    Undetermined,
}

impl Strand {
    /// vsearch's symbol for the strand, in `--tabbedout`.
    pub fn symbol(self) -> char {
        match self {
            Strand::Forward => '+',
            Strand::Reverse => '-',
            Strand::Undetermined => '?',
        }
    }
}

/// How many times more database sequences a word must be held by on one
/// strand than on the other to vote for it.
const HITS_FACTOR: u32 = 8;
/// How many times more votes one strand needs than the other.
const VOTES_FACTOR: u32 = 4;

/// Orients queries against the word counts of a database.
///
/// # Example
/// ```rust
/// # use rsearch::orient::{Orienter, Strand};
/// # use rsearch::words::word_counts;
/// let (counts, _) = word_counts(3, [&b"AACC"[..]]);
/// let mut orienter = Orienter::new(3, counts);
/// assert_eq!(orienter.orient(b"AAC"), (Strand::Forward, 1, 0));
/// assert_eq!(orienter.orient(b"GGTT"), (Strand::Reverse, 0, 2));
/// assert_eq!(orienter.orient(b"aac"), (Strand::Undetermined, 0, 0));
/// ```
pub struct Orienter {
    word_len: u32,
    counts: Vec<u32>,
    uniquer: Uniquer,
}

impl Orienter {
    /// `counts` holds, for each word of `word_len` residues, the number of
    /// database sequences holding it.
    pub fn new(word_len: u32, counts: Vec<u32>) -> Orienter {
        assert_eq!(counts.len(), 1 << (2 * word_len));
        Orienter {
            word_len,
            counts,
            uniquer: Uniquer::new(word_len),
        }
    }

    /// the strand of `seq`, and the votes for each strand.
    pub fn orient(&mut self, seq: &[u8]) -> (Strand, u32, u32) {
        let (mut forward, mut reverse) = (0u32, 0u32);
        for &word in self.uniquer.words(seq) {
            let hits = self.counts[word as usize];
            let hits_rev = self.counts[reverse_complement(word, self.word_len) as usize];
            // vsearch's unsigned arithmetic
            if hits > HITS_FACTOR.wrapping_mul(hits_rev) {
                forward += 1;
            } else if hits_rev > HITS_FACTOR.wrapping_mul(hits) {
                reverse += 1;
            }
        }
        let strand = if forward >= 1 && forward >= VOTES_FACTOR.wrapping_mul(reverse) {
            Strand::Forward
        } else if reverse >= 1 && reverse >= VOTES_FACTOR.wrapping_mul(forward) {
            Strand::Reverse
        } else {
            Strand::Undetermined
        };
        (strand, forward, reverse)
    }
}
