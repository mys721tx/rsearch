//! Words: the distinct k-mers of a sequence, as vsearch's `Uniquer`, and an
//! index of the sequences holding each word.
//!
//! A word is `k` residues of A, C, G, T or U in upper case; a word holding
//! anything else, lower case (masked) residues included, is skipped. A word
//! is coded on two bits a residue, A 0, C 1, G 2, T and U 3, the first
//! residue in the highest bits.

/// The shortest and the longest word lengths vsearch takes.
pub const MIN_WORD_LEN: u32 = 3;
pub const MAX_WORD_LEN: u32 = 15;

/// The 2-bit code of each byte, plus 4 for the bytes no word may hold.
const CODE: [u8; 256] = {
    let mut code = [4u8; 256];
    code[b'A' as usize] = 0;
    code[b'C' as usize] = 1;
    code[b'G' as usize] = 2;
    code[b'T' as usize] = 3;
    code[b'U' as usize] = 3;
    code
};

/// The word lengths up to which the words seen are marked in a table of
/// every word, rather than sorted out. A table of 4^12 entries is 64 MiB.
const MAX_TABLE_WORD_LEN: u32 = 12;

/// Finds the distinct words of sequences, in the order they first occur.
///
/// # Example
/// ```rust
/// # use rsearch::words::Uniquer;
/// let mut uniquer = Uniquer::new(3);
/// // AAA twice, and no word across the lower case c or the N
/// assert_eq!(uniquer.words(b"AAAAcGTUNACG"), &[0b000000, 0b101111, 0b000110]);
/// assert!(uniquer.words(b"AC").is_empty());
/// ```
pub struct Uniquer {
    word_len: u32,
    /// The words found in the last sequence.
    list: Vec<u32>,
    /// For each word, the number of the last sequence it was found in,
    /// plus one; empty above `MAX_TABLE_WORD_LEN`.
    stamps: Vec<u32>,
    /// The words seen in the last sequence, above `MAX_TABLE_WORD_LEN`.
    seen: std::collections::HashSet<u32>,
    /// The number of sequences read, plus one.
    stamp: u32,
}

impl Uniquer {
    pub fn new(word_len: u32) -> Uniquer {
        assert!((MIN_WORD_LEN..=MAX_WORD_LEN).contains(&word_len));
        let stamps = if word_len <= MAX_TABLE_WORD_LEN {
            vec![0; 1 << (2 * word_len)]
        } else {
            Vec::new()
        };
        Uniquer {
            word_len,
            list: Vec::new(),
            stamps,
            seen: std::collections::HashSet::new(),
            stamp: 0,
        }
    }

    pub fn word_len(&self) -> u32 {
        self.word_len
    }

    /// the distinct words of `seq`, in the order of their first
    /// occurrence, as vsearch lists them.
    pub fn words(&mut self, seq: &[u8]) -> &[u32] {
        self.list.clear();
        let k = self.word_len as usize;
        let mask = (1u32 << (2 * k)) - 1;
        let mut word = 0u32;
        // the residues read since the last one no word may hold
        let mut run = 0usize;
        if self.stamps.is_empty() {
            for &residue in seq {
                let code = CODE[usize::from(residue)];
                if code > 3 {
                    run = 0;
                    continue;
                }
                word = ((word << 2) | u32::from(code)) & mask;
                run += 1;
                if run >= k && self.seen.insert(word) {
                    self.list.push(word);
                }
            }
            self.seen.clear();
            return &self.list;
        }
        self.stamp = self.stamp.wrapping_add(1);
        if self.stamp == 0 {
            self.stamps.fill(0);
            self.stamp = 1;
        }
        let stamp = self.stamp;
        for &residue in seq {
            let code = CODE[usize::from(residue)];
            if code > 3 {
                run = 0;
                continue;
            }
            word = ((word << 2) | u32::from(code)) & mask;
            run += 1;
            if run >= k {
                let seen = &mut self.stamps[word as usize];
                if *seen != stamp {
                    *seen = stamp;
                    self.list.push(word);
                }
            }
        }
        &self.list
    }
}

/// The index of a database: for each word, the numbers of the sequences
/// holding it, in increasing order.
pub struct Index {
    word_len: u32,
    /// For each word, where its sequences start in `sequences`; one more
    /// entry ends the last.
    starts: Vec<u64>,
    sequences: Vec<u32>,
    /// The sequences holding no word.
    without_words: u64,
}

impl Index {
    /// indexes `seqs`, numbered from 0 in order.
    ///
    /// # Example
    /// ```rust
    /// # use rsearch::words::Index;
    /// let index = Index::new(3, [&b"ACGT"[..], b"GTA", b"AC"]);
    /// // ACG (6) in the first; GTA (44) in the second
    /// assert_eq!(index.sequences_of(6), &[0]);
    /// assert_eq!(index.sequences_of(44), &[1]);
    /// assert_eq!(index.counts().iter().sum::<u32>(), 3);
    /// assert_eq!(index.without_words(), 1);
    /// ```
    pub fn new<'a>(word_len: u32, seqs: impl IntoIterator<Item = &'a [u8]>) -> Index {
        let mut uniquer = Uniquer::new(word_len);
        let words = 1usize << (2 * word_len);
        // the words of every sequence, each list ended by its length
        let mut found: Vec<u32> = Vec::new();
        let mut counts = vec![0u32; words];
        let mut seq_words: Vec<u32> = Vec::new();
        let mut without_words = 0;
        for seq in seqs {
            let list = uniquer.words(seq);
            if list.is_empty() {
                without_words += 1;
            }
            for &word in list {
                counts[word as usize] += 1;
            }
            found.extend_from_slice(list);
            seq_words.push(list.len() as u32);
        }
        let mut starts = Vec::with_capacity(words + 1);
        let mut sum = 0u64;
        for &count in &counts {
            starts.push(sum);
            sum += u64::from(count);
        }
        starts.push(sum);
        // filled in sequence order, which keeps each list increasing
        let mut next: Vec<u64> = starts[..words].to_vec();
        let mut sequences = vec![0u32; sum as usize];
        let mut at = 0;
        for (number, &n) in seq_words.iter().enumerate() {
            for &word in &found[at..at + n as usize] {
                let slot = &mut next[word as usize];
                sequences[*slot as usize] = number as u32;
                *slot += 1;
            }
            at += n as usize;
        }
        Index {
            word_len,
            starts,
            sequences,
            without_words,
        }
    }

    /// rebuilds the index of a UDB file from its counts and its word list.
    pub fn from_parts(word_len: u32, counts: &[u32], sequences: Vec<u32>) -> Index {
        let mut starts = Vec::with_capacity(counts.len() + 1);
        let mut sum = 0u64;
        for &count in counts {
            starts.push(sum);
            sum += u64::from(count);
        }
        starts.push(sum);
        assert_eq!(sum, sequences.len() as u64);
        Index {
            word_len,
            starts,
            sequences,
            without_words: 0,
        }
    }

    pub fn word_len(&self) -> u32 {
        self.word_len
    }

    /// the number of sequences holding each word.
    pub fn counts(&self) -> Vec<u32> {
        self.starts
            .windows(2)
            .map(|pair| (pair[1] - pair[0]) as u32)
            .collect()
    }

    /// the numbers of the sequences holding `word`.
    pub fn sequences_of(&self, word: u32) -> &[u32] {
        let word = word as usize;
        &self.sequences[self.starts[word] as usize..self.starts[word + 1] as usize]
    }

    /// the sequence numbers of every word, in word order.
    pub fn all_sequences(&self) -> &[u32] {
        &self.sequences
    }

    /// the number of sequences holding no word.
    pub fn without_words(&self) -> u64 {
        self.without_words
    }

    /// the start of each word's sequences, one more ending the last, and
    /// the sequences.
    pub fn into_parts(self) -> (Vec<u64>, Vec<u32>) {
        (self.starts, self.sequences)
    }
}

/// counts, for each word, the sequences of `seqs` holding it, as vsearch's
/// database index does; and the sequences holding no word.
///
/// # Example
/// ```rust
/// # use rsearch::words::word_counts;
/// let (counts, without) = word_counts(3, [&b"ACGTACG"[..], b"ACG", b"AC"]);
/// // ACG (6) in two sequences, CGT (27) in one
/// assert_eq!((counts[6], counts[27], without), (2, 1, 1));
/// ```
pub fn word_counts<'a>(word_len: u32, seqs: impl IntoIterator<Item = &'a [u8]>) -> (Vec<u32>, u64) {
    let mut uniquer = Uniquer::new(word_len);
    let mut counts = vec![0u32; 1 << (2 * word_len)];
    let mut without_words = 0;
    for seq in seqs {
        let list = uniquer.words(seq);
        if list.is_empty() {
            without_words += 1;
        }
        for &word in list {
            counts[word as usize] += 1;
        }
    }
    (counts, without_words)
}

/// the reverse complement of `word`, of `word_len` residues.
///
/// # Example
/// ```rust
/// # use rsearch::words::reverse_complement;
/// // AAC and GTT
/// assert_eq!(reverse_complement(0b000001, 3), 0b101111);
/// // ACGT is its own
/// assert_eq!(reverse_complement(0b00011011, 4), 0b00011011);
/// ```
pub fn reverse_complement(word: u32, word_len: u32) -> u32 {
    let mut word = word;
    let mut reverse = 0;
    for _ in 0..word_len {
        reverse = (reverse << 2) | (3 - (word & 3));
        word >>= 2;
    }
    reverse
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A uniquer makes its table of every word on its first sequence. A
    /// worker that gets no sequence then holds no table. dbornt's table of
    /// 4^12 words is 64 MiB, and a thousand idle workers ran Windows out of
    /// memory.
    #[test]
    fn makes_the_table_on_the_first_sequence() {
        let mut uniquer = Uniquer::new(MAX_TABLE_WORD_LEN);
        assert_eq!(uniquer.stamps.capacity(), 0);
        assert_eq!(uniquer.words(b"ACGTACGTACGTA"), &[0x1b1b1b, 0x6c6c6c]);
        assert_eq!(uniquer.stamps.len(), 1 << (2 * MAX_TABLE_WORD_LEN));
        let mut longer = Uniquer::new(MAX_TABLE_WORD_LEN + 1);
        assert_eq!(longer.words(b"ACGTACGTACGTA"), &[0x6c6c6c]);
        assert_eq!(longer.stamps.capacity(), 0);
    }
}
