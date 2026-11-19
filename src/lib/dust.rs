//! DUST low-complexity masking, as vsearch's. Windows of 64 residues
//! stepped by 32 are each searched for their lowest-complexity region of
//! words of three. The region is masked when it scores above 20.

/// The residues searched at a time; windows overlap by half of it.
const WINDOW: usize = 64;
const HALF_WINDOW: usize = WINDOW / 2;
/// The score a region must beat to be masked.
const LEVEL: i32 = 20;
/// The length of a word.
const WORD: usize = 3;
/// The number of distinct words.
const WORDS: usize = 1 << (2 * WORD);

/// The 2-bit code DUST reads each residue as: A, C, G and T or U, any other
/// residue counting as A.
const CODE: [u8; 256] = {
    let mut code = [0u8; 256];
    code[b'C' as usize] = 1;
    code[b'c' as usize] = 1;
    code[b'G' as usize] = 2;
    code[b'g' as usize] = 2;
    code[b'T' as usize] = 3;
    code[b't' as usize] = 3;
    code[b'U' as usize] = 3;
    code[b'u' as usize] = 3;
    code
};

/// How masked residues are marked.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Style {
    /// The sequence is upper-cased, then each region lower-cased.
    Soft,
    /// Each region becomes N; the rest is left as it was.
    Hard,
}

/// The lowest-complexity region of a window: its score, and its first and
/// last positions in the window.
#[derive(Default)]
struct Region {
    score: i32,
    begin: usize,
    end: usize,
}

/// Finds the lowest-complexity region of `window`, at most `WINDOW` long.
///
/// A region from `a` to `b` scores `10 * pairs / (b - a - 1)` (integer
/// division), where `pairs` counts the pairs of equal words in it. The
/// first region of the best score wins, by start then end. Starts that
/// cannot beat `LEVEL` are skipped by vsearch's bound. By that bound, a
/// region's pairs are at most the sum, over its words but the first, of how
/// often each occurs earlier in the window.
fn worst_region(window: &[u8]) -> Region {
    let length = window.len();
    // the shortest region holds eight residues
    if length < WORD + 5 {
        return Region::default();
    }
    let starts = length - WORD + 1 - 5;
    let mut words = [0u8; WINDOW];
    let mut counts = [0u8; WORDS];
    // the excess of each word over the cost of a position, then the most a
    // run of them from there can reach
    let mut reach = [0i32; WINDOW];
    let mut word = 0usize;
    for (position, &residue) in window.iter().enumerate() {
        word = ((word << 2) | usize::from(CODE[usize::from(residue)])) & (WORDS - 1);
        words[position] = word as u8;
        if position >= WORD - 1 {
            reach[position] = 10 * i32::from(counts[word]) - (LEVEL + 1);
            counts[word] += 1;
        }
    }
    // Kadane, right to left
    for position in (WORD - 1..length - 1).rev() {
        let next = reach[position + 1];
        if next > 0 {
            reach[position] += next;
        }
    }

    let mut best = Region::default();
    for start in 0..starts {
        if reach[start + WORD] < 2 * (LEVEL + 1) {
            continue;
        }
        counts = [0; WORDS];
        let mut sum = 0i32;
        for offset in WORD - 1..length - start {
            let word = usize::from(words[start + offset]);
            let repeats = counts[word];
            if repeats != 0 {
                sum += i32::from(repeats);
                let score = 10 * sum / offset as i32;
                if score > best.score {
                    best = Region {
                        score,
                        begin: start,
                        end: start + offset,
                    };
                }
            }
            counts[word] = repeats + 1;
        }
    }
    best
}

/// Masks the low-complexity regions of `seq` in place with DUST, as
/// vsearch's `--qmask dust`. Soft masking upper-cases the sequence first,
/// then lower-cases each region as it was read. Hard masking writes N over
/// each region and leaves the rest alone.
///
/// # Example
/// ```rust
/// # use rsearch::dust::{dust, Style};
/// let mut seq = b"acgtACGTacgtNNNNNNNNNNNNNNNNNNNNacgt".to_vec();
/// dust(&mut seq, Style::Soft);
/// assert_eq!(seq, b"ACGTACGTACGTnnnnnnnnnnnnnnnnnnnnaCGT");
/// ```
pub fn dust(seq: &mut [u8], style: Style) {
    let original = seq.to_vec();
    if style == Style::Soft {
        seq.make_ascii_uppercase();
    }
    let mut i = 0;
    while i < seq.len() {
        let length = (seq.len() - i).min(WINDOW);
        let worst = worst_region(&original[i..i + length]);
        if worst.score > LEVEL {
            let region = i + worst.begin..=i + worst.end;
            match style {
                Style::Hard => seq[region].fill(b'N'),
                Style::Soft => {
                    for position in region {
                        seq[position] = original[position] | 0x20;
                    }
                }
            }
            // a region ending early in the window rewinds the next one
            if worst.end < HALF_WINDOW {
                i += HALF_WINDOW - worst.end;
            }
        }
        i += HALF_WINDOW;
    }
}

/// Writes N over every lowercase residue of `seq`, as vsearch's `--qmask
/// soft --hardmask`: every byte with the case bit set.
///
/// # Example
/// ```rust
/// # use rsearch::dust::hard_mask_lowercase;
/// let mut seq = b"acGTnN".to_vec();
/// hard_mask_lowercase(&mut seq);
/// assert_eq!(seq, b"NNGTNN");
/// ```
pub fn hard_mask_lowercase(seq: &mut [u8]) {
    for residue in seq {
        if *residue & 0x20 != 0 {
            *residue = b'N';
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_short_windows() {
        // seven residues never hold a region; eight can
        let mut seq = b"AAAAAAA".to_vec();
        dust(&mut seq, Style::Soft);
        assert_eq!(seq, b"AAAAAAA");
        let mut seq = b"AAAAAAAA".to_vec();
        dust(&mut seq, Style::Hard);
        assert_eq!(seq, b"NNNNNNNN");
        let mut seq = Vec::new();
        dust(&mut seq, Style::Soft);
        assert!(seq.is_empty());
    }

    #[test]
    fn test_hard_keeps_case() {
        let mut seq = b"acgtacgtGGGGGGGGGGGGGGGGtgca".to_vec();
        dust(&mut seq, Style::Hard);
        assert_eq!(seq, b"acgtacgtNNNNNNNNNNNNNNNNtgca");
    }
}
