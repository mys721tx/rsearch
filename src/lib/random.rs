//! Random numbers as vsearch draws them, so that a seed gives vsearch's
//! output. The module ports `std::mt19937_64` and SplitMix64, Lemire's
//! bounded draw, and vsearch's Fisher-Yates shuffle.

use std::hash::{BuildHasher, Hasher};
use std::time::{SystemTime, UNIX_EPOCH};

/// A source of 64-bit random numbers.
pub trait Generator {
    fn next_u64(&mut self) -> u64;
}

/// The 64-bit Mersenne Twister of C++'s `std::mt19937_64`.
///
/// # Example
/// ```rust
/// # use rsearch::random::{Generator, Mt19937_64};
/// // the check the C++ standard gives: the 10000th number of the default
/// // seed
/// let mut mt = Mt19937_64::new(5489);
/// let tenth_thousandth = (0..10_000).map(|_| mt.next_u64()).last().unwrap();
/// assert_eq!(tenth_thousandth, 9_981_545_732_273_789_042);
/// ```
#[derive(Clone)]
pub struct Mt19937_64 {
    state: [u64; Mt19937_64::N],
    index: usize,
}

impl Mt19937_64 {
    const N: usize = 312;
    const M: usize = 156;
    const MATRIX: u64 = 0xB502_6F5A_A966_19E9;
    const UPPER: u64 = 0xFFFF_FFFF_8000_0000;
    const LOWER: u64 = 0x0000_0000_7FFF_FFFF;

    /// seeds the generator as `std::mt19937_64(seed)` does.
    pub fn new(seed: u64) -> Mt19937_64 {
        let mut state = [0u64; Self::N];
        state[0] = seed;
        for i in 1..Self::N {
            let previous = state[i - 1];
            state[i] = 6_364_136_223_846_793_005u64
                .wrapping_mul(previous ^ (previous >> 62))
                .wrapping_add(i as u64);
        }
        Mt19937_64 {
            state,
            index: Self::N,
        }
    }

    /// draws the next block of numbers.
    fn twist(&mut self) {
        for i in 0..Self::N {
            let x = (self.state[i] & Self::UPPER) | (self.state[(i + 1) % Self::N] & Self::LOWER);
            let mut next = x >> 1;
            if x & 1 != 0 {
                next ^= Self::MATRIX;
            }
            self.state[i] = self.state[(i + Self::M) % Self::N] ^ next;
        }
        self.index = 0;
    }
}

impl Generator for Mt19937_64 {
    fn next_u64(&mut self) -> u64 {
        if self.index >= Self::N {
            self.twist();
        }
        let mut x = self.state[self.index];
        self.index += 1;
        x ^= (x >> 29) & 0x5555_5555_5555_5555;
        x ^= (x << 17) & 0x71D6_7FFF_EDA6_0000;
        x ^= (x << 37) & 0xFFF7_EEE0_0000_0000;
        x ^ (x >> 43)
    }
}

/// SplitMix64, as vsearch's: a counter stepped by the golden ratio, then
/// mixed.
///
/// # Example
/// ```rust
/// # use rsearch::random::{Generator, SplitMix64};
/// let mut splitmix = SplitMix64::new(0);
/// assert_eq!(splitmix.next_u64(), 0xE220_A839_7B1D_CDAF);
/// ```
#[derive(Clone)]
pub struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    pub fn new(seed: u64) -> SplitMix64 {
        SplitMix64 { state: seed }
    }
}

impl Generator for SplitMix64 {
    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
}

/// returns an unbiased number in `[0, range)`, as vsearch's
/// `random_bounded`: Lemire's multiply-shift with rejection. `range` must
/// not be 0.
pub fn bounded(generator: &mut impl Generator, range: u64) -> u64 {
    assert!(range > 0, "a bounded draw needs a range");
    let mut product = u128::from(generator.next_u64()) * u128::from(range);
    let mut low = product as u64;
    if low < range {
        let threshold = range.wrapping_neg() % range;
        while low < threshold {
            product = u128::from(generator.next_u64()) * u128::from(range);
            low = product as u64;
        }
    }
    (product >> 64) as u64
}

/// shuffles `data` in place, as vsearch's `random_shuffle`: Fisher-Yates
/// from the last element down.
///
/// # Example
/// ```rust
/// # use rsearch::random::{Mt19937_64, shuffle};
/// let mut deck: Vec<u32> = (0..8).collect();
/// shuffle(&mut deck, &mut Mt19937_64::new(1));
/// // vsearch --shuffle --randseed 1 on eight records
/// assert_eq!(deck, [3, 4, 5, 7, 6, 2, 0, 1]);
/// ```
pub fn shuffle<T>(data: &mut [T], generator: &mut impl Generator) {
    for i in (2..=data.len()).rev() {
        let j = bounded(generator, i as u64) as usize;
        data.swap(i - 1, j);
    }
}

/// draws `n` of the reads that `masses` count per record, as vsearch's
/// `--fastx_subsample`. It uses selection sampling over the reads in order,
/// each kept when a draw from the reads left falls below the number still
/// wanted. Returns how many reads of each record were drawn. `n` must not
/// exceed the total; when it equals it, every read is drawn without a draw.
///
/// # Example
/// ```rust
/// # use rsearch::random::{Mt19937_64, draw};
/// // vsearch --fastx_subsample --sizein --randseed 1 --sample_size 3
/// let drawn = draw(&[5, 1, 1, 2, 10, 1], 3, &mut Mt19937_64::new(1));
/// assert_eq!(drawn, [2, 0, 0, 1, 0, 0]);
/// ```
pub fn draw(masses: &[u64], n: u64, generator: &mut impl Generator) -> Vec<u64> {
    let total: u64 = masses.iter().sum();
    assert!(n <= total, "cannot draw more reads than there are");
    if n == total {
        return masses.to_vec();
    }
    let mut drawn = vec![0; masses.len()];
    let (mut left, mut checked) = (n, 0);
    let (mut record, mut in_record) = (0, 0);
    while left > 0 {
        if bounded(generator, total - checked) < left {
            drawn[record] += 1;
            left -= 1;
        }
        checked += 1;
        in_record += 1;
        if in_record >= masses[record] {
            record += 1;
            in_record = 0;
        }
    }
    drawn
}

/// returns a well-separated seed for item `index` of a run seeded with
/// `seed`, as vsearch's `RandomSeed::substream`, so that each item's
/// numbers depend only on the seed and the item's number.
///
/// # Example
/// ```rust
/// # use rsearch::random::{SplitMix64, Generator, substream};
/// let mut mixer = SplitMix64::new(1 ^ 0x9E37_79B9_7F4A_7C15);
/// assert_eq!(substream(1, 1), mixer.next_u64());
/// ```
pub fn substream(seed: u64, index: u64) -> u64 {
    SplitMix64::new(seed ^ index.wrapping_mul(0x9E37_79B9_7F4A_7C15)).next_u64()
}

/// returns the seed of a run: the one given, used as its unsigned bit
/// pattern as vsearch does, or for none or 0 one from the system.
pub fn seed(given: Option<i64>) -> u64 {
    match given {
        Some(seed) if seed != 0 => seed as u64,
        _ => system_seed(),
    }
}

/// a seed from the system: the random keys the standard library draws for
/// hashing, mixed with the time.
fn system_seed() -> u64 {
    let mut hasher = std::collections::hash_map::RandomState::new().build_hasher();
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |d| d.as_nanos());
    hasher.write_u128(nanos);
    match hasher.finish() {
        // 0 asks for a system seed and is never handed out
        0 => 1,
        seed => seed,
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_seed() {
        assert_eq!(seed(Some(42)), 42);
        assert_eq!(seed(Some(-1)), u64::MAX);
        assert_ne!(seed(None), 0);
        assert_ne!(seed(Some(0)), 0);
    }

    #[test]
    fn test_bounded() {
        let mut mt = Mt19937_64::new(7);
        for range in [1, 2, 3, 10, 1 << 40, u64::MAX] {
            assert!(bounded(&mut mt, range) < range);
        }
    }
}
