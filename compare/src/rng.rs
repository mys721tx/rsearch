//! A small seeded random generator (SplitMix64) with the helpers the input
//! generators use, named after Python's `random` module they replace.
//! A seed always gives the same stream, on every platform.

/// A seeded stream of random numbers.
#[derive(Debug, Clone)]
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Rng {
        // mix the seed so that nearby seeds start far apart
        let mut r = Rng(seed ^ 0x6a09_e667_f3bc_c909);
        r.next_u64();
        r
    }

    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    /// A number below `n` (`n` > 0), without modulo bias.
    pub fn below(&mut self, n: u64) -> u64 {
        assert!(n > 0, "below(0)");
        let zone = u64::MAX - u64::MAX % n;
        loop {
            let x = self.next_u64();
            if x < zone {
                return x % n;
            }
        }
    }

    /// A float in [0, 1).
    pub fn random(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// An integer from `a` to `b`, both included (Python's `randint`).
    pub fn randint(&mut self, a: i64, b: i64) -> i64 {
        assert!(a <= b, "randint({a}, {b})");
        let span = (b as i128 - a as i128 + 1) as u128;
        if span > u64::MAX as u128 {
            return self.next_u64() as i64;
        }
        (a as i128 + self.below(span as u64) as i128) as i64
    }

    /// An index below `n` (Python's `randrange(n)`).
    pub fn randrange(&mut self, n: usize) -> usize {
        self.below(n as u64) as usize
    }

    /// A float from `a` to `b` (Python's `uniform`).
    pub fn uniform(&mut self, a: f64, b: f64) -> f64 {
        a + (b - a) * self.random()
    }

    /// `k` random bits.
    pub fn getrandbits(&mut self, k: u32) -> u64 {
        if k >= 64 {
            self.next_u64()
        } else {
            self.next_u64() & ((1u64 << k) - 1)
        }
    }

    /// A random element.
    pub fn choice<T: Clone>(&mut self, items: &[T]) -> T {
        items[self.randrange(items.len())].clone()
    }

    /// A random byte of an alphabet.
    pub fn ch(&mut self, alphabet: &str) -> u8 {
        let b = alphabet.as_bytes();
        b[self.randrange(b.len())]
    }

    /// `n` random bytes of an alphabet, as text.
    pub fn chars(&mut self, alphabet: &str, n: usize) -> String {
        (0..n).map(|_| self.ch(alphabet) as char).collect()
    }

    /// A random element by weight (Python's `choices(items, weights)[0]`).
    pub fn weighted<T: Clone>(&mut self, items: &[T], weights: &[f64]) -> T {
        let total: f64 = weights.iter().sum();
        let mut x = self.random() * total;
        for (item, w) in items.iter().zip(weights) {
            if x < *w {
                return item.clone();
            }
            x -= w;
        }
        items[items.len() - 1].clone()
    }

    /// `k` distinct indices below `n`, in random order (Python's
    /// `sample(range(n), k)`).
    pub fn sample(&mut self, n: usize, k: usize) -> Vec<usize> {
        let mut pool: Vec<usize> = (0..n).collect();
        let k = k.min(n);
        for i in 0..k {
            let j = i + self.randrange(n - i);
            pool.swap(i, j);
        }
        pool.truncate(k);
        pool
    }

    /// Shuffles in place.
    pub fn shuffle<T>(&mut self, items: &mut [T]) {
        for i in (1..items.len()).rev() {
            let j = self.randrange(i + 1);
            items.swap(i, j);
        }
    }

    /// True with probability `p`.
    pub fn chance(&mut self, p: f64) -> bool {
        self.random() < p
    }
}

#[cfg(test)]
mod tests {
    use super::Rng;

    #[test]
    fn deterministic() {
        let (mut a, mut b) = (Rng::new(7), Rng::new(7));
        assert!((0..100).all(|_| a.next_u64() == b.next_u64()));
        assert_ne!(Rng::new(1).next_u64(), Rng::new(2).next_u64());
    }

    #[test]
    fn ranges() {
        let mut r = Rng::new(3);
        let mut seen = [false; 5];
        for _ in 0..1000 {
            let x = r.randint(-2, 2);
            assert!((-2..=2).contains(&x));
            seen[(x + 2) as usize] = true;
            let f = r.random();
            assert!((0.0..1.0).contains(&f));
            assert!(r.randrange(3) < 3);
        }
        assert!(seen.iter().all(|&s| s));
        r.randint(i64::MIN, i64::MAX);
        let s = r.sample(10, 4);
        assert_eq!(s.len(), 4);
        let mut d = s.clone();
        d.sort();
        d.dedup();
        assert_eq!(d.len(), 4);
        assert_eq!(r.weighted(&[1, 2], &[0.0, 1.0]), 2);
    }
}
