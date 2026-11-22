//! Lanes of 16-bit cells with the saturating operations of vsearch's vector
//! intrinsics, for the aligner. There are eight in plain Rust, SSE2, or
//! NEON, and sixteen in AVX2. [`crate::simd::level`] chooses which.

/// The most lanes of any implementation.
pub const MAX_LANES: usize = 16;

/// Lanes of 16-bit cells.
///
/// The AVX2 implementation's methods use AVX2 instructions. Call them only
/// from code that runs after [`crate::simd::level`] chose AVX2, inlined
/// into a function compiled with the feature, as `align16` does.
pub trait Lanes: Copy {
    /// The lanes.
    const N: usize;
    /// The two bits a lane of a comparison, lane `c` at bits `2c` and
    /// `2c + 1`.
    type Bits: Copy + Default;

    /// the bits of `bits` as a word.
    fn bits(bits: Self::Bits) -> u32;
    /// the first `N` cells of `cells`.
    fn load(cells: &[i16; MAX_LANES]) -> Self;
    /// the cells, in the first `N` places.
    fn store(self) -> [i16; MAX_LANES];
    fn dup(x: i16) -> Self;
    /// saturating addition and subtraction.
    fn add(self, b: Self) -> Self;
    fn sub(self, b: Self) -> Self;
    /// saturating subtraction of unsigned cells.
    fn sub_unsigned(self, b: Self) -> Self;
    fn max(self, b: Self) -> Self;
    fn min(self, b: Self) -> Self;
    fn and(self, b: Self) -> Self;
    fn xor(self, b: Self) -> Self;
    /// both bits of each lane where `self` is above `b`.
    fn mask_gt(self, b: Self) -> Self::Bits;

    #[inline(always)]
    fn zero() -> Self {
        Self::dup(0)
    }

    /// `self` with lane `c` set to `x`.
    #[inline(always)]
    fn with(self, c: usize, x: i16) -> Self {
        let mut cells = self.store();
        cells[c] = x;
        Self::load(&cells)
    }

    /// all ones in lane `c`, zero in the others.
    #[inline(always)]
    fn lane(c: usize) -> Self {
        Self::zero().with(c, -1)
    }
}

/// Eight lanes in plain Rust, for any machine.
#[derive(Clone, Copy)]
pub struct Scalar([i16; 8]);

impl Lanes for Scalar {
    const N: usize = 8;
    type Bits = u16;

    #[inline(always)]
    fn bits(bits: u16) -> u32 {
        u32::from(bits)
    }
    #[inline(always)]
    fn load(cells: &[i16; MAX_LANES]) -> Self {
        Scalar(std::array::from_fn(|i| cells[i]))
    }
    #[inline(always)]
    fn store(self) -> [i16; MAX_LANES] {
        std::array::from_fn(|i| if i < 8 { self.0[i] } else { 0 })
    }
    #[inline(always)]
    fn dup(x: i16) -> Self {
        Scalar([x; 8])
    }
    #[inline(always)]
    fn add(self, b: Self) -> Self {
        Scalar(std::array::from_fn(|i| self.0[i].saturating_add(b.0[i])))
    }
    #[inline(always)]
    fn sub(self, b: Self) -> Self {
        Scalar(std::array::from_fn(|i| self.0[i].saturating_sub(b.0[i])))
    }
    #[inline(always)]
    fn sub_unsigned(self, b: Self) -> Self {
        Scalar(std::array::from_fn(|i| {
            (self.0[i] as u16).saturating_sub(b.0[i] as u16) as i16
        }))
    }
    #[inline(always)]
    fn max(self, b: Self) -> Self {
        Scalar(std::array::from_fn(|i| self.0[i].max(b.0[i])))
    }
    #[inline(always)]
    fn min(self, b: Self) -> Self {
        Scalar(std::array::from_fn(|i| self.0[i].min(b.0[i])))
    }
    #[inline(always)]
    fn and(self, b: Self) -> Self {
        Scalar(std::array::from_fn(|i| self.0[i] & b.0[i]))
    }
    #[inline(always)]
    fn xor(self, b: Self) -> Self {
        Scalar(std::array::from_fn(|i| self.0[i] ^ b.0[i]))
    }
    #[inline(always)]
    fn mask_gt(self, b: Self) -> u16 {
        let mut mask = 0u16;
        for i in 0..8 {
            if self.0[i] > b.0[i] {
                mask |= 3 << (2 * i);
            }
        }
        mask
    }
}

#[cfg(target_arch = "x86_64")]
pub use x86::{Avx2, Sse2};

#[cfg(target_arch = "x86_64")]
mod x86 {
    use super::{Lanes, MAX_LANES};
    use std::arch::x86_64::*;

    /// Eight lanes of SSE2, part of every x86-64 processor.
    #[derive(Clone, Copy)]
    pub struct Sse2(__m128i);

    // SAFETY (for every unsafe block below): SSE2 is part of the x86-64
    // baseline, and the loads and stores read or write 16 bytes of a
    // 32-byte array.
    impl Lanes for Sse2 {
        const N: usize = 8;
        type Bits = u16;

        #[inline(always)]
        fn bits(bits: u16) -> u32 {
            u32::from(bits)
        }
        #[inline(always)]
        fn load(cells: &[i16; MAX_LANES]) -> Self {
            Sse2(unsafe { _mm_loadu_si128(cells.as_ptr().cast()) })
        }
        #[inline(always)]
        fn store(self) -> [i16; MAX_LANES] {
            let mut cells = [0i16; MAX_LANES];
            unsafe { _mm_storeu_si128(cells.as_mut_ptr().cast(), self.0) };
            cells
        }
        #[inline(always)]
        fn dup(x: i16) -> Self {
            Sse2(unsafe { _mm_set1_epi16(x) })
        }
        #[inline(always)]
        fn add(self, b: Self) -> Self {
            Sse2(unsafe { _mm_adds_epi16(self.0, b.0) })
        }
        #[inline(always)]
        fn sub(self, b: Self) -> Self {
            Sse2(unsafe { _mm_subs_epi16(self.0, b.0) })
        }
        #[inline(always)]
        fn sub_unsigned(self, b: Self) -> Self {
            Sse2(unsafe { _mm_subs_epu16(self.0, b.0) })
        }
        #[inline(always)]
        fn max(self, b: Self) -> Self {
            Sse2(unsafe { _mm_max_epi16(self.0, b.0) })
        }
        #[inline(always)]
        fn min(self, b: Self) -> Self {
            Sse2(unsafe { _mm_min_epi16(self.0, b.0) })
        }
        #[inline(always)]
        fn and(self, b: Self) -> Self {
            Sse2(unsafe { _mm_and_si128(self.0, b.0) })
        }
        #[inline(always)]
        fn xor(self, b: Self) -> Self {
            Sse2(unsafe { _mm_xor_si128(self.0, b.0) })
        }
        #[inline(always)]
        fn mask_gt(self, b: Self) -> u16 {
            unsafe { _mm_movemask_epi8(_mm_cmpgt_epi16(self.0, b.0)) as u16 }
        }
    }

    /// Sixteen lanes of AVX2.
    #[derive(Clone, Copy)]
    pub struct Avx2(__m256i);

    // SAFETY (for every unsafe block below): these run only once
    // `simd::level` found AVX2 on the processor, inlined into a function
    // compiled with it; the loads and stores read or write 32 bytes of a
    // 32-byte array.
    impl Lanes for Avx2 {
        const N: usize = 16;
        type Bits = u32;

        #[inline(always)]
        fn bits(bits: u32) -> u32 {
            bits
        }
        #[inline(always)]
        fn load(cells: &[i16; MAX_LANES]) -> Self {
            Avx2(unsafe { _mm256_loadu_si256(cells.as_ptr().cast()) })
        }
        #[inline(always)]
        fn store(self) -> [i16; MAX_LANES] {
            let mut cells = [0i16; MAX_LANES];
            unsafe { _mm256_storeu_si256(cells.as_mut_ptr().cast(), self.0) };
            cells
        }
        #[inline(always)]
        fn dup(x: i16) -> Self {
            Avx2(unsafe { _mm256_set1_epi16(x) })
        }
        #[inline(always)]
        fn add(self, b: Self) -> Self {
            Avx2(unsafe { _mm256_adds_epi16(self.0, b.0) })
        }
        #[inline(always)]
        fn sub(self, b: Self) -> Self {
            Avx2(unsafe { _mm256_subs_epi16(self.0, b.0) })
        }
        #[inline(always)]
        fn sub_unsigned(self, b: Self) -> Self {
            Avx2(unsafe { _mm256_subs_epu16(self.0, b.0) })
        }
        #[inline(always)]
        fn max(self, b: Self) -> Self {
            Avx2(unsafe { _mm256_max_epi16(self.0, b.0) })
        }
        #[inline(always)]
        fn min(self, b: Self) -> Self {
            Avx2(unsafe { _mm256_min_epi16(self.0, b.0) })
        }
        #[inline(always)]
        fn and(self, b: Self) -> Self {
            Avx2(unsafe { _mm256_and_si256(self.0, b.0) })
        }
        #[inline(always)]
        fn xor(self, b: Self) -> Self {
            Avx2(unsafe { _mm256_xor_si256(self.0, b.0) })
        }
        #[inline(always)]
        fn mask_gt(self, b: Self) -> u32 {
            unsafe { _mm256_movemask_epi8(_mm256_cmpgt_epi16(self.0, b.0)) as u32 }
        }
    }
}

#[cfg(target_arch = "aarch64")]
pub use arm::Neon;

#[cfg(target_arch = "aarch64")]
mod arm {
    use super::{Lanes, MAX_LANES};
    use std::arch::aarch64::*;

    /// Eight lanes of NEON, part of every aarch64 processor.
    #[derive(Clone, Copy)]
    pub struct Neon(int16x8_t);

    /// The two bits of each lane, for `mask_gt`.
    const MASK: [u16; 8] = [
        0x0003, 0x000c, 0x0030, 0x00c0, 0x0300, 0x0c00, 0x3000, 0xc000,
    ];

    // SAFETY (for every unsafe block below): NEON is part of the aarch64
    // baseline, and the loads and stores read or write 16 bytes of a
    // 32-byte array.
    impl Lanes for Neon {
        const N: usize = 8;
        type Bits = u16;

        #[inline(always)]
        fn bits(bits: u16) -> u32 {
            u32::from(bits)
        }
        #[inline(always)]
        fn load(cells: &[i16; MAX_LANES]) -> Self {
            Neon(unsafe { vld1q_s16(cells.as_ptr()) })
        }
        #[inline(always)]
        fn store(self) -> [i16; MAX_LANES] {
            let mut cells = [0i16; MAX_LANES];
            unsafe { vst1q_s16(cells.as_mut_ptr(), self.0) };
            cells
        }
        #[inline(always)]
        fn dup(x: i16) -> Self {
            Neon(unsafe { vdupq_n_s16(x) })
        }
        #[inline(always)]
        fn add(self, b: Self) -> Self {
            Neon(unsafe { vqaddq_s16(self.0, b.0) })
        }
        #[inline(always)]
        fn sub(self, b: Self) -> Self {
            Neon(unsafe { vqsubq_s16(self.0, b.0) })
        }
        #[inline(always)]
        fn sub_unsigned(self, b: Self) -> Self {
            Neon(unsafe {
                vreinterpretq_s16_u16(vqsubq_u16(
                    vreinterpretq_u16_s16(self.0),
                    vreinterpretq_u16_s16(b.0),
                ))
            })
        }
        #[inline(always)]
        fn max(self, b: Self) -> Self {
            Neon(unsafe { vmaxq_s16(self.0, b.0) })
        }
        #[inline(always)]
        fn min(self, b: Self) -> Self {
            Neon(unsafe { vminq_s16(self.0, b.0) })
        }
        #[inline(always)]
        fn and(self, b: Self) -> Self {
            Neon(unsafe { vandq_s16(self.0, b.0) })
        }
        #[inline(always)]
        fn xor(self, b: Self) -> Self {
            Neon(unsafe { veorq_s16(self.0, b.0) })
        }
        #[inline(always)]
        fn mask_gt(self, b: Self) -> u16 {
            unsafe { vaddvq_u16(vandq_u16(vcgtq_s16(self.0, b.0), vld1q_u16(MASK.as_ptr()))) }
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    /// the operations of `L` on fixed cells, as vsearch's intrinsics give
    /// them, in the first eight lanes and, for sixteen, the next eight.
    fn saturates_and_masks<L: Lanes>() {
        let mut a = [0i16; MAX_LANES];
        let mut b = [0i16; MAX_LANES];
        a[..8].copy_from_slice(&[i16::MAX, i16::MIN, 5, -5, 0, 1, 2, 3]);
        b[..8].copy_from_slice(&[1, 1, 7, -7, 0, 0, 3, 2]);
        if L::N == 16 {
            a[8..].copy_from_slice(&[0, 0, 0, 0, i16::MAX, 9, 9, 9]);
            b[8..].copy_from_slice(&[1, 1, 1, 1, 1, 9, 8, 10]);
        }
        let (va, vb) = (L::load(&a), L::load(&b));
        assert_eq!(
            va.add(vb).store()[..8],
            [i16::MAX, -32767, 12, -12, 0, 1, 5, 5]
        );
        assert_eq!(
            va.sub(vb).store()[..8],
            [32766, i16::MIN, -2, 2, 0, 1, -1, 1]
        );
        assert_eq!(vb.sub_unsigned(L::dup(-1)).store(), [0; MAX_LANES]);
        assert_eq!(va.max(vb).store()[..4], [i16::MAX, 1, 7, -5]);
        assert_eq!(va.min(vb).store()[..4], [1, i16::MIN, 5, -7]);
        // lanes 0, 3, 5 and 7 are above
        let mut expected = 0xccc3u32;
        if L::N == 16 {
            // and lanes 12 and 14
            expected |= 0x3300 << 16;
            assert_eq!(va.add(vb).store()[12], i16::MAX);
        }
        assert_eq!(L::bits(va.mask_gt(vb)), expected);
        assert_eq!(L::lane(3).store()[3], -1);
        assert_eq!(L::lane(3).store().iter().filter(|&&x| x != 0).count(), 1);
        assert_eq!(L::zero().with(L::N - 1, 7).store()[L::N - 1], 7);
        assert_eq!(va.and(L::dup(0)).store(), [0; MAX_LANES]);
        assert_eq!(va.xor(va).store(), [0; MAX_LANES]);
    }

    #[test]
    fn scalar_saturates_and_masks_as_vsearch() {
        saturates_and_masks::<Scalar>();
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn sse2_saturates_and_masks_as_vsearch() {
        saturates_and_masks::<Sse2>();
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn avx2_saturates_and_masks_as_vsearch() {
        if std::arch::is_x86_feature_detected!("avx2") {
            // SAFETY: the processor has AVX2
            #[target_feature(enable = "avx2")]
            fn run() {
                saturates_and_masks::<Avx2>();
            }
            unsafe { run() };
        }
    }

    #[cfg(target_arch = "aarch64")]
    #[test]
    fn neon_saturates_and_masks_as_vsearch() {
        saturates_and_masks::<Neon>();
    }
}
