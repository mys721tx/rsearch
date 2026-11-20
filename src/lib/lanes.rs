//! Eight 16-bit cells with the saturating operations of vsearch's vector
//! intrinsics: NEON on aarch64, SSE2 on x86-64, and an array elsewhere.

#[cfg(target_arch = "aarch64")]
mod imp {
    use std::arch::aarch64::*;

    pub type V = int16x8_t;

    #[inline(always)]
    pub fn from_array(a: [i16; 8]) -> V {
        unsafe { vld1q_s16(a.as_ptr()) }
    }
    #[inline(always)]
    pub fn to_array(v: V) -> [i16; 8] {
        let mut a = [0i16; 8];
        unsafe { vst1q_s16(a.as_mut_ptr(), v) };
        a
    }
    #[inline(always)]
    pub fn add(a: V, b: V) -> V {
        unsafe { vqaddq_s16(a, b) }
    }
    #[inline(always)]
    pub fn sub(a: V, b: V) -> V {
        unsafe { vqsubq_s16(a, b) }
    }
    #[inline(always)]
    pub fn sub_unsigned(a: V, b: V) -> V {
        unsafe {
            vreinterpretq_s16_u16(vqsubq_u16(
                vreinterpretq_u16_s16(a),
                vreinterpretq_u16_s16(b),
            ))
        }
    }
    #[inline(always)]
    pub fn max(a: V, b: V) -> V {
        unsafe { vmaxq_s16(a, b) }
    }
    #[inline(always)]
    pub fn min(a: V, b: V) -> V {
        unsafe { vminq_s16(a, b) }
    }
    #[inline(always)]
    pub fn and(a: V, b: V) -> V {
        unsafe { vandq_s16(a, b) }
    }
    #[inline(always)]
    pub fn xor(a: V, b: V) -> V {
        unsafe { veorq_s16(a, b) }
    }
    #[inline(always)]
    pub fn dup(x: i16) -> V {
        unsafe { vdupq_n_s16(x) }
    }
    #[inline(always)]
    pub fn shift_left(a: V) -> V {
        unsafe { vextq_s16::<7>(vdupq_n_s16(0), a) }
    }
    #[inline(always)]
    pub fn mask_gt(a: V, b: V) -> u16 {
        const MASK: [u16; 8] = [
            0x0003, 0x000c, 0x0030, 0x00c0, 0x0300, 0x0c00, 0x3000, 0xc000,
        ];
        unsafe { vaddvq_u16(vandq_u16(vcgtq_s16(a, b), vld1q_u16(MASK.as_ptr()))) }
    }
}

#[cfg(target_arch = "x86_64")]
mod imp {
    use std::arch::x86_64::*;

    pub type V = __m128i;

    #[inline(always)]
    pub fn from_array(a: [i16; 8]) -> V {
        unsafe { _mm_loadu_si128(a.as_ptr().cast()) }
    }
    #[inline(always)]
    pub fn to_array(v: V) -> [i16; 8] {
        let mut a = [0i16; 8];
        unsafe { _mm_storeu_si128(a.as_mut_ptr().cast(), v) };
        a
    }
    #[inline(always)]
    pub fn add(a: V, b: V) -> V {
        unsafe { _mm_adds_epi16(a, b) }
    }
    #[inline(always)]
    pub fn sub(a: V, b: V) -> V {
        unsafe { _mm_subs_epi16(a, b) }
    }
    #[inline(always)]
    pub fn sub_unsigned(a: V, b: V) -> V {
        unsafe { _mm_subs_epu16(a, b) }
    }
    #[inline(always)]
    pub fn max(a: V, b: V) -> V {
        unsafe { _mm_max_epi16(a, b) }
    }
    #[inline(always)]
    pub fn min(a: V, b: V) -> V {
        unsafe { _mm_min_epi16(a, b) }
    }
    #[inline(always)]
    pub fn and(a: V, b: V) -> V {
        unsafe { _mm_and_si128(a, b) }
    }
    #[inline(always)]
    pub fn xor(a: V, b: V) -> V {
        unsafe { _mm_xor_si128(a, b) }
    }
    #[inline(always)]
    pub fn dup(x: i16) -> V {
        unsafe { _mm_set1_epi16(x) }
    }
    #[inline(always)]
    pub fn shift_left(a: V) -> V {
        unsafe { _mm_slli_si128::<2>(a) }
    }
    #[inline(always)]
    pub fn mask_gt(a: V, b: V) -> u16 {
        unsafe { _mm_movemask_epi8(_mm_cmpgt_epi16(a, b)) as u16 }
    }
}

#[cfg(not(any(target_arch = "aarch64", target_arch = "x86_64")))]
mod imp {
    pub type V = [i16; 8];

    pub fn from_array(a: [i16; 8]) -> V {
        a
    }
    pub fn to_array(v: V) -> [i16; 8] {
        v
    }
    pub fn add(a: V, b: V) -> V {
        std::array::from_fn(|i| a[i].saturating_add(b[i]))
    }
    pub fn sub(a: V, b: V) -> V {
        std::array::from_fn(|i| a[i].saturating_sub(b[i]))
    }
    pub fn sub_unsigned(a: V, b: V) -> V {
        std::array::from_fn(|i| (a[i] as u16).saturating_sub(b[i] as u16) as i16)
    }
    pub fn max(a: V, b: V) -> V {
        std::array::from_fn(|i| a[i].max(b[i]))
    }
    pub fn min(a: V, b: V) -> V {
        std::array::from_fn(|i| a[i].min(b[i]))
    }
    pub fn and(a: V, b: V) -> V {
        std::array::from_fn(|i| a[i] & b[i])
    }
    pub fn xor(a: V, b: V) -> V {
        std::array::from_fn(|i| a[i] ^ b[i])
    }
    pub fn dup(x: i16) -> V {
        [x; 8]
    }
    pub fn shift_left(a: V) -> V {
        std::array::from_fn(|i| if i == 0 { 0 } else { a[i - 1] })
    }
    pub fn mask_gt(a: V, b: V) -> u16 {
        let mut mask = 0u16;
        for i in 0..8 {
            if a[i] > b[i] {
                mask |= 3 << (2 * i);
            }
        }
        mask
    }
}

pub use imp::*;

/// all lanes zero.
#[inline(always)]
pub fn zero() -> V {
    dup(0)
}

/// sets lane `c` of `v` to `x`.
#[inline]
pub fn set(v: &mut V, c: usize, x: i16) {
    let mut a = to_array(*v);
    a[c] = x;
    *v = from_array(a);
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn saturates_and_masks_as_vsearch() {
        let a = from_array([i16::MAX, i16::MIN, 5, -5, 0, 1, 2, 3]);
        let b = from_array([1, 1, 7, -7, 0, 0, 3, 2]);
        assert_eq!(to_array(add(a, b)), [i16::MAX, -32767, 12, -12, 0, 1, 5, 5]);
        assert_eq!(to_array(sub(a, b)), [32766, i16::MIN, -2, 2, 0, 1, -1, 1]);
        assert_eq!(to_array(sub_unsigned(b, dup(-1))), [0; 8]);
        assert_eq!(to_array(shift_left(a))[..3], [0, i16::MAX, i16::MIN]);
        // lanes 0, 3, 5 and 7 are above
        assert_eq!(mask_gt(a, b), 0xccc3);
    }
}
