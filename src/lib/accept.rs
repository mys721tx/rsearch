//! Hit acceptance: vsearch's criteria on a query and a target before they
//! are aligned (`search_acceptable_unaligned`).

/// How `value` compares with `ratio * reference`, exactly, as vsearch's
/// `abundance_ratio_cmp`: -1, 0 or 1.
///
/// # Example
/// ```rust
/// # use rsearch::accept::abundance_ratio_cmp;
/// assert_eq!(abundance_ratio_cmp(3, 0.5, 6), 0);
/// assert_eq!(abundance_ratio_cmp(3, 0.5, 9), -1);
/// assert_eq!(abundance_ratio_cmp(1, f64::INFINITY, 2), -1);
/// assert_eq!(abundance_ratio_cmp(1, 0.0, 2), 1);
/// ```
pub fn abundance_ratio_cmp(value: i64, ratio: f64, reference: i64) -> i32 {
    if reference <= 0 || ratio <= 0.0 {
        return i32::from(value > 0);
    }
    if !ratio.is_finite() {
        return -1;
    }
    const EXACT: i64 = 1 << 53;
    if value < EXACT && reference < EXACT {
        let product = ratio * reference as f64;
        let value = value as f64;
        return if value < product {
            -1
        } else if value > product {
            1
        } else {
            0
        };
    }
    // ratio = mantissa * 2^exponent, exactly; compare value * 2^-exponent
    // with mantissa * reference in 128 bits
    let bits = ratio.to_bits();
    let exponent = ((bits >> 52) & 0x7ff) as i32;
    let (mantissa, exponent) = if exponent == 0 {
        (bits & ((1 << 52) - 1), -1074)
    } else {
        ((bits & ((1 << 52) - 1)) | (1 << 52), exponent - 1075)
    };
    let left = i128::from(value);
    let right = i128::from(mantissa) * i128::from(reference);
    let (left, right) = if exponent >= 0 {
        (
            left,
            right.checked_shl(exponent as u32).unwrap_or(i128::MAX),
        )
    } else {
        let shift = (-exponent) as u32;
        if shift >= 127 {
            (left, 0)
        } else {
            (left.checked_shl(shift).unwrap_or(i128::MAX), right)
        }
    };
    match left.cmp(&right) {
        std::cmp::Ordering::Less => -1,
        std::cmp::Ordering::Equal => 0,
        std::cmp::Ordering::Greater => 1,
    }
}

/// The abundance criteria of a hit, as vsearch's `--maxqsize`,
/// `--mintsize`, `--minsizeratio` and `--maxsizeratio`.
#[derive(Clone, Copy, Debug)]
pub struct SizeFilters {
    pub max_query_size: i64,
    pub min_target_size: i64,
    pub min_size_ratio: f64,
    pub max_size_ratio: f64,
}

impl Default for SizeFilters {
    fn default() -> SizeFilters {
        SizeFilters {
            max_query_size: i64::MAX,
            min_target_size: 0,
            min_size_ratio: 0.0,
            max_size_ratio: f64::MAX,
        }
    }
}

impl SizeFilters {
    /// whether a query of abundance `query` may hit a target of `target`.
    pub fn accepts(&self, query: i64, target: i64) -> bool {
        query <= self.max_query_size
            && target >= self.min_target_size
            && abundance_ratio_cmp(query, self.min_size_ratio, target) >= 0
            && abundance_ratio_cmp(query, self.max_size_ratio, target) <= 0
    }
}
