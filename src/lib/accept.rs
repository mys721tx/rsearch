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

use crate::alignment::{Stats, parse_cigar};
use crate::scoring::{FOUR_BIT, INFINITE, Scoring};

/// What a hit is, after vsearch's criteria on its alignment.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Verdict {
    Accepted,
    /// Above `--weak_id` but below `--id`.
    Weak,
    Rejected,
}

/// vsearch's acceptance criteria of hits, with its defaults.
#[derive(Clone, Debug)]
pub struct Criteria {
    pub id: f64,
    pub weak_id: f64,
    pub iddef: usize,
    pub max_diffs: i64,
    pub max_subs: i64,
    pub max_gaps: i64,
    pub min_cols: i64,
    pub left_just: bool,
    pub right_just: bool,
    pub query_cov: f64,
    pub target_cov: f64,
    pub max_id: f64,
    pub mid: f64,
    pub id_prefix: usize,
    pub id_suffix: usize,
    pub max_qt: f64,
    pub min_qt: f64,
    pub min_sl: f64,
    pub max_sl: f64,
    pub sizes: SizeFilters,
    pub self_label: bool,
    pub self_id: bool,
}

impl Default for Criteria {
    fn default() -> Criteria {
        Criteria {
            id: -1.0,
            weak_id: 10.0,
            iddef: 2,
            max_diffs: i64::from(i32::MAX),
            max_subs: i64::from(i32::MAX),
            max_gaps: i64::from(i32::MAX),
            min_cols: 0,
            left_just: false,
            right_just: false,
            query_cov: 0.0,
            target_cov: 0.0,
            max_id: 1.0,
            mid: 0.0,
            id_prefix: 0,
            id_suffix: 0,
            max_qt: f64::MAX,
            min_qt: 0.0,
            min_sl: 0.0,
            max_sl: f64::MAX,
            sizes: SizeFilters::default(),
            self_label: false,
            self_id: false,
        }
    }
}

/// whether two sequences hold the same 4-bit codes, as vsearch's `seqcmp`.
fn same_codes(a: &[u8], b: &[u8]) -> bool {
    a.iter()
        .zip(b)
        .all(|(&x, &y)| FOUR_BIT[usize::from(x)] == FOUR_BIT[usize::from(y)])
}

impl Criteria {
    /// vsearch's `search_acceptable_unaligned`: whether a query may hit a
    /// target before they are aligned.
    #[allow(clippy::too_many_arguments)]
    pub fn unaligned(
        &self,
        query_label: &str,
        query: &[u8],
        query_size: i64,
        target_label: &str,
        target: &[u8],
        target_size: i64,
    ) -> bool {
        let qlen = query.len() as i64;
        let dlen = target.len() as i64;
        let (q, d) = (qlen as f64, dlen as f64);
        let shorter = if qlen < dlen {
            q >= self.min_sl * d
        } else {
            d >= self.min_sl * q
        };
        let longer = if qlen < dlen {
            q <= self.max_sl * d
        } else {
            d <= self.max_sl * q
        };
        self.sizes.accepts(query_size, target_size)
            && q >= self.min_qt * d
            && q <= self.max_qt * d
            && shorter
            && longer
            && query.len() >= self.id_prefix
            && target.len() >= self.id_prefix
            && same_codes(&query[..self.id_prefix], &target[..self.id_prefix])
            && query.len() >= self.id_suffix
            && target.len() >= self.id_suffix
            && same_codes(
                &query[query.len() - self.id_suffix..],
                &target[target.len() - self.id_suffix..],
            )
            && !(self.self_label && query_label == target_label)
            && !(self.self_id && qlen == dlen && same_codes(query, target))
    }

    /// vsearch's `search_acceptable_aligned`: what a hit is, given its
    /// alignment, the CIGAR, and the gap penalties.
    pub fn aligned(
        &self,
        stats: &Stats,
        cigar: &str,
        qlen: usize,
        tlen: usize,
        scoring: &Scoring,
    ) -> Verdict {
        let id = stats.id[self.iddef];
        let pairs = stats.matches + stats.mismatches;
        let mid = if pairs > 0 {
            100.0 * stats.matches as f64 / pairs as f64
        } else {
            0.0
        };
        let passes = id >= 100.0 * self.weak_id
            && stats.mismatches <= self.max_subs
            && stats.internal_gaps <= self.max_gaps
            && !uses_forbidden_gap(cigar, scoring)
            && stats.internal_length >= self.min_cols
            && (!self.left_just || stats.trim_q_left + stats.trim_t_left == 0)
            && (!self.right_just || stats.trim_q_right + stats.trim_t_right == 0)
            && (pairs as f64) >= self.query_cov * (qlen as i32) as f64
            && (pairs as f64) >= self.target_cov * tlen as f64
            && id <= 100.0 * self.max_id
            && mid >= self.mid
            && stats.mismatches + stats.internal_indels <= self.max_diffs;
        if !passes {
            Verdict::Rejected
        } else if id >= 100.0 * self.id {
            Verdict::Accepted
        } else {
            Verdict::Weak
        }
    }
}

/// whether an alignment opens a gap of infinite opening penalty, or
/// extends one of infinite extension penalty, as vsearch's
/// `alignment_uses_forbidden_gap`.
fn uses_forbidden_gap(cigar: &str, scoring: &Scoring) -> bool {
    let infinite = scoring
        .open_query
        .iter()
        .chain(&scoring.open_target)
        .chain(&scoring.extend_query)
        .chain(&scoring.extend_target)
        .any(|&p| p == INFINITE);
    if !infinite {
        return false;
    }
    let Ok(ops) = parse_cigar(cigar) else {
        return false;
    };
    let last = ops.len().saturating_sub(1);
    for (k, &(op, run)) in ops.iter().enumerate() {
        if op == b'M' {
            continue;
        }
        let end = if k == 0 {
            0
        } else if k == last {
            2
        } else {
            1
        };
        let (open, extend) = if op == b'I' {
            (scoring.open_query[end], scoring.extend_query[end])
        } else {
            (scoring.open_target[end], scoring.extend_target[end])
        };
        if open == INFINITE || (extend == INFINITE && run > 1) {
            return true;
        }
    }
    false
}
