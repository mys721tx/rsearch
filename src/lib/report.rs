//! Number formats of vsearch's reports that Rust's formatting lacks.

/// formats `value` as C's `%.{precision}g`: `precision` significant
/// digits, in plain notation when the exponent is from -4 to below
/// `precision` and in exponent notation (at least two exponent digits)
/// otherwise, without trailing zeros or a trailing decimal point.
///
/// # Example
/// ```rust
/// # use rsearch::report::general;
/// assert_eq!(general(0.0001, 2), "0.0001");
/// assert_eq!(general(0.036700000000000004, 2), "0.037");
/// assert_eq!(general(1.0, 2), "1");
/// assert_eq!(general(6.309573444801929e-5, 2), "6.3e-05");
/// assert_eq!(general(123.0, 2), "1.2e+02");
/// assert_eq!(general(-1.0, 2), "-1");
/// ```
pub fn general(value: f64, precision: usize) -> String {
    let precision = precision.max(1);
    if value.is_nan() {
        return "nan".to_string();
    }
    let sign = if value.is_sign_negative() { "-" } else { "" };
    let magnitude = value.abs();
    if magnitude.is_infinite() {
        return format!("{sign}inf");
    }
    // the exponent after rounding to the digits kept decides the notation
    let scientific = format!("{:.*e}", precision - 1, magnitude);
    let (mantissa, exponent) = scientific.split_once('e').unwrap();
    let exponent: i32 = exponent.parse().unwrap();
    if exponent < -4 || exponent >= precision as i32 {
        let mantissa = strip_zeros(mantissa);
        let exponent_sign = if exponent < 0 { '-' } else { '+' };
        format!("{sign}{mantissa}e{exponent_sign}{:02}", exponent.abs())
    } else {
        let decimals = (precision as i32 - 1 - exponent) as usize;
        let plain = format!("{magnitude:.decimals$}");
        format!("{sign}{}", strip_zeros(&plain))
    }
}

/// removes the trailing zeros of a number's fraction, and the decimal
/// point if nothing is left after it.
fn strip_zeros(number: &str) -> &str {
    if number.contains('.') {
        number.trim_end_matches('0').trim_end_matches('.')
    } else {
        number
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_general() {
        // printf("%.2g") of each
        for (value, expected) in [
            (0.0, "0"),
            (0.5, "0.5"),
            (0.25, "0.25"),
            (0.033_366_666, "0.033"),
            (0.000_123_4, "0.00012"),
            (0.000_099_99, "0.0001"),
            (9.96, "10"),
            (99.6, "1e+02"),
            (1.0e-9, "1e-09"),
            (6.3e-100, "6.3e-100"),
        ] {
            assert_eq!(general(value, 2), expected, "{value}");
        }
    }
}
