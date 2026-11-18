use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{Read, Write};
use std::str::FromStr;

use bio::io::fastq::{Reader, Writer};
use clap::{Arg, Command};

/// parses an `Option<&str>` to `Option<T>`.
///
/// # Example
/// ```rust
/// let foo = Some("123");
///
/// let bar = parse_option<u64>(foo);
/// assert_eq!(bar, Some(123));
/// ```
fn parse_option<T>(o: Option<&str>) -> Option<T>
where
    T: FromStr,
{
    o?.parse::<T>().ok()
}

fn main() {
    let matches =
        Command::new("QualF")
            .about("Filter FASTQ sequences by quality scores.")
            .arg(
                Arg::new("in")
                    .short('i')
                    .long("in")
                    .value_name("INPUT")
                    .help("Input FASTQ file, default to stdin"),
            )
            .arg(
                Arg::new("out")
                    .short('o')
                    .long("out")
                    .value_name("OUTPUT")
                    .help("Output FASTQ file, default to stdout"),
            )
            .arg(
                Arg::new("max")
                    .long("max")
                    .value_name("MAX_LEN")
                    .help("Maximum value for length filter, disabled by default"),
            )
            .arg(
                Arg::new("min")
                    .long("min")
                    .value_name("MIN_LEN")
                    .help("Minimum value for length filter, disabled by default"),
            )
            .arg(
                Arg::new("prob").long("prob").value_name("PROB").help(
                    "Probability of incorrect base call for quality filter, disabled by default",
                ),
            )
            .arg(
                Arg::new("bad").long("bad").value_name("BAD").help(
                    "Number of bad base call allowed for quality filter, disabled by default",
                ),
            )
            .get_matches();

    let fin: Box<dyn Read> = match matches.get_one::<String>("in").map(String::as_str) {
        Some(f) => Box::new(File::open(f).unwrap()),
        None => Box::new(io::stdin()),
    };

    let reader = Reader::new(fin);

    let fout: Box<dyn Write> = match matches.get_one::<String>("out").map(String::as_str) {
        Some(f) => Box::new(File::create(f).unwrap()),
        None => Box::new(io::stdout()),
    };

    let mut writer = Writer::new(fout);

    let min_len = parse_option(matches.get_one::<String>("min").map(String::as_str));

    let max_len = parse_option(matches.get_one::<String>("max").map(String::as_str));

    let prob = parse_option(matches.get_one::<String>("prob").map(String::as_str));

    let bad_call = parse_option(matches.get_one::<String>("bad").map(String::as_str));

    for record in reader.records() {
        // Since FASTQ uses '@' for both the description header and a quality symbol, there is no
        // reliable way to recover from a malformed record.
        // Panics if a FASTQ record cannot be read.
        let r = record.unwrap();

        if r.check().is_err() {
            continue;
        }

        let l = r.seq().len();

        if min_len.is_some_and(|m| l < m) || max_len.is_some_and(|m| l > m) {
            continue;
        }

        // Calculates probability of an mistake using Sanger encoding.
        // Sanger encoding starts from '!', or 33u8.
        let phred: f64 = r
            .qual()
            .iter()
            .map(|x| 10f64.powf(-(f64::from(*x - 33u8) / 10f64)))
            .sum();

        if prob.is_some_and(|p| phred > p) {
            continue;
        }

        // Counts the number of bad base call 'N'.
        let bad = r
            .seq()
            .iter()
            .filter(|x| x.eq_ignore_ascii_case(&b'N'))
            .count();

        if bad_call.is_some_and(|b| bad > b) {
            continue;
        }

        writer.write_record(&r).unwrap();
    }
}
