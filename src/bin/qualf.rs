extern crate jemallocator;
extern crate bio;
extern crate clap;
extern crate rsearch;

#[global_allocator]
static ALLOC: jemallocator::Jemalloc = jemallocator::Jemalloc;

use std::fs::File;
use std::io;
use std::io::{Read, Write};
use std::str::FromStr;

use bio::io::fastq::{Reader, Writer};
use clap::{App, Arg};

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
    let matches = App::new("QualF")
        .about("Filter FASTQ sequences by quality scores.")
        .arg(
            Arg::with_name("in")
                .short("i")
                .long("in")
                .value_name("INPUT")
                .help("Input FASTQ file, default to stdin")
                .takes_value(true),
        ).arg(
            Arg::with_name("out")
                .short("o")
                .long("out")
                .value_name("OUTPUT")
                .help("Output FASTQ file, default to stdout")
                .takes_value(true),
        ).arg(
            Arg::with_name("max")
                .long("max")
                .value_name("MAX_LEN")
                .help("Maximum value for length filter, disabled by default")
                .takes_value(true),
        ).arg(
            Arg::with_name("min")
                .long("min")
                .value_name("MIN_LEN")
                .help("Minimum value for length filter, disabled by default")
                .takes_value(true),
        ).arg(
            Arg::with_name("prob")
                .long("prob")
                .value_name("PROB")
                .help("Probability of incorrect base call for quality filter, disabled by default")
                .takes_value(true),
        ).arg(
            Arg::with_name("bad")
                .long("bad")
                .value_name("BAD")
                .help("Number of bad base call allowed for quality filter, disabled by default")
                .takes_value(true),
        ).get_matches();

    let fin: Box<Read> = match matches.value_of("in") {
        Some(f) => Box::new(File::open(f).unwrap()),
        None => Box::new(io::stdin()),
    };

    let reader = Reader::new(fin);

    let fout: Box<Write> = match matches.value_of("out") {
        Some(f) => Box::new(File::create(f).unwrap()),
        None => Box::new(io::stdout()),
    };

    let mut writer = Writer::new(fout);

    let min_len = parse_option(matches.value_of("min"));

    let max_len = parse_option(matches.value_of("max"));

    let prob = parse_option(matches.value_of("prob"));

    let bad_call = parse_option(matches.value_of("bad_call"));

    for record in reader.records() {
        // Since FASTQ uses '@' for both the description header and a quality symbol, there is no
        // reliable way to recover from a malformed record.
        // Panics if a FASTQ record cannot be read.
        let r = record.unwrap();

        if r.check().is_err() {
            continue;
        }

        let l = r.seq().len();

        if min_len.map_or(false, |m| l < m) || max_len.map_or(false, |m| l > m) {
            continue;
        }

        // Calculates probability of an mistake using Sanger encoding.
        // Sanger encoding starts from '!', or 33u8.
        let phred: f64 = r
            .qual()
            .iter()
            .map(|x| 10f64.powf(-(f64::from(*x - 33u8) / 10f64)))
            .sum();

        if prob.map_or(false, |p| phred > p) {
            continue;
        }

        // Counts the number of bad base call 'N'.
        let bad = r
            .seq()
            .iter()
            .filter(|x| x.to_ascii_uppercase() == b'N')
            .count();

        if bad_call.map_or(false, |b| bad > b) {
            continue;
        }

        writer.write_record(&r).unwrap();
    }
}
