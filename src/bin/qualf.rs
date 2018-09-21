extern crate bio;
extern crate clap;
extern crate rsearch;

use std::fs::File;
use std::io;
use std::io::{Read, Write};

use bio::io::fastq::{Reader, Writer};
use clap::{App, Arg};

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

    for record in reader.records() {
        let r = record.unwrap();

        if r.check().is_err() {
            continue;
        }

        let l = r.seq().len();

        if l < 225 || l > 275{
            continue;
        }

        let phred: f64 = r
            .qual()
            .iter()
            .map(|x| 10f64.powf(-(f64::from(*x - 33u8) / 10f64)))
            .sum();

        if phred > 0.5f64 {
            continue;
        }

        let bad = r
            .seq()
            .iter()
            .filter(|x| x.to_ascii_uppercase() == b'N')
            .count();

        if bad > 0 {
            continue;
        }

        writer.write_record(&r).unwrap();
    }
}
