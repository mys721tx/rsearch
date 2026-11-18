use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{BufReader, Read, Write};

use bio::io::{fasta, fastq};
use clap::{Arg, ArgAction, ArgMatches, Command, value_parser};

use rsearch::fastq::Reader;
use rsearch::filter::{Filter, abundance, annotate_ee};

/// parses a non-negative number of expected errors.
fn non_negative(s: &str) -> Result<f64, String> {
    match s.parse::<f64>() {
        Ok(x) if x >= 0.0 => Ok(x),
        _ => Err(format!("`{s}` is not a non-negative number")),
    }
}

/// parses a positive number of expected errors.
fn positive(s: &str) -> Result<f64, String> {
    match s.parse::<f64>() {
        Ok(x) if x > 0.0 => Ok(x),
        _ => Err(format!("`{s}` is not a positive number")),
    }
}

/// Writes records as FASTA or FASTQ.
enum Writer {
    Fasta(fasta::Writer<Box<dyn Write>>),
    Fastq(fastq::Writer<Box<dyn Write>>),
}

impl Writer {
    fn new(fout: Box<dyn Write>, as_fasta: bool) -> Writer {
        if as_fasta {
            Writer::Fasta(fasta::Writer::new(fout))
        } else {
            Writer::Fastq(fastq::Writer::new(fout))
        }
    }

    fn write(&mut self, id: &str, desc: Option<&str>, seq: &[u8], qual: &[u8]) {
        match self {
            Writer::Fasta(w) => w.write(id, desc, seq).unwrap(),
            Writer::Fastq(w) => w.write(id, desc, seq, qual).unwrap(),
        }
    }

    fn flush(&mut self) {
        match self {
            Writer::Fasta(w) => w.flush().unwrap(),
            Writer::Fastq(w) => w.flush().unwrap(),
        }
    }
}

fn filter_from(matches: &ArgMatches) -> Filter {
    let length = |name| matches.get_one::<u64>(name).map(|&x| x as usize);
    Filter {
        strip_left: length("strip-left").unwrap_or(0),
        strip_right: length("strip-right").unwrap_or(0),
        trunc_len: length("trunc-len"),
        trunc_keep: length("trunc-keep"),
        trunc_qual: matches.get_one::<u8>("trunc-qual").copied(),
        trunc_ee: matches.get_one::<f64>("trunc-prob").copied(),
        trunc_rate: matches.get_one::<f64>("trunc-rate").copied(),
        min_len: length("min").unwrap_or(1),
        max_len: length("max"),
        max_ns: length("bad"),
        max_ee: matches.get_one::<f64>("prob").copied(),
        max_rate: matches.get_one::<f64>("prob-rate").copied(),
        min_qual: matches.get_one::<u8>("min-qual").copied().unwrap_or(0),
        min_size: matches.get_one::<u64>("min-size").copied().unwrap_or(1),
        max_size: matches.get_one::<u64>("max-size").copied(),
    }
}

fn main() {
    let matches = Command::new("QualF")
        .about("Filter FASTQ sequences by quality scores.")
        .long_about(
            "Filter FASTQ sequences by quality scores.\n\n\
             Reads are trimmed first (strip, then truncate), and the filters apply \
             to the trimmed reads. Quality scores use the Phred+33 encoding.",
        )
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
            Arg::new("discard")
                .long("discard")
                .value_name("DISCARD")
                .help("Output file for discarded sequences, disabled by default"),
        )
        .arg(
            Arg::new("fasta")
                .long("fasta")
                .action(ArgAction::SetTrue)
                .help("Write output in FASTA format"),
        )
        .arg(
            Arg::new("ee")
                .long("ee")
                .action(ArgAction::SetTrue)
                .help("Annotate output labels with expected errors (;ee=)"),
        )
        .arg(
            Arg::new("strip-left")
                .long("strip-left")
                .value_name("LEN")
                .value_parser(value_parser!(u64))
                .help("Number of bases to remove from the 5' end"),
        )
        .arg(
            Arg::new("strip-right")
                .long("strip-right")
                .value_name("LEN")
                .value_parser(value_parser!(u64))
                .help("Number of bases to remove from the 3' end"),
        )
        .arg(
            Arg::new("trunc-len")
                .long("trunc-len")
                .value_name("LEN")
                .value_parser(value_parser!(u64).range(1..))
                .help("Truncate to this length, discarding shorter sequences"),
        )
        .arg(
            Arg::new("trunc-keep")
                .long("trunc-keep")
                .value_name("LEN")
                .value_parser(value_parser!(u64).range(1..))
                .help("Truncate to this length, keeping shorter sequences"),
        )
        .arg(
            Arg::new("trunc-qual")
                .long("trunc-qual")
                .value_name("QUAL")
                .value_parser(value_parser!(u8).range(0..=93))
                .help("Truncate at the first base with quality at or below this"),
        )
        .arg(
            Arg::new("trunc-prob")
                .long("trunc-prob")
                .value_name("PROB")
                .value_parser(non_negative)
                .help("Truncate before the expected errors exceed this"),
        )
        .arg(
            Arg::new("trunc-rate")
                .long("trunc-rate")
                .value_name("RATE")
                .value_parser(non_negative)
                .help("Truncate before the expected errors per base exceed this"),
        )
        .arg(
            Arg::new("max")
                .long("max")
                .value_name("MAX_LEN")
                .value_parser(value_parser!(u64).range(1..))
                .help("Maximum value for length filter, disabled by default"),
        )
        .arg(
            Arg::new("min")
                .long("min")
                .value_name("MIN_LEN")
                .value_parser(value_parser!(u64).range(1..))
                .help("Minimum value for length filter, default to 1"),
        )
        .arg(
            Arg::new("prob")
                .long("prob")
                .value_name("PROB")
                .value_parser(positive)
                .help("Maximum expected errors (sum of error probabilities), disabled by default"),
        )
        .arg(
            Arg::new("prob-rate")
                .long("prob-rate")
                .value_name("RATE")
                .value_parser(non_negative)
                .help("Maximum expected errors per base, disabled by default"),
        )
        .arg(
            Arg::new("min-qual")
                .long("min-qual")
                .value_name("QUAL")
                .value_parser(value_parser!(u8).range(0..=93))
                .help("Minimum quality of every base, disabled by default"),
        )
        .arg(
            Arg::new("bad")
                .long("bad")
                .value_name("BAD")
                .value_parser(value_parser!(u64))
                .help("Number of bad base call 'N' allowed, disabled by default"),
        )
        .arg(
            Arg::new("min-size")
                .long("min-size")
                .value_name("MIN_SIZE")
                .value_parser(value_parser!(u64).range(1..))
                .help("Minimum abundance (;size=), disabled by default"),
        )
        .arg(
            Arg::new("max-size")
                .long("max-size")
                .value_name("MAX_SIZE")
                .value_parser(value_parser!(u64).range(1..))
                .help("Maximum abundance (;size=), disabled by default"),
        )
        .get_matches();

    let fin: Box<dyn Read> = match matches.get_one::<String>("in").map(String::as_str) {
        Some(f) => Box::new(File::open(f).unwrap()),
        None => Box::new(io::stdin()),
    };

    let reader = Reader::new(BufReader::new(fin));

    let fout: Box<dyn Write> = match matches.get_one::<String>("out").map(String::as_str) {
        Some(f) => Box::new(File::create(f).unwrap()),
        None => Box::new(io::stdout()),
    };

    let as_fasta = matches.get_flag("fasta");

    let mut writer = Writer::new(fout, as_fasta);

    let mut discard = matches
        .get_one::<String>("discard")
        .map(|f| Writer::new(Box::new(File::create(f).unwrap()), as_fasta));

    let annotate = matches.get_flag("ee");

    let filter = filter_from(&matches);

    for record in reader.records() {
        // Since FASTQ uses '@' for both the description header and a quality symbol, there is no
        // reliable way to recover from a malformed record.
        // Panics if a FASTQ record cannot be read.
        let r = record.unwrap();

        if r.check().is_err() {
            continue;
        }

        let header = match r.desc() {
            Some(desc) => format!("{} {}", r.id(), desc),
            None => r.id().to_string(),
        };

        let outcome = filter.apply(r.seq(), r.qual(), abundance(&header));

        let out = if outcome.discarded {
            match discard.as_mut() {
                Some(w) => w,
                None => continue,
            }
        } else {
            &mut writer
        };

        let (seq, qual) = (outcome.trim(r.seq()), outcome.trim(r.qual()));

        if annotate {
            out.write(&annotate_ee(&header, outcome.ee), None, seq, qual);
        } else {
            out.write(r.id(), r.desc(), seq, qual);
        }
    }

    writer.flush();
    if let Some(w) = discard.as_mut() {
        w.flush();
    }
}
