use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{BufReader, Read, Write};

use clap::parser::ValueSource;
use clap::{Arg, ArgAction, ArgMatches, Command, value_parser};

use rsearch::fastx::{Format, Reader, Record};
use rsearch::filter::{Filter, Outcome};
use rsearch::label::{Label, abundance};
use rsearch::options;
use rsearch::writer::Writer;

/// Options that need quality scores, which cannot be used with FASTA input.
const QUALITY_OPTIONS: [&str; 10] = [
    "ascii",
    "qmin",
    "qmax",
    "prob",
    "prob-rate",
    "min-qual",
    "trunc-qual",
    "trunc-prob",
    "trunc-rate",
    "ee",
];

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

fn filter_from(matches: &ArgMatches) -> Filter {
    let length = |name| matches.get_one::<u64>(name).map(|&x| x as usize);
    Filter {
        strip_left: length("strip-left").unwrap_or(0),
        strip_right: length("strip-right").unwrap_or(0),
        trunc_len: length("trunc-len"),
        trunc_keep: length("trunc-keep"),
        trunc_qual: matches.get_one::<u8>("trunc-qual").map(|&q| i16::from(q)),
        trunc_ee: matches.get_one::<f64>("trunc-prob").copied(),
        trunc_rate: matches.get_one::<f64>("trunc-rate").copied(),
        min_len: length("min").unwrap_or(1),
        max_len: length("max"),
        max_ns: length("bad"),
        max_ee: matches.get_one::<f64>("prob").copied(),
        max_rate: matches.get_one::<f64>("prob-rate").copied(),
        min_qual: matches
            .get_one::<u8>("min-qual")
            .map_or(0, |&q| i16::from(q)),
        min_size: matches.get_one::<u64>("min-size").copied().unwrap_or(1),
        max_size: matches.get_one::<u64>("max-size").copied(),
        offset: ascii(matches),
        min_score: matches.get_one::<i16>("qmin").copied().unwrap_or(0),
        max_score: max_score(matches),
    }
}

/// the ASCII offset of quality symbols.
fn ascii(matches: &ArgMatches) -> u8 {
    matches
        .get_one::<String>("ascii")
        .map_or(33, |a| a.parse().unwrap())
}

/// the highest accepted quality score: by default the highest the offset
/// can encode in printable ASCII.
fn max_score(matches: &ArgMatches) -> i16 {
    matches
        .get_one::<i16>("qmax")
        .copied()
        .unwrap_or(126 - i16::from(ascii(matches)))
}

fn label_from(matches: &ArgMatches) -> Label {
    Label {
        sizeout: matches.get_flag("sizeout"),
        eeout: matches.get_flag("ee"),
        ..options::label_from(matches)
    }
}

/// prints an error and exits.
fn fail(message: &str) -> ! {
    eprintln!("Error: {message}");
    std::process::exit(1);
}

fn open(path: &str) -> Box<dyn Read> {
    Box::new(File::open(path).unwrap())
}

fn create(path: &str) -> Box<dyn Write> {
    Box::new(File::create(path).unwrap())
}

/// Trims and filters records, and writes them with their new labels.
struct Qualf {
    filter: Filter,
    label: Label,
}

impl Qualf {
    /// trims and filters a record; returns the outcome and the abundance.
    fn apply(&self, r: &Record) -> (Outcome, u64) {
        let size = if self.filter.uses_size() || self.label.sizeout {
            abundance(r.header())
        } else {
            1
        };
        let qual = r.is_fastq().then(|| r.qual());
        let outcome = self
            .filter
            .apply(r.seq(), qual, size)
            .unwrap_or_else(|e| fail(&format!("{e} in read {}", r.id())));
        (outcome, size)
    }

    /// writes the trimmed record, as the `ordinal`th of its output.
    fn write(
        &self,
        out: Option<&mut Writer<Box<dyn Write>>>,
        r: &Record,
        outcome: &Outcome,
        size: u64,
        ordinal: u64,
    ) {
        let Some(out) = out else {
            return;
        };
        let seq = outcome.trim(r.seq());
        let qual = if r.is_fastq() {
            outcome.trim(r.qual())
        } else {
            &[]
        };
        if self.label.is_identity() {
            out.write(r.id(), r.desc(), seq, qual).unwrap();
        } else {
            let header = self
                .label
                .format(r.header(), seq, size, ordinal, outcome.ee);
            out.write(&header, None, seq, qual).unwrap();
        }
    }
}

fn main() {
    let mut command = Command::new("QualF")
        .about("Filter FASTQ sequences by quality scores.")
        .long_about(
            "Filter FASTQ sequences by quality scores.\n\n\
             Reads are trimmed first (strip, then truncate), and the filters apply \
             to the trimmed reads. Quality scores use the Phred+33 encoding. FASTA \
             input can be filtered by length, N bases and abundance.",
        )
        .arg(
            Arg::new("in")
                .short('i')
                .long("in")
                .value_name("INPUT")
                .help("Input FASTQ or FASTA file, default to stdin"),
        )
        .arg(
            Arg::new("out")
                .short('o')
                .long("out")
                .value_name("OUTPUT")
                .help("Output file in the input format, default to stdout"),
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
        .arg(
            Arg::new("reverse")
                .long("reverse")
                .value_name("REVERSE")
                .help("Input file of reverse reads, filtered together with their pairs"),
        )
        .arg(
            Arg::new("rev-out")
                .long("rev-out")
                .value_name("OUTPUT")
                .requires("reverse")
                .help("Output file for kept reverse reads, disabled by default"),
        )
        .arg(
            Arg::new("rev-discard")
                .long("rev-discard")
                .value_name("DISCARD")
                .requires("reverse")
                .help("Output file for discarded reverse reads, disabled by default"),
        );
    command = options::label_args(command)
        .arg(
            Arg::new("sizeout")
                .long("sizeout")
                .action(ArgAction::SetTrue)
                .help("Annotate output labels with abundance (;size=)"),
        )
        .arg(
            Arg::new("ascii")
                .long("ascii")
                .value_name("OFFSET")
                .value_parser(["33", "64"])
                .help("ASCII offset of quality scores, default to 33"),
        )
        .arg(
            Arg::new("qmin")
                .long("qmin")
                .value_name("QUAL")
                .allow_negative_numbers(true)
                .value_parser(value_parser!(i16))
                .help("Minimum quality score accepted in the input, default to 0"),
        )
        .arg(
            Arg::new("qmax")
                .long("qmax")
                .value_name("QUAL")
                .allow_negative_numbers(true)
                .value_parser(value_parser!(i16))
                .help("Maximum quality score accepted in the input, default to 126 - ascii"),
        )
        .arg(
            Arg::new("fasta-width")
                .long("fasta-width")
                .value_name("WIDTH")
                .value_parser(value_parser!(usize))
                .help("Wrap FASTA sequence lines at this width, default to 0 (no wrapping)"),
        )
        .arg(
            Arg::new("quiet")
                .short('q')
                .long("quiet")
                .action(ArgAction::SetTrue)
                .help("Do not print the summary of kept and discarded sequences"),
        );

    let matches = command.get_matches_mut();

    let (offset, min_score, max_score) = (
        i16::from(ascii(&matches)),
        matches.get_one::<i16>("qmin").copied().unwrap_or(0),
        max_score(&matches),
    );
    let range_error = if min_score > max_score {
        Some("--qmin cannot be greater than --qmax")
    } else if offset + min_score < 33 {
        Some("the sum of --ascii and --qmin must be at least 33")
    } else if offset + max_score > 126 {
        Some("the sum of --ascii and --qmax must be at most 126")
    } else {
        None
    };
    if let Some(message) = range_error {
        command
            .error(clap::error::ErrorKind::ValueValidation, message)
            .exit();
    }

    let fin: Box<dyn Read> = match matches.get_one::<String>("in").map(String::as_str) {
        Some(f) => open(f),
        None => Box::new(io::stdin()),
    };

    let mut reader = Reader::new(BufReader::new(fin));

    let format = reader.format().unwrap();

    let mut reverse = matches
        .get_one::<String>("reverse")
        .map(|f| Reader::new(BufReader::new(open(f))));

    if let Some(rev) = reverse.as_mut() {
        let rev_format = rev.format().unwrap();
        if format.is_some() && rev_format.is_some() && rev_format != format {
            fail("forward and reverse reads must be in the same format");
        }
    }

    let is_fasta = format == Some(Format::Fasta);

    if is_fasta
        && let Some(option) = QUALITY_OPTIONS
            .iter()
            .find(|&&id| matches.value_source(id) == Some(ValueSource::CommandLine))
    {
        command
            .error(
                clap::error::ErrorKind::ArgumentConflict,
                format!("--{option} needs quality scores, which FASTA input lacks"),
            )
            .exit();
    }

    // FASTA input can only be written as FASTA
    let as_fasta = matches.get_flag("fasta") || is_fasta;

    let fout: Box<dyn Write> = match matches.get_one::<String>("out").map(String::as_str) {
        Some(f) => create(f),
        None => Box::new(io::stdout()),
    };

    let width = matches
        .get_one::<usize>("fasta-width")
        .copied()
        .unwrap_or(0);

    let mut writer = Writer::new(fout, as_fasta, width);

    let optional = |name: &str| {
        matches
            .get_one::<String>(name)
            .map(|f| Writer::new(create(f), as_fasta, width))
    };

    let mut discard = optional("discard");
    let mut rev_writer = optional("rev-out");
    let mut rev_discard = optional("rev-discard");

    let qualf = Qualf {
        filter: filter_from(&matches),
        label: label_from(&matches),
    };

    let (mut r, mut r2) = (Record::new(), Record::new());

    // ordinals of kept and discarded sequences, for --relabel
    let (mut kept, mut discarded) = (0, 0);

    // kept sequences (or pairs) with bases removed
    let mut truncated = 0;

    loop {
        // Since FASTQ uses '@' for both the description header and a quality symbol, there is no
        // reliable way to recover from a malformed record.
        // Panics if a record cannot be read.
        reader.read(&mut r).unwrap();

        if let Some(rev) = reverse.as_mut() {
            rev.read(&mut r2).unwrap();
            match (r.is_empty(), r2.is_empty()) {
                (true, false) => fail("more reverse reads than forward reads"),
                (false, true) => fail("more forward reads than reverse reads"),
                _ => {}
            }
        }

        if r.is_empty() {
            break;
        }

        // a pair is skipped if either read is invalid
        if r.check().is_err() || (reverse.is_some() && r2.check().is_err()) {
            continue;
        }

        let (outcome, size) = qualf.apply(&r);
        let rev = reverse.as_ref().map(|_| qualf.apply(&r2));

        // a pair is discarded if either read is
        if outcome.discarded || rev.as_ref().is_some_and(|(o, _)| o.discarded) {
            discarded += 1;
            qualf.write(discard.as_mut(), &r, &outcome, size, discarded);
            if let Some((o, s)) = &rev {
                qualf.write(rev_discard.as_mut(), &r2, o, *s, discarded);
            }
        } else {
            kept += 1;
            if outcome.truncated || rev.as_ref().is_some_and(|(o, _)| o.truncated) {
                truncated += 1;
            }
            qualf.write(Some(&mut writer), &r, &outcome, size, kept);
            if let Some((o, s)) = &rev {
                qualf.write(rev_writer.as_mut(), &r2, o, *s, kept);
            }
        }
    }

    for w in [
        Some(&mut writer),
        discard.as_mut(),
        rev_writer.as_mut(),
        rev_discard.as_mut(),
    ]
    .into_iter()
    .flatten()
    {
        w.flush().unwrap();
    }

    if !matches.get_flag("quiet") {
        eprintln!(
            "{kept} sequences kept (of which {truncated} truncated), {discarded} sequences discarded."
        );
    }
}
