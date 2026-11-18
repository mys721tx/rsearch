use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{BufReader, Read, Write};

use bio::io::{fasta, fastq};
use clap::parser::ValueSource;
use clap::{Arg, ArgAction, ArgGroup, ArgMatches, Command, value_parser};

use rsearch::fastx::{Format, Reader, Record};
use rsearch::filter::{Filter, Outcome};
use rsearch::label::{Label, Relabel, abundance};

/// Options that need quality scores, which cannot be used with FASTA input.
const QUALITY_OPTIONS: [&str; 7] = [
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

fn label_from(matches: &ArgMatches) -> Label {
    let relabel = if let Some(prefix) = matches.get_one::<String>("relabel") {
        Relabel::Prefix(prefix.clone())
    } else if matches.get_flag("relabel-self") {
        Relabel::Sequence
    } else if matches.get_flag("relabel-sha1") {
        Relabel::Sha1
    } else if matches.get_flag("relabel-md5") {
        Relabel::Md5
    } else {
        Relabel::None
    };
    Label {
        relabel,
        keep: matches.get_flag("relabel-keep"),
        sizeout: matches.get_flag("sizeout"),
        xsize: matches.get_flag("xsize"),
        eeout: matches.get_flag("ee"),
        xee: matches.get_flag("xee"),
        lengthout: matches.get_flag("lengthout"),
        xlength: matches.get_flag("xlength"),
        sample: matches.get_one::<String>("sample").cloned(),
        suffix: matches.get_one::<String>("label-suffix").cloned(),
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
        (self.filter.apply(r.seq(), qual, size), size)
    }

    /// writes the trimmed record, as the `ordinal`th of its output.
    fn write(
        &self,
        out: Option<&mut Writer>,
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
            out.write(r.id(), r.desc(), seq, qual);
        } else {
            let header = self
                .label
                .format(r.header(), seq, size, ordinal, outcome.ee);
            out.write(&header, None, seq, qual);
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
        )
        .arg(
            Arg::new("relabel")
                .long("relabel")
                .value_name("PREFIX")
                .help("Relabel sequences as PREFIX1, PREFIX2, ..."),
        )
        .arg(
            Arg::new("relabel-self")
                .long("relabel-self")
                .action(ArgAction::SetTrue)
                .help("Relabel sequences with the sequence itself"),
        )
        .arg(
            Arg::new("relabel-sha1")
                .long("relabel-sha1")
                .action(ArgAction::SetTrue)
                .help("Relabel sequences with the SHA1 digest of the sequence"),
        )
        .arg(
            Arg::new("relabel-md5")
                .long("relabel-md5")
                .action(ArgAction::SetTrue)
                .help("Relabel sequences with the MD5 digest of the sequence"),
        )
        .group(ArgGroup::new("relabeling").args([
            "relabel",
            "relabel-self",
            "relabel-sha1",
            "relabel-md5",
        ]))
        .arg(
            Arg::new("relabel-keep")
                .long("relabel-keep")
                .action(ArgAction::SetTrue)
                .help("Append the old label after the new one when relabelling"),
        )
        .arg(
            Arg::new("sizeout")
                .long("sizeout")
                .action(ArgAction::SetTrue)
                .help("Annotate output labels with abundance (;size=)"),
        )
        .arg(
            Arg::new("xsize")
                .long("xsize")
                .action(ArgAction::SetTrue)
                .help("Remove abundance annotations (;size=)"),
        )
        .arg(
            Arg::new("xee")
                .long("xee")
                .action(ArgAction::SetTrue)
                .help("Remove expected error annotations (;ee=)"),
        )
        .arg(
            Arg::new("lengthout")
                .long("lengthout")
                .action(ArgAction::SetTrue)
                .help("Annotate output labels with sequence length (;length=)"),
        )
        .arg(
            Arg::new("xlength")
                .long("xlength")
                .action(ArgAction::SetTrue)
                .help("Remove length annotations (;length=)"),
        )
        .arg(
            Arg::new("sample")
                .long("sample")
                .value_name("SAMPLE")
                .help("Annotate output labels with a sample name (;sample=)"),
        )
        .arg(
            Arg::new("label-suffix")
                .long("label-suffix")
                .value_name("SUFFIX")
                .help("Append this text to output labels"),
        );

    let matches = command.get_matches_mut();

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

    let mut writer = Writer::new(fout, as_fasta);

    let optional = |name: &str| {
        matches
            .get_one::<String>(name)
            .map(|f| Writer::new(create(f), as_fasta))
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
        w.flush();
    }
}
