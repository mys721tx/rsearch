use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{BufReader, Read, Write};

use clap::parser::ValueSource;
use clap::{Arg, ArgAction, ArgMatches, Command, value_parser};

use rsearch::fastx::{Format, Reader, Record};
use rsearch::label::{Label, abundance};
use rsearch::options;
use rsearch::quality::Conversion;
use rsearch::writer::Writer;

/// Options about input scores, so cannot be used with FASTA input.
const QUALITY_OPTIONS: [&str; 5] = ["ascii", "qmin", "qmax", "solexa", "qmin-out"];

/// The score written at every base of FASTA input, by default.
const FASTA_QUALITY: i16 = 41;

/// The lowest Solexa score.
const SOLEXA_MIN: i16 = -5;

/// prints an error and exits.
fn fail(message: &str) -> ! {
    eprintln!("Error: {message}");
    std::process::exit(1);
}

/// the ASCII offset of an option that takes 33 or 64, default to 33.
fn offset(matches: &ArgMatches, name: &str) -> u8 {
    matches
        .get_one::<String>(name)
        .map_or(33, |a| a.parse().unwrap())
}

/// the conversion the options ask for; `fasta` changes the default
/// highest output score to `FASTA_QUALITY`.
fn conversion_from(matches: &ArgMatches, fasta: bool) -> Conversion {
    let offset_in = offset(matches, "ascii");
    let offset_out = offset(matches, "ascii-out");
    let solexa = matches.get_flag("solexa");
    let score = |name: &str| matches.get_one::<i16>(name).copied();
    Conversion {
        offset: offset_in,
        min: score("qmin").unwrap_or(if solexa { SOLEXA_MIN } else { 0 }),
        max: score("qmax").unwrap_or(126 - i16::from(offset_in)),
        solexa,
        offset_out,
        min_out: score("qmin-out").unwrap_or(0),
        max_out: score("qmax-out").unwrap_or(if fasta {
            FASTA_QUALITY
        } else {
            126 - i16::from(offset_out)
        }),
    }
}

/// the first problem with the encodings, if any.
fn range_error(c: &Conversion) -> Option<&'static str> {
    let (ascii, ascii_out) = (i16::from(c.offset), i16::from(c.offset_out));
    if c.solexa && c.offset != 64 {
        Some("--solexa needs --ascii 64")
    } else if c.min > c.max {
        Some("--qmin cannot be greater than --qmax")
    } else if ascii + c.min < 33 {
        Some("the sum of --ascii and --qmin must be at least 33")
    } else if ascii + c.max > 126 {
        Some("the sum of --ascii and --qmax must be at most 126")
    } else if c.min_out > c.max_out {
        Some("--qmin-out cannot be greater than --qmax-out")
    } else if ascii_out + c.min_out < 33 {
        Some("the sum of --ascii-out and --qmin-out must be at least 33")
    } else if ascii_out + c.max_out > 126 {
        Some("the sum of --ascii-out and --qmax-out must be at most 126")
    } else {
        None
    }
}

fn main() {
    let score = |name: &'static str, help: &'static str| {
        Arg::new(name)
            .long(name)
            .value_name("QUAL")
            .allow_negative_numbers(true)
            .value_parser(value_parser!(i16))
            .help(help)
    };
    let command = Command::new("Conv")
        .about("Convert the quality encoding of FASTQ sequences, or FASTA to FASTQ.")
        .long_about(
            "Convert the quality encoding of FASTQ sequences, or FASTA to FASTQ.\n\n\
             FASTQ scores are read with the input encoding, converted from the Solexa \
             scale if asked, limited to the output range and written with the output \
             offset. FASTA sequences get the highest output score at every base.",
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
                .help("Output FASTQ file, default to stdout"),
        )
        .arg(
            Arg::new("ascii")
                .long("ascii")
                .value_name("OFFSET")
                .value_parser(["33", "64"])
                .help("ASCII offset of input quality scores, default to 33"),
        )
        .arg(score(
            "qmin",
            "Minimum quality score accepted in the input, default to 0 (-5 with --solexa)",
        ))
        .arg(score(
            "qmax",
            "Maximum quality score accepted in the input, default to 126 - ascii",
        ))
        .arg(
            Arg::new("solexa")
                .long("solexa")
                .action(ArgAction::SetTrue)
                .help(
                    "Input scores are on the Solexa scale, converted to Phred (needs --ascii 64)",
                ),
        )
        .arg(
            Arg::new("ascii-out")
                .long("ascii-out")
                .value_name("OFFSET")
                .value_parser(["33", "64"])
                .help("ASCII offset of output quality scores, default to 33"),
        )
        .arg(score(
            "qmin-out",
            "Lowest quality score written, default to 0",
        ))
        .arg(score(
            "qmax-out",
            "Highest quality score written, default to 126 - ascii-out, or 41 for FASTA input",
        ))
        .arg(
            Arg::new("sizeout")
                .long("sizeout")
                .action(ArgAction::SetTrue)
                .help("Annotate output labels with abundance (;size=)"),
        );
    let mut command = options::label_args(command);
    let matches = command.get_matches_mut();

    let fin: Box<dyn Read> = match matches.get_one::<String>("in").map(String::as_str) {
        Some(f) => Box::new(File::open(f).unwrap()),
        None => Box::new(io::stdin()),
    };

    let mut reader = Reader::new(BufReader::new(fin));

    let is_fasta = reader.format().unwrap() == Some(Format::Fasta);

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

    let conversion = conversion_from(&matches, is_fasta);
    if let Some(message) = range_error(&conversion) {
        command
            .error(clap::error::ErrorKind::ValueValidation, message)
            .exit();
    }
    let converter = conversion.converter();

    // the score written at every base of FASTA input
    let fasta_symbol = (i16::from(conversion.offset_out) + conversion.max_out) as u8;

    let label = Label {
        sizeout: matches.get_flag("sizeout"),
        ..options::label_from(&matches)
    };

    let fout: Box<dyn Write> = match matches.get_one::<String>("out").map(String::as_str) {
        Some(f) => Box::new(File::create(f).unwrap()),
        None => Box::new(io::stdout()),
    };
    let mut writer = Writer::new(fout, false, 0);

    let mut record = Record::new();
    let mut qual = Vec::new();
    let mut ordinal = 0;

    loop {
        reader.read(&mut record).unwrap();

        if record.is_empty() {
            break;
        }

        ordinal += 1;

        if is_fasta {
            qual.clear();
            qual.resize(record.seq().len(), fasta_symbol);
        } else if let Err(e) = converter.convert(record.qual(), &mut qual) {
            // the reads before are written, as vsearch does
            writer.flush().unwrap();
            fail(&format!("{e} in read {}", record.id()));
        }

        let header = record.header();
        if label.is_identity() {
            writer.write(header, None, record.seq(), &qual).unwrap();
        } else {
            let size = if label.sizeout { abundance(header) } else { 1 };
            let header = label.format(header, record.seq(), size, ordinal, None);
            writer.write(&header, None, record.seq(), &qual).unwrap();
        }
    }

    writer.flush().unwrap();
}
