use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader, Read, Write};

use clap::parser::ValueSource;
use clap::{Arg, ArgAction, ArgMatches, Command, value_parser};

use rsearch::fastx::{Format, Reader, Record};
use rsearch::label::{Label, abundance};
use rsearch::options;
use rsearch::quality::Conversion;
use rsearch::sff;
use rsearch::writer::Writer;

/// Options about input scores, which cannot be used with SFF input, as it
/// holds Phred scores as numbers.
const SFF_REFUSED: [&str; 4] = ["ascii", "qmin", "qmax", "solexa"];

/// Options about input scores or the lowest score written, which cannot be
/// used with FASTA input.
const FASTA_REFUSED: [&str; 5] = ["ascii", "qmin", "qmax", "solexa", "qmin-out"];

/// The score written at every base of FASTA input, by default.
const FASTA_QUALITY: i16 = 41;

/// The lowest Solexa score.
const SOLEXA_MIN: i16 = -5;

/// prints a warning, which may hold bytes that are not UTF-8.
fn warn(message: &[u8]) {
    let mut stderr = io::stderr().lock();
    let _ = stderr
        .write_all(b"Warning: ")
        .and_then(|()| stderr.write_all(message))
        .and_then(|()| stderr.write_all(b"\n"));
}

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
        .about("Convert the quality encoding of FASTQ sequences, or FASTA or SFF to FASTQ.")
        .long_about(
            "Convert the quality encoding of FASTQ sequences, or FASTA or SFF to FASTQ.\n\n\
             FASTQ scores are read with the input encoding, converted from the Solexa \
             scale if asked, limited to the output range and written with the output \
             offset. FASTA sequences get the highest output score at every base. SFF \
             reads keep their Phred scores, written with the output encoding.",
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
            Arg::new("clip")
                .long("clip")
                .action(ArgAction::SetTrue)
                .help("Remove the bases outside the clip points of SFF reads instead of writing them in lower case"),
        )
        .arg(
            Arg::new("sizeout")
                .long("sizeout")
                .action(ArgAction::SetTrue)
                .help("Annotate output labels with abundance (;size=)"),
        );
    let mut command = options::label_args(command);
    let matches = command.get_matches_mut();

    let fin: Box<dyn Read> = match matches.get_one::<String>("in").map(String::as_str) {
        Some(f) => {
            Box::new(File::open(f).unwrap_or_else(|e| fail(&format!("cannot open {f}: {e}"))))
        }
        None => Box::new(io::stdin()),
    };
    let mut fin = BufReader::new(fin);
    let is_sff = fin
        .fill_buf()
        .unwrap_or_else(|e| fail(&format!("cannot read the input: {e}")))
        .starts_with(sff::MAGIC);

    let mut input = if is_sff {
        Input::Sff(fin)
    } else {
        Input::Fastx(Reader::new(fin))
    };
    let is_fasta = match &mut input {
        Input::Fastx(reader) => {
            let format = reader
                .format()
                .unwrap_or_else(|_| fail("input is not SFF, FASTA or FASTQ"));
            format == Some(Format::Fasta)
        }
        Input::Sff(_) => false,
    };

    // options about input scores
    let refused: &[&str] = if is_sff {
        &SFF_REFUSED
    } else if is_fasta {
        &FASTA_REFUSED
    } else {
        &[]
    };
    let kind = if is_sff { "SFF" } else { "FASTA" };
    if let Some(option) = refused
        .iter()
        .find(|&&id| matches.value_source(id) == Some(ValueSource::CommandLine))
    {
        command
            .error(
                clap::error::ErrorKind::ArgumentConflict,
                format!("--{option} cannot be used with {kind} input"),
            )
            .exit();
    }
    if !is_sff && matches.get_flag("clip") {
        command
            .error(
                clap::error::ErrorKind::ArgumentConflict,
                "--clip needs SFF input",
            )
            .exit();
    }

    let conversion = conversion_from(&matches, is_fasta);
    if let Some(message) = range_error(&conversion) {
        command
            .error(clap::error::ErrorKind::ValueValidation, message)
            .exit();
    }

    let label = Label {
        sizeout: matches.get_flag("sizeout"),
        ..options::label_from(&matches)
    };

    let fout: Box<dyn Write> = match matches.get_one::<String>("out").map(String::as_str) {
        Some(f) => {
            Box::new(File::create(f).unwrap_or_else(|e| fail(&format!("cannot create {f}: {e}"))))
        }
        None => Box::new(io::stdout()),
    };
    let mut out = Output {
        writer: Writer::new(fout, false, 0),
        label,
    };

    match input {
        Input::Sff(fin) => convert_sff(fin, &mut out, &conversion, matches.get_flag("clip")),
        Input::Fastx(reader) => convert_fastx(reader, &mut out, &conversion, is_fasta),
    }

    out.writer.flush().unwrap();
}

/// The input: SFF, or FASTA or FASTQ.
enum Input<R: BufRead> {
    Sff(R),
    Fastx(Reader<R>),
}

/// Writes records with their new labels.
struct Output {
    writer: Writer<Box<dyn Write>>,
    label: Label,
}

impl Output {
    /// writes a record of abundance `size`, the `ordinal`th of the input.
    fn write(&mut self, header: &str, seq: &[u8], qual: &[u8], size: u64, ordinal: u64) {
        if self.label.is_identity() {
            self.writer.write(header, None, seq, qual).unwrap();
        } else {
            let header = self.label.format(header, seq, size, ordinal, None);
            self.writer.write(&header, None, seq, qual).unwrap();
        }
    }

    /// writes what is buffered and fails with `message`, so that the
    /// records before an error are written, as vsearch does.
    fn fail(&mut self, message: &str) -> ! {
        self.writer.flush().unwrap();
        fail(message);
    }
}

/// converts FASTQ scores, or gives FASTA sequences the highest score.
fn convert_fastx<R: BufRead>(
    mut reader: Reader<R>,
    out: &mut Output,
    conversion: &Conversion,
    is_fasta: bool,
) {
    let converter = conversion.converter();
    let fasta_symbol = (i16::from(conversion.offset_out) + conversion.max_out) as u8;

    let mut record = Record::new();
    let mut qual = Vec::new();
    let mut ordinal = 0;

    loop {
        if let Err(e) = reader.read(&mut record) {
            out.fail(&e.to_string());
        }

        if record.is_empty() {
            if let Some(warning) = reader.stripped_warning() {
                warn(&warning);
            }
            break;
        }

        ordinal += 1;

        if is_fasta {
            qual.clear();
            qual.resize(record.seq().len(), fasta_symbol);
        } else if let Err(e) = converter.convert(record.qual(), &mut qual) {
            out.fail(&format!("{e} in read {}", record.id()));
        }

        let header = record.header();
        let size = if out.label.sizeout {
            abundance(header)
        } else {
            1
        };
        out.write(header, record.seq(), &qual, size, ordinal);
    }
}

/// converts SFF reads: bases outside the clip points are lower case, or
/// with `clip` removed, and Phred scores are written with the output
/// encoding.
fn convert_sff<R: BufRead>(input: R, out: &mut Output, conversion: &Conversion, clip: bool) {
    // the scores are stored as numbers, not symbols
    let converter = Conversion {
        offset: 0,
        min: 0,
        max: 255,
        solexa: false,
        ..*conversion
    }
    .converter();

    let mut reader = sff::Reader::new(input).unwrap_or_else(|e| out.fail(&e.to_string()));
    let mut read = sff::SffRecord::default();
    let mut qual = Vec::new();
    let mut ordinal = 0;

    loop {
        match reader.read_record(&mut read) {
            Ok(true) => {}
            Ok(false) => break,
            Err(e) => out.fail(&e.to_string()),
        }
        ordinal += 1;

        let (start, end) = (read.clip_start, read.clip_end);
        // with inverted clip points, the clipped ends meet
        let kept = start.min(end)..end;
        read.bases.make_ascii_lowercase();
        read.bases[kept.clone()].make_ascii_uppercase();
        converter.convert(&read.quals, &mut qual).unwrap();

        let kept = if clip {
            if start > end {
                out.fail(
                    "Invalid SFF file. Clipping region is empty (clip start beyond clip end).",
                );
            }
            kept
        } else {
            0..read.bases.len()
        };
        out.write(
            &read.name,
            &read.bases[kept.clone()],
            &qual[kept],
            1,
            ordinal,
        );
    }
}
