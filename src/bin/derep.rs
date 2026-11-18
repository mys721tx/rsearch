use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{BufReader, BufWriter, Read, Write};

use clap::parser::ValueSource;
use clap::{Arg, ArgAction, ArgMatches, Command, value_parser};

use rsearch::cluster::{
    Cluster, Clusters, PrefixRead, Quality, dereplicate_prefixes, sort_by_size,
};
use rsearch::fastx::{Format, Reader, Record};
use rsearch::filter::QualityError;
use rsearch::label::{AbundanceError, Label, abundance, checked_abundance};
use rsearch::writer::Writer;

/// Options that need quality scores, which cannot be used with FASTA input.
const QUALITY_OPTIONS: [&str; 7] = [
    "ascii",
    "qmin",
    "qmax",
    "ascii-out",
    "qmin-out",
    "qmax-out",
    "best-qual",
];

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

fn create(path: &str) -> BufWriter<File> {
    BufWriter::new(File::create(path).unwrap())
}

/// the ASCII offset of an option that takes 33 or 64, default to 33.
fn offset(matches: &ArgMatches, name: &str) -> u8 {
    matches
        .get_one::<String>(name)
        .map_or(33, |a| a.parse().unwrap())
}

/// The accepted range of input quality scores.
struct ScoreRange {
    offset: u8,
    min: i16,
    max: i16,
}

impl ScoreRange {
    /// checks that every score is in range. As vsearch, reports the lowest
    /// score if it is below the minimum, and the highest otherwise.
    fn check(&self, qual: &[u8]) -> Result<(), QualityError> {
        // the validated ranges keep both within a byte. A symbol is in
        // range when its distance above the lowest is at most the width,
        // and a symbol below the lowest wraps around to a large distance.
        // Folding rather than stopping early lets this vectorize.
        let lowest = (i16::from(self.offset) + self.min) as u8;
        let width = (self.max - self.min) as u8;
        if !qual
            .iter()
            .fold(false, |out, &q| out | (q.wrapping_sub(lowest) > width))
        {
            return Ok(());
        }
        let lowest = i16::from(lowest);
        let error = |symbol: u8| QualityError {
            score: i16::from(symbol) - i16::from(self.offset),
            min_score: self.min,
            max_score: self.max,
        };
        let low = *qual.iter().min().unwrap();
        if i16::from(low) < lowest {
            Err(error(low))
        } else {
            Err(error(*qual.iter().max().unwrap()))
        }
    }
}

/// writes the clusters in vsearch's UCLUST format: a line per cluster and
/// per member, then a summary line per cluster.
fn write_uc(out: &mut impl Write, clusters: &[Cluster]) -> io::Result<()> {
    for (i, c) in clusters.iter().enumerate() {
        let head = &c.merged[0];
        writeln!(out, "S\t{i}\t{}\t*\t*\t*\t*\t*\t{head}\t*", c.seq.len())?;
        for j in 1..c.merged.len() {
            let strand = if c.reverse[j] { '-' } else { '+' };
            writeln!(
                out,
                "H\t{i}\t{}\t100.0\t{strand}\t0\t0\t*\t{}\t{head}",
                c.lengths[j], c.merged[j]
            )?;
        }
    }
    for (i, c) in clusters.iter().enumerate() {
        writeln!(out, "C\t{i}\t{}\t*\t*\t*\t*\t*\t{}\t*", c.size, c.merged[0])?;
    }
    out.flush()
}

/// writes vsearch's tab-separated cluster table: for each member, its
/// label, the cluster's label, the cluster's index, the member's index in
/// it, the number of members and the label of the first. As derep never
/// relabels, the cluster's label is that of its first member.
fn write_tabbed(out: &mut impl Write, clusters: &[Cluster]) -> io::Result<()> {
    for (i, c) in clusters.iter().enumerate() {
        let head = &c.merged[0];
        let count = c.merged.len();
        for (j, member) in c.merged.iter().enumerate() {
            writeln!(out, "{member}\t{head}\t{i}\t{j}\t{count}\t{head}")?;
        }
    }
    out.flush()
}

fn main() {
    let mut command = Command::new("DeRep")
        .about("Remove duplications in FASTA or FASTQ sequences.")
        .arg(
            Arg::new("in")
                .short('i')
                .long("in")
                .value_name("INPUT")
                .help("Input FASTA or FASTQ file, default to stdin"),
        )
        .arg(
            Arg::new("out")
                .short('o')
                .long("out")
                .value_name("OUTPUT")
                .help("Output file in the input format, default to stdout"),
        )
        .arg(
            Arg::new("strand")
                .short('s')
                .long("strand")
                .action(ArgAction::SetTrue)
                .help("Also merge reverse-complemented duplicates"),
        )
        .arg(
            Arg::new("same-label")
                .long("same-label")
                .action(ArgAction::SetTrue)
                .help("Merge only duplicates that also have identical labels"),
        )
        .arg(
            Arg::new("prefix")
                .long("prefix")
                .action(ArgAction::SetTrue)
                .conflicts_with_all(["strand", "same-label"])
                .help(
                    "Merge sequences into the longest sequence they are a prefix of; \
                     FASTA only",
                ),
        )
        .arg(
            Arg::new("unsorted")
                .long("unsorted")
                .action(ArgAction::SetTrue)
                .help("Output in order of first appearance instead of decreasing size"),
        )
        .arg(
            Arg::new("uc")
                .long("uc")
                .value_name("UC")
                .help("Output file of every cluster in UCLUST format"),
        )
        .arg(
            Arg::new("tabbedout")
                .long("tabbedout")
                .value_name("TABLE")
                .help("Output file of every cluster's members, tab-separated"),
        )
        .arg(
            Arg::new("no-sizein")
                .long("no-sizein")
                .action(ArgAction::SetTrue)
                .help("Count each input sequence once, ignoring abundances (;size=)"),
        )
        .arg(
            Arg::new("no-sizeout")
                .long("no-sizeout")
                .action(ArgAction::SetTrue)
                .help("Do not annotate output labels with abundance (;size=)"),
        )
        .arg(
            Arg::new("ascii")
                .long("ascii")
                .value_name("OFFSET")
                .value_parser(["33", "64"])
                .help("ASCII offset of input quality scores, default to 33"),
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
            Arg::new("ascii-out")
                .long("ascii-out")
                .value_name("OFFSET")
                .value_parser(["33", "64"])
                .help("ASCII offset of output quality scores, default to 33"),
        )
        .arg(
            Arg::new("qmin-out")
                .long("qmin-out")
                .value_name("QUAL")
                .value_parser(value_parser!(u8))
                .help("Lowest quality score written, default to 0"),
        )
        .arg(
            Arg::new("qmax-out")
                .long("qmax-out")
                .value_name("QUAL")
                .value_parser(value_parser!(u8))
                .help("Highest quality score written, default to 126 - ascii-out"),
        )
        .arg(
            Arg::new("best-qual")
                .long("best-qual")
                .action(ArgAction::SetTrue)
                .help(
                    "Keep the best quality score of merged sequences instead of averaging \
                     their error probabilities",
                ),
        );
    let matches = command.get_matches_mut();

    let range = ScoreRange {
        offset: offset(&matches, "ascii"),
        min: matches.get_one::<i16>("qmin").copied().unwrap_or(0),
        max: matches
            .get_one::<i16>("qmax")
            .copied()
            .unwrap_or(126 - i16::from(offset(&matches, "ascii"))),
    };
    let offset_out = offset(&matches, "ascii-out");
    let quality = Quality {
        offset: range.offset,
        offset_out,
        min_out: matches.get_one::<u8>("qmin-out").copied().unwrap_or(0),
        max_out: matches
            .get_one::<u8>("qmax-out")
            .copied()
            .unwrap_or(126 - offset_out),
        keep_best: matches.get_flag("best-qual"),
    };
    let (ascii, ascii_out) = (i16::from(range.offset), i16::from(offset_out));
    let range_error = if range.min > range.max {
        Some("--qmin cannot be greater than --qmax")
    } else if ascii + range.min < 33 {
        Some("the sum of --ascii and --qmin must be at least 33")
    } else if ascii + range.max > 126 {
        Some("the sum of --ascii and --qmax must be at most 126")
    } else if quality.min_out > quality.max_out {
        Some("--qmin-out cannot be greater than --qmax-out")
    } else if ascii_out + i16::from(quality.max_out) > 126 {
        Some("the sum of --ascii-out and --qmax-out must be at most 126")
    } else {
        None
    };
    if let Some(message) = range_error {
        command
            .error(clap::error::ErrorKind::ValueValidation, message)
            .exit();
    }

    let fin: Box<dyn Read> = match matches.get_one::<String>("in").map(String::as_str) {
        Some(f) => Box::new(File::open(f).unwrap()),
        None => Box::new(io::stdin()),
    };

    let mut reader = Reader::new(BufReader::new(fin));

    let format = reader.format().unwrap_or_else(|e| fail(&e.to_string()));
    let is_fastq = format == Some(Format::Fastq);

    if format == Some(Format::Fasta)
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

    let strand = matches.get_flag("strand");
    let prefix = matches.get_flag("prefix");
    if prefix && is_fastq {
        // vsearch --derep_prefix writes only FASTA, and changing the format
        // is conv's job
        command
            .error(
                clap::error::ErrorKind::ArgumentConflict,
                "--prefix cannot be used with FASTQ input",
            )
            .exit();
    }

    let sizein = !matches.get_flag("no-sizein");

    let mut seqs = if !is_fastq {
        Clusters::without_qualities()
    } else {
        Clusters::default().with_quality(quality)
    };
    if matches.get_flag("same-label") {
        seqs = seqs.by_label();
    }

    // with --prefix, every sequence is read before any is merged
    let mut prefix_reads = Vec::new();

    let mut record = Record::new();

    loop {
        if let Err(e) = reader.read(&mut record) {
            fail(&e.to_string());
        }

        if record.is_empty() {
            break;
        }

        // whole headers, as vsearch --notrunclabels
        let label = record.header();

        let size = if sizein {
            checked_abundance(label)
                .unwrap_or_else(|e: AbundanceError| fail(&format!("{e} in read {label}")))
        } else {
            1
        };

        if is_fastq && let Err(e) = range.check(record.qual()) {
            fail(&format!("{e} in read {}", record.id()));
        }

        if prefix {
            prefix_reads.push(PrefixRead {
                label: label.to_string(),
                size,
                annotated: abundance(label),
                seq: record.seq().to_vec(),
            });
        } else {
            seqs.add_read(label, size, &record, strand);
        }
    }

    if let Some(warning) = reader.stripped_warning() {
        warn(&warning);
    }

    let mut clusters = if prefix {
        dereplicate_prefixes(prefix_reads)
    } else {
        seqs.into_clusters()
    };

    if !matches.get_flag("unsorted") {
        sort_by_size(&mut clusters);
    }

    // the written label: that of the first merged sequence, with its size
    let written = Label {
        sizeout: !matches.get_flag("no-sizeout"),
        ..Label::default()
    };

    let fout: Box<dyn Write> = match matches.get_one::<String>("out").map(String::as_str) {
        Some(f) => Box::new(File::create(f).unwrap()),
        None => Box::new(io::stdout()),
    };

    let mut writer = Writer::new(fout, !is_fastq, 0);
    let table = quality.output_table();
    let mut qual = Vec::new();
    for c in &clusters {
        qual.clear();
        if let Some(q) = &c.qual {
            qual.extend(q.iter().map(|&q| table[usize::from(q)]));
        }
        let label = written.format(&c.merged[0], &c.seq, c.size, 0, None);
        writer.write(&label, None, &c.seq, &qual).unwrap();
    }
    writer.flush().unwrap();

    if let Some(f) = matches.get_one::<String>("uc") {
        write_uc(&mut create(f), &clusters).unwrap();
    }
    if let Some(f) = matches.get_one::<String>("tabbedout") {
        write_tabbed(&mut create(f), &clusters).unwrap();
    }
}
