use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{BufReader, Read, Write};

use clap::{Arg, ArgAction, Command, value_parser};

use rsearch::cluster::{Cluster, Clusters, sort_by_size};
use rsearch::fastx::{Format, Reader, Record};
use rsearch::label::{Label, size_annotation};
use rsearch::options;
use rsearch::writer::Writer;

fn main() {
    let command = Command::new("DeRep")
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
            Arg::new("min")
                .long("min")
                .value_name("MIN_SIZE")
                .value_parser(value_parser!(u64))
                .help("Minimum size of output sequences, disabled by default"),
        )
        .arg(
            Arg::new("max")
                .long("max")
                .value_name("MAX_SIZE")
                .value_parser(value_parser!(u64))
                .help("Maximum size of output sequences, disabled by default"),
        )
        .arg(
            Arg::new("top")
                .long("top")
                .value_name("TOP")
                .value_parser(value_parser!(usize))
                .help("Number of sequences to output, disabled by default"),
        )
        .arg(
            Arg::new("unsorted")
                .long("unsorted")
                .action(ArgAction::SetTrue)
                .help("Output in order of first appearance instead of decreasing size"),
        )
        .arg(
            Arg::new("map")
                .long("map")
                .value_name("MAP")
                .help("Output file mapping each merged sequence to its output sequence"),
        )
        .arg(
            Arg::new("fasta-width")
                .long("fasta-width")
                .value_name("WIDTH")
                .value_parser(value_parser!(usize))
                .help("Wrap FASTA sequence lines at this width, default to 0 (no wrapping)"),
        )
        .arg(
            Arg::new("fasta")
                .long("fasta")
                .action(ArgAction::SetTrue)
                .help("Write output in FASTA format"),
        )
        .arg(
            Arg::new("min-len")
                .long("min-len")
                .value_name("MIN_LEN")
                .value_parser(value_parser!(usize))
                .help("Minimum length of input sequences, default to 1"),
        )
        .arg(
            Arg::new("max-len")
                .long("max-len")
                .value_name("MAX_LEN")
                .value_parser(value_parser!(usize))
                .help("Maximum length of input sequences, disabled by default"),
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
            Arg::new("full-labels")
                .long("full-labels")
                .action(ArgAction::SetTrue)
                .help("Keep whole headers as labels instead of cutting them at the first space"),
        );
    let matches = options::label_args(command).get_matches();

    let fin: Box<dyn Read> = match matches.get_one::<String>("in").map(String::as_str) {
        Some(f) => Box::new(File::open(f).unwrap()),
        None => Box::new(io::stdin()),
    };

    let mut reader = Reader::new(BufReader::new(fin));

    let is_fastq = reader.format().unwrap() == Some(Format::Fastq);

    let strand = matches.get_flag("strand");

    let sizein = !matches.get_flag("no-sizein");
    let full_labels = matches.get_flag("full-labels");
    let min_len = matches.get_one::<usize>("min-len").copied().unwrap_or(1);
    let max_len = matches
        .get_one::<usize>("max-len")
        .copied()
        .unwrap_or(usize::MAX);

    let as_fasta = !is_fastq || matches.get_flag("fasta");

    let mut seqs = if as_fasta {
        Clusters::without_qualities()
    } else {
        Clusters::default()
    };

    let mut record = Record::new();

    loop {
        reader.read(&mut record).unwrap();

        if record.is_empty() {
            break;
        }

        let label = if full_labels {
            record.header()
        } else {
            record.id()
        };

        let size = if sizein {
            match size_annotation(label) {
                Some(0) => {
                    eprintln!("Error: invalid (zero) abundance annotation in read {label}");
                    std::process::exit(1);
                }
                Some(size) => size,
                None => 1,
            }
        } else {
            1
        };

        let len = record.seq().len();
        if len >= min_len && len <= max_len {
            seqs.add_read(label, size, &record, strand);
        }
    }

    let mut clusters = seqs.into_clusters();

    if !matches.get_flag("unsorted") {
        sort_by_size(&mut clusters);
    }

    let min = matches.get_one::<u64>("min").copied().unwrap_or(0);
    let max = matches.get_one::<u64>("max").copied().unwrap_or(u64::MAX);
    let top = matches
        .get_one::<usize>("top")
        .copied()
        .unwrap_or(usize::MAX);

    let selected = clusters
        .iter()
        .filter(|c| c.size >= min && c.size <= max)
        .take(top);

    // the written label: that of the first merged sequence, rewritten
    let written = Label {
        sizeout: !matches.get_flag("no-sizeout"),
        ..options::label_from(&matches)
    };
    let label =
        |i: usize, c: &Cluster| written.format(&c.merged[0], &c.seq, c.size, i as u64 + 1, None);

    let fout: Box<dyn Write> = match matches.get_one::<String>("out").map(String::as_str) {
        Some(f) => Box::new(File::create(f).unwrap()),
        None => Box::new(io::stdout()),
    };

    let width = matches
        .get_one::<usize>("fasta-width")
        .copied()
        .unwrap_or(0);
    let mut writer = Writer::new(fout, as_fasta, width);
    for (i, c) in selected.clone().enumerate() {
        let qual = c.qual.as_deref().unwrap_or_default();
        writer.write(&label(i, c), None, &c.seq, qual).unwrap();
    }
    writer.flush().unwrap();

    if let Some(f) = matches.get_one::<String>("map") {
        // labels without their size annotations
        let unsized_member = Label {
            xsize: true,
            ..Label::default()
        };
        let unsized_written = Label {
            sizeout: false,
            xsize: true,
            keep: false,
            ..written.clone()
        };
        let mut map = io::BufWriter::new(File::create(f).unwrap());
        for (i, c) in selected.enumerate() {
            let name = unsized_written.format(&c.merged[0], &c.seq, c.size, i as u64 + 1, None);
            for member in &c.merged {
                let member = unsized_member.format(member, &c.seq, c.size, 0, None);
                writeln!(map, "{member}\t{name}").unwrap();
            }
        }
        map.flush().unwrap();
    }
}
