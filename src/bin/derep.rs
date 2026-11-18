use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{BufReader, Read, Write};

use bio::io::{fasta, fastq};
use clap::{Arg, ArgAction, Command, value_parser};

use rsearch::cluster::{Cluster, Clusters, sort_by_size};
use rsearch::fastx::{Format, Reader, Record};
use rsearch::label::{Label, Relabel, size_annotation};

fn main() {
    let matches = Command::new("DeRep")
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
            Arg::new("relabel")
                .long("relabel")
                .value_name("PREFIX")
                .help("Relabel output sequences as PREFIX1, PREFIX2, ..."),
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
        .get_matches();

    let fin: Box<dyn Read> = match matches.get_one::<String>("in").map(String::as_str) {
        Some(f) => Box::new(File::open(f).unwrap()),
        None => Box::new(io::stdin()),
    };

    let mut reader = Reader::new(BufReader::new(fin));

    let is_fastq = reader.format().unwrap() == Some(Format::Fastq);

    let strand = matches.get_flag("strand");

    let mut seqs = Clusters::default();

    let mut record = Record::new();

    loop {
        reader.read(&mut record).unwrap();

        if record.is_empty() {
            break;
        }

        if size_annotation(record.id()) == Some(0) {
            eprintln!(
                "Error: invalid (zero) abundance annotation in read {}",
                record.id()
            );
            std::process::exit(1);
        }

        if !record.seq().is_empty() {
            seqs.add_fastx(&record, strand);
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

    let relabel = match matches.get_one::<String>("relabel") {
        Some(prefix) => Relabel::Prefix(prefix.clone()),
        None => Relabel::None,
    };

    // the written label: the first merged label with its size replaced,
    // or the new label, as vsearch --sizeout writes it
    let written = Label {
        relabel: relabel.clone(),
        sizeout: true,
        ..Label::default()
    };
    let label =
        |i: usize, c: &Cluster| written.format(&c.merged[0], &c.seq, c.size, i as u64 + 1, None);

    let fout: Box<dyn Write> = match matches.get_one::<String>("out").map(String::as_str) {
        Some(f) => Box::new(File::create(f).unwrap()),
        None => Box::new(io::stdout()),
    };

    if is_fastq {
        let mut writer = fastq::Writer::new(fout);
        for (i, c) in selected.clone().enumerate() {
            let qual = c.qual.as_deref().unwrap_or_default();
            writer.write(&label(i, c), None, &c.seq, qual).unwrap();
        }
        writer.flush().unwrap();
    } else {
        let width = matches
            .get_one::<usize>("fasta-width")
            .copied()
            .unwrap_or(0);
        let mut writer = fasta::Writer::new(fout);
        writer.set_linewrap((width > 0).then_some(width));
        for (i, c) in selected.clone().enumerate() {
            writer.write(&label(i, c), None, &c.seq).unwrap();
        }
        writer.flush().unwrap();
    }

    if let Some(f) = matches.get_one::<String>("map") {
        // labels without their size annotations
        let unsized_member = Label {
            xsize: true,
            ..Label::default()
        };
        let unsized_written = Label {
            relabel,
            xsize: true,
            ..Label::default()
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
