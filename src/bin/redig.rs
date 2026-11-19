use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{BufReader, BufWriter, Read, Write};

use clap::{Arg, Command};

use rsearch::cluster::reverse_complement_keeping_case;
use rsearch::fastx::{Format, Reader, Record};
use rsearch::label::checked_abundance;
use rsearch::writer::Writer;

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

fn create(path: &str) -> Box<dyn Write> {
    Box::new(BufWriter::new(
        File::create(path).unwrap_or_else(|e| fail(&format!("cannot create {path}: {e}"))),
    ))
}

/// The 4-bit IUPAC code of each byte, a bit per base (A, C, G, T or U);
/// 0 for the bytes that are not IUPAC codes.
const CODE: [u8; 256] = {
    let mut code = [0u8; 256];
    let codes: [(u8, u8); 16] = [
        (b'A', 1),
        (b'B', 14),
        (b'C', 2),
        (b'D', 13),
        (b'G', 4),
        (b'H', 11),
        (b'K', 12),
        (b'M', 3),
        (b'N', 15),
        (b'R', 5),
        (b'S', 6),
        (b'T', 8),
        (b'U', 8),
        (b'V', 7),
        (b'W', 9),
        (b'Y', 10),
    ];
    let mut i = 0;
    while i < codes.len() {
        let (letter, value) = codes[i];
        code[letter as usize] = value;
        code[letter.to_ascii_lowercase() as usize] = value;
        i += 1;
    }
    code
};

/// A restriction pattern: its IUPAC codes, and where it cuts the forward
/// and the reverse strand, from its start.
#[derive(Clone, Debug)]
struct Pattern {
    codes: Vec<u8>,
    cut_forward: i64,
    cut_reverse: i64,
}

/// parses a pattern of IUPAC codes with one ^ and one _, as vsearch's
/// `--cut_pattern`.
fn pattern(text: &str) -> Result<Pattern, String> {
    for (mark, missing) in [
        ('^', "no forward cut site (^) in the pattern"),
        ('_', "no reverse cut site (_) in the pattern"),
    ] {
        match text.matches(mark).count() {
            0 => return Err(missing.to_string()),
            1 => {}
            _ => return Err("more than one cut site in the pattern".to_string()),
        }
    }
    let bases: Vec<u8> = text.bytes().filter(|&c| c != b'^' && c != b'_').collect();
    if bases.is_empty() {
        return Err("no bases in the pattern".to_string());
    }
    if bases.iter().any(|&c| CODE[usize::from(c)] == 0) {
        return Err("expected IUPAC codes besides ^ and _".to_string());
    }
    let position = |mark: char, other: char| text.replace(other, "").find(mark).unwrap() as i64;
    Ok(Pattern {
        codes: bases.iter().map(|&c| CODE[usize::from(c)]).collect(),
        cut_forward: position('^', '_'),
        cut_reverse: position('_', '^'),
    })
}

/// The outputs, all but the forward fragments optional.
struct Outputs {
    forward: Writer<Box<dyn Write>>,
    reverse: Option<Writer<Box<dyn Write>>>,
    uncut: Option<Writer<Box<dyn Write>>>,
    rev_uncut: Option<Writer<Box<dyn Write>>>,
}

impl Outputs {
    fn all(&mut self) -> impl Iterator<Item = &mut Writer<Box<dyn Write>>> {
        std::iter::once(&mut self.forward).chain(
            [&mut self.reverse, &mut self.uncut, &mut self.rev_uncut]
                .into_iter()
                .flatten(),
        )
    }
}

/// writes a fragment to an output, if there is one.
fn write(writer: Option<&mut Writer<Box<dyn Write>>>, header: &str, seq: &[u8]) {
    if let Some(writer) = writer {
        writer
            .write(header, None, seq, &[])
            .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
    }
}

/// cuts a sequence at every match of the pattern, as vsearch's `--cut`:
/// matches may overlap, and a fragment that would be empty is skipped.
fn cut(
    record: &Record,
    pattern: &Pattern,
    outputs: &mut Outputs,
    coded: &mut Vec<u8>,
    rc: &mut Vec<u8>,
) {
    let (header, seq) = (record.header(), record.seq());
    if outputs.reverse.is_some() || outputs.rev_uncut.is_some() {
        reverse_complement_keeping_case(rc, seq);
    }
    coded.clear();
    coded.extend(seq.iter().map(|&c| CODE[usize::from(c)]));

    let length = seq.len() as i64;
    // the start of the next forward fragment, and the end of the next
    // reverse one
    let (mut start, mut rc_start) = (0i64, length);
    let mut matches = 0;
    let span = |from: i64, to: i64| from as usize..to as usize;
    if seq.len() >= pattern.codes.len() {
        for (i, window) in coded.windows(pattern.codes.len()).enumerate() {
            // codes match when they share a base
            if !window.iter().zip(&pattern.codes).all(|(&s, &p)| s & p != 0) {
                continue;
            }
            matches += 1;
            let i = i as i64;
            let fragment = i + pattern.cut_forward - start;
            let rc_fragment = rc_start - (length - (i + pattern.cut_reverse));
            rc_start -= rc_fragment;
            if fragment > 0 {
                write(
                    Some(&mut outputs.forward),
                    header,
                    &seq[span(start, start + fragment)],
                );
            }
            if rc_fragment > 0 && outputs.reverse.is_some() {
                write(
                    outputs.reverse.as_mut(),
                    header,
                    &rc[span(rc_start, rc_start + rc_fragment)],
                );
            }
            start += fragment;
        }
    }

    if matches > 0 {
        if length > start {
            write(
                Some(&mut outputs.forward),
                header,
                &seq[span(start, length)],
            );
        }
        if rc_start > 0 && outputs.reverse.is_some() {
            write(outputs.reverse.as_mut(), header, &rc[span(0, rc_start)]);
        }
    } else {
        write(outputs.uncut.as_mut(), header, seq);
        if outputs.rev_uncut.is_some() {
            write(outputs.rev_uncut.as_mut(), header, rc);
        }
    }
}

fn main() {
    let file = |name: &'static str, help: &'static str| {
        Arg::new(name).long(name).value_name("FILE").help(help)
    };
    let matches = Command::new("redig")
        .about("Cut sequences at a restriction pattern.")
        .long_about(
            "Cut sequences at a restriction pattern.\n\n\
             The pattern is IUPAC codes with ^ where the forward strand is cut and _ \
             where the reverse strand is, as G^AATT_C for EcoRI. Codes match when they \
             share a base, and matches may overlap. The fragments of each strand are \
             written in order, the reverse strand's from its 3' end, with the header of \
             their sequence. Sequences never cut go to --uncut. FASTA in and out, \
             unwrapped.",
        )
        .arg(
            Arg::new("in")
                .short('i')
                .long("in")
                .value_name("INPUT")
                .help("Input FASTA file, default to stdin"),
        )
        .arg(
            Arg::new("out")
                .short('o')
                .long("out")
                .value_name("OUTPUT")
                .help("Output file of the fragments of the forward strand, default to stdout"),
        )
        .arg(
            Arg::new("pattern")
                .long("pattern")
                .value_name("PATTERN")
                .required(true)
                .value_parser(pattern)
                .help("Restriction pattern, with ^ and _ at the cut sites"),
        )
        .arg(file(
            "rev-out",
            "Output file of the fragments of the reverse strand",
        ))
        .arg(file("uncut", "Output file of the sequences never cut"))
        .arg(file(
            "rev-uncut",
            "Output file of the reverse complements of the sequences never cut",
        ))
        .get_matches();
    let pattern = matches.get_one::<Pattern>("pattern").unwrap();

    let fin: Box<dyn Read> = match matches.get_one::<String>("in").map(String::as_str) {
        Some(f) => {
            Box::new(File::open(f).unwrap_or_else(|e| fail(&format!("cannot open {f}: {e}"))))
        }
        None => Box::new(io::stdin()),
    };
    let mut reader = Reader::new(BufReader::new(fin));
    if reader.format().unwrap_or_else(|e| fail(&e.to_string())) == Some(Format::Fastq) {
        fail("FASTA file expected, FASTQ file found");
    }

    let writer = |out: Box<dyn Write>| Writer::new(out, true, 0);
    let output = |name: &str| matches.get_one::<String>(name).map(|f| writer(create(f)));
    let mut outputs = Outputs {
        forward: writer(match matches.get_one::<String>("out").map(String::as_str) {
            Some(f) => create(f),
            None => Box::new(BufWriter::new(io::stdout())),
        }),
        reverse: output("rev-out"),
        uncut: output("uncut"),
        rev_uncut: output("rev-uncut"),
    };
    let fail_after = |outputs: &mut Outputs, message: &str| -> ! {
        // the fragments before an error are written, as vsearch does
        for writer in outputs.all() {
            let _ = writer.flush();
        }
        fail(message);
    };

    let mut record = Record::new();
    let (mut coded, mut rc) = (Vec::new(), Vec::new());
    loop {
        if let Err(e) = reader.read(&mut record) {
            fail_after(&mut outputs, &e.to_string());
        }
        if record.is_empty() {
            break;
        }
        // vsearch reads, and checks, every abundance
        if let Err(e) = checked_abundance(record.header()) {
            fail_after(&mut outputs, &format!("{e} in read {}", record.header()));
        }
        cut(&record, pattern, &mut outputs, &mut coded, &mut rc);
    }

    for writer in outputs.all() {
        writer
            .flush()
            .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
    }
    if let Some(warning) = reader.stripped_warning() {
        warn(&warning);
    }
}
