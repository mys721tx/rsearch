use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fmt::Write as _;
use std::fs::File;
use std::io;
use std::io::{BufReader, Read, Write};

use clap::{Arg, Command, value_parser};

use rsearch::fastx::{Format, Reader, Record};
use rsearch::quality::Encoding;

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

/// the symbols with a count, in increasing order.
fn present(counts: &[u64; 256]) -> impl Iterator<Item = u8> + '_ {
    (0..=255u8).filter(|&c| counts[usize::from(c)] > 0)
}

/// Counts of residues and score symbols, as vsearch's `--fastq_chars`.
struct Counts {
    reads: u64,
    residues: u64,
    /// Per residue, upper-cased: count and longest run.
    letters: [[u64; 256]; 4],
    max_run: [u64; 256],
    /// Per score symbol: count and reads ending in a run of it.
    scores: [[u64; 256]; 4],
    tails: [u64; 256],
    /// Lowest and highest score symbol of an N.
    n_scores: (u8, u8),
}

impl Counts {
    fn new() -> Counts {
        Counts {
            reads: 0,
            residues: 0,
            letters: [[0; 256]; 4],
            max_run: [0; 256],
            scores: [[0; 256]; 4],
            tails: [0; 256],
            n_scores: (u8::MAX, 0),
        }
    }

    /// adds a read; `tail` is the shortest run of a score symbol ending a
    /// read that counts as a tail.
    fn add(&mut self, seq: &[u8], qual: &[u8], tail: usize) {
        self.reads += 1;
        self.residues += seq.len() as u64;
        let (mut n_low, mut n_high) = self.n_scores;
        // the base of the current run, upper-cased, and its repeats
        let (mut base_run, mut repeats) = (0u8, 0u64);
        for (i, (&base, &score)) in seq.iter().zip(qual).enumerate() {
            // four tables in turn, so that consecutive increments rarely
            // wait on the same counter
            self.letters[i & 3][usize::from(base)] += 1;
            self.scores[i & 3][usize::from(score)] += 1;
            let base = base.to_ascii_uppercase();
            let is_n = base == b'N';
            n_low = n_low.min(if is_n { score } else { u8::MAX });
            n_high = n_high.max(if is_n { score } else { 0 });
            // as vsearch, a run's length is its number of repeats; the
            // longest is updated at every base, without a branch
            repeats = if i > 0 && base == base_run {
                repeats + 1
            } else {
                0
            };
            base_run = base;
            let longest = &mut self.max_run[usize::from(base)];
            *longest = (*longest).max(repeats);
        }
        self.n_scores = (n_low, n_high);
        if qual.len() >= tail
            && let Some(&last) = qual.last()
            && qual[qual.len() - tail..].iter().all(|&s| s == last)
        {
            self.tails[usize::from(last)] += 1;
        }
    }

    /// the report, as vsearch prints it.
    fn report(&self) -> String {
        let sum = |tables: &[[u64; 256]; 4]| -> [u64; 256] {
            std::array::from_fn(|c| tables.iter().map(|t| t[c]).sum())
        };
        let scores = sum(&self.scores);
        // residues are counted as read, and reported upper-cased
        let mut letters = sum(&self.letters);
        for c in b'a'..=b'z' {
            letters[usize::from(c.to_ascii_uppercase())] += letters[usize::from(c)];
            letters[usize::from(c)] = 0;
        }
        let mut out = format!("Read {} sequences.\n", self.reads);
        if self.reads == 0 {
            return out;
        }
        let lowest = present(&scores).next().unwrap_or(0);
        let highest = present(&scores).last().unwrap_or(0);
        let encoding = Encoding::guess(lowest, highest);
        let offset = i16::from(encoding.offset());
        let (lowest_i, highest_i) = (i16::from(lowest), i16::from(highest));
        writeln!(
            out,
            "Qmin {lowest_i}, Qmax {highest_i}, Range {}",
            highest_i - lowest_i + 1
        )
        .unwrap();
        writeln!(
            out,
            "Guess: -fastq_qmin {} -fastq_qmax {} -fastq_ascii {offset}",
            lowest_i - offset,
            highest_i - offset
        )
        .unwrap();
        out.push_str(match encoding {
            Encoding::Sanger => "Guess: Original Sanger format (phred+33)\n",
            Encoding::Illumina18 => "Guess: Illumina 1.8+ format (phred+33)\n",
            Encoding::Solexa => {
                "Guess: Solexa format (solexa+64, not supported;\n       \
                 convert with --fastq_convert --fastq_solexa)\n"
            }
            Encoding::Illumina13 => "Guess: Illumina 1.3+ format (phred+64)\n",
            Encoding::Illumina15 => "Guess: Illumina 1.5+ format (phred+64)\n",
        });
        out.push_str("\nLetter          N   Freq MaxRun\n------ ---------- ------ ------\n");
        let percent = |count: u64| 100.0 * count as f64 / self.residues as f64;
        for c in present(&letters) {
            let count = letters[usize::from(c)];
            write!(
                out,
                "     {} {count:>10} {:>5.1}% {:>6}",
                char::from(c),
                percent(count),
                self.max_run[usize::from(c)]
            )
            .unwrap();
            if c == b'N' {
                let (low, high) = self.n_scores;
                if low < high {
                    write!(out, "  Q={}..{}", char::from(low), char::from(high)).unwrap();
                } else {
                    write!(out, "  Q={}", char::from(low)).unwrap();
                }
            }
            out.push('\n');
        }
        out.push_str("\nChar  ASCII    Freq       Tails\n----  -----  ------  ----------\n");
        for c in present(&scores) {
            writeln!(
                out,
                " '{}'  {c:>5}  {:>5.1}%  {:>10}",
                char::from(c),
                percent(scores[usize::from(c)]),
                self.tails[usize::from(c)]
            )
            .unwrap();
        }
        out
    }
}

fn main() {
    let matches = Command::new("fqchar")
        .about("Report the residues and score symbols of FASTQ reads, and guess their encoding.")
        .long_about(
            "Report the residues and score symbols of FASTQ reads, and guess their encoding.\n\n\
             As vsearch --fastq_chars: the report gives the number of reads, the lowest \
             and highest score symbols, and the encoding and score range they suggest; \
             then for each residue, upper-cased, its count, share and longest run (less \
             one, as vsearch counts it), with the score symbols of N; then for each score \
             symbol its count, share and the number of reads ending in a run of at least \
             --tail of it.",
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
                .help("Output file of the report, default to stdout"),
        )
        .arg(
            Arg::new("tail")
                .long("tail")
                .value_name("LEN")
                .value_parser(value_parser!(u64).range(1..))
                .default_value("4")
                .help("Shortest run of a score symbol ending a read counted as a tail"),
        )
        .get_matches();
    let tail = matches.get_one::<u64>("tail").copied().unwrap_or(4) as usize;

    let fin: Box<dyn Read> = match matches.get_one::<String>("in").map(String::as_str) {
        Some(f) => {
            Box::new(File::open(f).unwrap_or_else(|e| fail(&format!("cannot open {f}: {e}"))))
        }
        None => Box::new(io::stdin()),
    };
    let mut reader = Reader::new(BufReader::new(fin));
    let format = reader.format().unwrap_or_else(|e| fail(&e.to_string()));
    if format == Some(Format::Fasta) {
        fail("FASTQ file expected, FASTA file found");
    }

    let mut counts = Counts::new();
    let mut record = Record::new();
    loop {
        if let Err(e) = reader.read(&mut record) {
            fail(&e.to_string());
        }
        if record.is_empty() {
            break;
        }
        counts.add(record.seq(), record.qual(), tail);
    }
    if let Some(warning) = reader.stripped_warning() {
        warn(&warning);
    }

    let mut fout: Box<dyn Write> = match matches.get_one::<String>("out").map(String::as_str) {
        Some(f) => {
            Box::new(File::create(f).unwrap_or_else(|e| fail(&format!("cannot create {f}: {e}"))))
        }
        None => Box::new(io::stdout()),
    };
    fout.write_all(counts.report().as_bytes())
        .and_then(|()| fout.flush())
        .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
}
