use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fmt::Write as _;
use std::fs::File;
use std::io;
use std::io::{BufReader, Read, Write};

use clap::{Arg, Command};

use rsearch::udb::Udb;

/// prints an error and exits.
fn fail(message: &str) -> ! {
    eprintln!("Error: {message}");
    std::process::exit(1);
}

/// The words the report lists, and how many of their sequences it shows.
const ROWS: usize = 11;
const ENTRIES: usize = 8;

/// A word the report lists: its code, how many sequences hold it, and the
/// first of them.
struct Row {
    word: u32,
    count: u32,
    /// The number of its first entry in the word list.
    start: u64,
    entries: Vec<u32>,
}

/// the `ROWS` words held by the most sequences, the smaller code first
/// among equal counts.
fn top_rows(counts: &[u32]) -> Vec<Row> {
    let mut top: Vec<Row> = Vec::with_capacity(ROWS + 1);
    for (word, &count) in counts.iter().enumerate() {
        if top.len() < ROWS || count > top[top.len() - 1].count {
            // after every row of a higher or equal count
            let at = top.partition_point(|row| row.count >= count);
            top.insert(
                at,
                Row {
                    word: word as u32,
                    count,
                    start: 0,
                    entries: Vec::new(),
                },
            );
            top.truncate(ROWS);
        }
    }
    // where each row's sequences start in the word list
    let mut order: Vec<usize> = (0..top.len()).collect();
    order.sort_by_key(|&i| top[i].word);
    let (mut start, mut word) = (0u64, 0usize);
    for i in order {
        let until = top[i].word as usize;
        start += counts[word..until]
            .iter()
            .map(|&c| u64::from(c))
            .sum::<u64>();
        word = until;
        top[i].start = start;
    }
    top
}

/// the residues of `word`, of `word_len` residues.
fn residues(word_len: u32, word: u32) -> String {
    (0..word_len)
        .map(|i| b"ACGT"[((word >> (2 * (word_len - i - 1))) & 3) as usize] as char)
        .collect()
}

/// a number in a column of 12, in thousands from 10,000.
fn column(value: f64) -> String {
    if value >= 10000.0 {
        format!("  {:9.1}k", value * 0.001)
    } else {
        format!("  {value:10.1}")
    }
}

/// vsearch's `--udbstats` report.
fn report(udb: &Udb, top: &[Row]) -> String {
    let word_len = udb.header.word_len();
    let slots = udb.counts.len() as u64;
    let sequences = u64::from(udb.header.sequences());
    let residues_total = udb.sequences.len() as u64;
    let words = udb.words;
    let mut out = String::new();
    let shortest = udb.lengths.iter().min().unwrap();
    let longest = udb.lengths.iter().max().unwrap();
    let _ = writeln!(
        out,
        "{residues_total} nt in {sequences} seqs, min {shortest}, max {longest}, avg {:.0}\n",
        residues_total as f64 / sequences as f64
    );

    // the number of words held by each number of sequences
    let mut histogram = vec![0u64; sequences as usize + 1];
    for &count in &udb.counts {
        histogram[count as usize] += 1;
    }
    let count_at_rank = |rank: u64| -> u64 {
        let mut seen = 0;
        for (count, &n) in histogram.iter().enumerate() {
            seen += n;
            if seen >= rank {
                return count as u64;
            }
        }
        unreachable!("the ranks asked for are below the number of words")
    };
    let median = (count_at_rank(slots / 2) + count_at_rank(slots / 2 + 1)) / 2;

    let _ = write!(
        out,
        "      Alphabet  nt\n    Word width  {word_len}\n     Word ones  {word_len}\n        \
         Spaced  No\n        Hashed  No\n         Coded  No\n       Stepped  No\n         \
         Slots  {slots} ({:.1}k)\n       DBAccel  {}%\n\n",
        slots as f64 / 1000.0,
        udb.header.db_accel()
    );
    let _ = write!(
        out,
        "{residues_total:10}  DB size ({:.1}k)\n{words:10}  Words\n{median:10}  Median \
         size\n{:10.1}  Mean size\n\n",
        residues_total as f64 / 1000.0,
        words as f64 / slots as f64
    );
    out.push_str("     iWord         sWord         Cap        Size  Row\n");
    out.push_str("----------  ------------  ----------  ----------  ---\n");
    for row in top {
        let padding = " ".repeat(12usize.saturating_sub(word_len as usize));
        let _ = write!(
            out,
            "{:10}  {padding}{}  {:10}  {:10} ",
            row.word,
            residues(word_len, row.word),
            0,
            row.count
        );
        for entry in &row.entries {
            let _ = write!(out, " {entry}");
        }
        if row.count as usize > ENTRIES {
            out.push_str("...");
        }
        out.push('\n');
    }
    let max = &top[0];
    let _ = write!(
        out,
        "\n\nWord width  {word_len}\nSlots       {slots}\nWords       {words}\nMax size    {} \
         ({})\n\n",
        max.count,
        residues(word_len, max.word)
    );
    out.push_str("   Size lo     Size hi  Total size   Nr. Words     Pct  TotPct\n");
    out.push_str("----------  ----------  ----------  ----------  ------  ------\n");
    // buckets of counts doubling in width, up to the number of sequences,
    // which stops one bucket short of it, as vsearch's
    let (mut low, mut high, mut next) = (0u64, 0u64, 0u64);
    let mut total_pct = 0.0;
    while low < sequences {
        let (mut count, mut size) = (0u64, 0u64);
        for value in next..=high {
            count += histogram[value as usize];
            size += value * histogram[value as usize];
        }
        next = high + 1;
        let pct = 100.0 * count as f64 / slots as f64;
        total_pct += pct;
        if low < high {
            let _ = write!(out, "{low:10}");
        } else {
            out.push_str("          ");
        }
        let _ = write!(
            out,
            "  {high:10}{}{}  {pct:5.1}%  {total_pct:5.1}%",
            column(size as f64),
            column(count as f64)
        );
        let dots = (pct / 3.0).round() as usize;
        if dots > 0 {
            out.push_str("  ");
            out.push_str(&"*".repeat(dots));
        }
        out.push('\n');
        low = high + 1;
        high = if high > 0 { high * 2 } else { 1 }.min(sequences);
    }
    out.push_str("----------  ----------  ----------  ----------\n");
    let _ = write!(
        out,
        "                      {}{}\n\n",
        column(words as f64),
        column(slots as f64)
    );
    let _ = write!(
        out,
        "{residues_total:10}  Upper\n{:10}  Lower (0.0%)\n{residues_total:10}  \
         Total\n{words:10}  Indexed words\n",
        0
    );
    out
}

fn main() {
    let matches = Command::new("udbstat")
        .about("Report a UDB file's word statistics.")
        .long_about(
            "Report a UDB file's word statistics.\n\n\
             Writes vsearch's --udbstats report: the size of the database, how many \
             sequences hold each word, the words held by the most sequences, and the \
             distribution of the words over these counts. The whole file is read and \
             checked as vsearch does.",
        )
        .arg(
            Arg::new("in")
                .short('i')
                .long("in")
                .value_name("INPUT")
                .help("Input UDB file, default to stdin"),
        )
        .arg(
            Arg::new("out")
                .short('o')
                .long("out")
                .value_name("OUTPUT")
                .help("Output report, default to stdout"),
        )
        .get_matches();

    // a file is read as it goes, knowing its size; stdin is read whole
    let mut top: Option<Vec<Row>> = None;
    let mut entries = |counts: &[u32], first: u64, block: &[u32]| {
        let rows = top.get_or_insert_with(|| top_rows(counts));
        rows_entries(rows, first, block);
    };
    let udb = match matches.get_one::<String>("in").map(String::as_str) {
        Some(f) => {
            let file = File::open(f).unwrap_or_else(|e| fail(&format!("cannot open {f}: {e}")));
            let size = file
                .metadata()
                .unwrap_or_else(|e| fail(&format!("cannot read {f}: {e}")))
                .len();
            Udb::read(&mut BufReader::new(file), size, &mut entries)
        }
        None => {
            let mut bytes = Vec::new();
            io::stdin()
                .lock()
                .read_to_end(&mut bytes)
                .unwrap_or_else(|e| fail(&format!("cannot read: {e}")));
            Udb::read(&mut &bytes[..], bytes.len() as u64, &mut entries)
        }
    }
    .unwrap_or_else(|e| fail(&e.to_string()));
    let top = top.unwrap_or_else(|| top_rows(&udb.counts));

    let mut fout: Box<dyn Write> = match matches.get_one::<String>("out").map(String::as_str) {
        Some(f) => {
            Box::new(File::create(f).unwrap_or_else(|e| fail(&format!("cannot create {f}: {e}"))))
        }
        None => Box::new(io::stdout()),
    };
    fout.write_all(report(&udb, &top).as_bytes())
        .and_then(|()| fout.flush())
        .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
}

/// keeps, for each row, its first entries found in `block`, the entries
/// from number `first` of the word list.
fn rows_entries(rows: &mut [Row], first: u64, block: &[u32]) {
    let end = first + block.len() as u64;
    for row in rows.iter_mut() {
        let wanted = (row.count as usize).min(ENTRIES);
        if row.entries.len() == wanted {
            continue;
        }
        let start = row.start + row.entries.len() as u64;
        let stop = start + (wanted - row.entries.len()) as u64;
        if start < end && stop > first {
            let from = start.max(first);
            let to = stop.min(end);
            row.entries
                .extend_from_slice(&block[(from - first) as usize..(to - first) as usize]);
        }
    }
}
