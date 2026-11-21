//! Benchmarks of selecting records by header.
//!
//! `matches` measures `Selector::matches` on the headers of parsed reads;
//! `binary` runs the getseq program end to end on a FASTQ file. It selects
//! by one label, 10,000 labels from a file, a substring, a word, and words
//! of a field from a file. `hits` runs it on the hits exsrch finds for
//! 100,000 reads of 150 bases against 25,000 sequences. It selects the
//! queries with hits, those without, and the targets hit.
//! compare/benches/getseq times vsearch on the same files.

use std::hint::black_box;
use std::path::Path;
use std::process::Command;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

use rsearch::fastx::{Reader, Record};
use rsearch::select::Selector;

mod input;
use input::*;

fn parse(text: &str) -> Vec<Record> {
    Reader::new(text.as_bytes())
        .records()
        .map(Result::unwrap)
        .collect()
}

fn bench_matches(c: &mut Criterion) {
    if skip() {
        return;
    }
    let reads = parse(&text());
    let mut group = c.benchmark_group(data::group("getseq/matches"));
    group.throughput(Throughput::Elements(reads.len() as u64));
    let labels = labels();
    for (name, selector) in [
        ("labels", Selector::labels(&labels, false)),
        ("substring", Selector::labels(["SAMPLE=S7;"], true)),
        ("word", Selector::words(["s7"], None)),
        ("field", Selector::words(["s7", "s42"], Some("sample"))),
    ] {
        group.bench_function(name, |b| {
            b.iter(|| {
                for r in &reads {
                    black_box(selector.matches(r.header()));
                }
            })
        });
    }
    group.finish();
}

/// A benchmark case: getseq's arguments.
struct Case<'a> {
    name: &'static str,
    getseq: Vec<&'a str>,
}

fn cases(inputs: &Inputs) -> Vec<Case<'_>> {
    let labels = inputs.labels.to_str().unwrap();
    let words = inputs.words.to_str().unwrap();
    let label = "read500;sample=s0;size=1";
    vec![
        Case {
            name: "label",
            getseq: vec!["--label", label],
        },
        Case {
            name: "labels",
            getseq: vec!["--labels", labels],
        },
        Case {
            name: "substring",
            getseq: vec!["--label", "SAMPLE=S7;", "--substr"],
        },
        Case {
            name: "word",
            getseq: vec!["--word", "s7"],
        },
        Case {
            name: "field_words",
            getseq: vec!["--words", words, "--field", "sample"],
        },
    ]
}

fn run(program: &str, args: &[&str]) {
    let status = Command::new(program).args(args).status().unwrap();
    assert!(status.success(), "{program} {args:?} failed");
}

fn bench_binary(c: &mut Criterion) {
    if skip() {
        return;
    }
    let inputs = input_files();
    let getseq = env!("CARGO_BIN_EXE_getseq");
    let fq = inputs.fastq.to_str().unwrap();

    let mut group = c.benchmark_group(data::group("getseq/binary"));
    group.throughput(Throughput::Elements(data::count(&inputs.fastq, READS)));
    group.sample_size(10);
    for case in cases(&inputs) {
        let args = [&["-i", fq, "-o", "/dev/null"][..], &case.getseq].concat();
        group.bench_function(format!("getseq/{}", case.name), |b| {
            b.iter(|| run(getseq, &args))
        });
    }
    group.finish();
}

fn bench_hits(c: &mut Criterion) {
    if skip() {
        return;
    }
    let (fa, db) = hit_files();
    let (fa, db) = (fa.to_str().unwrap(), db.to_str().unwrap());
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let hits = data::output(dir.join("getseq-hits.uc"));
    let hits = hits.to_str().unwrap();
    run(
        env!("CARGO_BIN_EXE_exsrch"),
        &["-i", fa, "--db", db, "-o", hits],
    );
    let getseq = env!("CARGO_BIN_EXE_getseq");
    let mut group = c.benchmark_group(data::group("getseq/hits"));
    group.throughput(Throughput::Elements(data::count(fa, 100_000)));
    group.sample_size(10);
    for (name, input, args) in [
        ("matched", fa, vec![]),
        ("notmatched", fa, vec!["--invert"]),
        ("dbmatched", db, vec!["--targets"]),
    ] {
        let args = [&["--hits", hits, "-i", input, "-o", "/dev/null"][..], &args].concat();
        group.bench_function(format!("getseq/{name}"), |b| b.iter(|| run(getseq, &args)));
    }
    group.finish();
}

criterion_group!(benches, bench_matches, bench_binary, bench_hits);
criterion_main!(benches);
