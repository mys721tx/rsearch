//! Times `vsearch --fastx_getseq` and `--fastx_getseqs` and their rsearch
//! pipelines on the FASTQ file of rsearch's getseq benchmark. It also times
//! `vsearch --search_exact --threads 1 --qmask none --dbmask none` writing
//! `--matched`, `--notmatched` and `--dbmatched` on its hit files, its
//! exact search included. It first checks that getseq, on exsrch's hits,
//! selects the same queries (see `common/vsearch.rs`).

use std::path::Path;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/getseq/input.rs"]
mod input;
#[path = "../common/vsearch.rs"]
mod vsearch;
use input::*;

/// The vsearch commands and their options, by name.
fn cases(inputs: &Inputs) -> Vec<(&'static str, &'static str, Vec<&str>)> {
    let labels = inputs.labels.to_str().unwrap();
    let words = inputs.words.to_str().unwrap();
    let label = "read500;sample=s0;size=1";
    vec![
        ("label", "--fastx_getseq", vec!["--label", label]),
        ("labels", "--fastx_getseqs", vec!["--labels", labels]),
        (
            "substring",
            "--fastx_getseq",
            vec!["--label", "SAMPLE=S7;", "--label_substr_match"],
        ),
        ("word", "--fastx_getseqs", vec!["--label_word", "s7"]),
        (
            "field_words",
            "--fastx_getseqs",
            vec!["--label_words", words, "--label_field", "sample"],
        ),
    ]
}

fn bench_binary(c: &mut Criterion) {
    if skip() {
        return;
    }
    let Some(vsearch) = vsearch::vsearch() else {
        return;
    };
    let inputs = input_files();
    let fq = inputs.fastq.to_str().unwrap();

    let mut group = c.benchmark_group(data::group("getseq/binary"));
    group.throughput(Throughput::Elements(data::count(&inputs.fastq, READS)));
    group.sample_size(10);
    for (name, command, options) in cases(&inputs) {
        let args = [
            &[
                command,
                fq,
                "--quiet",
                "--notrunclabels",
                "--fastqout",
                "/dev/null",
            ][..],
            &options,
        ]
        .concat();
        vsearch::bench(&mut group, format!("vsearch/{name}"), &vsearch, &args);
    }
    group.finish();
}

fn bench_hits(c: &mut Criterion) {
    if skip() {
        return;
    }
    let Some(vsearch) = vsearch::vsearch() else {
        return;
    };
    let (fa, db) = hit_files();
    let (fa, db) = (fa.to_str().unwrap(), db.to_str().unwrap());
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let hits = data::output(dir.join("getseq-hits.uc"));
    let hits = hits.to_str().unwrap();
    vsearch::run(
        &vsearch::tool("exsrch"),
        &["-i", fa, "--db", db, "-o", hits],
    );
    let vbase = [
        "--search_exact",
        fa,
        "--db",
        db,
        "--quiet",
        "--notrunclabels",
        "--threads",
        "1",
        "--fasta_width",
        "0",
        "--qmask",
        "none",
        "--dbmask",
        "none",
    ];
    // both select the same queries
    let (r, v) = (
        data::output(dir.join("getseq-hits-r.fa")),
        data::output(dir.join("getseq-hits-v.fa")),
    );
    vsearch::run(
        &vsearch::tool("getseq"),
        &["--hits", hits, "-i", fa, "-o", r.to_str().unwrap()],
    );
    vsearch::run(
        &vsearch,
        &[&vbase[..], &["--matched", v.to_str().unwrap()]].concat(),
    );
    assert!(
        std::fs::read(&r).unwrap() == std::fs::read(&v).unwrap(),
        "getseq and vsearch selected different queries"
    );
    let mut group = c.benchmark_group(data::group("getseq/hits"));
    group.throughput(Throughput::Elements(data::count(fa, 100_000)));
    group.sample_size(10);
    for (name, output) in [
        ("matched", vec!["--matched", "/dev/null"]),
        ("notmatched", vec!["--notmatched", "/dev/null"]),
        ("dbmatched", vec!["--dbmatched", "/dev/null"]),
    ] {
        let vargs = [&vbase[..], &output].concat();
        vsearch::bench(&mut group, format!("vsearch/{name}"), &vsearch, &vargs);
    }
    group.finish();
}

criterion_group!(benches, bench_binary, bench_hits);
criterion_main!(benches);
