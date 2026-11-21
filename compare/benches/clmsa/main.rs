//! Times `vsearch --cluster_smallmem --usersort --msaout --consout
//! --profile --threads 1` and its rsearch pipeline (clust then clmsa three
//! times) on the file of rsearch's clmsa benchmark. It first checks that
//! clmsa writes the same from clust's clusters (see `common/vsearch.rs`).

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/clmsa/input.rs"]
mod input;
#[path = "../common/vsearch.rs"]
mod vsearch;
use input::*;

fn output(program: &str, args: &[&str]) -> Vec<u8> {
    let output = Command::new(program)
        .args(args)
        .stderr(Stdio::null())
        .output()
        .unwrap();
    assert!(output.status.success(), "{program} {args:?} failed");
    output.stdout
}

fn bench_binary(c: &mut Criterion) {
    let Some(vsearch) = vsearch::vsearch() else {
        return;
    };
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let input = input_file(dir);
    let fa = input.to_str().unwrap();
    let hits = dir.join("clmsa-bench.uc");
    let hits = hits.to_str().unwrap();
    vsearch::run(
        &vsearch::tool("clust"),
        &["-i", fa, "-o", "/dev/null", "--hits", hits, "--id", "0.9"],
    );
    let path = |name: &str| dir.join(format!("clmsa-bench-{name}"));
    let (msa, cons, prof) = (path("v.msa"), path("v.cons"), path("v.prof"));
    let s = |p: &PathBuf| p.to_str().unwrap().to_string();
    let (msa, cons, prof) = (s(&msa), s(&cons), s(&prof));
    let vargs = [
        "--cluster_smallmem",
        fa,
        "--usersort",
        "--id",
        "0.9",
        "--msaout",
        &msa,
        "--consout",
        &cons,
        "--profile",
        &prof,
        "--sizein",
        "--sizeout",
        "--quiet",
        "--notrunclabels",
        "--qmask",
        "soft",
        "--fasta_width",
        "0",
        "--threads",
        "1",
    ];
    // both write the same
    vsearch::run(&vsearch, &vargs);
    let formats: [&[&str]; 3] = [&[], &["--consensus"], &["--profile"]];
    let clmsa = vsearch::tool("clmsa");
    for (format, written) in formats.iter().zip([&msa, &cons, &prof]) {
        let args = [&["--hits", hits, "-i", fa][..], format].concat();
        assert!(
            output(&clmsa, &args) == std::fs::read(written).unwrap(),
            "clmsa and vsearch wrote different alignments"
        );
    }
    let mut group = c.benchmark_group("clmsa/binary");
    group.throughput(Throughput::Elements(SEQUENCES as u64));
    group.sample_size(10);
    vsearch::bench(&mut group, "vsearch/all", &vsearch, &vargs);
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
