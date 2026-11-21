//! Times `vsearch --chimeras_denovo --threads 1` and its rsearch pipeline
//! on the file of rsearch's mchime benchmark, after checking that mchime
//! writes the same table (see `common/vsearch.rs`).

use std::path::Path;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/mchime/input.rs"]
mod input;
#[path = "../common/vsearch.rs"]
mod vsearch;
use input::*;

fn bench_binary(c: &mut Criterion) {
    if skip() {
        return;
    }
    let Some(vsearch) = vsearch::vsearch() else {
        return;
    };
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let (all, _, _) = input_files(dir);
    let all = all.to_str().unwrap();
    let mut group = c.benchmark_group(data::group("mchime/binary"));
    group.throughput(Throughput::Elements(data::count(all, SEQUENCES)));
    group.sample_size(10);
    let common = [
        "--chimeras_denovo",
        all,
        "--quiet",
        "--notrunclabels",
        "--qmask",
        "soft",
        "--fasta_width",
        "0",
        "--threads",
        "1",
    ];
    // both write the same tables
    let (r, v) = (
        data::output(dir.join("mchime-bench-r.tab")),
        data::output(dir.join("mchime-bench-v.tab")),
    );
    let (r, v) = (r.to_str().unwrap(), v.to_str().unwrap());
    // vsearch sorts the input by abundance, ties by label, as the
    // translated pipeline does with sortby. The check does the same, which
    // matters for real data, where abundances tie
    let sorted = data::output(dir.join("mchime-bench-sorted.fa"));
    let sorted = sorted.to_str().unwrap();
    vsearch::run(
        &vsearch::tool("sortby"),
        &["--size", "-i", all, "-o", sorted],
    );
    vsearch::run(
        &vsearch::tool("mchime"),
        &["-i", sorted, "-o", "/dev/null", "--table", r],
    );
    vsearch::run(&vsearch, &[&common[..], &["--tabbedout", v]].concat());
    assert!(
        std::fs::read(r).unwrap() == std::fs::read(v).unwrap(),
        "mchime and vsearch wrote different tables"
    );
    let args = [&common[..], &["--nonchimeras", "/dev/null"]].concat();
    vsearch::bench(&mut group, "vsearch/denovo", &vsearch, &args);
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
