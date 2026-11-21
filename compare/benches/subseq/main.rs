//! Times `vsearch --fastx_getsubseq`, selecting every record, and its
//! rsearch pipelines on the files of rsearch's subseq benchmark. The cases
//! are a region in the middle, and from a start to the end (see
//! `common/vsearch.rs`).

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/subseq/input.rs"]
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
    let inputs = input_files();
    let fq = inputs.fastq.to_str().unwrap();
    let fa = inputs.fasta.to_str().unwrap();

    let mut group = c.benchmark_group(data::group("subseq/binary"));
    group.throughput(Throughput::Elements(data::count(&inputs.fastq, READS)));
    group.sample_size(10);
    let regions: [(&str, &[&str]); 2] = [
        ("middle", &["--subseq_start", "50", "--subseq_end", "200"]),
        ("to_end", &["--subseq_start", "20"]),
    ];
    for (format, file, output) in [
        ("fastq", fq, &["--fastqout", "/dev/null"][..]),
        (
            "fasta",
            fa,
            &["--fastaout", "/dev/null", "--fasta_width", "0"],
        ),
    ] {
        for (region, region_args) in regions {
            let args = [
                &[
                    "--fastx_getsubseq",
                    file,
                    "--quiet",
                    "--notrunclabels",
                    "--label_substr_match",
                    "--label",
                    "",
                ][..],
                output,
                region_args,
            ]
            .concat();
            vsearch::bench(
                &mut group,
                format!("vsearch/{region}_{format}"),
                &vsearch,
                &args,
            );
        }
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
