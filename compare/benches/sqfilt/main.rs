//! Times `vsearch --fastx_filter` (or, for the masked share, `--fastx_mask
//! --qmask soft`) with sqfilt's filters and its rsearch pipelines on the
//! files of rsearch's sqfilt benchmark. It writes the discarded reads too
//! (see `common/vsearch.rs`).

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/sqfilt/input.rs"]
mod input;
#[path = "../common/vsearch.rs"]
mod vsearch;
use input::*;

/// The vsearch commands, by name.
fn cases(inputs: &Inputs) -> Vec<(&'static str, Vec<&str>)> {
    let fq = inputs.fastq.to_str().unwrap();
    let rev = inputs.reverse.to_str().unwrap();
    let fa = inputs.fasta.to_str().unwrap();
    let null = "/dev/null";
    let case = |name, file, options: &[&'static str]| {
        let output: &[&str] = if file == fa {
            &[
                "--fastaout",
                null,
                "--fastaout_discarded",
                null,
                "--fasta_width",
                "0",
            ]
        } else {
            &["--fastqout", null, "--fastqout_discarded", null]
        };
        let args = [
            &["--fastx_filter", file, "--quiet", "--notrunclabels"],
            output,
            options,
        ]
        .concat();
        (name, args)
    };
    vec![
        case("none", fq, &[]),
        case("min_len", fq, &["--fastq_minlen", "240"]),
        case("max_ns", fq, &["--fastq_maxns", "10"]),
        case("max_ee", fq, &["--fastq_maxee", "1"]),
        case("max_ee_rate", fq, &["--fastq_maxee_rate", "0.005"]),
        case("min_qual", fq, &["--fastq_minqual", "10"]),
        case("size", fq, &["--minsize", "2", "--maxsize", "4"]),
        case(
            "all",
            fq,
            &[
                "--fastq_minlen",
                "100",
                "--fastq_maxns",
                "20",
                "--fastq_maxee",
                "3",
                "--fastq_minqual",
                "5",
                "--minsize",
                "2",
            ],
        ),
        (
            "pairs",
            vec![
                "--fastx_filter",
                fq,
                "--reverse",
                rev,
                "--quiet",
                "--notrunclabels",
                "--fastqout",
                null,
                "--fastqout_rev",
                null,
                "--fastqout_discarded",
                null,
                "--fastqout_discarded_rev",
                null,
                "--fastq_maxee",
                "1",
            ],
        ),
        (
            "upper_pct",
            vec![
                "--fastx_mask",
                fq,
                "--quiet",
                "--notrunclabels",
                "--qmask",
                "soft",
                "--fastqout",
                null,
                "--min_unmasked_pct",
                "90",
            ],
        ),
        case(
            "fasta",
            fa,
            &["--fastq_minlen", "100", "--fastq_maxns", "10"],
        ),
    ]
}

fn bench_binary(c: &mut Criterion) {
    let Some(vsearch) = vsearch::vsearch() else {
        return;
    };
    let inputs = input_files();

    let mut group = c.benchmark_group("sqfilt/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    for (name, args) in cases(&inputs) {
        vsearch::bench(&mut group, format!("vsearch/{name}"), &vsearch, &args);
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
