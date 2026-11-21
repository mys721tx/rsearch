//! Times `vsearch --cluster_smallmem --usersort --clusters --threads 1` and
//! its rsearch pipeline (clust then clspl) on the file of rsearch's clspl
//! benchmark. It first checks that clspl writes the same files from clust's
//! clusters (see `common/vsearch.rs`).

use std::path::Path;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/clspl/input.rs"]
mod input;
#[path = "../common/vsearch.rs"]
mod vsearch;
use input::*;

/// The files written with `prefix`, in cluster order.
fn files(prefix: &Path) -> Vec<Vec<u8>> {
    (0..)
        .map_while(|i| std::fs::read(format!("{}{i}", prefix.display())).ok())
        .collect()
}

fn bench_binary(c: &mut Criterion) {
    let Some(vsearch) = vsearch::vsearch() else {
        return;
    };
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let input = input_file(dir);
    let fa = input.to_str().unwrap();
    let hits = dir.join("clspl-bench.uc");
    let hits = hits.to_str().unwrap();
    vsearch::run(
        &vsearch::tool("clust"),
        &["-i", fa, "-o", "/dev/null", "--hits", hits, "--id", "0.97"],
    );
    let (rdir, vdir) = (dir.join("clspl-bench-r"), dir.join("clspl-bench-v"));
    for d in [&rdir, &vdir] {
        let _ = std::fs::remove_dir_all(d);
        std::fs::create_dir_all(d).unwrap();
    }
    let (rprefix, vprefix) = (rdir.join("c"), vdir.join("c"));
    let (rp, vp) = (rprefix.to_str().unwrap(), vprefix.to_str().unwrap());
    let vargs = [
        "--cluster_smallmem",
        fa,
        "--usersort",
        "--id",
        "0.97",
        "--clusters",
        vp,
        "--quiet",
        "--notrunclabels",
        "--qmask",
        "soft",
        "--fasta_width",
        "0",
        "--threads",
        "1",
    ];
    // both write the same files
    vsearch::run(
        &vsearch::tool("clspl"),
        &["--hits", hits, "-i", fa, "--prefix", rp],
    );
    vsearch::run(&vsearch, &vargs);
    assert!(
        files(&rprefix) == files(&vprefix),
        "clspl and vsearch wrote different files"
    );
    let mut group = c.benchmark_group("clspl/binary");
    group.throughput(Throughput::Elements(SEQUENCES as u64));
    group.sample_size(10);
    vsearch::bench(&mut group, "vsearch/cluster", &vsearch, &vargs);
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
