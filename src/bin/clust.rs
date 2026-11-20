use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use clap::Command;

use rsearch::clustering::{cluster_args, run};

fn main() {
    let command = Command::new("clust")
        .about("Cluster sequences around centroids, in input order; write centroids and hits.")
        .long_about(
            "Cluster sequences around centroids, in input order; write centroids and \
             hits.\n\n\
             As vsearch --cluster_smallmem --usersort: each sequence is searched against \
             the centroids so far, as by dbsrch (lower case masked), and joins the cluster \
             of its best accepted hit, or becomes a centroid. The centroids are written \
             with the abundance of their cluster (;size=), in the input format; the hits, \
             with --hits, in UC format: an S record for each centroid and an H record for \
             each member, in input order, then a C record for each cluster. Sorted first \
             by sortby, the input gives --cluster_fast and --cluster_size.",
        );
    let matches = cluster_args(command)
        .mut_arg("id", |a| a.required(true))
        .get_matches();
    run(&matches, None);
}
