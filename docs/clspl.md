# clspl

Split records into one file per cluster, from clustering hits.

## Synopsis

```
clspl --hits FILE --prefix PREFIX [-i FILE]
```

## Description

clspl takes the hits of [clust](clust.md) or [unoise](unoise.md) and the
input they clustered, in the order they read it. The S and H records give
the cluster of each record, matched in order by header. Each cluster's
records are written unchanged, in input order and in the input format, to
PREFIX followed by the cluster's number from 0. A record without its hit,
or a hit without its record, is an error. No file is written before every
record is matched.

## Options

- `--hits FILE`: UC file of the clusters. Required.
- `-i`, `--in FILE`: Input FASTA or FASTQ file that was clustered.
  Default: standard input.
- `--prefix PREFIX`: Prefix of the files written, followed by each
  cluster's number. Required.

## vsearch

clspl is the `--clusters` output of vsearch's clustering commands. vsearch
rewrites the abundances of the records with `--sizeout`; clspl writes them
as they are, and wraps nothing.

## Examples

```console
$ sortby --size -i amplicons.fa > sorted.fa
$ clust -i sorted.fa --id 0.9 --hits clusters.uc > centroids.fa
$ clspl --hits clusters.uc -i sorted.fa --prefix cluster
$ grep '>' cluster0 cluster1
cluster0:>a1;size=60
cluster0:>a2;size=40
cluster0:>a5;size=4
cluster0:>a4;size=2
cluster1:>a3;size=25
```

## Exit status

0 on success, 1 for invalid input, including records and hits that do not
match, 2 for invalid arguments.

## See also

[clust](clust.md), [unoise](unoise.md), [clmsa](clmsa.md),
[getseq](getseq.md)
