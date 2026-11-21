# unoise

Denoise amplicons with the UNOISE3 rule; write centroids and hits.

## Synopsis

```
unoise [-i FILE] [-o FILE] [--hits FILE] [--alpha ALPHA] [--id SHARE]
       [--strand plus|both] [--maxaccepts N] [--maxrejects N]
       [--word-len N] [--min-word-matches N] [--size-order] [--no-sizein]
       [--no-sizeout] [criteria] [scores]
```

## Description

unoise is [clust](clust.md) with UNOISE's acceptance rule, for input sorted
by decreasing abundance (`sortby --size`). Each sequence is searched
against the centroids so far. It joins the cluster of its best hit of at
least `--id` identity whose centroid is at least 2^(alpha d + 1) times as
abundant as it, for d mismatches, or has no mismatch. Otherwise it becomes
a centroid. The centroids and hits are written as by clust.

## Options

- `--alpha ALPHA`: UNOISE's alpha: how much more abundant a centroid must
  be for each mismatch. Default: 2.
- `-i`, `--in FILE`: Input FASTA or FASTQ file. Default: standard input.
- `-o`, `--out FILE`: Output file of the centroids, in the input format.
  Default: standard output.
- `--hits FILE`: Output UC file of the clusters.
- `--strand plus|both`: Search the plus strand of the sequences, or both.
  Default: plus.
- `--maxaccepts N`: Stop after N accepted centroids on a strand, 0 for no
  limit. Default: 1.
- `--maxrejects N`: Stop after N rejected centroids on a strand, 0 for no
  limit. Default: 32.
- `--word-len N`: Length of the words indexed, 3 to 15. Default: 8.
- `--min-word-matches N`: Words a centroid must share with the sequence to
  be aligned. Default: taken from `--word-len`.
- `--size-order`: Join the most abundant of the accepted centroids, rather
  than the closest.
- `--no-sizein`: Count each sequence once, rather than by its `;size=`.
- `--no-sizeout`: Do not append the cluster's abundance to the centroids.
- `--id SHARE`: Lowest identity of a hit, 0 to 1. Default: 0.9.
- `--weak-id SHARE`: Refused. unoise has no weak hits.
- `--iddef 0..4`, `--max-diffs N`, `--max-subs N`, `--max-gaps N`,
  `--min-cols N`, `--left-just`, `--right-just`, `--query-cov SHARE`,
  `--target-cov SHARE`, `--max-id SHARE`, `--mid PCT`, `--id-prefix N`,
  `--id-suffix N`, `--max-qt RATIO`, `--min-qt RATIO`, `--min-sl SHARE`,
  `--max-sl SHARE`, `--max-qsize N`, `--min-tsize N`,
  `--min-size-ratio RATIO`, `--max-size-ratio RATIO`, `--self`, `--self-id`:
  The acceptance criteria of [dbsrch](dbsrch.md).
- `--match SCORE`, `--mismatch SCORE`, `--gapopen PENALTIES`,
  `--gapext PENALTIES`, `--n-mismatch`: The scores of [dbsrch](dbsrch.md).
  Defaults: 2, -4, 20I/2E, 2I/1E.

## vsearch

unoise is `--cluster_unoise` on input sorted by `sortby --size`, with
`--unoise_alpha` as `--alpha` and its minimum identity, vsearch's
`--weak_id` there, as `--id` (vsearch ignores `--id` for this command).
vsearch's defaults `--minsize 8`, `--minseqlength 32` and `--maxseqlength
50000` are `sqfilt --min-size 8 --min-len 32 --max-len 50000` before it.

## Examples

```console
$ sortby --size -i amplicons.fa | unoise --hits denoised.uc | grep '>'
>a1;size=62
>a2;size=40
>a3;size=25
>a5;size=4
$ grep -v '^C' denoised.uc | cut -f 1-4,9-10
S	0	160	*	a1;size=60	*
S	1	160	*	a2;size=40	*
S	2	160	*	a3;size=25	*
S	3	160	*	a5;size=4	*
H	0	160	99.4	a4;size=2	a1;size=60
```

vsearch's `--cluster_unoise amplicons.fa --minsize 2 --sizein --sizeout
--centroids c.fa`:

```console
$ sqfilt --min-len 32 --max-len 50000 --min-size 2 -i amplicons.fa | lcmask | sortby --size | unoise | conv --fasta | grep '>'
>a1;size=62
>a2;size=40
>a3;size=25
>a5;size=4
```

## Exit status

0 on success, 1 for invalid input, 2 for invalid arguments.

## See also

[clust](clust.md), [sortby](sortby.md), [uchime](uchime.md),
[otutab](otutab.md)
