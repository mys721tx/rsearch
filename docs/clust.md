# clust

Cluster sequences around centroids, in input order; write centroids and
hits.

## Synopsis

```
clust --id SHARE [-i FILE] [-o FILE] [--hits FILE] [--strand plus|both]
      [--maxaccepts N] [--maxrejects N] [--word-len N]
      [--min-word-matches N] [--size-order] [--no-sizein] [--no-sizeout]
      [criteria] [scores] [-t N]
```

## Description

clust clusters greedily, in the order of its input, as vsearch
`--cluster_smallmem --usersort` does. It reads every sequence first. Each
is then searched against the centroids so far, as by [dbsrch](dbsrch.md)
(lower case masked). It joins the cluster of its best accepted hit, or of
the most abundant accepted centroid with `--size-order`. Without an
accepted hit it becomes a new centroid.

The centroids are written to `-o` in cluster order, in the input format,
FASTA unwrapped, with the abundance of their cluster appended as `;size=`.
The abundance is the sum of its members' `;size=`, 1 for a member without
one. The hits go
to `--hits` in UC format: an S record for each centroid and an H record for
each member, in input order, then a C record for each cluster with its
abundance. An H record's CIGAR is whole, `=` only for a perfect alignment.

Sorted first by [sortby](sortby.md), the input gives vsearch's other
clustering commands. `sortby --length` with `--maxrejects 8` gives
`--cluster_fast`, and `sortby --size` gives `--cluster_size`.

## Options

- `-i`, `--in FILE`: Input FASTA or FASTQ file. Default: standard input.
- `-o`, `--out FILE`: Output file of the centroids, in the input format.
  Default: standard output.
- `--hits FILE`: Output UC file of the clusters.
- `-t`, `--threads N`: Cluster on N threads, 0 for one per logical core.
  The clusters are the same, written in the same order, on any number.
  Default: 0.
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
- `--id SHARE`: Accept hits of at least this identity, 0 to 1. Required.
- `--iddef 0..4`, `--weak-id SHARE`, `--max-diffs N`, `--max-subs N`,
  `--max-gaps N`, `--min-cols N`, `--left-just`, `--right-just`,
  `--query-cov SHARE`, `--target-cov SHARE`, `--max-id SHARE`, `--mid PCT`,
  `--id-prefix N`, `--id-suffix N`, `--max-qt RATIO`, `--min-qt RATIO`,
  `--min-sl SHARE`, `--max-sl SHARE`, `--max-qsize N`, `--min-tsize N`,
  `--min-size-ratio RATIO`, `--max-size-ratio RATIO`, `--self`, `--self-id`:
  The acceptance criteria of [dbsrch](dbsrch.md), the sequence as query and
  the centroids as targets.
- `--match SCORE`, `--mismatch SCORE`, `--gapopen PENALTIES`,
  `--gapext PENALTIES`, `--n-mismatch`: The scores of [dbsrch](dbsrch.md).
  Defaults: 2, -4, 20I/2E, 2I/1E.

## vsearch

clust is `--cluster_smallmem --usersort`, and with [sortby](sortby.md)
before it `--cluster_fast` and `--cluster_size`. Like [derep](derep.md) it
reads and writes abundances unless told not to, which is vsearch's
`--sizein --sizeout`. The rest of what vsearch's clustering commands do is
other tools. `--minseqlength 32` and `--maxseqlength 50000` are `sqfilt
--min-len 32 --max-len 50000` before it. Its masking is
[lcmask](lcmask.md). The relabelling and sorting of centroids
(`--relabel`, `--clusterout_sort`, `--clusterout_id`) are
[relab](relab.md) and [sortby](sortby.md) after it. `--clusters`,
`--msaout`, `--consout` and `--profile` are [clspl](clspl.md) and
[clmsa](clmsa.md) on the hits. Where
`--cluster_fast` writes `=` for an alignment that is perfect but for
terminal gaps, clust writes the whole CIGAR, which clmsa needs.

vsearch clusters on several threads in rounds of one sequence a thread,
searched at once against the centroids so far. Each sequence of a round is
then matched against the round's new centroids too. clust searches batches
of two sequences a thread at once against the centroids so far, then decides
them in order. A sequence is searched again, against the centroids then,
when a centroid found earlier in its batch could have changed its search.
That is when the centroid shares enough words with it and would have come
before the last candidate the search took, or when the search ran out of
candidates. Any other centroid would have been a candidate the search never
reached, so the sequence gets the same hits it gets on one thread. Like
vsearch's, its clusters are the same on any number of threads.

## Examples

```console
$ sortby --size -i amplicons.fa > sorted.fa
$ clust -i sorted.fa --id 0.97 --hits clusters.uc | grep '>'
>a1;size=62
>a2;size=40
>a3;size=25
>a5;size=4
$ cut -f 1-4,9-10 clusters.uc
S	0	160	*	a1;size=60	*
S	1	160	*	a2;size=40	*
S	2	160	*	a3;size=25	*
S	3	160	*	a5;size=4	*
H	0	160	99.4	a4;size=2	a1;size=60
C	0	62	*	a1;size=60	*
C	1	40	*	a2;size=40	*
C	2	25	*	a3;size=25	*
C	3	4	*	a5;size=4	*
```

vsearch's `--cluster_fast amplicons.fa --id 0.97 --sizein --sizeout
--centroids c.fa`:

```console
$ sqfilt --min-len 32 --max-len 50000 -i amplicons.fa | lcmask | sortby --length | clust --id 0.97 --maxrejects 8 | conv --fasta | grep ">"
>a1;size=62
>a2;size=40
>a3;size=25
>a5;size=4
```

## Exit status

0 on success, 1 for invalid input (including an invalid `;size=`), 2 for
invalid arguments, such as a missing `--id`.

## See also

[sortby](sortby.md), [unoise](unoise.md), [clspl](clspl.md),
[clmsa](clmsa.md), [otutab](otutab.md), [dbsrch](dbsrch.md)
