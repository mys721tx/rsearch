# clmsa

Align each cluster's members to its centroid, from clustering hits.

## Synopsis

```
clmsa --hits FILE [-i FILE] [-o FILE] [--msa | --consensus | --profile]
      [--no-sizein] [--no-sizeout]
```

## Description

clmsa takes the hits of [clust](clust.md) or [unoise](unoise.md) and the
input they clustered, in the order they read it. The S and H records give
each record's cluster, matched in order by header, and each member's
alignment to its centroid (its CIGAR and strand). The records of each
cluster are aligned together, members of the minus strand
reverse-complemented. Each column's residues are then counted, weighted by
abundance (`;size=`). One of three outputs is written, FASTA unwrapped:

- `--msa` (the default): each cluster's multiple alignment, the centroid
  marked with `*`, followed by its consensus;
- `--consensus`: each cluster's consensus sequence, labelled with its
  centroid, number of sequences and abundance;
- `--profile`: each cluster's profile, a line for each column with its
  position, consensus residue and the counts of A, C, G, T or U, gaps and
  ambiguous symbols,
  after the cluster's header and an empty line.

Every record is matched and every CIGAR checked before anything is
written.

## Options

- `--hits FILE`: UC file of the clusters. Required.
- `-i`, `--in FILE`: Input FASTA or FASTQ file that was clustered.
  Default: standard input.
- `-o`, `--out FILE`: Output file. Default: standard output.
- `--msa`: Write each cluster's multiple alignment and consensus, as
  vsearch `--msaout`.
- `--consensus`: Write each cluster's consensus sequence, as vsearch
  `--consout`.
- `--profile`: Write each cluster's profile, as vsearch `--profile`.
- `--no-sizein`: Count each sequence once, rather than by its `;size=`.
- `--no-sizeout`: Do not append abundances (`;size=`) to the labels.

## vsearch

clmsa is the `--msaout`, `--consout` and `--profile` outputs of vsearch's
clustering commands, with their `--sizein --sizeout` by default. It needs
whole CIGARs, as clust writes. `--cluster_fast`'s `=` for an alignment with
terminal gaps cannot be aligned. clmsa fails on a `=` between sequences of
different lengths.

## Examples

```console
$ sortby --size -i amplicons.fa > sorted.fa
$ clust -i sorted.fa --id 0.9 --hits clusters.uc > centroids.fa
$ clmsa --hits clusters.uc -i sorted.fa --consensus | grep '>'
>centroid=a1;seqs=4;size=106
>centroid=a3;seqs=1;size=25
$ clmsa --hits clusters.uc -i sorted.fa --msa | grep '>'
>*a1;size=60
>a2;size=40
>a5;size=4
>a4;size=2
>consensus
>*a3;size=25
>consensus
$ clmsa --hits clusters.uc -i sorted.fa --profile | head -5
>centroid=a1;seqs=4;size=106

0	G	0	0	106	0	0	0
1	C	0	106	0	0	0	0
2	G	0	0	106	0	0	0
```

## Exit status

0 on success, 1 for invalid input, including records and hits that do not
match, 2 for invalid arguments.

## See also

[clust](clust.md), [unoise](unoise.md), [clspl](clspl.md)
