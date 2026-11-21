# sortby

Sort reads by decreasing length or abundance.

## Synopsis

```
sortby (--length | --size) [-i FILE] [-o FILE]
```

## Description

sortby reads and checks the whole input, then writes it sorted:

- `--length`: by decreasing length, then decreasing abundance (`;size=`),
  then header, then input order;
- `--size`: by decreasing abundance, then header, then input order.

A record without `;size=` has an abundance of 1. Records are written
unchanged, in the input format, FASTA unwrapped.

The tools work on their input in the order given. sortby therefore comes
before those that vsearch sorts for. `sortby --length` before
[clust](clust.md) gives `--cluster_fast`. `sortby --size` gives
`--cluster_size` and the order of [unoise](unoise.md), [uchime](uchime.md)
and [mchime](mchime.md) de novo.

## Options

- `--length`: Sort by decreasing length, as vsearch `--sortbylength`.
- `--size`: Sort by decreasing abundance (`;size=`), as vsearch
  `--sortbysize`.
- `-i`, `--in FILE`: Input FASTA or FASTQ file. Default: standard input.
- `-o`, `--out FILE`: Output file, in the input format. Default: standard
  output.

## vsearch

sortby is `--sortbylength` and `--sortbysize`, and the sorting step inside
clustering and de novo chimera detection, which use the same order. The
filters of those commands are [sqfilt](sqfilt.md), `--topn` is
[topn](topn.md), and relabelling is [relab](relab.md). `--sortbysize in
--minsize 2 --topn 5 --output out` is `sqfilt --min-size 2 | sortby --size |
topn 5 | conv --fasta`.

## Examples

```console
$ sortby --size -i seqs.fa
>s1;size=3
ACGTTGCAAGGT
>s3;size=2
TTGACCGTAGCA
>s2
acgttgcaaggt
>s4
ACCTTGCAAGGTAAC
>s5
GGGAAATTTCCC
$ sortby --length -i seqs.fa | grep '>'
>s4
>s1;size=3
>s3;size=2
>s2
>s5
```

## Exit status

0 on success, 1 for invalid input (including an invalid `;size=`, after
which nothing is written), 2 for invalid arguments.

## See also

[shufl](shufl.md), [topn](topn.md), [clust](clust.md),
[uchime](uchime.md)
