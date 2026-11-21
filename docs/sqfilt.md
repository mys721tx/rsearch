# sqfilt

Keep or discard whole reads, or pairs of reads.

## Synopsis

```
sqfilt [-i FILE] [-o FILE] [--discarded FILE]
       [--reverse FILE [--rev-out FILE] [--rev-discarded FILE]]
       [--min-len N] [--max-len N] [--max-ns N] [--max-ee EE]
       [--max-ee-rate RATE] [--min-qual Q] [--min-size N] [--max-size N]
       [--min-upper-pct PCT] [--max-upper-pct PCT]
       [--min-non-n-pct PCT] [--max-non-n-pct PCT]
       [--ascii 33|64] [--qmin Q] [--qmax Q]
```

## Description

sqfilt discards a read when it fails any of its filters: length, number of
N, expected errors and expected errors per base, lowest score, abundance
(`;size=`), and the percentages of upper-case residues and of residues
other than N that masking leaves. The kept reads go to `-o`, the discarded
ones to `--discarded`, unchanged, in the input format, FASTA unwrapped.

With `--reverse`, the input holds forward reads and `--reverse` their mates
at the same positions, in the same format. A pair is discarded when either
read fails. The mates go to `--rev-out` and `--rev-discarded`.

It is the only tool that discards reads. `--min-len` is 1 by default, as
vsearch's `--fastq_minlen` and `--minseqlength` are for most commands.
Empty reads are therefore discarded unless `--min-len 0` is given.

## Options

- `-i`, `--in FILE`: Input FASTQ or FASTA file (the forward reads with
  `--reverse`). Default: standard input.
- `-o`, `--out FILE`: Output file of the kept reads. Default: standard
  output.
- `--discarded FILE`: Output file of the discarded reads.
- `--reverse FILE`: Input file of the reverse reads, filtered in pairs with
  the input.
- `--rev-out FILE`: Output file of the kept reverse reads. Needs
  `--reverse`.
- `--rev-discarded FILE`: Output file of the discarded reverse reads. Needs
  `--reverse`.
- `--min-len N`: Discard reads shorter than N. Default: 1.
- `--max-len N`: Discard reads longer than N.
- `--max-ns N`: Discard reads with more than N bases N.
- `--max-ee EE`: Discard reads with more expected errors than EE.
- `--max-ee-rate RATE`: Discard reads with more expected errors per base
  than RATE.
- `--min-qual Q`: Discard reads with a score below Q.
- `--min-size N`: Discard reads with an abundance (`;size=`, 1 when absent)
  below N.
- `--max-size N`: Discard reads with an abundance above N.
- `--min-upper-pct PCT`: Discard reads with a lower percentage of
  upper-case residues (the residues soft masking leaves).
- `--max-upper-pct PCT`: Discard reads with a higher percentage of
  upper-case residues.
- `--min-non-n-pct PCT`: Discard reads with a lower percentage of residues
  other than N (the residues hard masking leaves).
- `--max-non-n-pct PCT`: Discard reads with a higher percentage of residues
  other than N.
- `--ascii 33|64`: ASCII offset of the quality scores. Default: 33.
- `--qmin Q`: Lowest score accepted in the input. Default: 0.
- `--qmax Q`: Highest score accepted in the input. Default: 126 minus
  `--ascii`.

The options that read scores (`--max-ee`, `--max-ee-rate`, `--min-qual`,
`--ascii`, `--qmin`, `--qmax`) are refused for FASTA input.

## vsearch

sqfilt is the filtering half of `--fastq_filter` and `--fastx_filter`
(`--fastq_minlen`, `--fastq_maxlen`, `--fastq_maxns`, `--fastq_maxee`,
`--fastq_maxee_rate`, `--fastq_minqual`, `--minsize`, `--maxsize`, the
discarded outputs and `--reverse`). It is also `--minseqlength`,
`--maxseqlength`, `--minsize`, `--maxsize`, `--minuniquesize` and
`--maxuniquesize` of every command. Last, it is `--min_unmasked_pct` and
`--max_unmasked_pct` of `--fastx_mask`, split in two. After soft masking
they are `--min-upper-pct` and `--max-upper-pct`, and after `--hardmask`
`--min-non-n-pct` and `--max-non-n-pct`. Several vsearch commands discard
sequences shorter than 32 by default (`--minseqlength 32` for clustering,
searching and UDB files). In rsearch that is an explicit `sqfilt --min-len
32` before them.

## Examples

```console
$ sqfilt --max-ee 0.1 --discarded bad.fq -i reads.fq | grep '^@'
@r2
@r4
$ grep '^@' bad.fq
@r1
@r3
$ sqfilt --min-size 2 -i seqs.fa
>s1;size=3
ACGTTGCAAGGT
>s3;size=2
TTGACCGTAGCA
$ lcmask -i repeats.fa
>m1
ACGTTGCAAGGTcacacacacacacacacacaGGTACCAT
>m2
ACGTTGCAAGGTAGGTCAGTTCAGGACCTTAGCATACCAT
$ lcmask -i repeats.fa | sqfilt --min-upper-pct 80
>m2
ACGTTGCAAGGTAGGTCAGTTCAGGACCTTAGCATACCAT
```

Pairs, as vsearch's `--fastq_filter fwd.fq --reverse rev.fq --fastq_maxns 0`:

```console
$ sqfilt -i fwd.fq --reverse rev.fq --rev-out kept_rev.fq --max-ns 0 | grep '^@'
@p1/1
@p3/1
$ grep '^@' kept_rev.fq
@p1/2
@p3/2
```

## Exit status

0 on success, 1 for invalid input (a malformed record, a score out of
range, an invalid `;size=`, reverse reads more or fewer than the forward
ones), 2 for invalid arguments.

```console
$ sqfilt -i fwd.fq --reverse seqs.fa; echo "exit $?"
Error: forward and reverse reads must be in the same format
exit 1
```

## See also

[sqtrim](sqtrim.md), [lcmask](lcmask.md), [derep](derep.md),
[mergpa](mergpa.md), [topn](topn.md)
