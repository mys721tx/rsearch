# relab

Rewrite sequence labels: relabel, add or remove annotations.

## Synopsis

```
relab [-i FILE] [-o FILE]
      [--relabel PREFIX | --relabel-self | --relabel-sha1 | --relabel-md5]
      [--relabel-keep] [--sizeout] [--ee [--ascii 33|64]] [--lengthout]
      [--xsize] [--xee] [--xlength] [--sample NAME] [--label-suffix TEXT]
```

## Description

relab is the only tool that rewrites headers (besides the `;size=` that
[derep](derep.md), [rerep](rerep.md) and the clustering tools write as their
result). For each record:

1. The label is replaced by PREFIX and the record's number from 1
   (`--relabel`), by the sequence itself (`--relabel-self`), or by the SHA1
   or MD5 digest of the sequence upper-cased with U as T (`--relabel-sha1`,
   `--relabel-md5`). `--relabel-keep` adds the old header after a space.
   Without them the header is kept, less the `;size=`, `;ee=` and
   `;length=` annotations that `--xsize`, `--xee` and `--xlength` remove.
2. Then the suffix, `;sample=`, `;size=`, `;ee=` and `;length=` are
   appended, in that order, as asked.

`--sizeout` writes the abundance read from the header, or 1. `--ee` writes
the expected errors of a FASTQ read, the sum of 10^(-Q/10) over its bases,
left to right. Records are written in the input format, FASTA unwrapped.

## Options

- `-i`, `--in FILE`: Input FASTA or FASTQ file. Default: standard input.
- `-o`, `--out FILE`: Output file, in the input format. Default: standard
  output.
- `--sizeout`: Append the abundance (`;size=`), 1 when the header has none.
- `--ee`: Append the expected errors (`;ee=`). FASTQ only.
- `--ascii 33|64`: ASCII offset of the quality scores for `--ee`.
  Default: 33. Needs `--ee`.
- `--relabel PREFIX`: Relabel records as PREFIX1, PREFIX2, ...
- `--relabel-self`: Relabel records with their sequence.
- `--relabel-sha1`: Relabel records with the SHA1 digest of their sequence.
- `--relabel-md5`: Relabel records with the MD5 digest of their sequence.
- `--relabel-keep`: Keep the old header after the new label, separated by a
  space.
- `--xsize`: Remove abundance annotations (`;size=`).
- `--xee`: Remove expected-error annotations (`;ee=`).
- `--lengthout`: Append the sequence length (`;length=`).
- `--xlength`: Remove length annotations (`;length=`).
- `--sample NAME`: Append a sample name (`;sample=NAME`).
- `--label-suffix TEXT`: Append TEXT to the label.

Only one of the four relabelling options may be given.

## vsearch

relab stands for the label options every vsearch command takes:
`--relabel`, `--relabel_self`, `--relabel_sha1`, `--relabel_md5`,
`--relabel_keep`, `--sizeout`, `--eeout` (and `--fastq_eeout`, as `--ee`),
`--xsize`, `--xee`, `--lengthout`, `--xlength`, `--sample` and
`--label_suffix`. In vsearch they rewrite the records a command writes; in
rsearch they are a step at the end of the pipeline. `--relabel` numbers the
records relab receives, which are those the vsearch command would write.
vsearch's `--sizein` has no counterpart. relab always reads `;size=`.

## Examples

```console
$ relab --relabel OTU --sizeout -i seqs.fa
>OTU1;size=3
ACGTTGCAAGGT
>OTU2;size=1
acgttgcaaggt
>OTU3;size=2
TTGACCGTAGCA
>OTU4;size=1
ACCTTGCAAGGTAAC
>OTU5;size=1
GGGAAATTTCCC
$ relab --relabel-md5 --relabel-keep -i seqs.fa | head -2
>450979cc66007812800a187f961672a6 s1;size=3
ACGTTGCAAGGT
$ relab --ee --xsize --sample S1 -i reads.fq | grep '^@'
@r1;sample=S1;ee=0.3039
@r2;sample=S1;ee=0.01487
@r3;sample=S1;ee=0.4278
@r4;sample=S1;ee=0.004601
```

## Exit status

0 on success, 1 for invalid input (including an invalid `;size=` with
`--sizeout`), 2 for invalid arguments, such as `--ee` on FASTA input.

```console
$ relab --ee -i seqs.fa; echo "exit $?"
error: --ee cannot be used with FASTA input

Usage: relab [OPTIONS]

For more information, try '--help'.
exit 2
```

## See also

[derep](derep.md), [sqfilt](sqfilt.md), [conv](conv.md),
[mergpa](mergpa.md), [clust](clust.md)
