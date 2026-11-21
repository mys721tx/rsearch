# subseq

Keep one region of every read.

## Synopsis

```
subseq [-i FILE] [-o FILE] [--start POS] [--end POS]
```

## Description

subseq writes the region of each read from `--start` to `--end`, both
counted from 1 and included, with its scores. A read that ends before
`--end` is cut at its end, and one that ends before `--start` is written
empty. Headers are written unchanged, in the input format, FASTA
unwrapped.

## Options

- `-i`, `--in FILE`: Input FASTA or FASTQ file. Default: standard input.
- `-o`, `--out FILE`: Output file, in the input format. Default: standard
  output.
- `--start POS`: First base kept, from 1. Default: 1.
- `--end POS`: Last base kept. Default: the end of the read.

## vsearch

subseq is the region half of `--fastx_getsubseq` (`--subseq_start`,
`--subseq_end`); the selection half is [getseq](getseq.md). vsearch's
`--fastx_getsubseq in --label s1 --subseq_start 3 --subseq_end 6` is
`getseq --label s1 | subseq --start 3 --end 6`.

## Examples

```console
$ subseq --start 3 --end 6 -i seqs.fa | head -4
>s1;size=3
GTTG
>s2
gttg
$ getseq --label s4 -i seqs.fa | subseq --start 13
>s4
AAC
```

## Exit status

0 on success, 1 for invalid input, 2 for invalid arguments.

## See also

[getseq](getseq.md), [sqtrim](sqtrim.md), [redig](redig.md)
