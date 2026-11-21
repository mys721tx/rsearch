# rerep

Write each read as many times as its abundance (`;size=`).

## Synopsis

```
rerep [-i FILE] [-o FILE] [--no-sizeout]
```

## Description

rerep undoes dereplication. A read with `;size=N` is written N times. Each
copy has `;size=1` in place of the read's abundance, or no abundance with
`--no-sizeout`. A read without `;size=` is written once.
FASTQ copies keep their scores. Records are written in the input format,
FASTA unwrapped. When some reads have no abundance, rerep warns that it
assumed 1.

## Options

- `-i`, `--in FILE`: Input FASTA or FASTQ file. Default: standard input.
- `-o`, `--out FILE`: Output file, in the input format. Default: standard
  output.
- `--no-sizeout`: Remove the abundance (`;size=`) instead of writing
  `;size=1`.

## vsearch

rerep is `--rereplicate`. It writes `;size=1` on each copy by default, as
vsearch does with `--sizeout`; `--no-sizeout` gives vsearch's output
without it. Unlike vsearch it also reads FASTQ. vsearch's
`--rereplicate in --sizeout --output out` is `rerep | conv --fasta`.

## Examples

```console
$ sqfilt --min-size 2 -i seqs.fa | rerep
>s1;size=1
ACGTTGCAAGGT
>s1;size=1
ACGTTGCAAGGT
>s1;size=1
ACGTTGCAAGGT
>s3;size=1
TTGACCGTAGCA
>s3;size=1
TTGACCGTAGCA
$ rerep --no-sizeout -i seqs.fa > copies.fa
Warning: Missing abundance information for some input sequences, assumed 1
$ grep -c "^>" copies.fa
8
```

## Exit status

0 on success, 1 for invalid input (including an invalid `;size=`), 2 for
invalid arguments.

## See also

[derep](derep.md), [relab](relab.md), [conv](conv.md)
