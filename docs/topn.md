# topn

Keep the first N records.

## Synopsis

```
topn [-i FILE] [-o FILE] N
```

## Description

topn writes the first N records of its input, unchanged, in the input
format, FASTA unwrapped. As vsearch does, it reads and checks the whole
input. Invalid input after the first N records therefore still fails.

## Options

- `N`: Number of records to keep. Required.
- `-i`, `--in FILE`: Input FASTA or FASTQ file. Default: standard input.
- `-o`, `--out FILE`: Output file, in the input format. Default: standard
  output.

## vsearch

topn is the `--topn` option of every vsearch command that has one, as a
step after the tool that writes the records. `--sortbysize in --topn 5` is
`sortby --size | topn 5`.

## Examples

```console
$ sortby --size -i seqs.fa | topn 2
>s1;size=3
ACGTTGCAAGGT
>s3;size=2
TTGACCGTAGCA
```

## Exit status

0 on success, 1 for invalid input anywhere in it, 2 for invalid arguments.

## See also

[sortby](sortby.md), [derep](derep.md), [getseq](getseq.md)
