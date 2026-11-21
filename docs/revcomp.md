# revcomp

Reverse-complement reads.

## Synopsis

```
revcomp [-i FILE] [-o FILE]
```

## Description

revcomp writes each read reverse-complemented. Each base is complemented
keeping its case, IUPAC codes included, and U as A; FASTQ scores are
reversed. Headers are written unchanged, in the input format, FASTA
unwrapped.

## Options

- `-i`, `--in FILE`: Input FASTA or FASTQ file. Default: standard input.
- `-o`, `--out FILE`: Output file, in the input format. Default: standard
  output.

## vsearch

revcomp is `--fastx_revcomp`. Its relabelling is [relab](relab.md) after it,
its FASTA output [conv](conv.md).

## Examples

```console
$ revcomp -i reads.fq | head -4
@r1
TCGATCGATCCTAGCTAGCTACGTACGT
+
-./0123456789:<=>?@ABCDEFGHI
$ printf '>x\nACGTUacgtuRYN\n' | revcomp
>x
NRYaacgtAACGT
```

## Exit status

0 on success, 1 for invalid input, 2 for invalid arguments.

## See also

[dbornt](dbornt.md), [relab](relab.md), [conv](conv.md)
