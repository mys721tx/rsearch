# shufl

Shuffle the order of records.

## Synopsis

```
shufl [-i FILE] [-o FILE] [--seed SEED]
```

## Description

shufl reads and checks the whole input, then writes its records in a random
order. With a seed, the order is the one vsearch `--shuffle --randseed`
gives. Without one, or with 0, the seed comes from the system, and the order
changes from run to run. Records are written unchanged, in the input
format, FASTA unwrapped.

## Options

- `-i`, `--in FILE`: Input FASTA or FASTQ file. Default: standard input.
- `-o`, `--out FILE`: Output file, in the input format. Default: standard
  output.
- `--seed SEED`: Seed of the random numbers. Default: one from the system.

## vsearch

shufl is `--shuffle`, with `--randseed` as `--seed`. Its `--topn` is
[topn](topn.md) after it, its relabelling [relab](relab.md).

## Examples

```console
$ shufl --seed 1 -i seqs.fa | grep '>'
>s3;size=2
>s4
>s2
>s5
>s1;size=3
```

## Exit status

0 on success, 1 for invalid input, 2 for invalid arguments.

## See also

[subsmp](subsmp.md), [scrmbl](scrmbl.md), [sortby](sortby.md)
