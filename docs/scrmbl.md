# scrmbl

Shuffle the bases within each read.

## Synopsis

```
scrmbl [-i FILE] [-o FILE] [--kmer K] [--seed SEED]
```

## Description

scrmbl rearranges the bases of each read at random. With `--kmer 1`, the
default, the bases are shuffled; with K from 2 to 9, they are rearranged
keeping the count of every word of up to K bases. With a seed, each read is
scrambled as vsearch `--scramble --randseed` scrambles it. Without one, or
with 0, the seed comes from the system, and the output changes from run to
run. Scores and headers are written unchanged, in the input format, FASTA
unwrapped.

## Options

- `-i`, `--in FILE`: Input FASTA or FASTQ file. Default: standard input.
- `-o`, `--out FILE`: Output file, in the input format. Default: standard
  output.
- `--kmer K`: Keep the counts of words of up to K bases, 1 to 9.
  Default: 1.
- `--seed SEED`: Seed of the random numbers. Default: one from the system.

## vsearch

scrmbl is `--scramble`, with `--randseed` as `--seed` and
`--scramble_kmer` (vsearch's name for the word length) as `--kmer`.

## Examples

```console
$ scrmbl --seed 1 -i seqs.fa | head -4
>s1;size=3
TCAATAGTCGGG
>s2
ggtaggctacta
$ scrmbl --seed 1 --kmer 2 -i seqs.fa | head -2
>s1;size=3
AGTTGGCAACGT
```

## Exit status

0 on success, 1 for invalid input, 2 for invalid arguments.

## See also

[shufl](shufl.md), [subsmp](subsmp.md)
