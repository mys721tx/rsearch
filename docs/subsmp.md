# subsmp

Draw a random subset of records, or of reads by abundance.

## Synopsis

```
subsmp (--count N | --pct PCT) [-i FILE] [-o FILE] [--discarded FILE]
       [--by-size] [--allow-fewer] [--seed SEED]
```

## Description

subsmp reads and checks the whole input, then draws `--count` records, or
`--pct` percent of them rounded down, and writes them in input order,
unchanged. With `--by-size`, each record counts as many reads as its
abundance (`;size=`), and the reads are drawn. Each record with reads drawn
is written with their number as its abundance. `--discarded` takes the
records, or the reads, not drawn.

With a seed, the draw is the one vsearch `--fastx_subsample --randseed`
gives. Without one, or with 0, the seed comes from the system, and the draw
changes from run to run. Asking for more than there is fails, unless
`--allow-fewer` is given, in which case everything is drawn.

## Options

- `-i`, `--in FILE`: Input FASTA or FASTQ file. Default: standard input.
- `-o`, `--out FILE`: Output file of the records drawn. Default: standard
  output.
- `--discarded FILE`: Output file of the records, or reads, not drawn.
- `--count N`: Draw N records, or N reads with `--by-size`.
- `--pct PCT`: Draw this percentage of the records, or of the reads with
  `--by-size`, rounded down.
- `--by-size`: Count each record as its abundance, and write the reads
  drawn as it.
- `--allow-fewer`: Draw everything when asked for more than there is,
  instead of failing.
- `--seed SEED`: Seed of the random numbers. Default: one from the system.

One of `--count` and `--pct` is required.

## vsearch

subsmp is `--fastx_subsample`, with `--sample_size` as `--count`,
`--sample_pct` as `--pct`, `--sizein` as `--by-size` (vsearch's
`--sizeout` comes with it) and `--randseed` as `--seed`. vsearch fails
when asked for more than there is. `--allow-fewer` draws everything
instead.

## Examples

```console
$ subsmp --count 2 --seed 1 -i seqs.fa
>s1;size=3
ACGTTGCAAGGT
>s2
acgttgcaaggt
$ subsmp --by-size --count 4 --seed 1 -i seqs.fa
>s1;size=2
ACGTTGCAAGGT
>s2;size=1
acgttgcaaggt
>s4;size=1
ACCTTGCAAGGTAAC
$ subsmp --count 9 -i seqs.fa; echo "exit $?"
Error: Cannot subsample 9 reads from a sample of 5 (use --allow-fewer to keep them all)
exit 1
```

## Exit status

0 on success, 1 for invalid input or a draw larger than the input without
`--allow-fewer`, 2 for invalid arguments.

## See also

[shufl](shufl.md), [rerep](rerep.md), [topn](topn.md)
