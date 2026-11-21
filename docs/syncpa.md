# syncpa

Pair the reads of two files by label.

## Synopsis

```
syncpa [-i FILE] --reverse FILE [-o FILE] [--rev-out FILE]
       [--orphans FILE] [--rev-orphans FILE] [--separators CHARS]
```

## Description

syncpa matches each forward read with its mate in `--reverse`, whatever
their order. Mates share their header up to its first space or tab, less a
final 1 or 2 after one of the `--separators`. The forward reads with a mate
are written to `-o` in their order, and their mates in the same order to
`--rev-out`. The reads without a mate go to `--orphans` and
`--rev-orphans`. Reads are written as read, in the input format, FASTA
unwrapped. No abundance is read.

## Options

- `-i`, `--in FILE`: Input file of the forward reads, FASTA or FASTQ.
  Default: standard input.
- `-o`, `--out FILE`: Output file of the forward reads with a mate.
  Default: standard output.
- `--reverse FILE`: Input file of the reverse reads. Required.
- `--rev-out FILE`: Output file of the mates, in the order of the forward
  reads.
- `--orphans FILE`: Output file of the forward reads without a mate.
- `--rev-orphans FILE`: Output file of the reverse reads without a mate.
- `--separators CHARS`: Characters that may come before the final 1 or 2.
  Default: `/`.

## vsearch

syncpa is `--fastx_syncpairs`, with `--read_separators` as `--separators`.
Its FASTA outputs are [conv](conv.md) after it.

## Examples

```console
$ grep "^@" mixed.fq
@p3/2
@x9/2
@p1/2
$ syncpa -i fwd.fq --reverse mixed.fq --rev-out mates.fq --orphans orphans.fq | grep '^@'
@p1/1
@p3/1
$ grep '^@' mates.fq
@p1/2
@p3/2
$ grep "^@" orphans.fq
@p2/1
```

## Exit status

0 on success, 1 for invalid input, 2 for invalid arguments.

## See also

[mergpa](mergpa.md), [joinpa](joinpa.md), [sqfilt](sqfilt.md)
