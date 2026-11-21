# mkudb

Build a UDB database: sequences and word index.

## Synopsis

```
mkudb [-i FILE] [-o FILE] [--word-len N]
```

## Description

mkudb reads every record of its input, FASTA or FASTQ, then writes
vsearch's UDB file of it. The file holds the sequences as given and their
headers whole. It also holds, for each word of `--word-len` residues, the
sequences holding it. A word holding a lower-case (masked) residue, or a
symbol other than A, C, G, T or U, is not indexed. Masking is therefore
[lcmask](lcmask.md) before it. mkudb warns when some sequences give no word
at all. An input without sequences fails, as a UDB file without them cannot
be read back. It leaves `-o` empty.

[dbsrch](dbsrch.md), [sintax](sintax.md), [dbornt](dbornt.md) and
[uchime](uchime.md) read UDB databases. [udbinfo](udbinfo.md) and
[udbstat](udbstat.md) report on them. [conv](conv.md) turns them back into
FASTA.

## Options

- `-i`, `--in FILE`: Input FASTA or FASTQ file. Default: standard input.
- `-o`, `--out FILE`: Output UDB file. Default: standard output.
- `--word-len N`: Length of the indexed words, 3 to 15. Default: 8.

## vsearch

mkudb is `--makeudb_usearch`, with `--wordlength` as `--word-len`. vsearch
upper-cases every sequence as it reads it. It drops sequences shorter than
32 or longer than 50000 and masks with DUST. In rsearch these steps are
explicit. `--makeudb_usearch in --output db.udb` is `sqfilt --min-len 32
--max-len 50000 | lcmask | mkudb -o db.udb`. `--dbmask none` (or `soft`) is
`lcmask --clear` instead of `lcmask`, and `--hardmask` is `lcmask --clear |
lcmask --hard`.

## Examples

```console
$ sqfilt --min-len 32 --max-len 50000 -i db.fa | lcmask | mkudb -o db.udb
$ udbinfo -i db.udb | head -4
           Seqs  3
     SeqIx bits  32
          Alpha  nt (4)
     Word width  8
$ printf '>a\nACGT\n' | mkudb -o tiny.udb
Warning: 1 of 1 sequences give no word to index: they are shorter than the word length (8), or each of their words holds a lower-case residue or a symbol other than A, C, G, T or U
```

## Exit status

0 on success, 1 for invalid input, including an input without sequences,
2 for invalid arguments.

```console
$ printf '' | mkudb -o empty.udb; echo "exit $?"
Error: no sequences: a UDB file needs at least one
exit 1
```

## See also

[udbinfo](udbinfo.md), [udbstat](udbstat.md), [conv](conv.md),
[lcmask](lcmask.md), [dbsrch](dbsrch.md)
