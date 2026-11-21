# udbstat

Report a UDB file's word statistics.

## Synopsis

```
udbstat [-i FILE] [-o FILE]
```

## Description

udbstat reads a whole UDB file, checking it as vsearch does. It writes
vsearch's `--udbstats` report: a summary line of the database (its size in
residues, the number of sequences and their lengths), the index's
parameters, how many sequences hold each word (median and mean), the words
held by the most sequences, and the distribution of the words over these
counts.

## Options

- `-i`, `--in FILE`: Input UDB file. Default: standard input.
- `-o`, `--out FILE`: Output file of the report. Default: standard output.

## vsearch

udbstat is `--udbstats`, which writes its report to `--log` and nothing
without it. udbstat writes it to `-o`, the database's summary line
included. It also reads standard input. Where a spoiled length makes
vsearch run out of memory and abort, udbstat reads on and reports the error
that follows.

## Examples

```console
$ mkudb -i db.fa -o db.udb
$ udbstat -i db.udb | head -16
480 nt in 3 seqs, min 160, max 160, avg 160

      Alphabet  nt
    Word width  8
     Word ones  8
        Spaced  No
        Hashed  No
         Coded  No
       Stepped  No
         Slots  65536 (65.5k)
       DBAccel  100%

       480  DB size (0.5k)
       459  Words
         0  Median size
       0.0  Mean size
```

## Exit status

0 on success, 1 for invalid input, 2 for invalid arguments.

## See also

[udbinfo](udbinfo.md), [mkudb](mkudb.md)
