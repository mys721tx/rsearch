# udbinfo

Report a UDB file's header.

## Synopsis

```
udbinfo [-i FILE] [-o FILE]
```

## Description

udbinfo writes vsearch's `--udbinfo` report of a UDB file: the number of
sequences, the bits of a sequence index, the alphabet, the word length and
the other fields of its header, as the file holds them. Like vsearch, it
checks only the header, and that the file is long enough for its
sequences.

## Options

- `-i`, `--in FILE`: Input UDB file. Default: standard input.
- `-o`, `--out FILE`: Output file of the report. Default: standard output.

## vsearch

udbinfo is `--udbinfo`. vsearch writes the report to standard error and
`--log`; udbinfo writes it to `-o`. It also reads standard input, which
vsearch refuses. vsearch needs the file's size to check the number of
sequences, which udbinfo counts as it reads.

## Examples

```console
$ mkudb -i db.fa -o db.udb
$ udbinfo -i db.udb
           Seqs  3
     SeqIx bits  32
          Alpha  nt (4)
     Word width  8
          Slots  0
      Dict size  65536 (65.5k)
         DBstep  1
        DBAccel  100%
```

## Exit status

0 on success, 1 for invalid input (a file that is not UDB, or too short),
2 for invalid arguments.

```console
$ udbinfo -i db.fa; echo "exit $?"
Error: invalid UDB file
exit 1
```

## See also

[udbstat](udbstat.md), [mkudb](mkudb.md), [conv](conv.md)
