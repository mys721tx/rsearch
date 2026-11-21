# joinpa

Join read pairs with a pad.

## Synopsis

```
joinpa [-i FILE] --reverse FILE [-o FILE] [--pad SEQ] [--pad-qual QUALS]
       [--ascii 33|64]
```

## Description

joinpa writes each forward read followed by the pad and the reverse
complement of the read at the same position in `--reverse`, keeping its
case. The scores are those of the forward read, the pad's and those of the
mate, reversed. The forward header is kept whole. Both inputs are FASTQ,
with as many reads each.

## Options

- `-i`, `--in FILE`: Input FASTQ file of the forward reads. Default:
  standard input.
- `-o`, `--out FILE`: Output FASTQ file of the joined reads. Default:
  standard output.
- `--reverse FILE`: Input FASTQ file of the reverse reads. Required.
- `--pad SEQ`: Sequence put between the mates. Default: NNNNNNNN.
- `--pad-qual QUALS`: Score symbols of the pad, as long as it. Default: the
  symbol of Q40 for each base of the pad.
- `--ascii 33|64`: ASCII offset of the default pad scores. Default: 33.

## vsearch

joinpa is `--fastq_join`, with `--join_padgap` as `--pad` and
`--join_padgapq` as `--pad-qual`. Its relabelling is [relab](relab.md) after
it, and its FASTA output is [conv](conv.md). vsearch's `--fastq_join in
--reverse rev --relabel J --fastaout out` is `joinpa -i in --reverse rev | relab
--relabel J | conv --fasta`.

## Examples

```console
$ joinpa -i fwd.fq --reverse rev.fq --pad NN --pad-qual II | head -4
@p1/1
GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTNNATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCG
+
IIIIHHHGGGFFFFEEEDDDCCCCBBBAAA@@@@???>>>II::;;;<<<===>>>???@@@AABBBCCCDDDEEEFFFGGG
```

## Exit status

0 on success, 1 for invalid input (a malformed record, files of different
lengths), 2 for invalid arguments.

## See also

[mergpa](mergpa.md), [syncpa](syncpa.md), [relab](relab.md),
[conv](conv.md)
