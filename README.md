# rsearch

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

rsearch is a port of [vsearch](https://github.com/torognes/vsearch) to Rust,
split into 40 small tools that each do one thing. Together they reproduce
every vsearch command. A command that does several jobs becomes a pipeline
of tools. A family of commands that do one job becomes one tool. The
results are vsearch's, checked against vsearch's output by the tests and by
the comparison framework in [`compare/`](compare/README.md).

rsearch matches vsearch's results, not its command line. The options have
new names in words (`--max-ee` for `--fastq_maxee`, `--min-size` for
`--minsize`). Each tool takes only the options of its own job. The
defaults that vsearch builds into a command are explicit steps of the
pipeline. Each tool's page in [`docs/`](docs) names the vsearch commands and
options it stands for.

## Build

```bash
cargo build --release
```

The tools are built in `target/release/`; put that directory on your `PATH`
or copy the tools you need. Every tool prints its options with `-h` and a
longer description with `--help`.

The tests run every tool against expectations taken from vsearch's output.
They also check that each tool's `--help` agrees with its page and that
every example in the pages gives the output shown:

```bash
cargo test
```

## Shared function

- **Input and output.** `-i`/`--in` reads a file, standard input by default;
  `-o`/`--out` writes a file, standard output by default. So tools join with
  pipes.
- **Formats.** The tools read FASTA or FASTQ, recognised from the input,
  and write records in the format they read, FASTA unwrapped. Changing the
  format or wrapping lines is [`conv`](docs/conv.md)'s job only. End a
  pipeline with `conv --fasta` for vsearch's FASTA, wrapped at 80.
- **Headers.** Headers pass through whole and unchanged. No tool cuts a
  label at its first space or tab, as vsearch does unless given
  `--notrunclabels`. Rewriting headers is [`relab`](docs/relab.md)'s job
  only, except the `;size=` that dereplication and clustering write as
  their result.
- **Lengths.** The tools keep reads of every length, empty ones included;
  discarding reads is [`sqfilt`](docs/sqfilt.md)'s job only.
- **Order.** The tools work on their input in the order given. Where
  vsearch sorts first (clustering, de novo chimera detection), the pipeline
  sorts with [`sortby`](docs/sortby.md).
- **Masking.** The tools that index or align take lower-case residues as
  masked. Masking is [`lcmask`](docs/lcmask.md)'s job only; `lcmask --clear`
  turns masking off.
- **Hits.** The tools that find hits write them in vsearch's UC format;
  other tools turn hits into alignments, tables, taxonomy or selections.
- **Errors.** Invalid input fails with vsearch's checks and exit code 1;
  invalid arguments fail with exit code 2, as clap reports them. Abundances
  (`;size=`) are checked wherever vsearch checks them.
- **Compressed input.** No tool reads compressed files. Decompress into the
  pipe: `gzip -dc reads.fq.gz | sqfilt ...` or `bzip2 -dc reads.fq.bz2 |
  ...`.
- **Threads.** The tools that search, align, classify, merge, mask, or
  orient (`dbsrch`, `allpr`, `clust`, `unoise`, `uchime` with `--db`,
  `sintax`, `mergpa`, `lcmask`, `dbornt`) run on `-t/--threads` threads, one
  per logical core by default. Their output is the same on any number. It is
  vsearch's with `--threads 1`, in input order where vsearch on several
  threads writes queries as they finish. The others run on one thread.
- **Vector instructions.** The aligner, the reader, and reverse complements
  use the best vector instructions the processor has, chosen when the tool
  starts. These are AVX2, or else SSE2, on x86-64 and NEON on aarch64. The
  aligner takes sixteen targets at once with AVX2, where vsearch takes
  eight. `RSEARCH_SIMD` chooses a lower level (`sse2`, `neon`, or `scalar`
  on any machine) and fails with exit code 2 for one the processor lacks.
  The output is the same at every level.

## Tools

### Preprocessing

Records:

| Tool | Its one job | vsearch |
|---|---|---|
| [`conv`](docs/conv.md) | Change the format of records and the encoding of their scores: FASTQ, FASTA, SFF or UDB in; FASTQ or FASTA out | `--fastq_convert`, `--fasta2fastq`, `--sff_convert`, `--udb2fasta`; `--fastaout`, `--fasta_width` |
| [`sqtrim`](docs/sqtrim.md) | Shorten reads: strip their ends and truncate them | the trimming half of `--fastq_filter`, `--fastx_filter` |
| [`sqfilt`](docs/sqfilt.md) | Keep or discard whole reads, or pairs of reads | the filtering half of `--fastq_filter`, `--fastx_filter`; `--minseqlength`, `--maxseqlength`, `--minsize`, `--maxsize`, `--minuniquesize`, `--maxuniquesize`; `--min_unmasked_pct`, `--max_unmasked_pct` |
| [`relab`](docs/relab.md) | Rewrite sequence labels: relabel, add or remove annotations | the label options of every command |
| [`derep`](docs/derep.md) | Group identical sequences and count them | `--fastx_uniques`, `--derep_fulllength`, `--derep_id`, `--derep_prefix`, `--derep_smallmem` |
| [`rerep`](docs/rerep.md) | Write each read as many times as its abundance | `--rereplicate` |
| [`revcomp`](docs/revcomp.md) | Reverse-complement reads | `--fastx_revcomp` |
| [`lcmask`](docs/lcmask.md) | Mask low-complexity regions with DUST, or clear masking | `--fastx_mask`, `--maskfasta`; `--qmask`, `--dbmask`, `--hardmask` |
| [`redig`](docs/redig.md) | Cut sequences at a restriction pattern | `--cut` |
| [`scrmbl`](docs/scrmbl.md) | Shuffle the bases within each read | `--scramble` |
| [`getseq`](docs/getseq.md) | Select records by their headers or by hits | `--fastx_getseq`, `--fastx_getseqs`, the selection half of `--fastx_getsubseq`; `--matched`, `--notmatched`, `--dbmatched`, `--dbnotmatched` |
| [`subseq`](docs/subseq.md) | Keep one region of every read | the region half of `--fastx_getsubseq` |
| [`sortby`](docs/sortby.md) | Sort reads by decreasing length or abundance | `--sortbylength`, `--sortbysize`; the sorting step of clustering and chimera detection |
| [`shufl`](docs/shufl.md) | Shuffle the order of records | `--shuffle` |
| [`subsmp`](docs/subsmp.md) | Draw a random subset of records, or of reads by abundance | `--fastx_subsample` |
| [`topn`](docs/topn.md) | Keep the first N records | `--topn` |

Read pairs:

| Tool | Its one job | vsearch |
|---|---|---|
| [`mergpa`](docs/mergpa.md) | Merge overlapping read pairs | `--fastq_mergepairs` |
| [`joinpa`](docs/joinpa.md) | Join read pairs with a pad | `--fastq_join` |
| [`syncpa`](docs/syncpa.md) | Pair the reads of two files by label | `--fastx_syncpairs` |

Quality reports:

| Tool | Its one job | vsearch |
|---|---|---|
| [`fqchar`](docs/fqchar.md) | Report the residues and score symbols of FASTQ reads, and guess their encoding | `--fastq_chars` |
| [`fqstat`](docs/fqstat.md) | Report the lengths and quality scores of FASTQ reads | `--fastq_stats` |
| [`eestat`](docs/eestat.md) | Report quality scores and expected errors at each position | `--fastq_eestats` |
| [`eegrid`](docs/eegrid.md) | Report the reads kept by length and expected-error cutoffs | `--fastq_eestats2` |

### Analysis

Searching:

| Tool | Its one job | vsearch |
|---|---|---|
| [`dbsrch`](docs/dbsrch.md) | Search queries against a database; write hits | `--usearch_global` |
| [`exsrch`](docs/exsrch.md) | Find exact full-length matches in a database by hashing; write hits | `--search_exact` |
| [`allpr`](docs/allpr.md) | Align every pair of sequences in one file; write hits | `--allpairs_global` |

Clustering:

| Tool | Its one job | vsearch |
|---|---|---|
| [`clust`](docs/clust.md) | Cluster sequences around centroids, in input order; write centroids and hits | `--cluster_fast`, `--cluster_size`, `--cluster_smallmem` |
| [`unoise`](docs/unoise.md) | Denoise amplicons with the UNOISE3 rule; write centroids and hits | `--cluster_unoise` |
| [`clspl`](docs/clspl.md) | Split records into one file per cluster, from clustering hits | `--clusters` |
| [`clmsa`](docs/clmsa.md) | Align each cluster's members to its centroid: alignment, consensus or profile | `--msaout`, `--consout`, `--profile` |

Hits:

| Tool | Its one job | vsearch |
|---|---|---|
| [`hitfmt`](docs/hitfmt.md) | Write hits as alignments, BLAST tabular, user fields, SAM, FASTA pairs or segments | `--alnout`, `--blast6out`, `--userout`, `--samout`, `--fastapairs`, `--qsegout`, `--tsegout` |
| [`otutab`](docs/otutab.md) | Count hits per sample and OTU: classic, mothur or BIOM table | `--otutabout`, `--mothur_shared_out`, `--biomout` |
| [`lcatax`](docs/lcatax.md) | Write the last common ancestor of each query's top hits | `--lcaout` |

Chimeras:

| Tool | Its one job | vsearch |
|---|---|---|
| [`uchime`](docs/uchime.md) | Detect chimeras with UCHIME, de novo or against a database | `--uchime_denovo`, `--uchime2_denovo`, `--uchime3_denovo`, `--uchime_ref` |
| [`mchime`](docs/mchime.md) | Detect chimeras of two or more parents in long exact sequences | `--chimeras_denovo` |

Taxonomy and orientation:

| Tool | Its one job | vsearch |
|---|---|---|
| [`sintax`](docs/sintax.md) | Classify sequences by taxonomy with SINTAX | `--sintax` |
| [`dbornt`](docs/dbornt.md) | Orient sequences by database words; write them forward | `--orient` |

UDB files:

| Tool | Its one job | vsearch |
|---|---|---|
| [`mkudb`](docs/mkudb.md) | Build a UDB database: sequences and word index | `--makeudb_usearch` |
| [`udbinfo`](docs/udbinfo.md) | Report a UDB file's header | `--udbinfo` |
| [`udbstat`](docs/udbstat.md) | Report a UDB file's word statistics | `--udbstats` |

## Pipelines

Each pipeline gives the vsearch command's output, labels kept whole as with
vsearch's `--notrunclabels`. vsearch's own defaults appear as steps: the
`sqfilt` length limits (`--minseqlength 32` and `--maxseqlength 50000` for
the commands that have them), the masking by `lcmask`, the sorting by
`sortby`, and `conv --fasta` for vsearch's wrapped FASTA. `q.fa` are
queries, `db.fa` a database and `hits.uc` the hits.

| vsearch | rsearch |
|---|---|
| `--fastq_filter in.fq --fastq_trunclen 200 --fastq_maxee 1 --relabel R --fastaout out.fa` | `sqtrim -i in.fq --len 200 \| sqfilt --min-len 200 --max-ee 1 \| relab --relabel R \| conv --fasta > out.fa` |
| `--fastq_mergepairs f.fq --reverse r.fq --fastq_maxee 1 --fastqout m.fq` | `mergpa -i f.fq --reverse r.fq \| sqfilt --max-ee 1 > m.fq` |
| `--fastx_uniques in.fa --minuniquesize 2 --relabel U --sizeout --fastaout out.fa` | `sqfilt -i in.fa --max-len 50000 \| derep --no-sizein \| sqfilt --min-size 2 \| relab --relabel U --sizeout \| conv --fasta > out.fa` |
| `--derep_smallmem in.fa --fastaout out.fa` | `derep -i in.fa --unsorted --no-sizein --no-sizeout \| conv --fasta > out.fa` |
| `--sortbysize in.fa --minsize 2 --topn 5 --output out.fa` | `sqfilt -i in.fa --min-size 2 \| sortby --size \| topn 5 \| conv --fasta > out.fa` |
| `--fastx_getsubseq in.fa --label s1 --subseq_start 3 --subseq_end 6 --fastaout out.fa` | `getseq -i in.fa --label s1 \| subseq --start 3 --end 6 \| conv --fasta > out.fa` |
| `--fastx_mask in.fa --min_unmasked_pct 50 --fastaout out.fa` | `lcmask -i in.fa \| sqfilt --min-upper-pct 50 > out.fa` |
| `--cluster_fast in.fa --id 0.97 --sizein --sizeout --centroids c.fa --uc u.uc` | `sqfilt -i in.fa --min-len 32 --max-len 50000 \| lcmask \| sortby --length \| clust --maxrejects 8 --id 0.97 --hits u.uc \| conv --fasta > c.fa` |
| `--cluster_size in.fa --id 0.97 --sizein --sizeout --centroids c.fa` | `sqfilt -i in.fa --min-len 32 --max-len 50000 \| lcmask \| sortby --size \| clust --id 0.97 \| conv --fasta > c.fa` |
| `--cluster_unoise in.fa --sizein --sizeout --centroids c.fa` | `sqfilt -i in.fa --min-len 32 --max-len 50000 --min-size 8 \| lcmask \| sortby --size \| unoise \| conv --fasta > c.fa` |
| `... --consout cons.fa --clusters cl` | `clmsa --hits u.uc -i sorted.fa --consensus > cons.fa` and `clspl --hits u.uc -i sorted.fa --prefix cl`, on the sorted input clustered |
| `--uchime3_denovo in.fa --nonchimeras n.fa` | `sqfilt -i in.fa --max-len 50000 \| lcmask \| sortby --size \| uchime --uchime2 --abskew 16 \| conv --fasta > n.fa` |
| `--usearch_global q.fa --db db.fa --id 0.97 --alnout a.txt` | `lcmask -i db.fa \| sqfilt --min-len 32 --max-len 50000 > dbm.fa`, `lcmask -i q.fa \| dbsrch --db dbm.fa --id 0.97 > hits.uc`, `hitfmt hits.uc -i q.fa --db dbm.fa --alignments > a.txt` |
| `--usearch_global ... --otutabout t.tsv --lcaout l.tsv` | then `otutab hits.uc --db dbm.fa > t.tsv` and `lcatax hits.uc --db dbm.fa > l.tsv` |
| `--usearch_global ... --matched m.fa --dbnotmatched n.fa` | then `getseq --hits hits.uc -i q.fa > m.fa` and `getseq --hits hits.uc --targets --invert -i dbm.fa > n.fa` |
| `--search_exact q.fa --db db.fa --uc u.uc` | `sqfilt -i db.fa --max-len 50000 > dbf.fa`, `exsrch -i q.fa --db dbf.fa --maxhits 1 > u.uc` |
| `--sintax q.fa --db db.fa --randseed 1 --tabbedout t.tsv` | `sqfilt -i db.fa --min-len 32 --max-len 50000 > dbf.fa`, `lcmask --clear -i q.fa \| sintax --db dbf.fa --seed 1 > t.tsv` |
| `--orient q.fa --db db.fa --fastaout o.fa` | `lcmask -i db.fa \| sqfilt --max-len 50000 > dbm.fa`, `dbornt -i q.fa --db dbm.fa \| conv --fasta > o.fa` |
| `--makeudb_usearch in.fa --output db.udb` | `sqfilt -i in.fa --min-len 32 --max-len 50000 \| lcmask \| mkudb -o db.udb` |
| `--udb2fasta db.udb --output out.fa` | `conv -i db.udb --fasta > out.fa` |

The comparison framework's translator prints the pipeline for any vsearch
command line it supports, with the caveats that apply:

```bash
cargo run -p rscompare --release -- translate -- --cluster_fast in.fa --id 0.97 --centroids c.fa
```

The tools run in a pipe, but they are meant to have their intermediate
files saved. Each step's output is a file you can inspect, count or feed
to another step. The counts vsearch sums up in a report, such as the
number of reads truncated by its filter, come from these files.

## Example

The files of [`docs/examples/`](docs/examples) are small inputs for the
examples of every page; run these in that directory. Amplicons with their
abundances are checked for chimeras. The rest are clustered at 97% into
OTUs. Reads from two samples are mapped to the OTUs, which are then
counted per sample:

```console
$ sortby --size -i amplicons.fa | uchime --chimeras chimeras.fa > good.fa
$ grep '>' chimeras.fa
>a5;size=4
$ clust -i good.fa --id 0.97 | relab --relabel OTU --sizeout > otus.fa
$ grep '>' otus.fa
>OTU1;size=62
>OTU2;size=40
>OTU3;size=25
$ dbsrch -i queries.fa --db otus.fa --id 0.97 --strand both > hits.uc
$ otutab hits.uc --db otus.fa
#OTU ID	S1	S2
OTU1	3	0
OTU2	0	2
OTU3	1	1
```

## Comparison

[`compare/README.md`](compare/README.md) describes the comparison
framework. It runs vsearch and the translated rsearch pipeline on random,
malformed and real inputs. It then compares their outputs, exit codes and
errors.

## License

[GNU General Public License, version
3](http://www.gnu.org/licenses/gpl-3.0.html)
