# Comparing rsearch with vsearch

Untracked helpers for checking that rsearch's tools give vsearch's
results. They form the Python package `rscompare`. It translates a
vsearch command line into the rsearch pipeline that should give the same
results and runs both. It then compares their outputs, exit codes, error
messages, and warnings. Standard library only (Python 3.9 or later).

## Build

From the repository root, build vsearch from `../vsearch` into
`compare/vsearch-build/vsearch`, then rsearch's tools into
`target/release`:

```sh
./compare/build-vsearch.sh
cargo build --release
```

`build-vsearch.sh` takes the vsearch source directory and an output
directory as optional arguments. Tested on macOS on Apple silicon.

## Compare

Run from `compare/`. `list` shows the sources, their default seeds, and
what they cover. `run` runs sources at their default seeds, or at those
of `--seeds`, and `summary` shows the outcomes, the coverage, and the
caveats of a results file. `translate` prints the pipeline of a vsearch
command line, and `test` runs the unit tests:

```sh
cd compare
python3 -m rscompare list
python3 -m rscompare run derep clust
python3 -m rscompare run all --out results/run.jsonl
python3 -m rscompare run dbsrch --seeds 1-50,77 -j 8
python3 -m rscompare summary results/run.jsonl
python3 -m rscompare translate -- --cluster_fast in.fa --id 0.97 --centroids c.fa
python3 -m rscompare test
```

`run` prints a `DIFF` line for each failing case (the first 200). It
then prints the tally per source of `pass`, `fail`, `skip`, and `error`
(a case the translator refused, or a harness failure). It exits 1 if any
case failed.

`--vsearch` and `--bin` name the vsearch binary and the rsearch build
(defaults `compare/vsearch-build/vsearch` and `target/release`).
`--jobs` sets the processes (each tool stays single-threaded). `--tmp`
sets the working directory (`$TMPDIR/rscompare`). `--keep` keeps the
files of passing cases too. The files of a failing case are always kept
in `$TMPDIR/rscompare/<source>-<seed>/c<index>/`. vsearch's outputs are
in `v/` and rsearch's in `r/`.

`all` runs every source but `mkudb-long`. Its UDB files reach 4 GB. It
runs on its own, one case at a time:

```sh
python3 -m rscompare run mkudb-long -j 1
```

### Results

`--out` appends one JSON object a case:

| key | what |
|---|---|
| `source`, `seed`, `index`, `case` | where the case comes from |
| `kind` | `translated`, or `probe` for a case that gives its rsearch steps itself |
| `vsearch` | vsearch's argv |
| `rsearch` | the rsearch pipeline, a line a chain of steps |
| `outcome` | `pass`, `fail`, `skip`, or `error`. The last two have a `reason` |
| `diffs` | what differed: exit code, error, a warning, or the first different line of an output |
| `exit`, `stderr` | both programs' exit codes, and their error messages on a failure |
| `failed_step` | the first rsearch step that failed |
| `coverage` | the case's parameters, vsearch's exit code and error, and counts from vsearch's outputs (UC record types, records, table rows) |
| `caveats` | what the translation does not reproduce (below) |
| `seconds`, `wall` | time in vsearch, in the rsearch steps, and the whole case |

`summary` tallies the outcomes, the coverage values, and the caveats of
one or more result files so that two runs can be compared.

## The translator

`rscompare/translate.py` can be used alone:

```python
from rscompare.translate import translate
p = translate(['--fastx_uniques', 'in.fq', '--sizeout', '--fastqout', 'o.fq'], 'work')
p.steps      # [Step(tool, args, stdin, stdout), ...] in the order they run
p.outputs    # [Output(option, vsearch_file, rsearch_file, comparator), ...]
p.caveats    # what the pipeline does not reproduce
print(p.render())
```

It holds, once, every rule the zsh runners encoded by hand. The rules
come from the pipelines of `README.md`, the tool pages in `docs/`, and
the option map of `docs/vsearch.md`. They are:

- option names;
- `sqtrim` then `sqfilt` for the filters;
- `sqfilt` for vsearch's length filters (`--minseqlength`, 32 or 1 by
  command, `--maxseqlength 50000`), `--minsize`, `--minuniquesize`, and
  the unmasked shares;
- `lcmask` for `--qmask`, `--dbmask`, and `--hardmask` by each command's
  rules (searches, `--makeudb_usearch`, `--orient`, `--sintax`);
- `sortby` before `--cluster_fast`, `--cluster_size`,
  `--cluster_unoise`, and de novo chimeras;
- `--maxrejects 8` for `--cluster_fast`;
- `--sizein`/`--sizeout` against rsearch's defaults;
- `topn` and `relab` for the label options, and `conv --fasta` at
  `--fasta_width`;
- `hitfmt`, `getseq --hits`, `otutab`, `lcatax`, `clspl`, and `clmsa` on
  the hits;
- a search without `--maxhits` for `--dbmatched`, one with `--maxhits 1`
  for `--uc` without `--uc_allhits`;
- which rsearch output stands for each vsearch output. `--log` is the
  report of `--fastq_chars`, `--fastq_stats`, `--fastq_mergepairs`,
  `--orient`, `--sintax`, `--udbinfo`, and `--udbstats`, and otherwise a
  log that is not compared.

The steps run one after the other through files. They make a shell
pipeline without the concurrency. As in a pipeline, a failed step's
output still streams into the next step on its standard input. A step
that takes it as a file argument (a database, a hit file) does not run.
The exit code and error are the first failing step's.

An option it cannot translate is a `TranslationError` naming it (for
example `--alnout` with clustering, `--fastq_minlen 0`, compressed
input, and `--clusters` at a width other than 0). A default behaviour of
vsearch that the pipeline does not reproduce is a caveat. Each case
records its caveats:

- `--minseqlength 0` keeps empty sequences. `sqfilt` cannot keep them,
  and `--maxseqlength` is then not applied.
- Without `--notrunclabels`, vsearch cuts labels at spaces. rsearch
  keeps them whole.
- Without `--threads 1`, the multithreaded commands may order their
  outputs otherwise.
- Without `--randseed`, random draws differ.
- `sqfilt` checks `;size=`. Dereplication reads it only with `--sizein`.
- `getseq --targets` counts weak hits, and vsearch's `--dbmatched`
  counts accepted hits only. `clspl` writes the abundances of
  `--clusters` with `--sizeout` as they are. vsearch's `--search_exact`
  writes masked sequences.

## Comparators

`rscompare/comparators.py` holds the normalisations. It names each once:

| comparator | for |
|---|---|
| `bytes` | everything else |
| `uc-fast` | `--cluster_fast`'s `=` for an alignment perfect but for terminal gaps. clust writes the CIGAR in full there |
| `biom` | BIOM without its id, program, and date lines. Those lines name otutab's output, rsearch, and the time |
| `log-report` | a report in vsearch's `--log` without the log's own lines |
| `orient-log`, `sintax-log`, `udbinfo-log`, `udbstats-log` | the parts of `--log` those tools write to `--report` or `-o` |
| `alnout`, `sam` | `--alnout` without its opening lines, SAM without `@PG` and the database's file name |
| `cluster-files` | the files of `--clusters` by number |
| `upper-seq`, `fastq-to-fasta` | vsearch's output upper-cased or in FASTA for two probes |

Errors are compared by kind: `exit` (exit code only), `exact` (vsearch's
`Fatal error: X` against rsearch's `Error: X`), or a mapping of messages
that differ in wording (`abundance`, `joinpa`, `mergpa`, `mkudb`,
`mkudb-exact`, `udb`, `invalid`). The warnings compared are `stripped`,
`missing-abundance`, `nowords`, and `wordlength-adjusted`. rsearch's
exit code 2 for invalid arguments stands for vsearch's 1 where a case
says so.

## Sources

`rscompare/sources/` makes the cases. The generators in `gen/` write the
inputs from a seed as before. Every case is reproducible from its seed.

- `preprocess.py` and `analysis.py`: one source a former zsh runner.
  Each keeps the runner's generator, seeds, and option sets. Where a
  runner applied an rsearch tool to vsearch's own hits (hitfmt, otutab,
  lcatax, getseq `--hits`, clspl), the case is now vsearch's command
  whole. The translated pipeline then runs the search too.
- `fuzz.py`: `fuzz` (malformed FASTA and FASTQ against about 40 vsearch
  commands, from `input.zsh`) and `fuzz-report` (the quality reports,
  from `report-input.zsh`). Where the translated pipeline starts with
  another tool than the one `input.zsh` fuzzed, a probe runs
  `input.zsh`'s own command too.
- `datasets.py`: the real datasets of `datasets.json`.
  `fetch_datasets.py` downloads them into `$RSCOMPARE_DATA` (default
  `$TMPDIR/rsearch-data/<id>/`). Seed N is the manifest's Nth entry. An
  entry's command lines are read from `"vsearch"`, `"commands"`, or
  `"vsearch_commands"` (strings or argv lists, without the program).
  `{dir}` in them stands for the dataset's directory, and relative input
  names are taken from it. Output options are moved to the case's output
  directory. `--threads 1`, `--quiet`, and `--notrunclabels` are added.
  An entry whose files are missing is skipped.

A source is a function `make(T, seed, env)` returning `Case`s
(`rscompare/cases.py`). A case holds vsearch's argv with outputs under
`@OUT@` and the error kind and warnings to compare. It also says whether
outputs are compared when vsearch fails. It may hold a probe or a skip rule.

## The zsh runners

The 45 zsh runners and `input.zsh`/`report-input.zsh` are retired to
`retired-zsh/`. Each source above replaces the runner of the same name
(`getseq-hits.zsh` is `getseq-hits`, `conv-udb.zsh` `conv-udb`,
`input.zsh` `fuzz`, and `report-input.zsh` `fuzz-report`). It passed the
seeds the runner passed (see below). `gen/merge.py` and
`gen/vsearch_log_report.py` served only the runners.

The port's first full run (2026-11-20, rsearch `26fdb8c`, vsearch
`e1551256`; `results/final-2026-11-20.jsonl` and its `.summary.txt`). It
ran every source at the runners' seeds. 71,570 cases passed and none
failed. The runners' own rules skipped 18 (14 `conv-udb` files whose
magic number is spoiled, 4 `udbstat` files on which vsearch aborts). A
case is a vsearch command line. A runner counted each compared file. The
counts therefore differ from the runners' `pass=` totals.

## Time two builds

```sh
python3 compare/time.py BASELINE_DIR target/release
```

Runs the two builds alternately on the benchmark inputs in `target/tmp`
(written by `cargo bench`) and prints the median time of each. Build the
baseline from another commit with, for example,
`git archive HEAD | tar -x -C /tmp/base` and
`CARGO_TARGET_DIR=/tmp/base-target cargo build --release` there.
