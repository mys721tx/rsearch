"""Compressed input: the generators' files gzipped or bzip2'd, against
vsearch reading them itself (it knows them by their magic numbers) and
the translated pipeline decompressing them with gzip -dc or bzip2 -dc.
Malformed records inside the compressed files, and files cut short, are
compared too."""

import bz2
import gzip
import os
import random

from ..cases import Case
from . import argv, gen, read, source


def _pack(path: str, kind: str) -> str:
    with open(path, 'rb') as f:
        data = f.read()
    if kind == 'gzip':
        out = path + '.gz'
        with gzip.open(out, 'wb') as f:
            f.write(data)
    else:
        out = path + '.bz2'
        with bz2.open(out, 'wb') as f:
            f.write(data)
    return out


@source('compressed', list(range(1, 101)), 'gzip and bzip2 input: streaming '
        'and whole-file tools, pairs, databases, malformed records, files cut '
        'short')
def compressed(T, seed, env):
    r = random.Random(seed)
    kind = lambda: r.choice(['gzip', 'bzip2'])
    dirs = {}
    for g in ('derep', 'mergpa', 'dbsrch', 'fuzz'):
        d = os.path.join(T, g)
        os.makedirs(d, exist_ok=True)
        gen(g, d, seed)
        dirs[g] = d
    D, M, S, F = dirs['derep'], dirs['mergpa'], dirs['dbsrch'], dirs['fuzz']
    q, a = _pack(f'{D}/q.fq', kind()), _pack(f'{D}/a.fa', kind())
    out = [
        Case('uniques fq', argv(f'--fastx_uniques {q} --quiet --sizein --sizeout '
                                '--minseqlength 0 --notrunclabels --fastqout {O}/o')),
        Case('filter fq', argv(f'--fastx_filter {q} --quiet --notrunclabels '
                               '--fastq_stripleft 2 --fastq_maxee 1 --fastqout {O}/k '
                               '--fastqout_discarded {O}/d'),
             outputs_on_error=False),
        Case('revcomp fa', argv(f'--fastx_revcomp {a} --quiet --notrunclabels '
                                '--fasta_width 0 --fastaout {O}/o')),
        Case('sortbylength fa', argv(f'--sortbylength {a} --quiet --notrunclabels '
                                     '--minseqlength 0 --fasta_width 0 '
                                     '--output {O}/o')),
        Case('fasta at 80', argv(f'--fastx_revcomp {q} --quiet --notrunclabels '
                                 '--fastaout {O}/o')),
        Case('chars', argv(f'--fastq_chars {q} --quiet --log {{O}}/log')),
    ]
    f, rv = _pack(f'{M}/f.fq', kind()), _pack(f'{M}/r.fq', kind())
    out.append(Case('mergepairs', argv(
        f'--fastq_mergepairs {f} --reverse {rv} --quiet --fastqout {{O}}/m '
        '--fastqout_notmerged_fwd {O}/n --eetabbedout {O}/e --log {O}/log'),
        errors='mergpa', outputs_on_error=False))
    sq, sdb = _pack(f'{S}/q.fa', kind()), _pack(f'{S}/db.fa', kind())
    out.append(Case('search', argv(
        f'--usearch_global {sq} --db {sdb} --id 0.7 --uc {{O}}/uc --uc_allhits '
        '--quiet --notrunclabels --qmask soft --dbmask soft --minseqlength 0 '
        '--maxseqlength 100000 --threads 1'), success_only=['uc']))
    # malformed records inside compressed files
    ff, fa = _pack(f'{F}/f.fq', kind()), _pack(f'{F}/f.fa', kind())
    out += [
        Case('fuzz filter fq', argv(f'--fastx_filter {ff} --quiet --notrunclabels '
                                    '--fastqout {O}/o'),
             errors='exact', warnings=['stripped']),
        Case('fuzz revcomp fa', argv(f'--fastx_revcomp {fa} --quiet --notrunclabels '
                                     '--fasta_width 0 --fastaout {O}/o'),
             errors='exact', warnings=['stripped']),
        Case('fuzz uniques fa', argv(f'--fastx_uniques {fa} --quiet --notrunclabels '
                                     '--minseqlength 0 --sizein --sizeout '
                                     '--fasta_width 0 --fastaout {O}/o'),
             errors='exact', warnings=['stripped']),
    ]
    # a compressed file cut short: the exit code only, as what comes out
    # before the error is the decompressor's (macOS gzip -dc writes nothing,
    # bzip2 -dc what it decoded)
    data = open(q, 'rb').read()
    cut = f'{q}.cut{os.path.splitext(q)[1]}'
    with open(cut, 'wb') as h:
        h.write(data[:r.randint(1, max(1, len(data) - 1))])
    out.append(Case('cut short', argv(f'--fastx_revcomp {cut} --quiet '
                                      '--notrunclabels --fastqout {O}/o'),
                    outputs_on_error=False))
    return out
