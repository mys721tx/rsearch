"""Malformed FASTA and FASTQ (compare/gen/fuzz.py) against the vsearch
commands, ported from compare/input.zsh and compare/report-input.zsh.

Each check compares the output, the exit code, the error message
(vsearch's "Fatal error:" with rsearch's "Error:") and the warning about
stripped characters. Every check runs the translated pipeline; where that
pipeline starts with another tool than the one input.zsh fuzzed (conv
--fasta, relab, topn, subseq, mkudb, dbornt, allpr), input.zsh's command
also runs, as a probe named with "probe"."""

import os

from ..cases import Case
from ..translate import Step
from . import argv, gen, seeds, source


def _files(T):
    with open(f'{T}/odb.fa', 'w') as f:
        f.write('>d1\nACGTACGTNACGTTGCAAGTCCGT\n>d2\nTTGACCAGTAGGCATTAGCAGGTA\n')
    with open(f'{T}/oq.fa', 'w') as f:
        f.write('>q1\nACGTTGCAAGTC\n>q2\nTACCTGCTAATG\n>q3\nGGGGGGG\n')
    with open(f'{T}/sdb.fa', 'w') as f:
        f.write('>t1;tax=d:A,g:B\nACGTACGTNACGTTGCAAGTCCGTAACCGGTTACGATCGATCGGCTA\n'
                '>t2;tax=d:C\nTTGACCAGTAGGCATTAGCAGGTAACGGATCCATGCAAGGCTTGCA\n')
    with open(f'{T}/sq.fa', 'w') as f:
        f.write('>q1\nACGTTGCAAGTCCGTAACCGGTTACGATCGATCG\n'
                '>q2\nTACCTGCTAATGCCTACTGGTCAAGGCAGGTA\n')
    with open(f'{T}/edb.fa', 'w') as f:
        f.write('>e1\nACGT\n>e2\nacgtn\n>e3\nA\n')
    with open(f'{T}/eq.fa', 'w') as f:
        f.write('>q1\nACGT\n>q2\nA\n>q3\nACGTN\n')
    with open(f'{T}/ddb.fa', 'w') as f:
        f.write('>d1\nACGTACGTTGCAACGTAGCT\n>d2\nacgtACGTTGCAACGTAGCA\n'
                '>d3\nTTGCAACG\n')
    with open(f'{T}/dq.fa', 'w') as f:
        f.write('>p1\nACGTACGTTGCAACGTAGCT\n>p2\nAGCTACGTTGCAACGTACGT\n>p3\n\n')
    with open(f'{T}/udb.fa', 'w') as f:
        f.write('>u1;size=9\nACGTACGTTGCAACGTAGCTAGGCTAACGTTAGCAT\n'
                '>u2;size=8\nTTGCAACGGTCATTGACCATGGTACCAGTTAACGGA\n')


def _probe(*steps, out='out'):
    """input.zsh's rsearch command: steps of (tool, args), the first on its
    -i, each into the next; the last one's output stands for `out`."""
    def make(rdir):
        made, prev = [], None
        for k, (tool, args) in enumerate(steps):
            dst = os.path.join(rdir, f'p{k}')
            made.append(Step(tool, list(args), prev, dst))
            prev = dst
        return made, {out: prev}
    return make


@source('fuzz', seeds(1, 400), 'input checks of every tool on malformed '
        'FASTA and FASTQ (compare/input.zsh)')
def fuzz(T, seed, env):
    gen('fuzz', T, seed)
    _files(T)
    out = []
    Q = '--quiet --notrunclabels'

    def chk(name, v, probe=None, ok_on_success=False, errors='exact',
            probe_errors=None, out_opt=None):
        a = argv(v, T=T)
        kw = dict(warnings=['stripped'], outputs_on_error=not ok_on_success)
        # where the pipeline streams into a tool other than the one fuzzed,
        # what it writes before an error is the probe's to check
        out.append(Case(name, a, errors=errors,
                        **dict(kw, outputs_on_error=kw['outputs_on_error']
                               and not probe)))
        if probe:
            out.append(Case(name + ' probe', a, errors=probe_errors or errors,
                            probe=probe, **kw))

    for f in ('fq', 'fa'):
        F = f'{{T}}/f.{f}'
        o = '--fastqout' if f == 'fq' else '--fasta_width 0 --fastaout'
        oo = o.split()[-1].lstrip('-')
        chk(f'conv-{f}', ('--fastq_convert ' if f == 'fq' else '--fasta2fastq ')
            + F + ' --quiet --fastqout {O}/out',
            probe=_probe(('relab', ['-i', f'{T}/f.{f}']),
                         *([('conv', [])] if f == 'fa' else []), out='fastqout'))
        out[-1].name = f'relab-{f} probe'
        chk(f'conv-fasta-{f}', f'--fastx_subsample {F} --quiet --sample_pct 100 '
            f'--minseqlength 0 --notrunclabels --fastaout {{O}}/out',
            probe=_probe(('conv', ['-i', f'{T}/f.{f}', '--fasta']), out='fastaout'),
            ok_on_success=True)
        chk(f'topn-{f}', f'--sortbysize {F} {Q} --minseqlength 0 --fasta_width 0 '
            '--topn 1000000 --output {O}/out')
        out.append(Case(f'topn-{f} probe', argv(
            f'--fastx_subsample {F} {Q} --sample_pct 100 --minseqlength 0 {o} '
            '{O}/out', T=T), errors='exact', warnings=['stripped'],
            outputs_on_error=False,
            probe=_probe(('topn', ['-i', f'{T}/f.{f}', '1000000']), out=oo)))
        chk(f'derep-{f}', f'--fastx_uniques {F} --quiet --minseqlength 0 --notrunclabels '
            f'--sizein --sizeout {o} {{O}}/out')
        chk(f'sqfilt-{f}', f'--fastx_filter {F} {Q} {o} {{O}}/out')
        chk(f'revcomp-{f}', f'--fastx_revcomp {F} --quiet {o} {{O}}/out')
        so = '--fastqout' if f == 'fq' else '--fasta_width 0 --fastaout'
        chk(f'syncpa-fwd-{f}', f'--fastx_syncpairs {F} --reverse /dev/null --quiet '
            f'{so}_orphans {{O}}/out')
        chk(f'syncpa-rev-{f}', f'--fastx_syncpairs /dev/null --reverse {F} --quiet '
            f'{so}_orphans_rev {{O}}/out')
        chk(f'lcmask-{f}', f'--fastx_mask {F} {Q} --minseqlength 0 {o} {{O}}/out',
            ok_on_success=True)
        chk(f'getseq-{f}', f'--fastx_getseq {F} {Q} --label_substr_match --label a '
            f'{o} {{O}}/out')
        chk(f'subseq-{f}', f"--fastx_getsubseq {F} {Q} --label_substr_match --label '' "
            f'--subseq_start 2 --subseq_end 9 {o} {{O}}/out',
            probe=_probe(('subseq', ['--start', '2', '--end', '9', '-i',
                                     f'{T}/f.{f}']),
                         ('getseq', ['--substr', '--label', '']), out=oo))
        chk(f'shufl-{f}', f'--shuffle {F} {Q} --minseqlength 0 --fasta_width 0 '
            '--randseed 7 --output {O}/out')
        chk(f'subsmp-{f}', f'--fastx_subsample {F} {Q} --minseqlength 0 --randseed 3 '
            f'--sample_pct 60 {o} {{O}}/out')
        chk(f'scrmbl-{f}', f'--scramble {F} {Q} --randseed 5 --scramble_kmer 3 '
            f'{o} {{O}}/out')
        # mkudb reads everything before writing; vsearch upper-cases, so
        # the probe compares the UDB file of the input through lcmask --clear
        mk = lambda rdir, f=f: ([
            Step('mkudb', ['-i', f'{T}/f.{f}', '--word-len', '3'], None,
                 os.path.join(rdir, 'p0')),
            Step('lcmask', ['--clear', '-i', f'{T}/f.{f}'], None,
                 os.path.join(rdir, 'p1')),
            Step('mkudb', ['--word-len', '3'], os.path.join(rdir, 'p1'),
                 os.path.join(rdir, 'p2'))], {'output': os.path.join(rdir, 'p2')})
        chk(f'mkudb-{f}', f'--makeudb_usearch {F} {Q} --minseqlength 0 --dbmask none '
            '--wordlength 3 --output {O}/out', probe=mk, ok_on_success=True,
            errors='mkudb-exact')
        chk(f'dbornt-{f}', f'--orient {F} --db {{T}}/odb.fa --wordlength 3 {Q} {o} '
            '{O}/out', probe=_probe(('dbornt', ['-i', f'{T}/f.{f}', '--db',
                                                f'{T}/odb.fa', '--word-len', '3']),
                                    out=oo))
        chk(f'dbornt-db-{f}', f'--orient {{T}}/oq.fa --db {F} --dbmask soft '
            f'--wordlength 3 {Q} --minseqlength 1 --tabbedout {{O}}/out',
            probe=lambda rdir, f=f: ([Step('dbornt', [
                '-i', f'{T}/oq.fa', '--db', f'{T}/f.{f}', '--word-len', '3',
                '--tabbedout', os.path.join(rdir, 'tab')], None,
                os.path.join(rdir, 'p0'))], {'tabbedout': os.path.join(rdir, 'tab')}))
        chk(f'sintax-{f}', f'--sintax {F} --db {{T}}/sdb.fa --wordlength 3 --randseed 3 '
            f'--threads 1 {Q} --minseqlength 0 --tabbedout {{O}}/out')
        chk(f'sintax-db-{f}', f'--sintax {{T}}/sq.fa --db {F} --dbmask soft '
            f'--wordlength 3 --randseed 3 --threads 1 {Q} --minseqlength 0 '
            '--sintax_random --tabbedout {O}/out')
        chk(f'exsrch-{f}', f'--search_exact {F} --db {{T}}/edb.fa --uc {{O}}/out '
            f'--uc_allhits --strand both {Q} --threads 1', ok_on_success=True)
        chk(f'exsrch-db-{f}', f'--search_exact {{T}}/eq.fa --db {F} --uc {{O}}/out '
            f'--uc_allhits {Q} --threads 1 --minseqlength 0')
        chk(f'allpr-{f}', f'--allpairs_global {F} --acceptall {Q} --qmask none '
            '--minseqlength 0 --threads 1 --uc {O}/out',
            probe=_probe(('allpr', ['-i', f'{T}/f.{f}', '--accept-all']), out='uc'))
        chk(f'dbsrch-{f}', f'--usearch_global {F} --db {{T}}/ddb.fa --id 0.5 '
            '--wordlength 3 --maxaccepts 0 --maxrejects 0 --strand both --uc {O}/out '
            f'--uc_allhits --qmask soft --dbmask soft {Q} --minseqlength 0 --threads 1',
            ok_on_success=True)
        chk(f'dbsrch-db-{f}', f'--usearch_global {{T}}/dq.fa --db {F} --id 0.5 '
            '--wordlength 3 --maxaccepts 4 --uc {O}/out --uc_allhits --qmask soft '
            f'--dbmask soft {Q} --minseqlength 0 --maxseqlength 1000000000 --threads 1')
        big = '--minseqlength 0 --maxseqlength 1000000000 --threads 1'
        chk(f'clust-{f}', f'--cluster_smallmem {F} --usersort --id 0.8 --wordlength 3 '
            f'--strand both --uc {{O}}/out --sizein {Q} --qmask soft {big}')
        chk(f'unoise-{f}', f'--cluster_unoise {F} --minsize 1 --wordlength 3 '
            f'--uc {{O}}/out --sizein {Q} --qmask soft {big}')
        chk(f'uchime-{f}', f'--uchime_denovo {F} --uchimeout {{O}}/out {Q} '
            f'--qmask soft {big}')
        chk(f'mchime-{f}', f'--chimeras_denovo {F} --tabbedout {{O}}/out {Q} '
            f'--qmask soft {big}')
    F = '{T}/f'
    chk('rerep-fa', f'--rereplicate {F}.fa {Q} --sizeout --fasta_width 0 --output {{O}}/out')
    chk('joinpa-fq', f'--fastq_join {F}.fq --reverse {F}.fq --quiet --fastqout {{O}}/out')
    chk('mergpa-fq', f'--fastq_mergepairs {F}.fq --reverse {F}.fq --quiet '
        '--fastqout {O}/out', ok_on_success=True)
    chk('redig-fa', f'--cut {F}.fa --cut_pattern G^AATT_C --quiet --fasta_width 0 '
        '--fastaout {O}/out')
    chk('sortby-fq', f'--sortbylength {F}.fq {Q} --minseqlength 0 --fasta_width 0 '
        '--output {O}/out')
    chk('sortby-fa', f'--sortbysize {F}.fa {Q} --minseqlength 0 --fasta_width 0 '
        '--output {O}/out')
    chk('uchime-ref-fa', f'--uchime_ref {F}.fa --db {{T}}/udb.fa --uchimeout {{O}}/out '
        f'{Q} --qmask soft --dbmask soft --minseqlength 0 --threads 1',
        ok_on_success=True)
    return out


@source('fuzz-report', seeds(1, 400), 'the quality reports on malformed FASTQ '
        '(compare/report-input.zsh)')
def fuzz_report(T, seed, env):
    gen('fuzz', T, seed)
    return [Case(cmd, argv(f'--{cmd} {T}/f.fq --quiet {o} {{O}}/out'),
                 errors='exact', success_only=[o.lstrip('-')])
            for cmd, o in (('fastq_chars', '--log'), ('fastq_stats', '--log'),
                           ('fastq_eestats', '--output'),
                           ('fastq_eestats2', '--output'))]
