"""Cases for the preprocessing tools, ported from the zsh runners
(compare/<tool>.zsh, now retired) with the same generators and seeds."""

import os
import random

from ..cases import OUT, Case
from ..translate import Step
from . import argv, gen, read, seeds, source


def fx(name: str):
    """vsearch's output option for a file's format: (records, discarded)."""
    if name.endswith('.fq'):
        return '--fastqout', '--fastqout_discarded', ''
    return '--fastaout', '--fastaout_discarded', '--fasta_width 0'


@source('derep', seeds(1, 5), 'qualities, encodings, uc, tabbedout, '
        '--same-label, --prefix')
def derep(T, seed, env):
    gen('derep', T, seed)
    C = '--quiet --minseqlength 0 --notrunclabels --fasta_width 0 --sizein --sizeout'
    out = []
    for name, f, v in [
            ('default', 'q.fq', ''), ('strand', 'q.fq', '--strand both'),
            ('ascii64', 'q64.fq', '--fastq_ascii 64'),
            ('ascii64-out33', 'q64.fq', '--fastq_ascii 64 --fastq_asciiout 33'),
            ('ascii33-out64', 'q.fq', '--fastq_asciiout 64'),
            ('qout-range', 'q.fq', '--fastq_qminout 5 --fastq_qmaxout 30'),
            ('qout-range-strand', 'q64.fq', '--fastq_ascii 64 --fastq_qminout 3 '
             '--fastq_qmaxout 35 --strand both'),
            ('best', 'q.fq', '--fastq_qout_max'),
            ('best-strand', 'q.fq', '--fastq_qout_max --strand both '
             '--fastq_qmaxout 39')]:
        out.append(Case(name, argv(f'--fastx_uniques {T}/{f} {C} --fastqout '
                                   '{O}/out ' + v)))
    for name, f, v in [('uc', 'a.fa', ''), ('uc-strand', 'a.fa', '--strand both'),
                       ('uc-fq', 'q.fq', '')]:
        o = fx(f)[0]
        out.append(Case(name, argv(f'--fastx_uniques {T}/{f} {C} {o} {{O}}/out '
                                   '--uc {O}/uc --tabbedout {O}/tab ' + v)))
    out += [
        Case('unsized', argv(f'--fastx_uniques {T}/a.fa --quiet --minseqlength 0 '
                             '--notrunclabels --fasta_width 0 --strand both '
                             '--fastaout {O}/out --uc {O}/uc --tabbedout {O}/tab')),
        Case('derep_id', argv(f'--derep_id {T}/a.fa {C} --output {{O}}/out '
                              '--uc {O}/uc')),
        Case('derep_id-strand', argv(f'--derep_id {T}/a.fa {C} --strand both '
                                     '--output {O}/out')),
        Case('prefix', argv(f'--derep_prefix {T}/p.fa {C} --output {{O}}/out '
                            '--uc {O}/uc')),
        Case('prefix-nosizein', argv(
            f'--derep_prefix {T}/p.fa --quiet --minseqlength 0 --notrunclabels '
            '--fasta_width 0 --sizeout --output {O}/out --uc {O}/uc')),
        Case('smallmem', argv(f'--derep_smallmem {T}/a.fa {C} --fastaout {{O}}/out')),
    ]
    return out


@source('conv', seeds(1, 3), 'FASTQ and FASTA: encodings, Solexa, ranges, '
        'labels, errors')
def conv(T, seed, env):
    gen('conv', T, seed)
    out = []
    for f, v in [('q33.fq', ''), ('q33.fq', '--fastq_asciiout 64 --fastq_qmaxout 62'),
                 ('q33.fq', '--fastq_qminout 3 --fastq_qmaxout 40'),
                 ('q33.fq', '--fastq_qmax 41'), ('q33.fq', '--fastq_qmin 2'),
                 ('q64.fq', '--fastq_ascii 64'),
                 ('q64.fq', '--fastq_ascii 64 --fastq_asciiout 64'),
                 ('q64.fq', '--fastq_ascii 64 --fastq_qminout 2 --fastq_qmaxout 41'),
                 ('sol.fq', '--fastq_ascii 64 --fastq_solexa'),
                 ('sol.fq', '--fastq_ascii 64 --fastq_solexa --fastq_asciiout 64'),
                 ('sol.fq', '--fastq_ascii 64 --fastq_solexa --fastq_qmin 0'),
                 ('sol.fq', '--fastq_ascii 64'), ('a.fa', ''),
                 ('a.fa', '--fastq_asciiout 64 --fastq_qmaxout 50'),
                 ('a.fa', '--fastq_qmaxout 93'), ('bad.fq', ''),
                 ('bad.fq', '--fastq_asciiout 64 --fastq_qmaxout 62'),
                 ('bad.fa', '')]:
        cmd = '--fasta2fastq' if f.endswith('.fa') else '--fastq_convert'
        out.append(Case(f'{f} {v}', argv(f'{cmd} {T}/{f} --quiet --fastqout '
                                         '{O}/out ' + v)))
    # FASTA output, against --fastx_subsample, which writes every read
    for f, v in [('q33.fq', ''), ('q33.fq', '--fasta_width 0'),
                 ('q33.fq', '--fasta_width 7'), ('q64.fq', ''), ('a.fa', ''),
                 ('a.fa', '--fasta_width 0'), ('a.fa', '--fasta_width 1')]:
        out.append(Case(f'fasta {f} {v}', argv(
            f'--fastx_subsample {T}/{f} --quiet --sample_pct 100 --minseqlength 0 '
            '--notrunclabels --fastaout {O}/out ' + v)))
    return out


@source('sff', seeds(1, 60), 'SFF: clipping, encodings, labels, truncated files')
def sff(T, seed, env):
    gen('sff', T, seed)
    out = []
    for v in ['', '--sff_clip', '--fastq_asciiout 64 --fastq_qmaxout 50 --fastq_qminout 3']:
        out.append(Case(f'sff {v}', argv(f'--sff_convert {T}/r.sff --quiet '
                                         '--fastqout {O}/out ' + v)))
    # FASTA: the bases --sff_convert writes, against conv --fasta
    for clip in (False, True):
        def probe(rdir, clip=clip):
            o = os.path.join(rdir, 'fasta')
            return ([Step('conv', ['-i', f'{T}/r.sff', '--fasta'] +
                          (['--clip'] if clip else []), None, o)],
                    {'fastqout': o})
        out.append(Case(f'fasta clip={clip}', argv(
            f'--sff_convert {T}/r.sff --quiet --fastqout {{O}}/out' +
            (' --sff_clip' if clip else '')), probe=probe,
            comparators={'fastqout': 'fastq-to-fasta'}))
    data = open(f'{T}/r.sff', 'rb').read()
    for cut in (1, 2, 3):
        n = random.Random(seed * 10 + cut).randint(0, len(data))
        with open(f'{T}/t{cut}.sff', 'wb') as f:
            f.write(data[:n])
        out.append(Case(f'cut at {n}', argv(f'--sff_convert {T}/t{cut}.sff '
                                            '--quiet --fastqout {O}/out')))
    return out


@source('relab', seeds(1, 30), 'every label option on awkward headers, '
        'abundance errors, ee')
def relab(T, seed, env):
    gen('relab', T, seed)
    opts = ["", "--relabel R", "--relabel '' --relabel_keep --sizeout",
            "--relabel_self --relabel_keep", "--relabel_sha1 --lengthout",
            "--relabel_md5 --sizeout --sample S", "--sizeout", "--xsize",
            "--xee", "--lengthout", "--xlength", "--xsize --xee --xlength",
            "--sample S1 --label_suffix ';'",
            "--label_suffix ';k=v;' --sizeout --lengthout",
            "--label_suffix _s --xsize --sizeout", "--relabel_keep --xsize"]
    out = []

    def c(f, v):
        o, cmd, w = '--fastqout', '--fastx_filter', ''
        if f.endswith('.fa'):
            o, w = '--fastaout', '--fasta_width 0'
        elif f == 'e.fq':
            cmd = '--fastq_convert'
        out.append(Case(f'{f} {v}', argv(f'{cmd} {T}/{f} --quiet {o} {{O}}/out '
                                         f'{w} ' + v.replace('{', '{{'))))
    for f in ('q33.fq', 'a.fa', 'e.fq', 'bad.fa'):
        for v in opts:
            c(f, v)
    for f, v in [('q33.fq', '--eeout'),
                 ('q33.fq', '--eeout --sizeout --lengthout --sample S'),
                 ('q33.fq', "--eeout --relabel R --label_suffix ';'"),
                 ('q33.fq', '--eeout --xee'), ('q64.fq', '--eeout --fastq_ascii 64'),
                 ('q64.fq', '--eeout --fastq_ascii 64 --xsize --lengthout'),
                 ('q33.fq', '--eeout --fastq_ascii 64')]:
        c(f, v)
    return out


def _kept_discarded(T, f, cmd, v, **kw):
    o, d, w = fx(f)
    return Case(f'{f} {v}', argv(f'{cmd} {T}/{f} --quiet --notrunclabels {o} '
                                 f'{{O}}/k {d} {{O}}/x {w} ' + v), **kw)


@source('sqtrim', seeds(1, 20), "--fastx_filter's kept and discarded reads")
def sqtrim(T, seed, env):
    gen('sqtrim', T, seed)
    out = []
    c = lambda f, v: out.append(_kept_discarded(
        T, f, '--fastx_filter', v, outputs_on_error=False))
    for f in ('q33.fq', 'a.fa'):
        for v in ['', '--fastq_stripleft 3', '--fastq_stripright 7',
                  '--fastq_stripleft 20 --fastq_stripright 30',
                  '--fastq_trunclen_keep 25',
                  '--fastq_stripleft 2 --fastq_trunclen_keep 10']:
            c(f, v)
    for q in (0, 2, 15, 25, 30):
        c('q33.fq', f'--fastq_truncqual {q}')
    for e in ('0', '0.01', '0.1', '0.5', '1', '2.5'):
        c('q33.fq', f'--fastq_truncee {e}')
    for e in ('0', '0.001', '0.005', '0.02', '0.1'):
        c('q33.fq', f'--fastq_truncee_rate {e}')
    c('q33.fq', '--fastq_stripleft 1 --fastq_stripright 2 --fastq_trunclen_keep 30 '
      '--fastq_truncqual 12 --fastq_truncee 0.8 --fastq_truncee_rate 0.03')
    c('q64.fq', '--fastq_ascii 64 --fastq_truncqual 15 --fastq_truncee 1')
    c('q64.fq', '--fastq_truncqual 15')
    c('hi.fq', '--fastq_qmax 41')
    c('hi.fq', '--fastq_qmax 41 --fastq_truncqual 20')
    c('hi.fq', '--fastq_qmax 41 --fastq_trunclen_keep 5')
    c('hi.fq', '--fastq_qmax 60 --fastq_truncee 0.05')
    c('hi.fq', '--fastq_qmin 5 --fastq_stripleft 4')
    c('bad.fq', '')
    c('bad.fq', '--fastq_stripleft 2 --fastq_truncqual 15')
    return out


@source('sqfilt', seeds(1, 20), 'every filter, pairs, and --fastx_mask '
        'unmasked shares')
def sqfilt(T, seed, env):
    gen('sqfilt', T, seed)
    out = []
    # rsearch exits 2 on invalid arguments, as clap does, where vsearch exits 1
    kw = dict(outputs_on_error=False, exit_map={(1, 2): 1})
    c = lambda f, v: out.append(_kept_discarded(T, f, '--fastx_filter', v, **kw))
    for f in ('q33.fq', 'a.fa'):
        for v in ['', '--fastq_minlen 20', '--fastq_maxlen 25', '--fastq_maxns 0',
                  '--fastq_maxns 2', '--minsize 3 --maxsize 7',
                  '--fastq_minlen 10 --fastq_maxlen 30 --fastq_maxns 1 --minsize 2']:
            c(f, v)
    for e in ('0.01', '0.1', '0.5', '1', '3'):
        c('q33.fq', f'--fastq_maxee {e}')
    for e in ('0', '0.001', '0.005', '0.02'):
        c('q33.fq', f'--fastq_maxee_rate {e}')
    for q in (0, 5, 15, 25):
        c('q33.fq', f'--fastq_minqual {q}')
    c('q64.fq', '--fastq_ascii 64 --fastq_maxee 0.5 --fastq_minqual 10')
    c('q64.fq', '--fastq_maxee 0.5')
    c('hi.fq', '--fastq_qmax 41')
    c('hi.fq', '--fastq_qmax 60 --fastq_maxee 0.01')
    c('hi.fq', '--fastq_qmin 3')
    for v in ['', '--fastq_maxee 0.5', '--fastq_minlen 15 --fastq_maxns 1',
              '--fastq_minqual 10 --minsize 3']:
        out.append(Case(f'pairs {v}', argv(
            f'--fastx_filter {T}/f.fq --reverse {T}/r.fq --quiet --notrunclabels '
            '--fastqout {O}/k --fastqout_discarded {O}/x --fastqout_rev {O}/k2 '
            '--fastqout_discarded_rev {O}/x2 ' + v), **kw))

    def m(f, v):
        o = '--fastqout' if f.endswith('.fq') else '--fasta_width 0 --fastaout'
        out.append(Case(f'mask {f} {v}', argv(
            f'--fastx_mask {T}/{f} --quiet --notrunclabels --qmask soft {o} '
            '{O}/k ' + v), **kw))
    for x in ('0', '33.3', '50', '80', '100'):
        m('a.fa', f'--min_unmasked_pct {x}')
        m('a.fa', f'--max_unmasked_pct {x}')
        m('q33.fq', f'--min_unmasked_pct {x} --max_unmasked_pct 95')
        m('h.fa', f'--hardmask --min_unmasked_pct {x}')
        m('h.fa', f'--hardmask --max_unmasked_pct {x}')
    return out


@source('rerep', seeds(1, 30), 'FASTA, --no-sizeout, and FASTQ through conv')
def rerep(T, seed, env):
    gen('rerep', T, seed)
    C = '--quiet --notrunclabels --fasta_width 0 --output {O}/out'
    w = ['missing-abundance']
    out = [Case('sizeout', argv(f'--rereplicate {T}/a.fa {C} --sizeout'), warnings=w),
           Case('no-sizeout', argv(f'--rereplicate {T}/a.fa {C}'), warnings=w)]
    # FASTQ, which vsearch does not read: rerep | conv --fasta against
    # --rereplicate on the same reads in FASTA
    import subprocess
    subprocess.run([os.path.join(env.bindir, 'conv'), '-i', f'{T}/q.fq',
                    '--fasta', '--fasta-width', '0', '-o', f'{T}/q.fa'],
                   stderr=subprocess.DEVNULL)

    def probe(rdir):
        t, o = os.path.join(rdir, 'tmp1'), os.path.join(rdir, 'out')
        return ([Step('rerep', ['-i', f'{T}/q.fq'], None, t),
                 Step('conv', ['--fasta', '--fasta-width', '0'], t, o)],
                {'output': o})
    out.append(Case('fastq', argv(f'--rereplicate {T}/q.fa {C} --sizeout'),
                    warnings=w, probe=probe))
    return out


@source('topn', seeds(1, 5), 'derep | topn N against --fastx_uniques --topn N')
def topn(T, seed, env):
    gen('derep', T, seed)
    C = '--quiet --sizein --sizeout --notrunclabels --minseqlength 0 --fasta_width 0'
    out = []
    for f in ('a.fa', 'q.fq'):
        o = fx(f)[0]
        for n in (1, 2, 3, 5, 8, 13, 50, 1000):
            for s in ('', '--strand both'):
                out.append(Case(f'{f} topn {n} {s}', argv(
                    f'--fastx_uniques {T}/{f} {C} {o} {{O}}/out --topn {n} {s}')))
    return out


@source('revcomp', seeds(1, 20), 'every nucleotide code in both cases')
def revcomp(T, seed, env):
    gen('revcomp', T, seed)
    return [Case(f, argv(f'--fastx_revcomp {T}/{f} --quiet --notrunclabels '
                         f'--fasta_width 0 {fx(f)[0]} {{O}}/out'))
            for f in ('a.fa', 'q.fq')]


@source('sortby', seeds(1, 20), 'ties of length, abundance and header; FASTQ')
def sortby(T, seed, env):
    gen('sortby', T, seed)
    return [Case(f'{k} {f}', argv(f'--sortby{k} {T}/{f} --quiet --notrunclabels '
                                  '--minseqlength 0 --fasta_width 0 --output {O}/out'))
            for k in ('length', 'size') for f in ('a.fa', 'q.fq')]


@source('getseq', seeds(1, 30), 'every selection mode, selected and discarded')
def getseq(T, seed, env):
    gen('getseq', T, seed)
    out = []

    def c(cmd, f, v):
        if f.endswith('.fq'):
            o, nm, w = '--fastqout', '--notmatchedfq', ''
        else:
            o, nm, w = '--fastaout', '--notmatched', '--fasta_width 0'
        out.append(Case(f'{f} {cmd} {v}', argv(
            f'{cmd} {T}/{f} --quiet --notrunclabels {w} {o} {{O}}/k {nm} {{O}}/x ')
            + v))
    for f in ('a.fa', 'q.fq'):
        for l in ['seq1', 'SEQ1', 's1', 'x', 'ok', 'a b', 'seq1;sample=s1']:
            c('--fastx_getseq', f, ['--label', l])
            c('--fastx_getseq', f, ['--label', l, '--label_substr_match'])
            c('--fastx_getseqs', f, ['--label_word', l])
            c('--fastx_getseqs', f, ['--label_word', l, '--label_field', 'sample'])
            c('--fastx_getseqs', f, ['--label_word', l, '--label_field', 'tag'])
        c('--fastx_getseqs', f, ['--labels', f'{T}/labels.txt'])
        c('--fastx_getseqs', f, ['--labels', f'{T}/labels.txt', '--label_substr_match'])
        c('--fastx_getseqs', f, ['--label_words', f'{T}/words.txt'])
        c('--fastx_getseqs', f, ['--label_words', f'{T}/words.txt',
                                 '--label_field', 'sample'])
    for l in ('seq1', 's1'):
        c('--fastx_getseq', 'bad.fa', ['--label', l])
        c('--fastx_getseqs', 'bad.fa', ['--label_word', l])
    return out


@source('subseq', seeds(1, 20), 'many regions, FASTA and FASTQ, abundances of zero')
def subseq(T, seed, env):
    gen('subseq', T, seed)
    r = random.Random(seed)
    regions = []
    for _ in range(8):
        s = r.randint(1, 45)
        regions.append((s, r.randint(s, 60)))
    out = []
    for f in ('a.fa', 'q.fq', 'bad.fa', 'bad.fq'):
        o, w = ('--fastqout', '') if f.endswith('.fq') else ('--fastaout', '--fasta_width 0')
        base = argv(f'--fastx_getsubseq {T}/{f} --quiet --notrunclabels '
                    f'--label_substr_match {w} {o} {{O}}/out') + ['--label', '']
        out.append(Case(f'{f}', base))
        for s, e in regions:
            for v in ([f'--subseq_start', str(s), '--subseq_end', str(e)],
                      ['--subseq_start', str(s)], ['--subseq_end', str(e)]):
                out.append(Case(f'{f} {" ".join(v)}', base + v))
    return out


@source('fqchar', seeds(1, 60), "reports against vsearch's --log, several tails")
def fqchar(T, seed, env):
    gen('fqchar', T, seed)
    return [Case(f'tail {t}', argv(f'--fastq_chars {T}/q.fq --quiet --fastq_tail {t} '
                                   '--log {O}/log')) for t in (1, 2, 4, 7)]


def _report_cases(T, cmd, out_opt, pairs):
    return [Case(f'{f} {v}', argv(f'{cmd} {T}/{f} --quiet {out_opt} {{O}}/out ' + v),
                 outputs_on_error=False) for f, v in pairs]


@source('fqstat', seeds(1, 200), "reports against --log at offsets 33 and 64")
def fqstat(T, seed, env):
    gen('fqstat', T, seed)
    return [Case(f'{f} {v}', argv(f'--fastq_stats {T}/{f} --quiet --log {{O}}/log ' + v),
                 outputs_on_error=False)
            for f, v in [('q.fq', ''), ('q.fq', '--fastq_qmax 41'),
                         ('s.fq', '--fastq_ascii 64 --fastq_qmin -5'),
                         ('s.fq', '--fastq_ascii 64')]]


@source('eestat', seeds(1, 200), "tables against --output at offsets 33 and 64")
def eestat(T, seed, env):
    gen('fqstat', T, seed)
    return _report_cases(T, '--fastq_eestats', '--output', [
        ('q.fq', ''), ('q.fq', '--fastq_qmax 41'),
        ('s.fq', '--fastq_ascii 64 --fastq_qmin -5'), ('s.fq', '--fastq_ascii 64')])


@source('eegrid', seeds(1, 200), "grids against --output, several cutoffs")
def eegrid(T, seed, env):
    gen('fqstat', T, seed)
    return _report_cases(T, '--fastq_eestats2', '--output', [
        ('q.fq', ''),
        ('q.fq', '--length_cutoffs 1,*,1 --ee_cutoffs 0.1,1,0.5,3'),
        ('q.fq', '--length_cutoffs 3,100,7 --ee_cutoffs 0.25'),
        ('q.fq', '--length_cutoffs 200,*,25 --fastq_qmax 41'),
        ('s.fq', '--length_cutoffs 2,*,3 --ee_cutoffs 1,2,5 --fastq_ascii 64 '
         '--fastq_qmin -5'),
        ('s.fq', '--fastq_ascii 64')])


@source('shufl', seeds(1, 30), 'random sizes and seeds, negative and extreme')
def shufl(T, seed, env):
    gen('shufl', T, seed)
    return [Case(f'{f} seed {s}', argv(
        f'--shuffle {T}/{f} --quiet --notrunclabels --minseqlength 0 --fasta_width 0 '
        f'--randseed {s} --output {{O}}/out'))
        for s in read(f'{T}/seeds').split() for f in ('a.fa', 'q.fq')]


@source('subsmp', seeds(1, 30), 'draws by count and share, by size or not')
def subsmp(T, seed, env):
    gen('subsmp', T, seed)
    out = []
    for line in read(f'{T}/draws').split('\n'):
        if not line:
            continue
        s, size = line.split('|')
        vsize = size.replace('--count', '--sample_size').replace('--pct', '--sample_pct')
        for f in ('a.fa', 'q.fq'):
            o, d, w = fx(f)
            for by in ('', '--sizein --sizeout'):
                out.append(Case(f'{f} seed {s} {size} {by}', argv(
                    f'--fastx_subsample {T}/{f} --quiet --notrunclabels --minseqlength 0 '
                    f'{w} {o} {{O}}/k {d} {{O}}/x --randseed {s} {vsize} {by}')))
    return out


@source('scrmbl', seeds(1, 20), 'every k from 1 to 9, random seeds')
def scrmbl(T, seed, env):
    gen('scrmbl', T, seed)
    return [Case(f'{f} seed {s} k {k}', argv(
        f'--scramble {T}/{f} --quiet --notrunclabels {fx(f)[2]} {fx(f)[0]} {{O}}/out '
        f'--randseed {s} --scramble_kmer {k}'))
        for s in read(f'{T}/seeds').split() for k in range(1, 10)
        for f in ('a.fa', 'q.fq')]


@source('lcmask', seeds(1, 200), 'DUST soft and hard, hard masking alone, clearing')
def lcmask(T, seed, env):
    gen('lcmask', T, seed)
    out = []
    for f in ('a.fa', 'q.fq'):
        o, _, w = fx(f)
        base = f'--fastx_mask {T}/{f} --quiet --notrunclabels --minseqlength 0 {w} {o} {{O}}/out '
        for v in ('--qmask dust', '--qmask dust --hardmask', '--qmask soft --hardmask'):
            out.append(Case(f'{f} {v}', argv(base + v)))

        # clearing: vsearch writes the input as it is, upper-cased here on
        # the sequence lines
        def probe(rdir, f=f, key=o.lstrip('-')):
            o_ = os.path.join(rdir, 'out')
            return [Step('lcmask', ['-i', f'{T}/{f}', '--clear'], None, o_)], \
                {key: o_}
        out.append(Case(f'{f} --clear', argv(
            f'--fastx_mask {T}/{f} --quiet --notrunclabels --minseqlength 0 '
            f'--fasta_width 0 --qmask none {o} {{O}}/out'), probe=probe,
            comparators={o.lstrip('-'): 'upper-seq'}))
    return out


@source('syncpa', seeds(1, 300), 'shuffled pairs, orphans, separator sets')
def syncpa(T, seed, env):
    gen('syncpa', T, seed)
    out = []
    for sep in ('', '_/', '._/'):
        for f, r, x in (('f.fa', 'r.fa', 'fasta'), ('f.fq', 'r.fq', 'fastq')):
            v = argv(f'--fastx_syncpairs {T}/{f} --reverse {T}/{r} --quiet '
                     f'--fasta_width 0 --{x}out {{O}}/1 --{x}out_rev {{O}}/2 '
                     f'--{x}out_orphans {{O}}/3 --{x}out_orphans_rev {{O}}/4')
            if sep:
                v += ['--read_separators', sep]
            out.append(Case(f'{f} sep[{sep}]', v, errors='exact'))
    return out


@source('joinpa', seeds(1, 300), 'pairs of any lengths, five pads')
def joinpa(T, seed, env):
    gen('joinpa', T, seed)
    base = argv(f'--fastq_join {T}/f.fq --reverse {T}/r.fq --quiet --fastqout {{O}}/out')
    return [Case(f'pad {k}', base + v, errors='joinpa') for k, v in enumerate([
        [], ['--join_padgap', '', '--join_padgapq', ''], ['--join_padgap', 'ACG'],
        ['--join_padgap', 'x y', '--join_padgapq', '!~#'],
        ['--join_padgap', 'NN', '--fastq_ascii', '64']])]


@source('mergpa', seeds(1, 700), 'pairs from random fragments, four option sets')
def mergpa(T, seed, env):
    gen('mergpa', T, seed)
    base = argv(f'--fastq_mergepairs {T}/f.fq --reverse {T}/r.fq --quiet '
                '--fastqout {O}/1 --fastqout_notmerged_fwd {O}/2 '
                '--fastqout_notmerged_rev {O}/3 --eetabbedout {O}/4 --log {O}/log')
    return [Case(f'opts {k}', base + argv(v), errors='mergpa', outputs_on_error=False)
            for k, v in enumerate([
                '', '--fastq_minovlen 5 --fastq_maxdiffs 3 --fastq_allowmergestagger',
                '--fastq_minovlen 7 --fastq_maxdiffpct 5 --fastq_qmaxout 60 '
                '--fastq_qminout 3',
                '--fastq_minovlen 20 --fastq_maxdiffs 0 --fastq_qmin 2 --fastq_qmax 44'])]


@source('redig', seeds(1, 1000), 'planted sites, patterns with ^ and _ anywhere')
def redig(T, seed, env):
    gen('redig', T, seed)
    p = read(f'{T}/pattern')
    return [Case(f'pattern {p}', argv(
        f'--cut {T}/in.fa --quiet --fasta_width 0 --fastaout {{O}}/1 --fastaout_rev '
        '{O}/2 --fastaout_discarded {O}/3 --fastaout_discarded_rev {O}/4')
        + ['--cut_pattern', p], errors='abundance')]
