"""Cases for the analysis tools, ported from the zsh runners with the
same generators and seeds.

Where a runner applied an rsearch tool to vsearch's own hits (hitfmt,
otutab, lcatax, getseq --hits, clspl), the case is now vsearch's command
whole and the translated pipeline runs the search too, which checks the
search and the formatting together."""

import os
import subprocess

from ..cases import Case
from . import argv, gen, read, seeds, source

VFLAGS = {'acceptall', 'leftjust', 'rightjust', 'self', 'selfid',
          'top_hits_only', 'n_mismatch', 'sizeorder'}


def params(T):
    return read(f'{T}/params').split()


def vopts(words):
    """A generator's parameter words (vsearch option names without dashes,
    values after them) as vsearch options."""
    out, i = [], 0
    while i < len(words):
        n = words[i]
        if n in VFLAGS:
            out.append('--' + n)
        else:
            out += ['--' + n, words[i + 1]]
            i += 1
        i += 1
    return out


def sh(env, tool, *args) -> int:
    """Runs an rsearch tool (to make a UDB database both programs read)."""
    return subprocess.run([os.path.join(env.bindir, tool)] + list(args),
                          stdout=subprocess.DEVNULL,
                          stderr=subprocess.DEVNULL).returncode


def pipe(env, cmds, src, dst) -> int:
    """Runs rsearch tools one into the next, from src into dst."""
    cur = src
    for k, (tool, args) in enumerate(cmds):
        out = dst if k == len(cmds) - 1 else f'{dst}.{k}'
        with open(out, 'wb') as f:
            code = subprocess.run([os.path.join(env.bindir, tool), '-i', cur]
                                  + args, stdout=f,
                                  stderr=subprocess.DEVNULL).returncode
        if code:
            return code
        cur = out
    return 0


def nonempty(path):
    return os.path.exists(path) and os.path.getsize(path) > 0


MASK_V = {'dust': ['--dbmask', 'dust'], 'soft': ['--dbmask', 'soft'],
          'none': ['--dbmask', 'none'],
          'dusthard': ['--dbmask', 'dust', '--hardmask'],
          'softhard': ['--dbmask', 'soft', '--hardmask'],
          'nonehard': ['--dbmask', 'none', '--hardmask']}


@source('mkudb', seeds(1, 600), 'UDB files through sqfilt and lcmask, six '
        'masking choices, word lengths 3 to 13')
def mkudb(T, seed, env, long=None):
    if long:
        gen('mkudb', T, long, long)
    else:
        gen('mkudb', T, seed)
    k, minlen, maxlen, mask, fmt = params(T)
    return [Case(f'k={k} mask={mask} fmt={fmt}', argv(
        f'--makeudb_usearch {T}/in.{fmt} --output {{O}}/udb --quiet --notrunclabels '
        f'--wordlength {k} --minseqlength {minlen} --maxseqlength {maxlen}')
        + MASK_V[mask], errors='mkudb', warnings=['nowords'],
        success_only=['output', 'nowords'],
        coverage={'k': k, 'mask': mask, 'fmt': fmt})]


@source('mkudb-long', [13, 14, 15], 'the long words: files of 4^k words '
        '(256 MB at 13, 4 GB at 15)')
def mkudb_long(T, seed, env):
    return mkudb(T, seed, env, long=seed)


@source('udbinfo', seeds(1, 500), 'headers valid or with one field wrong')
def udbinfo(T, seed, env):
    gen('udbinfo', T, seed)
    return [Case('udbinfo', argv(f'--udbinfo {T}/in.udb --quiet --log {{O}}/log'),
                 errors='invalid')]


def _udb_skip(T, magic=False):
    def skip(vr, rr):
        if vr.exit == 134:
            return 'vsearch aborted'
        if magic:
            with open(f'{T}/in.udb', 'rb') as f:
                head = f.read(4)
            if head != b'FBDU' and 'Invalid UDB file' in vr.stderr:
                return 'not a UDB file to conv'
        return None
    return skip


@source('udbstat', seeds(1, 700), "mkudb's files, half of them spoiled in one "
        'part; vsearch aborts skipped')
def udbstat(T, seed, env):
    gen('udbstat', T, seed, os.path.join(env.bindir, 'mkudb'))
    return [Case('udbstats', argv(f'--udbstats {T}/in.udb --quiet --log {{O}}/log'),
                 errors='udb', success_only=['log'], skip=_udb_skip(T),
                 coverage={'what': read(f'{T}/what').strip()})]


@source('conv-udb', seeds(1, 700), "udbstat's spoiled files against "
        '--udb2fasta, unwrapped and at 80')
def conv_udb(T, seed, env):
    gen('udbstat', T, seed, os.path.join(env.bindir, 'mkudb'))
    return [Case(f'width {w}', argv(f'--udb2fasta {T}/in.udb --output {{O}}/fa '
                                    f'--quiet --fasta_width {w}'),
                 errors='udb', skip=_udb_skip(T, magic=True)) for w in (0, 80)]


@source('dbornt', seeds(1, 600), 'queries cut from a database, six database '
        'masking choices, three query ones, FASTA or UDB')
def dbornt(T, seed, env):
    for f in ('q.fa', 'q.fq'):
        if os.path.exists(f'{T}/{f}'):
            os.remove(f'{T}/{f}')
    gen('dbornt', T, seed)
    k, dbmask, qmask, dbkind, fmt = params(T)
    db, vmask = f'{T}/db.fa', MASK_V[dbmask]
    if dbkind == 'udb':
        masking = {'dust': [('lcmask', [])], 'dusthard': [('lcmask', ['--hard'])],
                   'soft': [], 'softhard': [('lcmask', ['--hard', '--no-dust'])]}.get(
                       dbmask, [('lcmask', ['--clear'])])
        rdb = f'{T}/rdb.fa'
        if pipe(env, [('sqfilt', ['--max-len', '50000'])] + masking, db, rdb) == 0 \
                and nonempty(rdb) and sh(env, 'mkudb', '-i', rdb, '--word-len', k,
                                         '-o', f'{T}/db.udb') == 0:
            db, vmask = f'{T}/db.udb', []
    out = '--fastqout' if fmt == 'fq' else '--fastaout'
    records = '' if qmask == 'none' else f'{out} {{O}}/out --notmatched {{O}}/not'
    v = argv(f'--orient {T}/q.{fmt} --db {db} {records} --tabbedout {{O}}/tab --quiet '
             f'--notrunclabels --fasta_width 0 --wordlength {k} --qmask {qmask} '
             '--log {O}/log') + vmask
    return [Case(f'k={k} dbmask={dbmask} qmask={qmask} db={dbkind} fmt={fmt}', v,
                 errors='abundance', warnings=['nowords', 'wordlength-adjusted'],
                 success_only=['log', 'nowords', 'wordlength-adjusted'],
                 coverage={'dbmask': dbmask, 'qmask': qmask, 'db': dbkind},
                 counts={'fwd': ('tabbedout', 'col 2 +'),
                         'rev': ('tabbedout', 'col 2 -'),
                         'und': ('tabbedout', 'col 2 ?')})]


@source('sintax', seeds(1, 600), 'taxonomies on random trees, strands, cutoffs, '
        'random ties, seeds, masking, FASTA or UDB')
def sintax(T, seed, env):
    for f in ('q.fa', 'q.fq'):
        if os.path.exists(f'{T}/{f}'):
            os.remove(f'{T}/{f}')
    gen('sintax', T, seed)
    k, strand, cutoff, rand, rseed, dbmask, qmode, dbkind, fmt = params(T)
    db, q = f'{T}/db.fa', f'{T}/q.{fmt}'
    if dbkind == 'udb':
        rdb = f'{T}/rdb.fa'
        cmds = [('sqfilt', ['--min-len', '32', '--max-len', '50000'])]
        if dbmask == 'none':
            cmds.append(('lcmask', ['--clear']))
        if pipe(env, cmds, db, rdb) == 0 and nonempty(rdb) and sh(
                env, 'mkudb', '-i', rdb, '--word-len', k, '-o', f'{T}/db.udb') == 0:
            db = f'{T}/db.udb'
    if qmode == 'mask':
        # vsearch never masks queries: their lower case turned into N
        lines = read(q).split('\n')
        for i, l in enumerate(lines):
            if (i % 4 == 1) if fmt == 'fq' else not l.startswith('>'):
                lines[i] = ''.join('N' if c.islower() else c for c in l)
        q = f'{T}/v.q'
        with open(q, 'w') as f:
            f.write('\n'.join(lines))
    v = argv(f'--sintax {q} --db {db} --tabbedout {{O}}/tab --quiet --notrunclabels '
             f'--threads 1 --strand {strand} --randseed {rseed} --wordlength {k} '
             f'--dbmask {dbmask} --log {{O}}/log')
    if cutoff != '0':
        v += ['--sintax_cutoff', cutoff]
    if rand == '1':
        v.append('--sintax_random')
    return [Case(f'k={k} strand={strand} cutoff={cutoff} random={rand} '
                 f'dbmask={dbmask} q={qmode} db={dbkind} fmt={fmt}', v,
                 errors='abundance', success_only=['log'],
                 coverage={'strand': strand, 'db': dbkind, 'q': qmode},
                 counts={'lines': ('tabbedout', 'lines')})]


@source('exsrch', seeds(1, 600), 'few sequences copied in either case, strands, '
        '--maxhits, --self, abundance filters')
def exsrch(T, seed, env):
    for f in ('q.fa', 'q.fq'):
        if os.path.exists(f'{T}/{f}'):
            os.remove(f'{T}/{f}')
    gen('exsrch', T, seed)
    w = params(T)
    fmt, o, vo, i = w[0], w[1:], [], 0
    while i < len(o):
        n = o[i]
        if n == 'strand':
            vo += ['--strand', 'both']
        elif n == 'self':
            vo.append('--self')
        else:
            vo += ['--' + n, o[i + 1]]
            i += 1
        i += 1
    return [Case(f'fmt={fmt} {" ".join(o)}', argv(
        f'--search_exact {T}/q.{fmt} --db {T}/db.fa --uc {{O}}/uc --uc_allhits '
        '--quiet --notrunclabels --threads 1') + vo, errors='abundance',
        success_only=['uc'], counts={'uc': ('uc', 'uc')})]


HITFMT_FIELDS = ('query+target+evalue+id+pctpv+pctgaps+pairs+gaps+qlo+qhi+tlo+thi+'
                 'pv+ql+tl+qs+ts+alnlen+opens+exts+raw+bits+aln+caln+qstrand+'
                 'tstrand+qrow+trow+qframe+tframe+mism+ids+qcov+tcov+id0+id1+id2+'
                 'id3+id4+qilo+qihi+tilo+tihi+diffs+mid+qseq+tseq+qrowdots+'
                 'trowdots+qlor+qhir+tlor+thir')


@source('hitfmt', seeds(1, 600), "--usearch_global's hit outputs: every "
        'identity definition, scores, row lengths, DUST masking')
def hitfmt(T, seed, env):
    gen('hitfmt', T, seed)
    w = params(T)
    mask, rowlen, o = w[0], w[1], w[2:]
    v = argv(f'--usearch_global {T}/q.fa --db {T}/db.fa --quiet --notrunclabels '
             f'--threads 1 --qmask {mask} --dbmask {mask}') + o + argv(
        f'--rowlen {rowlen} --uc {{O}}/uc --uc_allhits --alnout {{O}}/aln '
        f'--blast6out {{O}}/b6 --userout {{O}}/user --userfields {HITFMT_FIELDS} '
        '--samout {O}/sam --samheader --fastapairs {O}/pairs --qsegout {O}/qseg '
        '--tsegout {O}/tseg --fasta_width 0')
    return [Case(f'mask={mask} rowlen={rowlen} {" ".join(o)}', v,
                 outputs_on_error=False, coverage={'mask': mask},
                 counts={'uc': ('uc', 'uc')})]


@source('otutab', seeds(1, 600), 'the tables of exact and global searches')
def otutab(T, seed, env):
    gen('otutab', T, seed)
    w = params(T)
    cmd = '--usearch_global' if w[0] == 'global' else '--search_exact'
    return [Case(f'{w[0]} {" ".join(w[1:])}', argv(
        f'{cmd} {T}/q.fa --db {T}/db.fa --quiet --notrunclabels --threads 1 '
        '--minseqlength 1') + w[1:] + argv(
        '--uc {O}/uc --uc_allhits --otutabout {O}/tsv --mothur_shared_out {O}/mo '
        '--biomout {O}/biom'), outputs_on_error=False,
        coverage={'mode': w[0]}, counts={'uc': ('uc', 'uc')})]


@source('lcatax', seeds(1, 600), 'taxonomies on random tree paths, many hits, '
        'cutoffs')
def lcatax(T, seed, env):
    gen('lcatax', T, seed)
    w = params(T)
    return [Case(f'cutoff={w[0]} {" ".join(w[1:])}', argv(
        f'--usearch_global {T}/q.fa --db {T}/db.fa --quiet --notrunclabels '
        '--threads 1 --minseqlength 1') + w[1:] + argv(
        f'--lca_cutoff {w[0]} --uc {{O}}/uc --uc_allhits --lcaout {{O}}/lca'),
        outputs_on_error=False, counts={'uc': ('uc', 'uc')})]


@source('getseq-hits', seeds(1, 600), "--usearch_global's matched and not "
        "matched queries and targets, from hitfmt's searches")
def getseq_hits(T, seed, env):
    gen('hitfmt', T, seed)
    w = params(T)
    mask, o = w[0], w[2:]
    return [Case(f'mask={mask} {" ".join(o)}', argv(
        f'--usearch_global {T}/q.fa --db {T}/db.fa --quiet --notrunclabels '
        f'--threads 1 --qmask {mask} --dbmask {mask}') + o + argv(
        '--uc {O}/uc --uc_allhits --matched {O}/m --notmatched {O}/n '
        '--dbmatched {O}/dm --dbnotmatched {O}/dn --fasta_width 0'),
        outputs_on_error=False,
        counts={'matched': ('matched', 'records'),
                'dbmatched': ('dbmatched', 'records')})]


@source('allpr', seeds(1, 1300), 'every criterion, identity definition, score '
        'and gap penalty')
def allpr(T, seed, env):
    gen('allpr', T, seed)
    vo = vopts(params(T))
    return [Case(' '.join(vo), argv(
        f'--allpairs_global {T}/in.fa --uc {{O}}/uc --quiet --notrunclabels '
        '--qmask none --minseqlength 0 --maxseqlength 100000 --threads 1') + vo,
        success_only=['uc'], counts={'uc': ('uc', 'uc')})]


@source('dbsrch', seeds(1, 1000), 'related sequences, UDB databases, accepts and '
        'rejects, words, strands, criteria, scores')
def dbsrch(T, seed, env):
    gen('dbsrch', T, seed)
    w = params(T)
    db = f'{T}/db.fa'
    if w[0] == 'udb':
        sh(env, 'mkudb', '-i', db, '-o', f'{T}/db.udb', '--word-len', w[1])
        db, w = f'{T}/db.udb', w[2:]
    vo = vopts(w)
    return [Case(' '.join(vo), argv(
        f'--usearch_global {T}/q.fa --db {db} --uc {{O}}/uc --uc_allhits --quiet '
        '--notrunclabels --qmask soft --dbmask soft --minseqlength 0 '
        '--maxseqlength 100000 --threads 1') + vo, success_only=['uc'],
        coverage={'db': db[-3:]}, counts={'uc': ('uc', 'uc')})]


CLUSTER = {'smallmem': '--cluster_smallmem {in_} --usersort',
           'fast': '--cluster_fast {in_}', 'size': '--cluster_size {in_}'}
COMMON = ('--quiet --notrunclabels --qmask soft --minseqlength 0 '
          '--maxseqlength 100000 --threads 1 --fasta_width 0')


def _cluster_params(T):
    w = params(T)
    mode, rest = w[0], w[1:]
    sizes = {'nosizein': False, 'nosizeout': False}
    words = []
    for x in rest:
        if x in sizes:
            sizes[x] = True
        else:
            words.append(x)
    return mode, words, sizes


@source('clust', seeds(1, 1000), 'input order or sorted by length or size; '
        'limits, words, strands, abundances, criteria, scores')
def clust(T, seed, env):
    gen('clust', T, seed)
    mode, words, sizes = _cluster_params(T)
    s = ([] if sizes['nosizein'] else ['--sizein']) + \
        ([] if sizes['nosizeout'] else ['--sizeout'])
    vo = vopts(words)
    return [Case(f'{mode} {" ".join(vo + s)}', argv(
        CLUSTER[mode].format(in_=f'{T}/in.fa') + ' --uc {O}/uc --centroids {O}/c '
        + COMMON) + s + vo, success_only=['uc', 'centroids'],
        coverage={'mode': mode}, counts={'uc': ('uc', 'uc')})]


@source('unoise', seeds(1, 1000), 'abundant templates and their variants; alpha, '
        'identity, limits, words, strands, abundances')
def unoise(T, seed, env):
    gen('unoise', T, seed)
    w = params(T)
    s = ([] if 'nosizein' in w else ['--sizein']) + \
        ([] if 'nosizeout' in w else ['--sizeout'])
    vo = vopts([x for x in w if x not in ('nosizein', 'nosizeout')])
    return [Case(' '.join(vo + s), argv(
        f'--cluster_unoise {T}/in.fa --uc {{O}}/uc --centroids {{O}}/c --quiet '
        '--notrunclabels --qmask soft --minseqlength 1 --maxseqlength 100000 '
        '--threads 1 --fasta_width 0') + s + vo, success_only=['uc', 'centroids'],
        counts={'uc': ('uc', 'uc')})]


def _id_strand_words(words):
    vo, i = [], 0
    while i < len(words):
        if words[i] in ('id', 'strand', 'wordlength'):
            vo += ['--' + words[i], words[i + 1]]
            i += 2
        elif words[i] in VFLAGS:
            i += 1
        else:
            i += 2
    return vo


@source('clspl', seeds(1, 600), "--clusters, from clust's random inputs")
def clspl(T, seed, env):
    gen('clust', T, seed)
    mode, words, _ = _cluster_params(T)
    vo = _id_strand_words(words)
    return [Case(f'{mode} {" ".join(vo)}', argv(
        CLUSTER[mode].format(in_=f'{T}/in.fa') + ' --uc {O}/uc --clusters '
        '{O}/cluster --sizein ' + COMMON) + vo, success_only=['uc', 'clusters'],
        coverage={'mode': mode})]


@source('clmsa', seeds(1, 600), "--msaout, --consout and --profile, from clust's "
        'random inputs')
def clmsa(T, seed, env):
    gen('clust', T, seed)
    mode, words, sizes = _cluster_params(T)
    vo = _id_strand_words(words)
    s = ([] if sizes['nosizein'] else ['--sizein']) + \
        ([] if sizes['nosizeout'] else ['--sizeout'])
    return [Case(f'{mode} {" ".join(vo + s)}', argv(
        CLUSTER[mode].format(in_=f'{T}/in.fa') + ' --msaout {O}/msa --consout '
        '{O}/cons --profile {O}/prof ' + COMMON) + s + vo,
        success_only=['msaout', 'consout', 'profile'], coverage={'mode': mode},
        counts={'clusters': ('consout', 'records')})]


@source('uchime', seeds(1, 1000), 'the four UCHIME commands, scores, thresholds, '
        'skew, alignment width; all five outputs')
def uchime(T, seed, env):
    gen('uchime', T, seed)
    w = params(T)
    mode, vo = w[0], ['--' + x if k % 2 == 0 else x for k, x in enumerate(w[1:])]
    db = f'--db {T}/db.fa --dbmask soft' if mode == 'uchime_ref' else ''
    return [Case(f'{mode} {" ".join(vo)}', argv(
        f'--{mode} {T}/in.fa {db} --nonchimeras {{O}}/non --chimeras {{O}}/chim '
        '--borderline {O}/bord --uchimeout {O}/tab --uchimealns {O}/aln '
        + COMMON) + vo,
        success_only=['nonchimeras', 'chimeras', 'borderline', 'uchimeout',
                      'uchimealns'],
        coverage={'mode': mode},
        counts={'Y': ('uchimeout', 'col 17 Y'), 'N': ('uchimeout', 'col 17 N'),
                'border': ('uchimeout', 'col 17 ?')})]


@source('mchime', seeds(1, 1000), 'chimeras of two to five long parents; parts, '
        'parents, segments, differences, skew, width')
def mchime(T, seed, env):
    gen('mchime', T, seed)
    w = params(T)
    vo = ['--' + x if k % 2 == 0 else x for k, x in enumerate(w)]
    return [Case(' '.join(vo), argv(
        f'--chimeras_denovo {T}/in.fa --nonchimeras {{O}}/non --chimeras {{O}}/chim '
        '--tabbedout {O}/tab --alnout {O}/aln ' + COMMON) + vo,
        success_only=['nonchimeras', 'chimeras', 'tabbedout', 'alnout'],
        counts={'Y': ('tabbedout', 'lines')})]
