"""The translator on command lines: the pipelines of the vsearch commands,
the option names, and the errors and caveats.

Run with `python3 -m rscompare test` (or `python3 -m unittest discover
rscompare/tests` from compare/). Files named here do not exist, so no
input is known to be FASTA and FASTA outputs always go through conv."""

import unittest

from rscompare.translate import translate, TranslationError
from rscompare.vsearch_options import Args

W = 'w'


def render(cmd):
    return translate(cmd.split(), W).render()


def outputs(cmd):
    return {o.option: (o.rsearch, o.comparator)
            for o in translate(cmd.split(), W).outputs}


def caveats(cmd):
    return translate(cmd.split(), W).caveats


class ToolsMdPreprocessing(unittest.TestCase):
    """The preprocessing pipelines, one by one."""

    def test_filter_trunclen(self):
        self.assertEqual(render(
            '--fastx_filter in --fastq_trunclen 200 --fastq_maxee 1 --relabel R '
            '--fastaout out'),
            'sqtrim -i in --len 200 | sqfilt --min-len 200 --max-ee 1 | '
            'relab --relabel R | conv --fasta > w/tmp4.fastaout')

    def test_filter_pairs(self):
        self.assertEqual(render(
            '--fastq_filter in --reverse rev --fastq_stripleft 3 --fastq_maxee 1 '
            '--fastqout k --fastqout_rev k2'),
            'sqtrim -i in --strip-left 3 > w/tmp1.sqtrim\n'
            'sqtrim -i rev --strip-left 3 > w/tmp2.sqtrim\n'
            'sqfilt --reverse w/tmp2.sqtrim --rev-out w/tmp4.rev_kept --max-ee 1 '
            '< w/tmp1.sqtrim > w/tmp3.kept')

    def test_uniques(self):
        self.assertEqual(render(
            '--fastx_uniques in --minseqlength 32 --minuniquesize 2 --topn 10 '
            '--relabel U --sizeout --fastaout o'),
            'sqfilt -i in --min-len 32 --max-len 50000 | derep --no-sizein | '
            'sqfilt --min-size 2 | topn 10 | relab --relabel U | conv --fasta '
            '> w/tmp6.fastaout')

    def test_derep_smallmem(self):
        self.assertEqual(render('--derep_smallmem in --fastaout out'),
                         'sqfilt -i in --max-len 50000 | derep --unsorted '
                         '--no-sizein --no-sizeout | conv --fasta > w/tmp3.fastaout')

    def test_rereplicate(self):
        self.assertEqual(render('--rereplicate in --sizeout --output out'),
                         'rerep -i in | conv --fasta > w/tmp2.output')
        self.assertEqual(render('--rereplicate in --output out'),
                         'rerep -i in --no-sizeout | conv --fasta > w/tmp2.output')

    def test_join(self):
        self.assertEqual(render('--fastq_join in --reverse rev --relabel J '
                                '--fastaout out'),
                         'joinpa -i in --reverse rev | relab --relabel J | '
                         'conv --fasta > w/tmp3.fastaout')

    def test_mergepairs_filters(self):
        self.assertEqual(render(
            '--fastq_mergepairs in --reverse rev --fastq_truncqual 2 '
            '--fastq_minlen 50 --fastq_maxee 1 --fastqout m'),
            'sqtrim -i in --qual 2 > w/tmp1.sqtrim\n'
            'sqtrim -i rev --qual 2 > w/tmp2.sqtrim\n'
            'sqfilt -i w/tmp1.sqtrim --reverse w/tmp2.sqtrim --rev-out w/tmp3.rev '
            '--min-len 50 > w/tmp4.sqfilt\n'
            'mergpa -i w/tmp4.sqfilt --reverse w/tmp3.rev | sqfilt --max-ee 1 '
            '> w/tmp6.sqfilt')

    def test_mergepairs_lengths(self):
        self.assertEqual(render(
            '--fastq_mergepairs in --reverse rev --fastq_minmergelen 100 '
            '--fastq_maxmergelen 300 --fastqout m'),
            'mergpa -i in --reverse rev | sqfilt --min-len 100 --max-len 300 '
            '> w/tmp2.sqfilt')

    def test_getsubseq(self):
        self.assertEqual(render('--fastx_getsubseq in --label s1 --subseq_start 3 '
                                '--subseq_end 6 --fastaout o'),
                         'getseq -i in --label s1 | subseq --start 3 --end 6 | '
                         'conv --fasta > w/tmp3.fastaout')

    def test_sortbysize(self):
        self.assertEqual(render('--sortbysize in --minsize 2 --topn 5 --output out'),
                         'sqfilt -i in --max-len 50000 --min-size 2 | sortby --size '
                         '| topn 5 | conv --fasta > w/tmp4.output')

    def test_fastq_convert(self):
        self.assertEqual(render('--fastq_convert in --fastq_ascii 64 --relabel R '
                                '--fastqout o'),
                         'conv -i in --ascii 64 | relab --relabel R > w/tmp2.relab')

    def test_cut(self):
        self.assertEqual(render('--cut in --cut_pattern G^AATT_C --fastaout f '
                                '--fastaout_discarded d'),
                         "redig -i in --pattern 'G^AATT_C' --rev-out w/tmp2.rev "
                         "--uncut w/tmp3.uncut --rev-uncut w/tmp4.rev_uncut | "
                         "conv --fasta > w/tmp5.fastaout\n"
                         "conv --fasta < w/tmp3.uncut > w/tmp6.fastaout_discarded")

    def test_mask(self):
        self.assertEqual(render('--fastx_mask in --fastaout out'),
                         'lcmask -i in | sqfilt --max-len 50000 | conv --fasta '
                         '> w/tmp3.fastaout')
        self.assertEqual(render('--maskfasta in --output out'),
                         'lcmask -i in | sqfilt --max-len 50000 | conv --fasta '
                         '> w/tmp3.output')
        self.assertIn('lcmask -i in --hard --no-dust |',
                      render('--fastx_mask in --qmask soft --hardmask --fastaout o'))
        self.assertIn('lcmask -i in | sqfilt --max-len 50000 --min-upper-pct 50',
                      render('--fastx_mask in --min_unmasked_pct 50 --fastaout o'))
        self.assertIn('lcmask -i in --hard | sqfilt --max-len 50000 '
                      '--min-non-n-pct 50',
                      render('--fastx_mask in --hardmask --min_unmasked_pct 50 '
                             '--fastaout o'))


class ToolsMdAnalysis(unittest.TestCase):
    """The analysis pipelines, one by one."""

    def test_cluster_fast(self):
        self.assertEqual(render('--cluster_fast in --id 0.97 --sizein --sizeout '
                                '--centroids c --uc u'),
                         'sqfilt -i in --min-len 32 --max-len 50000 | lcmask | '
                         'sortby --length | clust --maxrejects 8 --id 0.97 --hits '
                         'w/tmp4.uc | conv --fasta > w/tmp6.centroids')
        self.assertEqual(outputs('--cluster_fast in --id 0.97 --uc u')['uc'][1],
                         'uc-fast')

    def test_cluster_size(self):
        self.assertIn('sortby --size | clust --id 0.97 --hits',
                      render('--cluster_size in --id 0.97 --sizein --sizeout '
                             '--centroids c'))

    def test_cluster_smallmem(self):
        self.assertIn('lcmask | clust --id 0.97 --no-sizein --no-sizeout',
                      render('--cluster_smallmem in --usersort --id 0.97 '
                             '--centroids c'))

    def test_cluster_unoise(self):
        self.assertIn('sqfilt -i in --min-len 32 --max-len 50000 --min-size 8 | '
                      'lcmask | sortby --size | unoise --hits',
                      render('--cluster_unoise in --sizein --sizeout --centroids c'))
        r = render('--cluster_unoise in --minsize 2 --weak_id 0.8 --unoise_alpha 3 '
                   '--id 0.5 --centroids c')
        self.assertIn('unoise --id 0.8 --alpha 3 --no-sizein', r)
        self.assertNotIn('0.5', r)

    def test_cluster_outputs(self):
        self.assertEqual(render(
            '--cluster_fast in --id 0.97 --sizein --sizeout --consout cons '
            '--clusters cl --fasta_width 0'),
            'sqfilt -i in --min-len 32 --max-len 50000 | lcmask | sortby --length '
            '> w/tmp3.sortby\n'
            'clust --maxrejects 8 --id 0.97 --hits w/tmp4.uc < w/tmp3.sortby > '
            'w/tmp5.centroids\n'
            'clspl --hits w/tmp4.uc -i w/tmp3.sortby --prefix w/clusters > '
            'w/tmp6.clspl\n'
            'clmsa --hits w/tmp4.uc -i w/tmp3.sortby --consensus > w/tmp7.consout')

    def test_usearch_global_alnout(self):
        self.assertEqual(render('--usearch_global q --db db --id 0.97 --alnout a'),
                         'lcmask -i q > w/tmp1.lcmask\n'
                         'lcmask -i db | sqfilt --min-len 32 --max-len 50000 > '
                         'w/tmp3.db\n'
                         'dbsrch --db w/tmp3.db --id 0.97 < w/tmp1.lcmask > '
                         'w/tmp4.uc\n'
                         'hitfmt w/tmp4.uc -i w/tmp1.lcmask --db w/tmp3.db '
                         '--alignments > w/tmp5.alnout')

    def test_usearch_global_tables(self):
        r = render('--usearch_global q --db db --id 0.97 --otutabout t --lcaout l')
        self.assertIn('otutab w/tmp4.uc --db w/tmp3.db > w/tmp5.otutabout', r)
        self.assertIn('lcatax w/tmp4.uc --db w/tmp3.db > w/tmp6.lcaout', r)
        self.assertEqual(outputs('--usearch_global q --db db --id 1 --biomout b')
                         ['biomout'][1], 'biom')

    def test_usearch_global_matched(self):
        r = render('--usearch_global q --db db --id 0.97 --matched m '
                   '--dbnotmatched n')
        self.assertIn('getseq --hits w/tmp4.uc -i w/tmp1.lcmask | conv --fasta', r)
        self.assertIn('getseq --hits w/tmp4.uc -i w/tmp3.db --targets --invert', r)

    def test_usearch_global_qmask_none(self):
        r = render('--usearch_global q --db db --id 0.97 --qmask none --alnout a '
                   '--rowlen 30 --gapopen 20I/2E --n_mismatch')
        self.assertIn('lcmask -i q --clear', r)
        # hitfmt shows the queries as they are, with the search's scores
        self.assertIn('hitfmt w/tmp4.uc -i q --db w/tmp3.db --alignments '
                      '--rowlen 30 --gapopen 20I/2E --n-mismatch', r)

    def test_uc_top_hit_and_dbmatched(self):
        r = render('--usearch_global q --db db --id 0.97 --uc u')
        self.assertIn('--maxhits 1 < w/tmp1.lcmask > w/tmp4.uc', r)
        r = render('--usearch_global q --db db --id 0.97 --dbmatched m --maxhits 2 '
                   '--uc u --uc_allhits')
        # --dbmatched counts every hit: a search without --maxhits (decision 28)
        self.assertIn('--id 0.97 --maxhits 2 < w/tmp1.lcmask > w/tmp4.uc', r)
        self.assertIn('--id 0.97 < w/tmp1.lcmask > w/tmp5.uc', r)
        self.assertIn('getseq --hits w/tmp5.uc -i w/tmp3.db --targets', r)

    def test_orient(self):
        self.assertEqual(render('--orient q --db db --fastaout o'),
                         'sqfilt -i db --max-len 50000 | lcmask > w/tmp2.lcmask\n'
                         'dbornt -i q --db w/tmp2.lcmask --not-oriented '
                         'w/tmp3.notmatched | conv --fasta > w/tmp5.fastaout')
        self.assertNotIn('lcmask', render('--orient q --db db --dbmask soft '
                                          '--fastaout o'))

    def test_uchime(self):
        self.assertEqual(render('--uchime_denovo in --nonchimeras n'),
                         'sqfilt -i in --max-len 50000 | lcmask | sortby --size | '
                         'uchime | conv --fasta > w/tmp5.nonchimeras')
        self.assertIn('uchime --uchime2 --abskew 16 --table w/uchimeout',
                      render('--uchime3_denovo in --nonchimeras n --uchimeout t'))
        self.assertEqual(render('--uchime_ref in --db ref --nonchimeras n'),
                         'lcmask -i in > w/tmp1.lcmask\n'
                         'lcmask -i ref | sqfilt --max-len 50000 > w/tmp3.db\n'
                         'uchime --db w/tmp3.db < w/tmp1.lcmask | conv --fasta > '
                         'w/tmp5.nonchimeras')

    def test_makeudb(self):
        self.assertEqual(render('--makeudb_usearch in --output db.udb'),
                         'sqfilt -i in --min-len 32 --max-len 50000 | lcmask | '
                         'mkudb > w/tmp3.udb')
        self.assertIn('| lcmask --clear | mkudb',
                      render('--makeudb_usearch in --dbmask none --output o'))
        self.assertIn('| lcmask --clear | mkudb',
                      render('--makeudb_usearch in --dbmask soft --output o'))
        self.assertIn('| lcmask --clear | lcmask --hard | mkudb',
                      render('--makeudb_usearch in --hardmask --output o'))

    def test_sintax(self):
        self.assertEqual(render('--sintax q --db db --randseed 1 --tabbedout t'),
                         'sqfilt -i db --min-len 32 --max-len 50000 > w/tmp1.db\n'
                         'lcmask -i q --clear | sintax --db w/tmp1.db --seed=1 > '
                         'w/tmp3.sintax')

    def test_search_exact(self):
        self.assertEqual(render('--search_exact q --db db --uc u --uc_allhits'),
                         'sqfilt -i db --max-len 50000 > w/tmp1.db\n'
                         'exsrch -i q --db w/tmp1.db > w/tmp2.uc')
        self.assertIn('exsrch -i q --db w/tmp1.db --maxhits 1',
                      render('--search_exact q --db db --uc u'))

    def test_allpairs(self):
        self.assertEqual(render('--allpairs_global in --id 0.9 --uc u'),
                         'sqfilt -i in --max-len 50000 | lcmask | allpr --id 0.9 > '
                         'w/tmp3.uc')

    def test_udb2fasta(self):
        self.assertEqual(render('--udb2fasta db.udb --output out'),
                         'conv -i db.udb --fasta > w/tmp1.output')


class Options(unittest.TestCase):

    def test_names(self):
        r = render('--allpairs_global in --qmask soft --minseqlength 1 '
                   '--maxseqlength 100 --id 0.5 --weak_id 0.4 --maxdiffs 3 '
                   '--maxsubs 2 --maxgaps 1 --mincols 4 --query_cov 0.5 '
                   '--target_cov 0.6 --maxid 0.99 --mid 80 --idprefix 2 '
                   '--idsuffix 3 --maxqt 2 --minqt 0.5 --minsl 0.5 --maxsl 2 '
                   '--maxqsize 9 --mintsize 2 --minsizeratio 0.5 --maxsizeratio 3 '
                   '--leftjust --rightjust --self --selfid --n_mismatch --iddef 1 '
                   '--match 3 --mismatch -2 --gapopen 20 --gapext 2 --maxhits 4 '
                   '--top_hits_only --maxaccepts 3 --wordlength 9 --uc u')
        self.assertEqual(r, 'sqfilt -i in --max-len 100 | allpr --id 0.5 '
                         '--iddef 1 --weak-id 0.4 --max-diffs 3 --max-subs 2 '
                         '--max-gaps 1 --min-cols 4 --query-cov 0.5 --target-cov '
                         '0.6 --max-id 0.99 --mid 80 --id-prefix 2 --id-suffix 3 '
                         '--max-qt 2 --min-qt 0.5 --min-sl 0.5 --max-sl 2 '
                         '--max-qsize 9 --min-tsize 2 --min-size-ratio 0.5 '
                         '--max-size-ratio 3 --match 3 --mismatch -2 --gapopen 20 '
                         '--gapext 2 --left-just --right-just --self --self-id '
                         '--n-mismatch --maxhits 4 --top-hits-only > w/tmp2.uc')

    def test_words(self):
        self.assertIn('dbsrch -i q --db w/tmp1.db --id 1 --strand both '
                      '--maxaccepts 0 --maxrejects 0 --word-len 10 '
                      '--min-word-matches 3',
                      render('--usearch_global q --db db --qmask soft --dbmask soft '
                             '--id 1 --strand both --maxaccepts 0 --maxrejects 0 '
                             '--wordlength 10 --minwordmatches 3 --uc u'))

    def test_quality(self):
        self.assertEqual(render('--fastq_convert in --fastq_ascii 64 --fastq_qmin 2 '
                                '--fastq_qmax 40 --fastq_asciiout 33 '
                                '--fastq_qminout 3 --fastq_qmaxout 39 '
                                '--fastq_solexa --fastqout o'),
                         'conv -i in --ascii 64 --qmin 2 --qmax 40 --ascii-out 33 '
                         '--qmin-out 3 --qmax-out 39 --solexa > w/tmp1.conv')

    def test_labels(self):
        self.assertIn("relab --relabel R --sample S --label-suffix ';x' "
                      "--relabel-keep --sizeout --xsize --lengthout --xlength "
                      "--xee --ee --ascii 64",
                      render("--fastq_convert in --fastq_ascii 64 --relabel R "
                             "--sample S --label_suffix ;x --relabel_keep "
                             "--sizeout --xsize --lengthout --xlength --xee "
                             "--eeout --fastqout o"))

    def test_fasta_width(self):
        self.assertIn('conv --fasta --fasta-width 0',
                      render('--fastx_revcomp in --fasta_width 0 --fastaout o'))
        self.assertIn('conv --fasta --fasta-width 7',
                      render('--fastx_revcomp in --fasta_width 7 --fastaout o'))

    def test_seed(self):
        self.assertIn('--seed=-3', render('--shuffle in --randseed -3 --output o'))
        self.assertIn('no --randseed', ' '.join(caveats('--shuffle in --output o')))

    def test_subsample(self):
        self.assertIn('subsmp --count 5 --by-size --seed=1 --discarded',
                      render('--fastx_subsample in --sample_size 5 --sizein '
                             '--sizeout --randseed 1 --fastqout k '
                             '--fastqout_discarded d'))

    def test_eestats2(self):
        self.assertEqual(render('--fastq_eestats2 in --length_cutoffs 3,*,7 '
                                '--ee_cutoffs 0.5,1 --output o'),
                         'eegrid -i in --min-len 3 --len-step 7 --max-ee 0.5,1 > '
                         'w/tmp1.eegrid')

    def test_logs(self):
        self.assertEqual(outputs('--fastq_chars in --fastq_tail 3 --log l'),
                         {'log': ('w/log', 'log-report')})
        self.assertEqual(outputs('--fastq_stats in --log l')['log'][1],
                         'log-report')
        self.assertEqual(outputs('--udbinfo d --log l')['log'][1], 'udbinfo-log')
        self.assertEqual(outputs('--udbstats d --log l')['log'][1], 'udbstats-log')
        self.assertEqual(outputs('--orient q --db d --tabbedout t --log l')
                         ['log'][1], 'orient-log')
        self.assertEqual(outputs('--sintax q --db d --tabbedout t --log l')
                         ['log'][1], 'sintax-log')
        self.assertEqual(outputs('--fastq_mergepairs f --reverse r --fastqout m '
                                 '--log l')['log'], ('w/log', 'log-report'))
        # elsewhere --log is a log of the run, not compared
        self.assertEqual(outputs('--fastx_revcomp in --log l --fastqout o'),
                         {'fastqout': ('w/tmp1.revcomp', 'bytes')})

    def test_mergepairs_outputs(self):
        self.assertEqual(render('--fastq_mergepairs in --reverse rev '
                                '--fastq_minovlen 5 --fastqout m '
                                '--fastqout_notmerged_fwd n --eetabbedout e '
                                '--log l'),
                         'mergpa -i in --reverse rev --min-overlap 5 --not-merged '
                         'w/tmp1.notmerged_fwd --ee-table w/eetabbedout --report '
                         'w/log > w/tmp2.merged')

    def test_args(self):
        a = Args(['vsearch', '--fastx_filter=in', '--fastq_maxee=1', '-quiet'])
        self.assertEqual((a.command, a.input, a.opts),
                         ('fastx_filter', 'in', {'fastq_maxee': '1',
                                                 'quiet': True}))
        a = Args(['--shuffle', 'in', '--randseed', '-5'])
        self.assertEqual(a.opts['randseed'], '-5')


class EveryCommand(unittest.TestCase):

    MINIMAL = {
        'usearch_global': '--db d --id 0.9 --uc u', 'search_exact': '--db d --uc u',
        'orient': '--db d --tabbedout t', 'sintax': '--db d --tabbedout t',
        'uchime_ref': '--db d --nonchimeras n', 'allpairs_global': '--id 0.9 --uc u',
        'cluster_fast': '--id 0.9 --centroids c', 'cluster_size': '--id 0.9 --uc u',
        'cluster_smallmem': '--usersort --id 0.9 --uc u',
        'cluster_unoise': '--centroids c', 'cut': '--cut_pattern G^AATT_C --fastaout o',
        'fastq_join': '--reverse r --fastqout o',
        'fastq_mergepairs': '--reverse r --fastqout o',
        'fastx_syncpairs': '--reverse r --fastqout o',
        'fastx_getsubseq': '--label a --fastaout o', 'fastx_getseq': '--label a --fastaout o',
        'fastx_getseqs': '--label_word a --fastaout o',
        'fastx_subsample': '--sample_pct 50 --randseed 1 --fastaout o',
    }

    def test_every_command(self):
        from rscompare.vsearch_options import COMMANDS
        for cmd in sorted(COMMANDS - {'help', 'version', 'h', 'v'}):
            extra = self.MINIMAL.get(cmd, '--output o')
            argv = f'--{cmd} in {extra}'.split()
            try:
                p = translate(argv, W)
            except TranslationError:
                # commands whose outputs have other names
                for out in ('--fastaout o', '--fastqout o', '--nonchimeras n',
                            '--tabbedout t', '--log l'):
                    try:
                        p = translate(f'--{cmd} in {out}'.split(), W)
                        if p.outputs:
                            break
                    except TranslationError:
                        continue
                else:
                    self.fail(f'--{cmd} does not translate')
            self.assertTrue(p.steps, cmd)
            self.assertTrue(p.outputs, cmd)


class Errors(unittest.TestCase):

    def bad(self, cmd, text):
        with self.assertRaises(TranslationError) as e:
            translate(cmd.split(), W)
        self.assertIn(text, str(e.exception))

    def test_unknown(self):
        self.bad('--fastx_filter in --no_such_option 3', '--no_such_option')
        self.bad('--fastx_filter in extra', "unexpected argument 'extra'")
        self.bad('--fastq_maxee 1', 'no vsearch command')
        self.bad('--help', 'no rsearch translation')

    def test_untranslatable(self):
        self.bad('--fastx_filter in --gzip_decompress --fastqout o',
                 '--gzip_decompress')
        self.bad('--cluster_fast in --id 0.9 --alnout a', '--alnout with')
        self.bad('--cluster_fast in --id 0.9 --clusters c', '--fasta_width 80')
        self.bad('--fastx_filter in --fastq_minlen 0 --fastqout o', 'empty read')
        self.bad('--fastx_uniques in --minuniquesize 2 --fastaout o',
                 'without --sizeout')
        self.bad('--fastx_filter in --relabel R --fastqout o --fastqout_discarded d',
                 '--relabel')
        self.bad('--fastq_mergepairs in --reverse r --fastq_maxee 1 '
                 '--fastqout_notmerged_fwd n', 'not written as not merged')
        self.bad('--usearch_global q --db d --uc u', 'needs --id')
        self.bad('--usearch_global q --db d --id 1 --userout u', 'userfields')
        self.bad('--usearch_global q --db d --id 1 --fastapairs p', '--fasta_width 0')
        self.bad('--orient q --db d --qmask none --fastaout o', '--qmask none')
        self.bad('--sortbysize in --minseqlength 0 --minsize 2 --output o',
                 'empty sequences')
        self.bad('--fastx_subsample in --sizein --sample_pct 5 --fastaout o',
                 '--by-size')
        self.bad('--uchime_denovo in --fasta_score --nonchimeras n', '--fasta_score')
        self.bad('--fastx_revcomp in --fastq_ascii 64 --fastqout o', '--fastq_ascii')

    def test_caveats(self):
        c = caveats('--fastx_uniques in --minseqlength 0 --fastaout o')
        self.assertTrue(any('--maxseqlength 50000 not applied' in x for x in c))
        self.assertTrue(any('--notrunclabels' in x for x in c))
        self.assertFalse(caveats('--fastx_uniques in --minseqlength 1 --sizein '
                                 '--notrunclabels --fastaout o'))
        self.assertTrue(any('threads' in x for x in caveats(
            '--sintax q --db d --tabbedout t')))
        self.assertFalse(any('threads' in x for x in caveats(
            '--sintax q --db d --tabbedout t --threads 1 --randseed 1')))
        self.assertTrue(any('ignores' in x for x in caveats(
            '--usearch_global q --db d --id 1 --band 3 --uc u')))


if __name__ == '__main__':
    unittest.main()
