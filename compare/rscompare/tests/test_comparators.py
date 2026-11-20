"""The comparators' normalisations, on small files."""

import os
import tempfile
import unittest

from rscompare import comparators as C


class Files(unittest.TestCase):

    def setUp(self):
        self.d = tempfile.mkdtemp()

    def f(self, name, text):
        p = os.path.join(self.d, name)
        with open(p, 'wb') as h:
            h.write(text.encode() if isinstance(text, str) else text)
        return p

    def test_bytes_and_missing(self):
        self.assertIsNone(C.FILES['bytes'](self.f('a', 'x\n'), self.f('b', 'x\n')))
        self.assertIn('line 1', C.FILES['bytes'](self.f('a', 'x\n'),
                                                 self.f('b', 'y\n')))
        # a missing file reads as empty
        self.assertIsNone(C.FILES['bytes'](self.f('e', ''), '/nonexistent/x'))

    def test_uc_fast(self):
        v = self.f('v', 'S\t0\t5\t*\t*\t*\t*\t*\ta\t*\nH\t0\t4\t100.0\t+\t0\t0\t=\tb\ta\n')
        r = self.f('r', 'S\t0\t5\t*\t*\t*\t*\t*\ta\t*\nH\t0\t4\t100.0\t+\t0\t0\t4MD\tb\ta\n')
        self.assertIsNone(C.FILES['uc-fast'](v, r))
        self.assertIsNotNone(C.FILES['bytes'](v, r))
        r2 = self.f('r2', 'S\t0\t5\t*\t*\t*\t*\t*\ta\t*\nH\t0\t4\t100.0\t+\t0\t0\t2MI2M\tb\ta\n')
        self.assertIsNotNone(C.FILES['uc-fast'](v, r2))

    def test_biom(self):
        v = self.f('v', '{\n\t"id":"v.biom",\n\t"format":"x",\n\t"generated_by":"vsearch",\n\t"date":"1",\n}\n')
        r = self.f('r', '{\n\t"id":"-",\n\t"format":"x",\n\t"generated_by":"rsearch",\n\t"date":"2",\n}\n')
        self.assertIsNone(C.FILES['biom'](v, r))

    def test_log_report(self):
        log = ('vsearch v1\ncmd\nStarted x\nreport 1\nreport 2\n\n'
               'Finished y\nElapsed\nMax memory\n')
        self.assertEqual(C.log_report(log.encode()), b'report 1\nreport 2\n')
        self.assertIsNone(C.FILES['log-report'](self.f('v', log),
                                                self.f('r', 'report 1\nreport 2\n')))

    def test_udbstats_log(self):
        log = ''.join(f'{k}\n' for k in range(1, 13))
        self.assertEqual(C.udbstats_log(log.encode()), b'4\n5\n6\n7\n8\n')

    def test_logs_by_range(self):
        v = self.f('v', 'x\nForward oriented 3\nReverse 1\nTotal number 4\nFinished\n')
        r = self.f('r', 'Forward oriented 3\nReverse 1\nTotal number 4\n')
        self.assertIsNone(C.FILES['orient-log'](v, r))
        v = self.f('v2', 'a\n  Seqs 3\n  x\n  DBAccel 1\nb\n')
        self.assertIsNone(C.FILES['udbinfo-log'](v, self.f('r2', '  Seqs 3\n  x\n  DBAccel 1\n')))

    def test_alnout_sam(self):
        self.assertIsNone(C.FILES['alnout'](self.f('v', 'vsearch\ncmd\naln\n'),
                                            self.f('r', 'aln\n')))
        v = self.f('v', '@HD\n@SQ\tSN:t\tUR:file:/a/db.fa\n@PG\tvsearch\nq\t0\n')
        r = self.f('r', '@HD\n@SQ\tSN:t\tUR:file:/b/x\n@PG\trsearch\nq\t0\n')
        self.assertIsNone(C.FILES['sam'](v, r))

    def test_upper_and_fasta(self):
        self.assertEqual(C.upper_seq(b'>a b\nacg\n'), b'>a b\nACG\n')
        self.assertEqual(C.upper_seq(b'@a\nacg\n+\nII#\n'), b'@a\nACG\n+\nII#\n')
        fq = b'@a\n' + b'A' * 81 + b'\n+\n' + b'I' * 81 + b'\n@e\n\n+\n\n'
        self.assertEqual(C.fastq_to_fasta(80)(fq),
                         b'>a\n' + b'A' * 80 + b'\nA\n>e\n')
        self.assertEqual(C.fastq_to_fasta(0)(fq), b'>a\n' + b'A' * 81 + b'\n>e\n\n')

    def test_cluster_files(self):
        os.makedirs(os.path.join(self.d, 'v'))
        os.makedirs(os.path.join(self.d, 'r'))
        for side, prefix in (('v', 'c'), ('r', 'clusters')):
            for k in (0, 1):
                self.f(f'{side}/{prefix}{k}', f'>s{k}\nA\n')
        self.assertIsNone(C.cluster_files(os.path.join(self.d, 'v', 'c'),
                                          os.path.join(self.d, 'r', 'clusters')))
        self.f('r/clusters2', '>x\nA\n')
        self.assertIn('files differ', C.cluster_files(
            os.path.join(self.d, 'v', 'c'), os.path.join(self.d, 'r', 'clusters')))


class Errors(unittest.TestCase):

    def test_exact(self):
        self.assertEqual(C.error_kinds('exact', 'x\nFatal error: Bad thing\n',
                                       'Error: Bad thing\n'),
                         ('Bad thing', 'Bad thing'))

    def test_kinds(self):
        v, r = C.error_kinds('abundance',
                             'Fatal error: Invalid (zero) abundance annotation in FASTA file header',
                             'Error: invalid (zero) abundance annotation in read a;size=0')
        self.assertEqual(v, r)
        v, r = C.error_kinds('joinpa',
                             'Fatal error: Invalid (zero) abundance annotation in FASTA file header',
                             'Error: invalid (zero) abundance annotation in read a;size=0')
        self.assertEqual(v, r)
        v, r = C.error_kinds('udb', 'Fatal error: Unable to read from UDB file',
                             'Error: the UDB file is truncated')
        self.assertEqual((v, r), ('truncated', 'truncated'))
        v, r = C.error_kinds('invalid', 'Fatal error: Invalid UDB file header',
                             'Error: invalid header')
        self.assertEqual((v, r), ('invalid', 'invalid'))
        v, r = C.error_kinds('mkudb-exact',
                             'Fatal error: Cannot write a UDB file for an empty database',
                             'Error: no sequences: a UDB file needs at least one')
        self.assertEqual(v, r)

    def test_warnings(self):
        self.assertEqual(C.WARNINGS['nowords'](
            'WARNING: 2 of 5 sequences yielded no k-mer', 'Warning: 2 of 5 sequences give no word'),
            ('2 of 5', '2 of 5'))
        self.assertEqual(C.WARNINGS['missing-abundance'](
            'Missing abundance\nMissing abundance', 'Missing abundance x\nMissing abundance'),
            ('2', '2'))


if __name__ == '__main__':
    unittest.main()
