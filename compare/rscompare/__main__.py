"""python3 -m rscompare COMMAND ...  (run from compare/)

  run SOURCE... [--seeds 1-50,70] [--jobs N] [--out results.jsonl]
      runs the sources' cases (all: every source) and prints a DIFF line
      for each failure and the tally per source
  translate [--workdir DIR] -- VSEARCH ARGS...
      prints the rsearch pipeline for a vsearch command line
  list
      the sources, their default seeds and what they cover
  summary RESULTS.jsonl...
      outcomes, coverage and caveats of recorded runs
  test
      the unit tests
"""

import argparse
import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # compare/
REPO = os.path.dirname(HERE)


def parse_seeds(text):
    out = []
    for part in text.split(','):
        if '-' in part[1:]:
            a, b = part.split('-', 1)
            out += list(range(int(a), int(b) + 1))
        elif part:
            out.append(int(part))
    return out


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if not argv or argv[0] in ('-h', '--help'):
        print(__doc__)
        return 0
    cmd, rest = argv[0], argv[1:]
    if cmd == 'translate':
        from .translate import translate, TranslationError
        ap = argparse.ArgumentParser(prog='rscompare translate')
        ap.add_argument('--workdir', default='.')
        ap.add_argument('--json', action='store_true')
        if '--' in rest:
            k = rest.index('--')
            opts, vargs = rest[:k], rest[k + 1:]
        else:
            opts, vargs = [], rest
        a = ap.parse_args(opts)
        try:
            p = translate(vargs, a.workdir)
        except TranslationError as e:
            print(f'untranslatable: {e}', file=sys.stderr)
            return 2
        if a.json:
            import json
            print(json.dumps(p.to_json(), indent=1))
            return 0
        print(p.render())
        for o in p.outputs:
            print(f'# --{o.option} {o.vsearch} <- {o.rsearch} ({o.comparator})')
        for c in p.caveats:
            print(f'# caveat: {c}')
        return 0
    if cmd == 'list':
        from .sources import SOURCES
        for n, s in SOURCES.items():
            seeds = s.seeds
            if n == 'datasets':
                from .sources.datasets import default_seeds
                seeds = default_seeds()
            if not seeds:
                print(f'{n:14} {"-":10} {s.about}')
                continue
            rng = f'{seeds[0]}-{seeds[-1]}' if seeds == list(
                range(seeds[0], seeds[-1] + 1)) else ','.join(map(str, seeds))
            print(f'{n:14} {rng:10} {s.about}')
        return 0
    if cmd == 'summary':
        from .runner import summary
        summary(rest)
        return 0
    if cmd == 'test':
        import unittest
        suite = unittest.defaultTestLoader.discover(
            os.path.join(os.path.dirname(os.path.abspath(__file__)), 'tests'),
            top_level_dir=HERE)
        ok = unittest.TextTestRunner(verbosity=1).run(suite).wasSuccessful()
        return 0 if ok else 1
    if cmd == 'run':
        from .sources import SOURCES
        from .runner import run, print_tally
        ap = argparse.ArgumentParser(prog='rscompare run')
        ap.add_argument('sources', nargs='+')
        ap.add_argument('--seeds', type=parse_seeds)
        ap.add_argument('--jobs', '-j', type=int, default=os.cpu_count() or 4)
        ap.add_argument('--vsearch', default=os.path.join(HERE, 'vsearch-build',
                                                          'vsearch'))
        ap.add_argument('--bin', default=os.path.join(REPO, 'target', 'release'))
        ap.add_argument('--out', help='append JSON lines to this file')
        ap.add_argument('--tmp', default=os.path.join(
            os.environ.get('TMPDIR', '/tmp'), 'rscompare'))
        ap.add_argument('--keep', action='store_true',
                        help='keep the files of passing cases too')
        ap.add_argument('--quiet', '-q', action='store_true')
        a = ap.parse_args(rest)
        names = list(SOURCES) if a.sources == ['all'] else a.sources
        for n in names:
            if n not in SOURCES:
                print(f'unknown source {n}; see "list"', file=sys.stderr)
                return 2
        tally = run(names, a.seeds, a.jobs, a.tmp, os.path.abspath(a.vsearch),
                    os.path.abspath(a.bin), a.out, a.keep, a.quiet)
        print_tally(tally)
        bad = sum(t.get('fail', 0) + t.get('error', 0) for t in tally.values())
        return 1 if bad else 0
    print(__doc__)
    return 2


if __name__ == '__main__':
    sys.exit(main())
