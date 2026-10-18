#!/usr/bin/env python3
"""Standalone public numerical checks; no unpublished chapter source needed."""
import csv
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / 'scripts'


def need(ok, name):
    if not ok:
        raise SystemExit('BAYES_COMPANION_VERIFY_ERROR ' + name)


def run(name, *args):
    return subprocess.check_output([sys.executable, str(SCRIPTS / name), *args], text=True)


def same(observed, expected, location='record'):
    if isinstance(expected, dict):
        need(isinstance(observed, dict) and set(observed) == set(expected), location+' keys')
        for k in expected:
            same(observed[k], expected[k], location + '.' + k)
    elif isinstance(expected, list):
        need(isinstance(observed, list) and len(observed) == len(expected), location+' length')
        for n, (a, b) in enumerate(zip(observed, expected)):
            same(a, b, location + '[' + str(n) + ']')
    elif isinstance(expected, (float, int)) and not isinstance(expected, bool):
        need(isinstance(observed, (float, int)) and abs(observed - expected) <=
             1e-8 * max(1, abs(expected)), location+' value')
    else:
        need(observed == expected, location+' value')


def main():
    manifest = ROOT / 'SHA256SUMS'
    entries = {}
    for line in manifest.read_text().splitlines():
        match = re.fullmatch(r'([a-f0-9]{64})  (.+)', line)
        need(match is not None and match.group(2) not in entries, 'manifest format')
        entries[match.group(2)] = match.group(1)
    # as_posix(): the manifest uses forward slashes on every OS, Windows included.
    actual = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*')
              if p.is_file() and p != manifest and 'review' not in p.relative_to(ROOT).parts
              and '__pycache__' not in p.parts and p.suffix != '.pyc'}
    need(actual == set(entries), 'file inventory mismatch')
    for name, digest in entries.items():
        p = ROOT / name
        need(p.is_file() and not p.is_symlink() and
             hashlib.sha256(p.read_bytes()).hexdigest() == digest, 'SHA256 ' + name)
    with (ROOT/'data/unit02_help_seeking.csv').open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    need(len(rows) == 12 and sum(int(r['help_sought_within_7_days']) for r in rows) == 8,
         'Unit 2 raw counts')
    expected = json.loads((ROOT/'expected-results/unit02_exact.json').read_text())
    need(expected['posterior'] == [10, 6] and
         abs(expected['posterior_mean'] - float(Fraction(10, 16))) < 1e-12,
         'Unit 2 Beta update')
    text_cases = {
        'inspect_unit01.py': ('pairs=64',),
        'analyze_unit02.py': ('n=12 y=8', 'posterior=Beta(10,6)', 'posterior_mean=0.625'),
        'analyze_unit03.py': ('wilson_95=[0.391,0.862]', 'posterior=Beta(10,6)'),
        'analyze_unit05.py': ('n=72', 'OLS_slope_adjusted=2.766'),
        'analyze_unit06.py': ('n=80 groups=8', 'group_W01'),
        'analyze_unit07.py': ('n=36 desk=8/24 poster=0/12',),
        'analyze_unit08.py': ('UNIT08_ZERO poster=0/12', 'UNIT08_MATCHED'),
    }
    for script, anchors in text_cases.items():
        result = run(script)
        for anchor in anchors:
            need(anchor in result, script + ': ' + anchor)
    json_cases = {
        'analyze_unit04.py': 'unit04_exact.json',
        'analyze_unit04_checks.py': 'unit04_checks.json',
        'simulate_unit04_joint.py': 'unit04_joint_reference.json',
        'analyze_unit04_likelihood.py': 'unit04_likelihood_reference.json',
    }
    for script, name in json_cases.items():
        observed = json.loads(run(script))
        expected = json.loads((ROOT/'expected-results'/name).read_text())
        if script == 'simulate_unit04_joint.py':
            # This field records the software used to make each run. The
            # committed reference was made under a different NumPy version
            # from the pinned companion environment; its version string is
            # provenance, not a simulated statistic or a frozen requirement.
            need(observed['engine'].pop('numpy_version') == np.__version__,
                 'Unit 4 joint runtime metadata')
            need(isinstance(expected['engine'].pop('numpy_version'), str),
                 'Unit 4 joint reference metadata')
        same(observed, expected, script)
    with tempfile.TemporaryDirectory() as temp:
        run('build_unit09_handoff.py', '--outdir', temp)
        report = json.loads((Path(temp)/'unit09_case_evidence.json').read_text())
        need(report['counts'] == {'desk': {'reports': 8, 'n': 24},
                                  'poster': {'reports': 0, 'n': 12}},
             'Unit 9 handoff')
    print('BAYES_COMPANION_VERIFY_OK files=' + str(len(entries)) +
          ' analyses=' + str(len(text_cases)+len(json_cases)) +
          ' unit09=yes models=conditional no_mcmc=yes')


if __name__ == '__main__':
    main()
