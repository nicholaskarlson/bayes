#!/usr/bin/env python3
"""Audit the illustrative finite grid with exact rational arithmetic."""
import csv
from fractions import Fraction as F
import hashlib
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CSV = ROOT / 'expected-results/prelude_grid.csv'
DATA = ROOT / 'data/unit02_help_seeking.csv'
FIELDS = ('theta', 'prior', 'relative_fit', 'product', 'posterior')


def require(ok, label):
    if not ok:
        raise SystemExit('PRELUDE_GRID_VERIFY_ERROR ' + label)


def near(printed, exact):
    return abs(F(printed) - exact) <= F(1, 2_000_000)


def main():
    with DATA.open(newline='', encoding='utf-8') as stream:
        observations = list(csv.DictReader(stream))
    require(len(observations) == 12 and
            [row['participant_id'] for row in observations] ==
            [f'H{i:03d}' for i in range(1, 13)] and
            all(row['help_sought_within_7_days'] in ('0', '1') for row in observations) and
            sum(int(row['help_sought_within_7_days']) for row in observations) == 8,
            'source CSV n=12, y=8, IDs and binary values')
    expected_hash = next((line.split()[0] for line in
                          (ROOT / 'SHA256SUMS').read_text(encoding='utf-8').splitlines()
                          if line.endswith('  data/unit02_help_seeking.csv')), None)
    require(expected_hash == hashlib.sha256(DATA.read_bytes()).hexdigest(),
            'source CSV SHA256SUMS')
    with CSV.open(newline='', encoding='utf-8') as stream:
        reader = csv.DictReader(stream)
        require(tuple(reader.fieldnames or ()) == FIELDS, 'CSV columns')
        rows = list(reader)
    theta = [F(i, 10) for i in range(1, 10)]
    require(len(rows) == 9, 'nine candidates')
    raw_prior = [t * (1 - t) for t in theta]
    prior = [x / sum(raw_prior) for x in raw_prior]
    fit = [(t / F(7, 10)) ** 8 * ((1-t) / F(3, 10)) ** 4 for t in theta]
    product = [a * b for a, b in zip(prior, fit)]
    posterior = [x / sum(product) for x in product]
    expected = [theta, prior, fit, product, posterior]
    for i, row in enumerate(rows):
        for field, values in zip(FIELDS, expected):
            require(near(row[field], values[i]), f'row {i+1} {field}')
    require(near('0.432173', sum(product)), 'normalizer')
    require(near('0.745038', sum(p for t, p in zip(theta, posterior) if t > F(1, 2))),
            'grid tail')
    continuous_tail = F(sum(comb(15, k) for k in range(10)), 2**15)
    require(near('0.849121', continuous_tail), 'continuous beta tail')
    require(continuous_tail != sum(p for t, p in zip(theta, posterior) if t > F(1, 2)),
            'grid is not continuous beta')
    print('PRELUDE_GRID_INDEPENDENT_OK candidates=9 grid_gt_half=0.745038 '
          'continuous_gt_half=0.849121')


if __name__ == '__main__':
    main()
