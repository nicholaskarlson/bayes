#!/usr/bin/env python3
"""Author-side deterministic fictional workshop score generator; not inferential truth."""
import argparse
import csv
from pathlib import Path
from random import Random

ROOT = Path(__file__).resolve().parents[1]
SIZES = (4, 5, 6, 8, 10, 12, 15, 20)
# Generator-only latent centers, never a reader-facing answer or model prior.
CENTERS = (73, 48, 52, 53, 55, 57, 58, 60)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'data/unit06_workshop_confidence.csv')
    args = parser.parse_args()
    rng = Random(20260928)
    rows = []
    for group, (size, center) in enumerate(zip(SIZES, CENTERS), 1):
        for index in range(1, size + 1):
            value = round(rng.gauss(center, 7.5), 1)
            if not 0 <= value <= 100:
                raise ValueError('Generator made an out-of-range teaching score')
            rows.append((f'W{group:02d}-S{index:02d}', f'W{group:02d}', f'{value:.1f}'))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream, lineterminator='\n')
        writer.writerow(('student_id', 'workshop_id', 'confidence_score'))
        writer.writerows(rows)
    print(f'UNIT06_GENERATOR_OK rows={len(rows)} output={args.output}')

if __name__ == '__main__':
    main()
