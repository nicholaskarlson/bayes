#!/usr/bin/env python3
"""Generate a fictional paired-score teaching file; never use generator truth as evidence."""

import csv
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if len(sys.argv) > 2:
    raise SystemExit("usage: generate_unit04.py [OUTPUT.csv]")
OUT = Path(sys.argv[1]) if len(sys.argv) == 2 else ROOT / "data" / "unit04_reflection_attention_paired.csv"
RNG = random.Random(20260927)

rows = []
for number in range(1, 65):
    before = round(49.0 + RNG.gauss(0, 7.5), 1)
    change = RNG.gauss(2.6, 4.4)
    if number in (17, 42):
        change -= 8.0  # valid but influential fictional observations
    after = round(before + change, 1)
    assert 0 <= before <= 100 and 0 <= after <= 100
    rows.append((f"P{number:03d}", f"{before:.1f}", f"{after:.1f}"))

with OUT.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, lineterminator="\n")
    writer.writerow(("participant_id", "before_score", "after_score"))
    writer.writerows(rows)
print("BAYESIAN_UNIT04_GENERATED", "n=64", str(OUT))
