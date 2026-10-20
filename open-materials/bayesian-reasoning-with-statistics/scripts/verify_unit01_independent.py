#!/usr/bin/env python3
"""Audit the Unit 1 reader route using integer tenths, without importing it."""
import csv
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
with (ROOT / "data/unit04_reflection_attention_paired.csv").open(newline="", encoding="utf-8") as stream:
    rows = list(csv.reader(stream))
assert rows[0] == ["participant_id", "before_score", "after_score"]
assert len(rows) == 65
assert [r[0] for r in rows[1:]] == [f"P{i:03d}" for i in range(1, 65)]

def tenths(value):
    whole, dot, fraction = value.partition(".")
    assert dot == "." and len(fraction) == 1 and whole.isdigit()
    return int(whole) * 10 + int(fraction)

d = [tenths(row[2]) - tenths(row[1]) for row in rows[1:]]
assert (sum(d), sum(x < 0 for x in d), sum(x > 0 for x in d), d[0]) == (1405, 20, 44, -2)
result = subprocess.check_output([sys.executable, str(ROOT / "scripts/inspect_unit01.py")], text=True)
expected = {
    "pairs=64", "mean_after_minus_before=2.1953125", "negative_differences=20",
    "positive_differences=44", "P001_after_minus_before=-0.2",
    "scope=description_only; no causal effect or Bayesian model estimated",
}
assert set(result.strip().splitlines()) == expected
assert (ROOT / "figures/unit01_route.svg").is_file()
print("BAYESIAN_UNIT01_INDEPENDENT_OK pairs=64 mean=2.1953125 negative=20")
