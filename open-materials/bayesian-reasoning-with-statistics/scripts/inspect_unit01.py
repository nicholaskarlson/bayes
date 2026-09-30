#!/usr/bin/env python3
"""Read the unchanged Unit 4 teaching CSV for Unit 1's descriptive handoff."""
import csv
from decimal import Decimal
from pathlib import Path

RESOURCE = Path(__file__).resolve().parents[1]
DATA = RESOURCE / "data/unit04_reflection_attention_paired.csv"

with DATA.open(newline="", encoding="utf-8") as handle:
    reader = csv.DictReader(handle)
    if reader.fieldnames != ["participant_id", "before_score", "after_score"]:
        raise SystemExit("UNIT01_INPUT_ERROR unexpected CSV columns")
    rows = list(reader)

if len(rows) != 64 or [r["participant_id"] for r in rows] != [f"P{i:03d}" for i in range(1, 65)]:
    raise SystemExit("UNIT01_INPUT_ERROR expected 64 unique ordered pairs")
differences = []
for row in rows:
    before, after = Decimal(row["before_score"]), Decimal(row["after_score"])
    if not (0 <= before <= 100 and 0 <= after <= 100):
        raise SystemExit("UNIT01_INPUT_ERROR rating outside illustrative scale")
    differences.append(after - before)
mean = sum(differences) / len(differences)
print(f"pairs={len(rows)}")
print(f"mean_after_minus_before={mean:.7f}")
print(f"negative_differences={sum(d < 0 for d in differences)}")
print(f"positive_differences={sum(d > 0 for d in differences)}")
print(f"P001_after_minus_before={differences[0]:.1f}")
print("scope=description_only; no causal effect or Bayesian model estimated")
