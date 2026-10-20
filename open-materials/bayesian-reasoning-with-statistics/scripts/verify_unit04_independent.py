#!/usr/bin/env python3
"""Independent raw-data and exact-reference audit; never imports analyze_unit04."""

import csv
from decimal import Decimal
from hashlib import sha256
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile

from scipy.stats import t as student_t

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "unit04_reflection_attention_paired.csv"
SPEC = ROOT / "model-specs" / "unit04_model.json"
PARAMS = ROOT / "model-specs" / "unit04_parameters.csv"
EXPECTED = ROOT / "expected-results" / "unit04_exact.json"
ANCHORS = ROOT / "expected-results" / "unit04_anchors.csv"


def require(condition, message):
    if not condition:
        raise AssertionError("UNIT04_INDEPENDENT_ERROR " + message)


with tempfile.TemporaryDirectory(prefix="bayesian-unit04-regenerate-") as temp:
    paths = [Path(temp) / name for name in ("a.csv", "b.csv")]
    for p in paths:
        subprocess.run([sys.executable, str(ROOT / "scripts" / "generate_unit04.py"), str(p)],
                       check=True, stdout=subprocess.DEVNULL)
    require(paths[0].read_bytes() == paths[1].read_bytes() == DATA.read_bytes(),
            "deterministic generator does not reproduce committed CSV")
    svg = Path(temp) / "figure.svg"
    subprocess.run([sys.executable, str(ROOT / "scripts" / "make_unit04_figure.py"), str(svg)],
                   check=True, stdout=subprocess.DEVNULL)
    require(svg.read_bytes() == (ROOT / "figures" / "unit04_prior_posterior.svg").read_bytes(),
            "generated SVG differs from committed figure")

record = json.loads(EXPECTED.read_text(encoding="utf-8"))
spec = json.loads(SPEC.read_text(encoding="utf-8"))
dictionary = json.loads((ROOT / "dictionaries" / "unit04_dictionary.json").read_text(encoding="utf-8"))
require(sha256(DATA.read_bytes()).hexdigest() == record["source"]["csv_sha256"], "CSV digest")
require(sha256(SPEC.read_bytes()).hexdigest() == record["source"]["model_spec_sha256"], "model digest")
require(record["source"]["source_base_commit"] == "f7efe9c43211bd6ecb52bcda951d6a34a9e7b044", "base commit")
require(record["source"]["artifact_commit"] is None, "pre-commit record must not invent a commit")
require(record["sampler"]["kind"] == "none", "sampler mislabeled")
require(dictionary["rows"] == 64 and dictionary["synthetic"] is True, "dictionary")

with PARAMS.open(newline="", encoding="utf-8") as handle:
    params = {r["parameter"]: float(r["value"]) for r in csv.DictReader(handle)}
for key in ("m0", "kappa0", "a0", "b0"):
    require(params[key] == float(spec["prior"][key]), "R/Python prior parameter drift: " + key)
require(params["meaningful_change"] == spec["illustrative_meaningful_change_points"], "threshold drift")
require([params["rope_low"], params["rope_high"]] == spec["illustrative_negligible_change_interval_points"], "ROPE drift")

with DATA.open(newline="", encoding="utf-8") as handle:
    reader = csv.DictReader(handle)
    require(reader.fieldnames == ["participant_id", "before_score", "after_score"], "CSV schema")
    rows = list(reader)
require(len(rows) == 64, "number of pairs")
require([r["participant_id"] for r in rows] == [f"P{i:03d}" for i in range(1,65)], "IDs")
for row in rows:
    for field in ("before_score", "after_score"):
        value = Decimal(row[field])
        require(Decimal(0) <= value <= Decimal(100), "score range")
        require(value.as_tuple().exponent == -1, "one-decimal source")
differences = [float(Decimal(r["after_score"]) - Decimal(r["before_score"])) for r in rows]
n = len(differences)
average = math.fsum(differences)/n
ss = math.fsum((d-average)**2 for d in differences)
sd = math.sqrt(ss/(n-1))
se = sd/math.sqrt(n)
t_value = average/se
frequentist = {
    "mean_difference": average, "sd_difference": sd, "se": se, "t": t_value,
    "p_two_sided": 2*student_t.sf(abs(t_value),n-1),
    "ci95_lower": average+student_t.ppf(.025,n-1)*se,
    "ci95_upper": average+student_t.ppf(.975,n-1)*se,
    "paired_dz": average/sd
}
for key, value in frequentist.items():
    observed = record["frequentist"][key] if key in record["frequentist"] else record["frequentist"]["ci95"][0 if key.endswith("lower") else 1]
    require(math.isclose(value, observed, rel_tol=0, abs_tol=2e-12), "frequentist " + key)
require(record["frequentist"]["negative_count"] == sum(d<0 for d in differences), "negative count")


def independently_update(k0):
    m0, a0, b0 = params["m0"], params["a0"], params["b0"]
    precision_count = k0+n
    location = m0 + (n/precision_count)*(average-m0)
    shape = a0+n/2
    scale_parameter = b0 + (ss + (k0*n/precision_count)*(average-m0)**2)/2
    df = 2*shape
    location_scale = math.sqrt(scale_parameter/(shape*precision_count))
    predictive_scale = math.sqrt(scale_parameter/shape*(1+1/precision_count))
    return {
        "m_n": location, "a_n": shape, "b_n": scale_parameter,
        "posterior_mu_scale": location_scale,
        "mu_ci95": [location+location_scale*student_t.ppf(q,df) for q in (.025,.975)],
        "p_mu_gt_zero": student_t.sf((0-location)/location_scale,df),
        "p_mu_gt_two": student_t.sf((2-location)/location_scale,df),
        "p_mu_in_minus1_plus1": student_t.cdf((1-location)/location_scale,df)-student_t.cdf((-1-location)/location_scale,df),
        "predictive_new_ci95": [location+predictive_scale*student_t.ppf(q,df) for q in (.025,.975)]
    }


for expected in record["bayesian_exact"]["sensitivity"]:
    got = independently_update(float(expected["kappa0"]))
    for key, value in got.items():
        target = expected[key]
        if isinstance(value, list):
            require(all(math.isclose(x,y,rel_tol=0,abs_tol=2e-12) for x,y in zip(value,target)), "posterior "+key)
        else:
            require(math.isclose(value,target,rel_tol=0,abs_tol=2e-12), "posterior "+key)

with ANCHORS.open(newline="", encoding="utf-8") as handle:
    anchor_rows = {r["field"]: float(r["value"]) for r in csv.DictReader(handle)}
published = {
    "n": n,
    "mean_difference": frequentist["mean_difference"],
    "t": frequentist["t"],
    "ci95_lower": frequentist["ci95_lower"],
    "ci95_upper": frequentist["ci95_upper"],
    "m_n": independently_update(4)["m_n"],
    "posterior_ci95_lower": independently_update(4)["mu_ci95"][0],
    "posterior_ci95_upper": independently_update(4)["mu_ci95"][1],
    "p_mu_gt_two": independently_update(4)["p_mu_gt_two"]
}
require(set(anchor_rows)==set(published), "anchor field set")
for key,value in published.items():
    require(math.isclose(anchor_rows[key],value,rel_tol=0,abs_tol=2e-10), "anchor "+key)

require((ROOT/"figures"/"unit04_prior_posterior.svg").is_file(), "figure missing")
print("BAYESIAN_UNIT04_INDEPENDENT_OK", "n=64", "models=3")
