#!/usr/bin/env python3
"""Exact replicated-cohort checks and leave-one-out influence; no simulation/MCMC."""

import csv
import hashlib
import json
import math
import statistics
import sys
from pathlib import Path

from scipy.stats import f, t

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/unit04_reflection_attention_paired.csv"
MODEL = ROOT / "model-specs/unit04_model.json"
CHECKS = ROOT / "model-specs/unit04_check_plan.json"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_data():
    with DATA.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 64
    assert [r["participant_id"] for r in rows] == [f"P{i:03d}" for i in range(1, 65)]
    return [float(r["after_score"]) - float(r["before_score"]) for r in rows]


def update(values, prior):
    n = len(values)
    mean = statistics.fmean(values)
    ss = math.fsum((v - mean) ** 2 for v in values)
    m0, k0, a0, b0 = (float(prior[k]) for k in ("m0", "kappa0", "a0", "b0"))
    kn = k0 + n
    mn = (k0*m0+n*mean)/kn
    an = a0+n/2
    bn = b0+ss/2+k0*n*(mean-m0)**2/(2*kn)
    mu_scale = math.sqrt(bn/(an*kn))
    return dict(kn=kn, mn=mn, an=an, bn=bn,
                p_mu_gt_two=float(t.sf((2-mn)/mu_scale, 2*an)))


def replicated(n, kappa, a, b, location, observed_mean, observed_variance):
    """Integrate over shared mean and variance, not independent marginal draws."""
    df = 2*a
    mean_scale = math.sqrt((b/a)*(1/kappa+1/n))
    variance_scale = b/a
    return {
        "mean_95": [float(x) for x in t.ppf([.025, .975], df, loc=location, scale=mean_scale)],
        "mean_cdf_at_observed": float(t.cdf(observed_mean, df, loc=location, scale=mean_scale)),
        "sample_variance_95": [float(x*variance_scale) for x in f.ppf([.025, .975], n-1, df)],
        "sample_variance_cdf_at_observed": float(f.cdf(observed_variance/variance_scale, n-1, df)),
    }


def analyze():
    values = read_data()
    model = json.loads(MODEL.read_text(encoding="utf-8"))
    plan = json.loads(CHECKS.read_text(encoding="utf-8"))
    assert plan["schema_version"] == "unit04-check-plan-1"
    assert plan["model_spec"] == MODEL.name
    assert plan["replicated_cohort_size"] == len(values)
    assert plan["tail_rule"] == "abs(difference)>8"
    threshold = float(plan["absolute_tail_cutoff_points"])
    assert threshold == 8
    prior = model["prior"]
    n = len(values)
    mean = statistics.fmean(values)
    variance = statistics.variance(values)
    post = update(values, prior)
    leave_one_out = [(f"P{i+1:03d}", update(values[:i]+values[i+1:], prior)["p_mu_gt_two"])
                     for i in range(n)]
    tail_scale = math.sqrt(post["bn"]/post["an"]*(1+1/post["kn"]))
    tail_prob = float(t.cdf((-threshold-post["mn"])/tail_scale, 2*post["an"])
                      + t.sf((threshold-post["mn"])/tail_scale, 2*post["an"]))
    return {
        "schema_version": "unit04-exact-checks-1",
        "source": {"csv_sha256": sha(DATA), "model_spec_sha256": sha(MODEL),
                   "check_plan_sha256": sha(CHECKS), "base_commit": "b958f9e21194d9eb80e79608a5fee10f858e1aa4",
                   "artifact_commit": None},
        "observed": {"n": n, "mean": mean, "sample_variance": variance,
                     "absolute_tail_count": sum(abs(v)>threshold for v in values)},
        "prior_replicated_cohort": replicated(n, float(prior["kappa0"]), float(prior["a0"]),
                                               float(prior["b0"]), float(prior["m0"]), mean, variance),
        "posterior_replicated_cohort": replicated(n, post["kn"], post["an"], post["bn"],
                                                   post["mn"], mean, variance),
        "individual_predictive_tail": {"cutoff_points": threshold, "probability": tail_prob,
                                        "expected_count_for_64": n*tail_prob,
                                        "joint_tail_count_interval": None},
        "leave_one_out": {"full_p_mu_gt_two": post["p_mu_gt_two"],
                          "min": dict(zip(("omitted_id", "p_mu_gt_two"), min(leave_one_out, key=lambda x:x[1]))),
                          "max": dict(zip(("omitted_id", "p_mu_gt_two"), max(leave_one_out, key=lambda x:x[1]))),
                          "pairs_examined": n},
        "interpretation": "Model-dependent replicated-data checks, not ordinary p values; expected tail count is not a joint count check; omission is influence, not correction for selection or measurement."
    }


if __name__ == "__main__":
    result = json.dumps(analyze(), indent=2, sort_keys=True, allow_nan=False) + "\n"
    if len(sys.argv) == 1:
        print(result, end="")
    elif len(sys.argv) == 2:
        Path(sys.argv[1]).write_text(result, encoding="utf-8")
    else:
        raise SystemExit("usage: analyze_unit04_checks.py [OUTPUT.json]")
