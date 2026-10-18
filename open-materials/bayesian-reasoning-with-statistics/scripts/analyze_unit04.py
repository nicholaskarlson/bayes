#!/usr/bin/env python3
"""Primary exact paired frequentist/Bayesian analysis. No MCMC is used."""

import csv
import hashlib
import json
import math
import statistics
import sys
from pathlib import Path

from scipy.stats import t

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "unit04_reflection_attention_paired.csv"
SPEC = ROOT / "model-specs" / "unit04_model.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_differences():
    with DATA.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 64
    assert [r["participant_id"] for r in rows] == [f"P{i:03d}" for i in range(1, 65)]
    before = [float(r["before_score"]) for r in rows]
    after = [float(r["after_score"]) for r in rows]
    assert all(0 <= x <= 100 for x in before + after)
    return [a - b for a, b in zip(after, before)]


def posterior(values, prior, kappa_override=None):
    n = len(values)
    mean = statistics.fmean(values)
    ss = math.fsum((x - mean) ** 2 for x in values)
    m0, a0, b0 = (float(prior[key]) for key in ("m0", "a0", "b0"))
    k0 = float(prior["kappa0"] if kappa_override is None else kappa_override)
    kn = k0 + n
    mn = (k0 * m0 + n * mean) / kn
    an = a0 + n / 2
    bn = b0 + ss / 2 + k0 * n * (mean - m0) ** 2 / (2 * kn)
    df = 2 * an
    scale_mu = math.sqrt(bn / (an * kn))
    scale_new = math.sqrt(bn * (kn + 1) / (an * kn))
    rv = t(df, loc=mn, scale=scale_mu)
    return {
        "kappa0": k0, "kappa_n": kn, "m_n": mn, "a_n": an, "b_n": bn,
        "df": df, "posterior_mu_scale": scale_mu,
        "mu_ci95": [float(rv.ppf(.025)), float(rv.ppf(.975))],
        "p_mu_gt_zero": float(rv.sf(0)),
        "p_mu_gt_two": float(rv.sf(2)),
        "p_mu_in_minus1_plus1": float(rv.cdf(1) - rv.cdf(-1)),
        "predictive_new_scale": scale_new,
        "predictive_new_ci95": [float(x) for x in t.ppf([.025, .975], df, loc=mn, scale=scale_new)]
    }


def analyze():
    values = read_differences()
    spec = json.loads(SPEC.read_text(encoding="utf-8"))
    prior = spec["prior"]
    n = len(values)
    mean = statistics.fmean(values)
    sd = statistics.stdev(values)
    se = sd / math.sqrt(n)
    test = mean / se
    frequentist = {
        "n": n, "mean_difference": mean, "sd_difference": sd, "se": se,
        "t": test, "df": n-1, "p_two_sided": float(2 * t.sf(abs(test), n-1)),
        "ci95": [float(mean + q*se) for q in t.ppf([.025, .975], n-1)],
        "paired_dz": mean / sd,
        "negative_count": sum(x < 0 for x in values),
        "min_difference": min(values), "max_difference": max(values)
    }
    k0 = float(prior["kappa0"])
    a0, b0 = float(prior["a0"]), float(prior["b0"])
    prior_pred_scale = math.sqrt(b0 * (k0+1)/(a0*k0))
    record = {
        "schema_version": "unit04-exact-1",
        "source": {"csv_sha256": digest(DATA), "model_spec_sha256": digest(SPEC),
                   "source_base_commit": "f7efe9c43211bd6ecb52bcda951d6a34a9e7b044",
                   "artifact_commit": None, "record_stage": "pre-commit pilot"},
        "question": {"target": spec["target"], "contrast": spec["contrast"], "design": spec["design"]},
        "frequentist": frequentist,
        "bayesian_exact": {
            "likelihood": spec["likelihood"], "prior": prior,
            "prior_predictive_new_ci95": [float(x) for x in t.ppf([.025,.975], 2*a0, loc=float(prior["m0"]), scale=prior_pred_scale)],
            "primary": posterior(values, prior),
            "sensitivity": [posterior(values, prior, x) for x in spec["sensitivity_kappa0"]]
        },
        "sampler": {"kind": "none", "mcmc_diagnostics": "not applicable; exact conjugate calculation"},
        "claim_boundary": spec["causal_claim_boundary"],
        "reporting_note": "All probabilities conditional on paired Normal likelihood and stated Normal-Inverse-Gamma prior"
    }
    return record


if __name__ == "__main__":
    output = json.dumps(analyze(), indent=2, sort_keys=True, allow_nan=False) + "\n"
    if len(sys.argv) == 1:
        print(output, end="")
    elif len(sys.argv) == 2:
        Path(sys.argv[1]).write_text(output, encoding="utf-8")
    else:
        raise SystemExit("usage: analyze_unit04.py [OUTPUT.json]")
