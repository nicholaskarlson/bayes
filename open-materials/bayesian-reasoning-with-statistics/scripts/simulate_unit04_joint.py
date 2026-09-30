#!/usr/bin/env python3
"""Direct conjugate prior/posterior draws of entire paired-difference cohorts."""

import csv
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/unit04_reflection_attention_paired.csv"
MODEL = ROOT / "model-specs/unit04_model.json"
PLAN = ROOT / "model-specs/unit04_joint_plan.json"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def data():
    with DATA.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 64 and [r["participant_id"] for r in rows] == [f"P{i:03d}" for i in range(1, 65)]
    return np.array([float(r["after_score"]) - float(r["before_score"]) for r in rows])


def skew(values, axis):
    centers = values.mean(axis=axis, keepdims=True)
    centered = values - centers
    variance = np.mean(centered**2, axis=axis)
    third = np.mean(centered**3, axis=axis)
    return third / variance**1.5


def summary_stat(values, plan):
    return {
        "mean": float(np.mean(values)),
        "sample_variance": float(np.var(values, ddof=1)),
        "count_abs_gt_8": int(np.sum(np.abs(values) > plan["absolute_tail_cutoff_points"])),
        "count_negative": int(np.sum(values < 0)),
        "moment_skewness": float(skew(values, axis=0)),
    }


def probability(mask):
    p = float(np.mean(mask))
    return {"estimate": p, "mcse": math.sqrt(p*(1-p)/len(mask))}


def sample_cohorts(parameters, plan, seed, observed):
    rng = np.random.default_rng(seed)
    reps, n = plan["replicates_per_distribution"], plan["cohort_size"]
    m,k,a,b = (parameters[key] for key in ("m","k","a","b"))
    # Gamma(shape=a, scale=1/b) is precision 1/v for IG(shape=a, scale=b).
    variance = 1 / rng.gamma(shape=a, scale=1/b, size=reps)
    mean = rng.normal(loc=m, scale=np.sqrt(variance/k))
    replicated = rng.normal(loc=mean[:,None], scale=np.sqrt(variance[:,None]), size=(reps,n))
    centers = np.mean(replicated, axis=1)
    spread = np.var(replicated, axis=1, ddof=1)
    tails = np.sum(np.abs(replicated)>plan["absolute_tail_cutoff_points"], axis=1)
    negatives = np.sum(replicated<0, axis=1)
    asymmetry = skew(replicated, axis=1)
    large = np.any(np.abs(replicated)>plan["large_change_cutoff_points"], axis=1)
    return {
        "mean_of_replicated_means": float(np.mean(centers)),
        "se_mean_of_replicated_means": float(np.std(centers,ddof=1)/math.sqrt(reps)),
        "mean_of_replicated_sample_variances": float(np.mean(spread)),
        "se_mean_of_replicated_sample_variances": float(np.std(spread,ddof=1)/math.sqrt(reps)),
        "mean_tail_count": float(np.mean(tails)),
        "se_mean_tail_count": float(np.std(tails,ddof=1)/math.sqrt(reps)),
        "mean_negative_count": float(np.mean(negatives)),
        "se_mean_negative_count": float(np.std(negatives,ddof=1)/math.sqrt(reps)),
        "tail_count_95": [int(x) for x in np.quantile(tails,[.025,.975])],
        "moment_skewness_95": [float(x) for x in np.quantile(asymmetry,[.025,.975])],
        "p_tail_count_ge_observed": probability(tails >= observed["count_abs_gt_8"]),
        "p_abs_skew_ge_observed": probability(np.abs(asymmetry) >= abs(observed["moment_skewness"])),
        "p_negative_count_ge_observed": probability(negatives >= observed["count_negative"]),
        "p_any_abs_change_gt_15": probability(large),
    }


def analyze():
    values = data()
    spec = json.loads(MODEL.read_text(encoding="utf-8"))
    plan = json.loads(PLAN.read_text(encoding="utf-8"))
    assert plan["schema_version"] == "unit04-joint-predictive-1"
    assert len(values) == plan["cohort_size"] == 64
    assert plan["absolute_tail_cutoff_points"] == 8 and plan["large_change_cutoff_points"] == 15
    observed = summary_stat(values,plan)
    prior = spec["prior"]
    n = len(values)
    data_mean = float(np.mean(values))
    ss = float(np.sum((values-data_mean)**2))
    m0,k0,a0,b0 = (float(prior[key]) for key in ("m0","kappa0","a0","b0"))
    kn = k0+n;mn=(k0*m0+n*data_mean)/kn;an=a0+n/2
    bn=b0+ss/2+k0*n*(data_mean-m0)**2/(2*kn)
    prior_par=dict(m=m0,k=k0,a=a0,b=b0)
    post_par=dict(m=mn,k=kn,a=an,b=bn)
    result={
        "schema_version":"unit04-joint-reference-1",
        "source":{"csv_sha256":sha(DATA),"model_spec_sha256":sha(MODEL),"plan_sha256":sha(PLAN),
                  "source_base_commit":"d28810c0e33645d4991cbe575ae8de16a46b8947","artifact_commit":None},
        "engine":{"kind":"direct_conjugate_predictive_draws_not_mcmc","numpy_version":np.__version__,
                  "replicates_per_distribution":plan["replicates_per_distribution"],
                  "prior_seed":plan["prior_seed_python"],"posterior_seed":plan["posterior_seed_python"]},
        "observed":observed,
        "prior":sample_cohorts(prior_par,plan,plan["prior_seed_python"],observed),
        "posterior":sample_cohorts(post_par,plan,plan["posterior_seed_python"],observed),
        "analytic_controls":{"prior_replicated_mean":m0,"prior_replicated_variance_mean":b0/(a0-1),
                             "posterior_replicated_mean":mn,"posterior_replicated_variance_mean":bn/(an-1)},
        "interpretation":"Checks reuse the data and are model-dependent. Tail count and skewness are joint replicated statistics; Monte Carlo error is reported. No MCMC convergence or causal claim."
    }
    return result


if __name__ == "__main__":
    output=json.dumps(analyze(),indent=2,sort_keys=True,allow_nan=False)+"\n"
    if len(sys.argv)==1:print(output,end="")
    elif len(sys.argv)==2:Path(sys.argv[1]).write_text(output,encoding="utf-8")
    else:raise SystemExit("usage: simulate_unit04_joint.py [OUTPUT.json]")
