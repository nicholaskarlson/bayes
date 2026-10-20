#!/usr/bin/env python3
"""Independent exact check audit from decimal CSV, model card and finite F/t laws."""

import csv
from decimal import Decimal
from hashlib import sha256
import json
import math
from pathlib import Path
import subprocess
import sys

from scipy.stats import f, t

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/unit04_reflection_attention_paired.csv"
MODEL = ROOT / "model-specs/unit04_model.json"
PLAN = ROOT / "model-specs/unit04_check_plan.json"
EXPECTED = ROOT / "expected-results/unit04_checks.json"
ANCHORS = ROOT / "expected-results/unit04_checks_anchors.csv"


def require(ok, message):
    if not ok:
        raise AssertionError("UNIT04_CHECKS_ERROR " + message)


def near(actual, target, name, tolerance=2e-8):
    require(math.isfinite(actual) and abs(actual-target) <= tolerance, name)


record = json.loads(EXPECTED.read_text(encoding="utf-8"))
model = json.loads(MODEL.read_text(encoding="utf-8"))
plan = json.loads(PLAN.read_text(encoding="utf-8"))
for name, path in (("csv_sha256", DATA), ("model_spec_sha256", MODEL), ("check_plan_sha256", PLAN)):
    require(record["source"][name] == sha256(path.read_bytes()).hexdigest(), name)
require(record["source"]["base_commit"] == "b958f9e21194d9eb80e79608a5fee10f858e1aa4", "base")
require(record["source"]["artifact_commit"] is None, "no invented artifact commit")
require(plan["tail_rule"] == "abs(difference)>8" and plan["absolute_tail_cutoff_points"] == 8, "tail rule")
require(plan["model_spec"] == MODEL.name, "check/model link")
require(plan["replicated_cohort_size"] == 64, "replicate size")
with DATA.open(newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle))
require(len(rows)==64, "number of pairs")
require([r["participant_id"] for r in rows] == [f"P{i:03d}" for i in range(1,65)], "IDs")
d = [float(Decimal(r["after_score"])-Decimal(r["before_score"])) for r in rows]
n = len(d)
threshold = Decimal(str(plan["absolute_tail_cutoff_points"]))
count = sum(abs(Decimal(r["after_score"])-Decimal(r["before_score"]))>threshold for r in rows)
prior = model["prior"]
m0, k0, a0, b0 = (float(prior[k]) for k in ("m0","kappa0","a0","b0"))


def update(values):
    size = len(values)
    mean = math.fsum(values)/size
    sum_squares = math.fsum(x*x for x in values)-size*mean*mean
    k = k0+size
    m = m0+(size/k)*(mean-m0)
    a = a0+size/2
    b = b0+(sum_squares+k0*size/k*(mean-m0)**2)/2
    return m,k,a,b,mean,sum_squares/(size-1)


mn,kn,an,bn,observed_mean,observed_var = update(d)
near(record["observed"]["mean"], observed_mean, "observed mean",2e-12)
near(record["observed"]["sample_variance"], observed_var, "observed variance",2e-12)
require(record["observed"]["absolute_tail_count"] == count == 5, "tail count")


def predictive(m,k,a,b):
    # A new cohort shares latent mean and variance. Mean is t, sample variance is scaled F.
    mean_scale = math.sqrt((b/a)*(n+k)/(n*k))
    df = 2*a
    var_scale = b/a
    return {
        "mean_95": [m+mean_scale*t.ppf(q,df) for q in (.025,.975)],
        "mean_cdf_at_observed": t.cdf((observed_mean-m)/mean_scale,df),
        "sample_variance_95": [var_scale*f.ppf(q,n-1,df) for q in (.025,.975)],
        "sample_variance_cdf_at_observed": f.cdf(observed_var/var_scale,n-1,df),
    }


for name, params in (("prior_replicated_cohort",(m0,k0,a0,b0)),
                     ("posterior_replicated_cohort",(mn,kn,an,bn))):
    result = predictive(*params)
    for field,value in result.items():
        recorded = record[name][field]
        if isinstance(value,list):
            require(len(recorded)==2 and value[0]<value[1], name+" "+field+" shape")
            for i in range(2): near(value[i],recorded[i],name+" "+field+" "+str(i))
        else: near(value,recorded,name+" "+field)

scale = math.sqrt((bn/an)*(kn+1)/kn)
prob = t.cdf((-8-mn)/scale,2*an)+t.sf((8-mn)/scale,2*an)
tail = record["individual_predictive_tail"]
near(tail["probability"],prob,"individual tail")
near(tail["expected_count_for_64"],n*prob,"expected count")
require(tail["joint_tail_count_interval"] is None,"joint tail count must remain unclaimed")

loo = []
for i in range(n):
    m,k,a,b,_,_ = update(d[:i]+d[i+1:])
    probability=t.sf((2-m)/math.sqrt(b/(a*k)),2*a)
    loo.append((i+1,probability))
lo,hi=min(loo,key=lambda z:z[1]),max(loo,key=lambda z:z[1])
influence=record["leave_one_out"]
require(influence["pairs_examined"]==64,"LOO coverage")
near(influence["full_p_mu_gt_two"],t.sf((2-mn)/math.sqrt(bn/(an*kn)),2*an),"full probability")
for label,(idx,value) in (("min",lo),("max",hi)):
    require(influence[label]["omitted_id"]==f"P{idx:03d}","LOO identity "+label)
    near(influence[label]["p_mu_gt_two"],value,"LOO value "+label)

with ANCHORS.open(newline="",encoding="utf-8") as handle:
    anchors={r["field"]:float(r["value"]) for r in csv.DictReader(handle)}
expected_anchors={
    "observed_mean":observed_mean,"observed_sample_variance":observed_var,
    "observed_tail_count":count,"prior_mean_lower":predictive(m0,k0,a0,b0)["mean_95"][0],
    "prior_mean_upper":predictive(m0,k0,a0,b0)["mean_95"][1],
    "posterior_mean_lower":predictive(mn,kn,an,bn)["mean_95"][0],
    "posterior_mean_upper":predictive(mn,kn,an,bn)["mean_95"][1],
    "posterior_var_lower":predictive(mn,kn,an,bn)["sample_variance_95"][0],
    "posterior_var_upper":predictive(mn,kn,an,bn)["sample_variance_95"][1],
    "posterior_mean_cdf":predictive(mn,kn,an,bn)["mean_cdf_at_observed"],
    "posterior_var_cdf":predictive(mn,kn,an,bn)["sample_variance_cdf_at_observed"],
    "individual_tail_probability":prob,"expected_tail_count":n*prob,
    "loo_min_id_number":lo[0],"loo_min_probability":lo[1],
    "loo_max_id_number":hi[0],"loo_max_probability":hi[1]
}
require(set(anchors)==set(expected_anchors),"anchor field set")
for key,value in expected_anchors.items(): near(anchors[key],value,"anchor "+key)

primary=json.loads(subprocess.check_output([sys.executable,str(ROOT/"scripts/analyze_unit04_checks.py")],text=True))


def compare_records(got, reference, label):
    if isinstance(got, dict) and isinstance(reference, dict):
        require(set(got)==set(reference),label+" keys")
        for key in got: compare_records(got[key],reference[key],label+"/"+key)
    elif isinstance(got,list) and isinstance(reference,list):
        require(len(got)==len(reference),label+" length")
        for i,(a,b) in enumerate(zip(got,reference)):
            compare_records(a,b,label+"/"+str(i))
    elif isinstance(got,(int,float)) and not isinstance(got,bool) and isinstance(reference,(int,float)) and not isinstance(reference,bool):
        near(float(got),float(reference),label)
    else:
        require(got==reference,label+" value")


compare_records(primary,record,"primary record")
note=(ROOT/"technical-notes/unit04_deep_dive.md").read_text(encoding="utf-8")
for text in ("replicated cohort", "leave-one-out", "not an ordinary p value", "P059", "P042"):
    require(text in note,"missing technical audit safeguard "+text)
print("BAYESIAN_UNIT04_CHECKS_INDEPENDENT_OK n=64 exact_replications=mean,variance influence=64 tail=expected_only")
