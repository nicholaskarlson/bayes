#!/usr/bin/env python3
"""Independent raw-data and simulation audit using Python's standard-library RNG."""

import csv
from decimal import Decimal
from hashlib import sha256
import json
import math
from pathlib import Path
import random
import statistics
import subprocess
import sys
import tempfile

from scipy.stats import t

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/"data/unit04_reflection_attention_paired.csv"
MODEL = ROOT/"model-specs/unit04_model.json"
PLAN = ROOT/"model-specs/unit04_joint_plan.json"
PARAMS = ROOT/"model-specs/unit04_joint_parameters.csv"
REFERENCE = ROOT/"expected-results/unit04_joint_reference.json"
ANCHORS = ROOT/"expected-results/unit04_joint_anchors.csv"
EXACT_CHECKS = ROOT/"expected-results/unit04_checks.json"


def require(ok, message):
    if not ok: raise AssertionError("UNIT04_JOINT_ERROR "+message)


def near(got, expected, name, tol=2e-8):
    require(math.isfinite(got) and abs(got-expected) <= tol,name)


reference=json.loads(REFERENCE.read_text(encoding="utf-8"))
plan=json.loads(PLAN.read_text(encoding="utf-8"))
model=json.loads(MODEL.read_text(encoding="utf-8"))
exact=json.loads(EXACT_CHECKS.read_text(encoding="utf-8"))
with PARAMS.open(newline="",encoding="utf-8") as handle:
    numerical={r["parameter"]:float(r["value"]) for r in csv.DictReader(handle)}
with ANCHORS.open(newline="",encoding="utf-8") as handle:
    anchors={r["field"]:float(r["value"]) for r in csv.DictReader(handle)}
require(numerical=={k:float(plan[v]) for k,v in {
    "replicates_per_distribution":"replicates_per_distribution",
    "absolute_tail_cutoff_points":"absolute_tail_cutoff_points",
    "large_change_cutoff_points":"large_change_cutoff_points",
    "r_seed":"r_seed"}.items()},"R/Python plan drift")
for name,path in (("csv_sha256",DATA),("model_spec_sha256",MODEL),("plan_sha256",PLAN)):
    require(reference["source"][name]==sha256(path.read_bytes()).hexdigest(),"source "+name)
require(reference["source"]["source_base_commit"]=="d28810c0e33645d4991cbe575ae8de16a46b8947","base commit")
require(reference["source"]["artifact_commit"] is None,"invented artifact commit")
require(reference["engine"]["kind"]=="direct_conjugate_predictive_draws_not_mcmc","sampler label")
require(reference["engine"]["replicates_per_distribution"]==30000,"replicate count")
with DATA.open(newline="",encoding="utf-8") as handle:
    rows=list(csv.DictReader(handle))
require(len(rows)==64 and [r["participant_id"] for r in rows]==[f"P{i:03d}" for i in range(1,65)],"data schema")
d=[float(Decimal(r["after_score"])-Decimal(r["before_score"])) for r in rows]
n=len(d);mean=math.fsum(d)/n;ss=math.fsum((v-mean)**2 for v in d)
center=[v-mean for v in d]
skew=(math.fsum(v**3 for v in center)/n)/((math.fsum(v*v for v in center)/n)**1.5)
observed={"mean":mean,"sample_variance":ss/(n-1),
          "count_abs_gt_8":sum(abs(v)>8 for v in d),"count_negative":sum(v<0 for v in d),
          "moment_skewness":skew}
for key,value in observed.items():near(float(reference["observed"][key]),float(value),"observed "+key,2e-12)
prior=model["prior"]
m0,k0,a0,b0=(float(prior[k]) for k in ("m0","kappa0","a0","b0"))
kn=k0+n;mn=m0+n/kn*(mean-m0);an=a0+n/2
bn=b0+(ss+k0*n/kn*(mean-m0)**2)/2
controls={"prior_replicated_mean":m0,"prior_replicated_variance_mean":b0/(a0-1),
          "posterior_replicated_mean":mn,"posterior_replicated_variance_mean":bn/(an-1)}
for key,value in controls.items():near(value,reference["analytic_controls"][key],"control "+key)
near(exact["individual_predictive_tail"]["expected_count_for_64"],
     n*(t.cdf((-8-mn)/math.sqrt(bn/an*(1+1/kn)),2*an)+
        t.sf((8-mn)/math.sqrt(bn/an*(1+1/kn)),2*an)),"exact tail control")


def independent_draws(params,seed,reps=6000):
    rng=random.Random(seed)
    m,k,a,b=params
    means=[];variances=[];tails=[];negatives=[];skews=[];large=[]
    for _ in range(reps):
        precision=rng.gammavariate(a,1/b)
        variance=1/precision
        latent_mean=rng.gauss(m,math.sqrt(variance/k))
        draws=[rng.gauss(latent_mean,math.sqrt(variance)) for _ in range(n)]
        xbar=math.fsum(draws)/n
        residuals=[x-xbar for x in draws]
        m2=math.fsum(x*x for x in residuals)/n
        means.append(xbar);variances.append(m2*n/(n-1))
        tails.append(sum(abs(x)>8 for x in draws))
        negatives.append(sum(x<0 for x in draws))
        skews.append((math.fsum(x**3 for x in residuals)/n)/m2**1.5)
        large.append(any(abs(x)>15 for x in draws))
    return dict(mean_of_replicated_means=(statistics.fmean(means),statistics.stdev(means)/math.sqrt(reps)),
                mean_of_replicated_sample_variances=(statistics.fmean(variances),statistics.stdev(variances)/math.sqrt(reps)),
                mean_tail_count=(statistics.fmean(tails),statistics.stdev(tails)/math.sqrt(reps)),
                mean_negative_count=(statistics.fmean(negatives),statistics.stdev(negatives)/math.sqrt(reps)),
                p_tail_count_ge_observed=(statistics.fmean(x>=observed["count_abs_gt_8"] for x in tails),None),
                p_abs_skew_ge_observed=(statistics.fmean(abs(x)>=abs(skew) for x in skews),None),
                p_negative_count_ge_observed=(statistics.fmean(x>=observed["count_negative"] for x in negatives),None),
                p_any_abs_change_gt_15=(statistics.fmean(large),None))


for label,params,seed in (("prior",(m0,k0,a0,b0),plan["independent_python_seed"]),
                          ("posterior",(mn,kn,an,bn),plan["independent_python_seed"]+1)):
    independent=independent_draws(params,seed)
    committed=reference[label]
    for name,(value,se) in independent.items():
        if name.startswith("p_"):
            target=committed[name]
            near(target["mcse"],math.sqrt(target["estimate"]*(1-target["estimate"])/30000),label+" "+name+" MCSE")
            combined=math.sqrt(target["mcse"]**2+value*(1-value)/6000)
            require(abs(value-target["estimate"])<=6*combined+0.01,label+" "+name)
        else:
            target=committed[name]
            require(abs(value-target)<=6*math.hypot(se,committed["se_"+name])+0.05,label+" "+name)
    mean_target=controls[label+"_replicated_mean"]
    var_target=controls[label+"_replicated_variance_mean"]
    require(abs(committed["mean_of_replicated_means"]-mean_target)<=6*committed["se_mean_of_replicated_means"]+0.005,label+" exact mean")
    require(abs(committed["mean_of_replicated_sample_variances"]-var_target)<=6*committed["se_mean_of_replicated_sample_variances"]+0.02,label+" exact variance")
near(reference["posterior"]["mean_tail_count"],exact["individual_predictive_tail"]["expected_count_for_64"],
     "joint tail expectation",6*reference["posterior"]["se_mean_tail_count"]+0.02)
negative=t.cdf((0-mn)/math.sqrt(bn/an*(1+1/kn)),2*an)
near(reference["posterior"]["mean_negative_count"],n*negative,"negative expectation",
     6*reference["posterior"]["se_mean_negative_count"]+0.02)

anchor_expected={}
for label in ("prior","posterior"):
    item=reference[label]
    for name in ("mean_of_replicated_means","se_mean_of_replicated_means",
                 "mean_of_replicated_sample_variances","se_mean_of_replicated_sample_variances",
                 "mean_tail_count","se_mean_tail_count","mean_negative_count","se_mean_negative_count"):
        anchor_expected[label+"_"+name]=item[name]
    for name in ("p_tail_count_ge_observed","p_abs_skew_ge_observed","p_negative_count_ge_observed","p_any_abs_change_gt_15"):
        anchor_expected[label+"_"+name]=item[name]["estimate"]
        anchor_expected[label+"_"+name+"_mcse"]=item[name]["mcse"]
require(set(anchor_expected)==set(anchors),"R anchors field set")
for key,value in anchor_expected.items():near(anchors[key],value,"R anchor "+key)

current=json.loads(subprocess.check_output([sys.executable,str(ROOT/"scripts/simulate_unit04_joint.py")],text=True))
require(current["engine"]["kind"]==reference["engine"]["kind"],"primary engine label")
for label in ("prior","posterior"):
    for key in ("p_tail_count_ge_observed","p_abs_skew_ge_observed","p_any_abs_change_gt_15"):
        p=current[label][key]["estimate"];q=reference[label][key]["estimate"]
        require(abs(p-q)<=6*math.hypot(current[label][key]["mcse"],reference[label][key]["mcse"])+.005,"current primary "+label+key)
note=(ROOT/"technical-notes/unit04_deep_dive.md").read_text(encoding="utf-8")
for safeguard in ("shared", "joint simulation", "skewness", "direct conjugate", "not MCMC", "figures/unit04_joint_checks.svg"):
    require(safeguard in note,"technical audit "+safeguard)
with tempfile.TemporaryDirectory() as directory:
    generated=Path(directory)/"joint.svg"
    subprocess.check_call([sys.executable,str(ROOT/"scripts/make_unit04_joint_figure.py"),str(generated)])
    figure=ROOT/"figures/unit04_joint_checks.svg"
    require(generated.read_bytes()==figure.read_bytes(),"joint SVG drift")
    require('aria-labelledby="title desc"' in figure.read_text(encoding="utf-8"),"figure description")
print("BAYESIAN_UNIT04_JOINT_INDEPENDENT_OK cohorts=30000+30000 standard_library_audit=6000+6000")
