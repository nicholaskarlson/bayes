#!/usr/bin/env python3
"""Build a deterministic, auditable teaching report from the Unit 7 case."""
import argparse
import csv
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/unit07_help_seeking_routes.csv'
REFERENCE=ROOT/'expected-results/unit07_proportions.json'
MODEL=ROOT/'model-specs/unit07_proportions.json'
PLAN=ROOT/'model-specs/unit09_handoff.json'

def digest(path):return sha256(path.read_bytes()).hexdigest()
def need(ok,message):
    if not ok:raise SystemExit('UNIT09_HANDOFF_ERROR '+message)
def probabilities(a,b,n):
    prior=Fraction(1)
    for i in range(n):prior*=Fraction(b+i,a+b+i)
    return float(prior)
def build(outdir):
    with DATA.open(newline='',encoding='utf-8') as h:rows=list(csv.DictReader(h))
    need(len(rows)==36 and list(rows[0])==['volunteer_id','recruitment_route','reported_contact_7d'],'CSV schema')
    need(len({r['volunteer_id'] for r in rows})==36,'unique IDs')
    need(all(r['recruitment_route'] in ('desk','poster') and r['reported_contact_7d'] in ('0','1') for r in rows),'values')
    counts={route:(len(subset),sum(int(r['reported_contact_7d']) for r in subset))
            for route in ('desk','poster') for subset in [[r for r in rows if r['recruitment_route']==route]]}
    need(counts=={'desk':(24,8),'poster':(12,0)},'counts')
    ref=json.loads(REFERENCE.read_text());model=json.loads(MODEL.read_text());plan=json.loads(PLAN.read_text())
    need(ref['spec_id']==model['spec_id']==plan['case'],'spec identity')
    need(ref['data_sha256']==digest(DATA),'data identity')
    need((ref['n'],ref['desk_n'],ref['desk_y'],ref['poster_n'],ref['poster_y'])==(36,24,8,12,0),'reference counts')
    for label,prior in [('primary',1),('sensitivity',2)]:
        v=ref[label]
        need(v['prior']==[prior,prior] and v['desk_posterior']==[prior+8,prior+16] and
             v['poster_posterior']==[prior,prior+12],'posterior shapes '+label)
        need(abs(v['desk_mean']-float(Fraction(prior+8,2*prior+24)))<1e-13,'desk mean '+label)
        need(abs(v['poster_mean']-float(Fraction(prior,2*prior+12)))<1e-13,'poster mean '+label)
        need(abs(v['prior_predictive_poster_zero']-probabilities(prior,prior,12))<1e-13,'prior predictive '+label)
        need(abs(v['posterior_predictive_poster_zero']-probabilities(prior,prior+12,12))<1e-13,'posterior predictive '+label)
    p=ref['primary'];s=ref['sensitivity']
    evidence={'schema':'unit09-handoff-evidence-v1','sources':{
        'data_sha256':digest(DATA),'unit07_reference_sha256':digest(REFERENCE),
        'unit07_model_sha256':digest(MODEL),'unit09_plan_sha256':digest(PLAN)},
        'counts':{'desk':{'reports':8,'n':24},'poster':{'reports':0,'n':12}},
        'raw_desk_minus_poster':float(Fraction(8,24)),
        'primary':{key:p[key] for key in ('prior','desk_posterior','poster_posterior','desk_mean','poster_mean','difference_mean','difference_95','P_difference_gt_zero','P_difference_gt_0_10','prior_predictive_poster_zero','posterior_predictive_poster_zero')},
        'sensitivity':{key:s[key] for key in ('prior','desk_posterior','poster_posterior','difference_mean','difference_95','P_difference_gt_zero','P_difference_gt_0_10','prior_predictive_poster_zero','posterior_predictive_poster_zero')},
        'scope':plan['scope'],'computation':'exact, deterministic Unit 7 reference independently audited in Python and base R; no MCMC; Unit 9 does not reintegrate the difference tail'}
    report=f'''# Fictional route case — reproducible research handoff

**Status.** Teaching report generated from unchanged Unit 7 records and previously independently audited exact results. This is neither a preregistered real study nor an intervention recommendation. The accompanying evidence JSON records full-precision values and SHA-256 source identities.

## Question and study card

Among 36 fictional convenience volunteers, how do seven-day self-reported contact frequencies differ by recruitment route? One row is one volunteer. Route was not assigned, there is no campus sampling frame, and the report is not a validated service log. The modeled target is θ_D−θ_P for comparable volunteers under this regime; it is not the effect of changing route.

## Data and familiar description

Desk: **8/24** reported contact (0.333). Poster: **0/12** (0.000). The recorded desk-minus-poster gap is **33.3 percentage points**. The relative change using a zero poster fraction as denominator is undefined. These counts describe the file without a prior or a population claim.

## Model, prior and conditional result

Within each route, the working likelihood treats binary reports as independent Bernoulli observations with one stable route probability. Separate Beta(1,1) priors update to Beta(9,17) for desk and Beta(1,13) for poster. The posterior means are **{p['desk_mean']:.3f}** and **{p['poster_mean']:.3f}**. The posterior mean of θ_D−θ_P is **{p['difference_mean']:.3f}**; its central 95% credible interval is **[{p['difference_95'][0]:.3f}, {p['difference_95'][1]:.3f}]**. P(θ_D−θ_P>0 | model, records) is **{p['P_difference_gt_zero']:.3f}**. At an *illustrative*, unvalidated 10-percentage-point magnitude, P(θ_D−θ_P>0.10) is **{p['P_difference_gt_0_10']:.3f}**. These conditional probabilities concern association between routes under a specified model.

## Checks and sensitivity

Before this count, the Beta(1,1) model assigned probability **{p['prior_predictive_poster_zero']:.3f}** to zero reports among twelve poster volunteers. After these twelve zeros, its posterior predictive probability of zero in another comparable twelve is **{p['posterior_predictive_poster_zero']:.3f}**. Neither number assesses held-out generalization or validates self-report and selection. Under alternative independent Beta(2,2) priors, mean θ_D−θ_P is **{s['difference_mean']:.3f}**, its central 95% interval is **[{s['difference_95'][0]:.3f}, {s['difference_95'][1]:.3f}]**, P(θ_D−θ_P>0) is **{s['P_difference_gt_zero']:.3f}**, and P(θ_D−θ_P>0.10) is **{s['P_difference_gt_0_10']:.3f}**. This prior assigned probability **{s['prior_predictive_poster_zero']:.3f}** to the zero count before the case; it deserves closer scrutiny. Do not select a prior solely for a desired interval.

## Reportable interpretation and next evidence

“Among these 36 fictional convenience volunteers, eight of 24 desk recruits and none of twelve poster recruits reported contact within seven days. Under separate Beta(1,1) priors and a within-route Bernoulli model, the posterior mean desk-minus-poster reporting difference is {p['difference_mean']:.3f} with a 95% credible interval [{p['difference_95'][0]:.3f}, {p['difference_95'][1]:.3f}]; the Beta(2,2) sensitivity interval includes zero.” This is a conditional association, not an effect or a campus rate. Verify recruitment, missingness, measurement and prior information; if a route *effect* is sought, define the intervention, target population and assignment/comparison design before collecting new data. No policy decision follows from this exercise.

## Reproduction and provenance

From the repository root run `python3 open-materials/bayesian-reasoning-with-statistics/scripts/analyze_unit07.py`, `Rscript open-materials/bayesian-reasoning-with-statistics/scripts/analyze_unit07.R`, `make verify`, and `Rscript open-materials/bayesian-reasoning-with-statistics/scripts/verify_unit07_in_R.R`. Generate this report with `python3 open-materials/bayesian-reasoning-with-statistics/scripts/build_unit09_handoff.py`. The report uses the committed Unit 7 exact reference, checked by the public Python and base-R routes for its interval and tail probabilities; this builder audits identities and simple arithmetic but does not reintegrate those tails. See `technical-notes/unit09_handoff_audit.md` for the chain of evidence. Computation agreement cannot repair recruitment, outcome measurement or a causal design.
'''
    outdir.mkdir(parents=True,exist_ok=True)
    (outdir/'unit09_case_evidence.json').write_text(json.dumps(evidence,sort_keys=True,indent=2)+'\n',encoding='utf-8')
    (outdir/'unit09_case_report.md').write_text(report,encoding='utf-8')
    print('UNIT09_HANDOFF_BUILD_OK report='+str(outdir/'unit09_case_report.md')+' evidence='+str(outdir/'unit09_case_evidence.json'))
    return evidence,report

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir',type=Path,default=ROOT/'review')
    args=parser.parse_args()
    build(args.outdir)
