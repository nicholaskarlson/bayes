#!/usr/bin/env python3
"""Independent standard-library scalar audit of the eight-group Normal hierarchy."""
import csv
import hashlib
import json
from math import sqrt
from pathlib import Path
from statistics import NormalDist
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[1]
CSV=ROOT/'data/unit06_workshop_confidence.csv'
REFERENCE=ROOT/'expected-results/unit06_partial_pooling.json'
ANCHORS=ROOT/'expected-results/unit06_anchors.csv'
SPEC=ROOT/'model-specs/unit06_partial_pooling.json'
NORM=NormalDist()

def near(x,y,tol=1e-8):
    assert abs(x-y)<tol,(x,y)

records=list(csv.DictReader(CSV.open(newline='',encoding='utf-8')))
labels=[f'W{k:02d}' for k in range(1,9)]
counts=[4,5,6,8,10,12,15,20]
assert len(records)==80 and list(records[0])==['student_id','workshop_id','confidence_score']
assert [(r['student_id'],r['workshop_id']) for r in records]==[(f'{g}-S{i:02d}',g) for g,n in zip(labels,counts) for i in range(1,n+1)]
y=[[float(r['confidence_score']) for r in records if r['workshop_id']==g] for g in labels]
assert all(len(z)==n and all(0<=v<=100 for v in z) for n,z in zip(counts,y))
means=[sum(z)/len(z) for z in y]
sse=sum(sum((v-m)**2 for v in z) for z,m in zip(y,means))
spec=json.loads(SPEC.read_text(encoding='utf-8'));ref=json.loads(REFERENCE.read_text(encoding='utf-8'))
assert spec['spec_id']==ref['spec_id']=='unit06-workshop-partial-pooling-v1'
assert 'tau=5' in spec['group_prior'] and 'sigma=8' in spec['likelihood']
assert ref['data_sha256']==hashlib.sha256(CSV.read_bytes()).hexdigest()
with tempfile.TemporaryDirectory() as tmp:
    for k in range(2):
        path=Path(tmp)/f'copy{k}.csv'
        subprocess.run([sys.executable,str(ROOT/'scripts/generate_unit06.py'),'--output',str(path)],check=True,capture_output=True)
        assert path.read_bytes()==CSV.read_bytes()
near(ref['pooled_mean'],sum(map(sum,y))/80)
near(ref['within_sse'],sse);near(ref['within_sd'],sqrt(sse/72))
near(ref['replicated_expected_within_sse'],72*64)
assert ref['within_df']==72 and ref['n']==80 and ref['sizes']==dict(zip(labels,counts))
for g,m in zip(labels,means):near(ref['raw_means'][g],m)
assert ref['score_range']==[min(map(min,y)),max(map(max,y))]
prior_sd=sqrt(100+25+64/4)
for x,p in zip(ref['prior_predictive_W01_group_mean_95'],[.025,.975]):near(x,55+NORM.inv_cdf(p)*prior_sd)
score_sd=sqrt(100+25+64)
near(ref['prior_predictive_score_outside_0_100'],NORM.cdf(-55/score_sd)+NORM.cdf(-45/score_sd))
near(ref['conditional_prior_icc'],25/89)

with ANCHORS.open(newline='', encoding='utf-8') as stream: anchors=list(csv.DictReader(stream))
assert len(anchors)==24
for case,tau in [('primary_tau5',5),('sensitivity_tau2',2),('sensitivity_tau10',10)]:
    # Marginalize each theta_j in ybar_j | mu, then integrate mu analytically.
    noise=[tau*tau+64/n for n in counts]
    mu_var=1/(.01+sum(1/v for v in noise))
    mu=mu_var*(.55+sum(m/v for m,v in zip(means,noise)))
    weight=[tau*tau/v for v in noise]
    theta=[w*m+(1-w)*mu for w,m in zip(weight,means)]
    conditional=[1/(1/tau**2+n/64) for n in counts]
    vtheta=[v+(1-w)**2*mu_var for v,w in zip(conditional,weight)]
    cov_AH=(1-weight[0])*(1-weight[7])*mu_var
    sd_diff=sqrt(vtheta[0]+vtheta[7]-2*cov_AH)
    diff=theta[0]-theta[7]
    p=NORM.cdf(diff/sd_diff)
    target=ref[case]
    near(target['mu_mean'],mu);near(target['mu_sd'],sqrt(mu_var))
    near(target['A_minus_H_mean'],diff);near(target['A_minus_H_sd'],sd_diff)
    near(target['P_A_gt_H'],p)
    case_anchors=[a for a in anchors if a['case']==case]
    assert [a['group'] for a in case_anchors]==labels
    for j,g in enumerate(labels):
        lo=theta[j]+NORM.inv_cdf(.025)*sqrt(vtheta[j])
        hi=theta[j]+NORM.inv_cdf(.975)*sqrt(vtheta[j])
        for value,actual in [(target['theta_mean'][j],theta[j]),(target['theta_sd'][j],sqrt(vtheta[j])),
                             (target['theta_95'][j][0],lo),(target['theta_95'][j][1],hi),
                             (target['weights'][j],weight[j])]:near(value,actual)
        row=case_anchors[j]
        assert int(row['n'])==counts[j]
        for field,actual in [('raw',means[j]),('posterior',theta[j]),('lower',lo),('upper',hi),
                             ('weight',weight[j]),('mu',mu),('p_A_gt_H',p)]:near(float(row[field]),actual)
assert ref['computation']=='exact joint Normal integration with fixed sigma and tau; no MCMC'
assert (ROOT/'figures/unit06_pooling.svg').exists()
output=subprocess.check_output([sys.executable,str(ROOT/'scripts/analyze_unit06.py')],text=True)
for term in ('n=80 groups=8','posterior_mu=55.032','group_W01 n=4 raw=63.625 posterior=60.272',
             'group_W08 n=20 raw=59.345 posterior=58.856','tau=2 P_A_gt_H=0.401',
             'tau=5 P_A_gt_H=0.653','observed_SSE=4540.640 replicated_expected_SSE=4608.000',
             'exact no MCMC'):
    assert term in output,term
print('BAYESIAN_UNIT06_INDEPENDENT_OK n=80 groups=8 exact_hierarchy=yes')
