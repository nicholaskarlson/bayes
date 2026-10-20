#!/usr/bin/env python3
"""Independent standard-library audit: raw rows, binomial-sum beta CDF, Wilson."""
import csv
from fractions import Fraction
import json
from math import comb,sqrt
from pathlib import Path
from statistics import NormalDist
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
with (ROOT/'data/unit02_help_seeking.csv').open(newline='',encoding='utf-8') as stream:
    rows=list(csv.reader(stream))
assert rows[0]==['participant_id','help_sought_within_7_days']
assert [row[0] for row in rows[1:]]==[f'H{i:03d}' for i in range(1,13)]
assert all(row[1] in ('0','1') for row in rows[1:])
n=len(rows)-1;y=sum(int(row[1]) for row in rows[1:])
assert (n,y)==(12,8)
spec=json.loads((ROOT/'model-specs/unit03_intervals_decision.json').read_text(encoding='utf-8'))
assert spec['spec_id']=='unit03-intervals-decision-v1'
assert spec['source']=='data/unit02_help_seeking.csv'
assert spec['source_model']=='model-specs/unit02_model.json'
assert spec['bayesian']['prior']==[2,2]
assert spec['bayesian']['sensitivity_priors']==[[1,1],[8,8]]
assert spec['frequentist']['procedure']=='Wilson score interval'
assert spec['frequentist']['confidence_level']==.95
assert spec['illustrative_decision']['threshold_theta']==.6
assert [spec['illustrative_decision'][key] for key in ('utility_if_above','utility_if_not_above','utility_no_action')]==[1,-1,0]
source=json.loads((ROOT/spec['source_model']).read_text(encoding='utf-8'))
assert source['spec_id']=='unit02-beta-binomial-v1' and source['primary_prior']['alpha']==2 and source['primary_prior']['beta']==2
record=json.loads((ROOT/'expected-results/unit03_intervals_decision.json').read_text(encoding='utf-8'))
assert record['spec_id']==spec['spec_id'] and (record['n'],record['y'])==(n,y)
assert record['computation']=='exact; no MCMC'

def beta_cdf(x,a,b):
    m=a+b-1
    return sum(comb(m,j)*x**j*(1-x)**(m-j) for j in range(a,m+1))
def quantile(q,a,b):
    low,high=0.,1.
    for _ in range(90):
        middle=(low+high)/2
        if beta_cdf(middle,a,b)<q:low=middle
        else:high=middle
    return (low+high)/2

observed=Fraction(y,n)
assert abs(record['observed_proportion']-float(observed))<1e-14
z=NormalDist().inv_cdf(.975)
p=float(observed)
denom=1+z*z/n
mid=(p+z*z/(2*n))/denom
half=z*sqrt(p*(1-p)/n+z*z/(4*n*n))/denom
for actual,expected in zip(record['wilson_95'],[mid-half,mid+half]):
    assert abs(actual-expected)<1e-12
assert record['posterior']==[10,6]
assert record['posterior_mean']==.625
for actual,expected in zip(record['posterior_equal_tail_95'],[quantile(.025,10,6),quantile(.975,10,6)]):
    assert abs(actual-expected)<1e-9
for item,prior in zip([record['primary_decision']]+record['sensitivity'],[[2,2],[1,1],[8,8]]):
    a,b=prior
    assert item['prior']==prior and item['posterior']==[a+y,b+n-y]
    exact=1-beta_cdf(Fraction(3,5),a+y,b+n-y)
    mass=float(exact)
    assert abs(item['probability_above_0_6']-mass)<1e-12
    assert abs(item['expected_action_utility']-(2*mass-1))<1e-12
    assert item['higher_expected_utility']==('action' if mass>.5 else 'no_action')
output=subprocess.check_output([sys.executable,str(ROOT/'scripts/analyze_unit03.py')],text=True)
for anchor in ('n=12 y=8 observed=0.667','wilson_95=[0.391,0.862]',
 'posterior=Beta(10,6) mean=0.625 credible_95=[0.384,0.837]',
 'Beta(2,2) P(theta>0.6)=0.597 toy_expected_utility=0.194 toy_choice=action',
 'Beta(8,8) P(theta>0.6)=0.387 toy_expected_utility=-0.225 toy_choice=no_action',
 'scope=fictional conditional model'):
    assert anchor in output,anchor
assert (ROOT/'figures/unit03_intervals.svg').is_file()
print('BAYESIAN_UNIT03_INDEPENDENT_OK n=12 y=8 wilson=yes posterior=yes sensitivity=yes')
