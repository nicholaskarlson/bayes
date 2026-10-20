#!/usr/bin/env python3
"""Integer/Fraction audit independent of scipy analyzer and its imports."""
import csv
from fractions import Fraction
import json
from math import comb, factorial
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
with (ROOT/'data/unit02_help_seeking.csv').open(newline='', encoding='utf-8') as stream:
    rows=list(csv.reader(stream))
assert rows[0]==['participant_id','help_sought_within_7_days']
assert len(rows)==13 and [r[0] for r in rows[1:]]==[f'H{i:03d}' for i in range(1,13)]
assert all(r[1] in ('0','1') for r in rows[1:])
y=sum(int(r[1]) for r in rows[1:]); n=12
assert y==8
spec=json.loads((ROOT/'model-specs/unit02_model.json').read_text(encoding='utf-8'))
assert spec['spec_id']=='unit02-beta-binomial-v1'
assert spec['primary_prior']=={'distribution':'Beta','alpha':2,'beta':2}
assert spec['sensitivity_priors']==[{'alpha':1,'beta':1},{'alpha':8,'beta':8}]

def B(a,b):
    return Fraction(factorial(a-1)*factorial(b-1),factorial(a+b-1))
def pmf(k,n,a,b):
    return comb(n,k)*B(a+k,b+n-k)/B(a,b)
def cdf(x,a,b):
    m=a+b-1
    return sum(comb(m,j)*x**j*(1-x)**(m-j) for j in range(a,m+1))
def quantile(q,a,b):
    lo,hi=0.,1.
    for _ in range(90):
        mid=(lo+hi)/2
        if cdf(mid,a,b)<q: lo=mid
        else: hi=mid
    return (lo+hi)/2
reference=json.loads((ROOT/'expected-results/unit02_exact.json').read_text(encoding='utf-8'))
assert reference['spec_id']==spec['spec_id'] and reference['n']==n and reference['y']==y
assert reference['prior']==[2,2] and reference['posterior']==[10,6]
assert reference['computation']=='exact; no MCMC'
assert reference['sample_proportion']==y/n
checks={
 'prior_predictive_ge_10_of_12':sum(pmf(k,n,2,2) for k in range(10,13)),
 'prior_predictive_outer_0_to_2_or_10_to_12':sum(pmf(k,n,2,2) for k in list(range(3))+list(range(10,13))),
 'posterior_mean':Fraction(10,16),
 'posterior_probability_above_half':1-cdf(.5,10,6),
 'posterior_predictive_ge_8_of_12':sum(pmf(k,n,10,6) for k in range(8,13)),
}
for name,number in checks.items():
    assert abs(reference[name]-float(number))<1e-11,(name,number)
for printed,calculated in zip(reference['posterior_equal_tail_95'],(quantile(.025,10,6),quantile(.975,10,6))):
    assert abs(printed-calculated)<1e-9
for item,(a,b) in zip(reference['sensitivity'],[(1,1),(8,8)]):
    assert item['prior']==[a,b] and item['posterior']==[a+8,b+4]
    assert abs(item['mean']-Fraction(a+8,a+b+12))<1e-11
    assert abs(item['probability_above_half']-(1-cdf(.5,a+8,b+4)))<1e-11
output=subprocess.check_output([sys.executable,str(ROOT/'scripts/analyze_unit02.py')],text=True)
for anchor in ('n=12 y=8 observed=0.667','posterior=Beta(10,6)',
               'prior_P(Y>=10_of_12)=0.154','posterior_mean=0.625',
               'P(theta>0.5)=0.849','posterior_predictive_P(Yrep>=8_of_12)=0.519',
               'sensitivity_Beta(8,8)_mean=0.571 P(theta>0.5)=0.779'):
    assert anchor in output,anchor
assert (ROOT/'figures/unit02_update.svg').is_file()
print('BAYESIAN_UNIT02_INDEPENDENT_OK n=12 y=8 posterior=Beta(10,6)')
