#!/usr/bin/env python3
"""Recompute the Primer anchors without using its generation script or SciPy."""
import csv
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
r=json.loads((ROOT/'expected-results/appendix_primer.json').read_text(encoding='utf-8'))


def need(ok, label):
    if not ok:
        raise SystemExit('PRIMER_INDEPENDENT_ERROR '+label)


def near(x,y,tol=1e-12):
    return abs(x-y)<=tol

s=r['screening']
p=F(8,100); sens=F(9,10); fpr=F(15,100)
flag_b=p*sens
flag_not=(1-p)*fpr
first=flag_b/(flag_b+flag_not)
second=first*sens/(first*sens+(1-first)*fpr)
for key,value in [('base_rate',p),('sensitivity',sens),('false_positive_rate',fpr),
                  ('flagged_benefit_share',flag_b),('flagged_nonbenefit_share',flag_not),
                  ('flagged_total_share',flag_b+flag_not),('benefit_given_first_flag',first),
                  ('benefit_given_two_flags',second)]:
    need(near(s[key],float(value)),key)
need(s['benefit_given_first_flag_fraction']==str(first),'first fraction')
need(s['benefit_given_two_flags_fraction']==str(second),'second fraction')
need('conditional independence' in s['second_flag_assumption'],'conditional independence')

c=r['count_model']
data=ROOT/'data/unit02_help_seeking.csv'
with data.open(newline='', encoding='utf-8') as f:
    rows=list(csv.DictReader(f))
need(len(rows)==12 and [row['participant_id'] for row in rows]==
     [f'H{i:03d}' for i in range(1,13)] and
     all(row['help_sought_within_7_days'] in ('0','1') for row in rows),
     'source observations and identifiers')
n=len(rows); y=sum(int(row['help_sought_within_7_days']) for row in rows)
listed=next((line.split()[0] for line in (ROOT/'SHA256SUMS').read_text(encoding='utf-8').splitlines()
             if line.endswith('  data/unit02_help_seeking.csv')),None)
need(listed==hashlib.sha256(data.read_bytes()).hexdigest(),'source CSV SHA256SUMS')
need(c['n']==n==12 and c['y']==y==8 and c['prior_shapes']==[2,2] and
     c['posterior_shapes']==[2+y,2+n-y], 'observed count and shapes')
with (ROOT/'expected-results/prelude_grid.csv').open(newline='', encoding='utf-8') as f:
    grid=list(csv.DictReader(f))
need(len(grid)==9,'grid inventory')
need(near(c['grid_mass_at_half'],float(grid[4]['posterior']),1e-6),'grid at half')
need(near(c['grid_tail_above_half'],sum(float(a['posterior']) for a in grid[5:]),2e-6),'grid tail rounded')
# For integer Beta shapes, the tail at 1/2 is a finite binomial sum.
continuous=F(sum(math.comb(15,k) for k in range(10)),2**15)
need(near(c['continuous_tail_above_half'],float(continuous)),'continuous tail')
for q,key in ((F(5,8),'beta_density_at_0_625'),(F(1,2),'beta_density_at_0_5')):
    value=F(math.factorial(15), math.factorial(9)*math.factorial(5))*q**9*(1-q)**5
    need(near(c[key],float(value)),'density '+key)
need(near(c['posterior_mean'],10/16),'posterior mean')
need(near(c['posterior_mode'],9/14),'posterior mode')
need(near(c['posterior_sd'],math.sqrt(60/(256*17))),'posterior sd')
need(near(c['mle'],8/12) and near(c['plug_in_proportion_se'],math.sqrt((8/12)*(4/12)/12)),'MLE and SE')

w=r['workshop_primary_tau5']
unit6=json.loads((ROOT/'expected-results/unit06_partial_pooling.json').read_text(encoding='utf-8'))
for key,ref in [('center',unit6['primary_tau5']['mu_mean']),
                ('w01_raw',unit6['raw_means']['W01']),('w01_posterior',unit6['primary_tau5']['theta_mean'][0]),
                ('w08_raw',unit6['raw_means']['W08']),('w08_posterior',unit6['primary_tau5']['theta_mean'][7])]:
    need(near(w[key],ref),key)
for f in ('appendix_primer_grid_density.svg','appendix_primer_sampling_posterior.svg','appendix_primer_shrinkage.svg'):
    need((ROOT/'figures'/f).is_file(),'figure '+f)
print('PRIMER_INDEPENDENT_OK screening=exact grid=distinct beta=exact workshops=existing')
