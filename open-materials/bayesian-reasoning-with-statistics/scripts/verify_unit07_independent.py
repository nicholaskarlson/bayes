#!/usr/bin/env python3
"""Standard-library audit using integer-Beta polynomial CDF and Simpson integration."""
import csv
from fractions import Fraction
import hashlib
import json
from math import comb,exp,lgamma
from pathlib import Path
import subprocess
import sys
import tempfile
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/unit07_help_seeking_routes.csv'
REFERENCE=ROOT/'expected-results/unit07_proportions.json'
ANCHORS=ROOT/'expected-results/unit07_anchors.csv'
SPEC=ROOT/'model-specs/unit07_proportions.json'
def need(ok,message):
    if not ok:raise AssertionError('UNIT07_INDEPENDENT_ERROR '+message)
def near(actual,expected,label,tol=2e-7):
    need(abs(actual-expected)<tol,f'{label}: {actual} vs {expected}')
def pdf(x,a,b):
    if not 0<x<1:return 0.
    return exp((a-1)*__import__('math').log(x)+(b-1)*__import__('math').log1p(-x)-lgamma(a)-lgamma(b)+lgamma(a+b))
def cdf(x,a,b):
    if x<=0:return 0.
    if x>=1:return 1.
    n=a+b-1
    return sum(comb(n,k)*x**k*(1-x)**(n-k) for k in range(a,n+1))
def simpson(func,lo,hi,N=4096):
    if hi<=lo:return 0.
    step=(hi-lo)/N
    return (func(lo)+func(hi)+4*sum(func(lo+i*step) for i in range(1,N,2))+2*sum(func(lo+i*step) for i in range(2,N,2)))*step/3
def tail(a,b,c,d,h):
    if h<=-1:return 1.
    if h>=1:return 0.
    lo=max(0.,h);hi=min(1.,1+h)
    main=simpson(lambda x:pdf(x,a,b)*cdf(x-h,c,d),lo,hi)
    above=1-cdf(1+h,a,b) if h<0 else 0.
    return main+above
def quantile(a,b,c,d,p):
    lo,hi=-1.,1.
    for _ in range(42):
        mid=(lo+hi)/2
        if 1-tail(a,b,c,d,mid)<p:lo=mid
        else:hi=mid
    return (lo+hi)/2
with DATA.open(newline='',encoding='utf-8') as f:rows=list(csv.DictReader(f))
need(len(rows)==36 and list(rows[0])==['volunteer_id','recruitment_route','reported_contact_7d'],'CSV schema')
need([(r['volunteer_id'],r['recruitment_route']) for r in rows]==[(f'D{i:03d}','desk') for i in range(1,25)]+[(f'P{i:03d}','poster') for i in range(1,13)],'identities and order')
need(all(r['reported_contact_7d'] in ('0','1') for r in rows),'binary values')
need((sum(int(r['reported_contact_7d']) for r in rows[:24]),sum(int(r['reported_contact_7d']) for r in rows[24:]))==(8,0),'counts')
with tempfile.TemporaryDirectory() as tmp:
    output=Path(tmp)/'generated.csv'
    subprocess.run([sys.executable,str(ROOT/'scripts/generate_unit07.py'),'--output',str(output)],check=True,capture_output=True)
    need(output.read_bytes()==DATA.read_bytes(),'generator byte identity')
ref=json.loads(REFERENCE.read_text(encoding='utf-8'));spec=json.loads(SPEC.read_text(encoding='utf-8'))
need(ref['spec_id']==spec['spec_id']=='unit07-help-seeking-routes-v1','model identity')
need(ref['data_sha256']==hashlib.sha256(DATA.read_bytes()).hexdigest(),'data digest')
need((ref['n'],ref['desk_n'],ref['desk_y'],ref['poster_n'],ref['poster_y'])==(36,24,8,12,0),'reference counts')
with ANCHORS.open(newline='', encoding='utf-8') as f:anchors=list(csv.DictReader(f))
need([v['case'] for v in anchors]==['primary','sensitivity'],'anchor cases')
for label,prior in [('primary',1),('sensitivity',2)]:
    a,b,c,d=prior+8,prior+16,prior,prior+12
    v=ref[label];r=next(x for x in anchors if x['case']==label)
    need(v['prior']==[prior,prior] and v['desk_posterior']==[a,b] and v['poster_posterior']==[c,d],'shapes')
    for key,value in [('a_desk',a),('b_desk',b),('a_poster',c),('b_poster',d)]:near(float(r[key]),value,key,1e-12)
    for key,value in [('desk_mean',Fraction(a,a+b)),('poster_mean',Fraction(c,c+d)),
                      ('difference_mean',Fraction(a,a+b)-Fraction(c,c+d))]:
        near(v[key],value,key,1e-12);near(float(r[key]),value,key+' anchor',1e-12)
    # At h=0, expand the integer-Beta CDF into a finite binomial sum,
    # then integrate each polynomial term exactly as a beta ratio.
    n=c+d-1
    p0=sum(comb(n,k)*exp(lgamma(a+k)+lgamma(b+n-k)-lgamma(a+b+n)-lgamma(a)-lgamma(b)+lgamma(a+b)) for k in range(c,n+1))
    near(v['P_difference_gt_zero'],p0,label+' exact finite sum',2e-12)
    near(float(r['P_diff_gt_zero']),p0,label+' p0 anchor',2e-12)
    p10=tail(a,b,c,d,.1)
    near(v['P_difference_gt_0_10'],p10,label+' Simpson tail')
    near(float(r['P_diff_gt_0_10']),p10,label+' p10 anchor')
    for p,x,field in [(0.025,v['difference_95'][0],'difference_lower'),(.975,v['difference_95'][1],'difference_upper')]:
        near(x,quantile(a,b,c,d,p),label+' difference quantile',4e-7)
        near(float(r[field]),x,label+' '+field,1e-12)
    for shape in [('desk',a,b),('poster',c,d)]:
        key,aa,bb=shape
        lo,hi=v[key+'_95']
        near(cdf(lo,aa,bb),.025,key+' lower CDF',1e-9)
        near(cdf(hi,aa,bb),.975,key+' upper CDF',1e-9)
    # Exact ratio for P(zero in 12) under integer Beta shapes.
    prior_zero=Fraction(1,1)
    replicated_zero=Fraction(1,1)
    for i in range(12):
        prior_zero*=Fraction(prior+i,2*prior+i)
        replicated_zero*=Fraction(d+i,c+d+i)
    for key,value,col in [('prior_predictive_poster_zero',prior_zero,'prior_zero'),('posterior_predictive_poster_zero',replicated_zero,'replicated_zero')]:
        near(v[key],value,key,1e-12);near(float(r[col]),value,col,1e-12)
need(ref['scope']=='fictional convenience volunteers; route association only; no causal effect or campus rate; exact no MCMC','scope')
need((ROOT/'figures/unit07_proportions.svg').is_file(),'figure')
output=subprocess.check_output([sys.executable,str(ROOT/'scripts/analyze_unit07.py')],text=True)
for term in ('n=36 desk=8/24 poster=0/12','difference=0.275 difference_95=[0.042,0.490]',
             'difference=0.232 difference_95=[-0.017,0.456]','P_diff_gt_zero=0.987',
             'prior_zero=0.077 replicated_zero=0.520','exact no MCMC'):
    need(term in output,'Python output '+term)
print('BAYESIAN_UNIT07_INDEPENDENT_OK n=36 sparse_zero=yes exact_difference=yes')
