#!/usr/bin/env python3
"""Standard-library matrix/Normal audit from Unit 5 CSV, no primary imports."""
import csv
import hashlib
import json
from math import gamma,pi,sqrt
from pathlib import Path
from statistics import NormalDist
import subprocess
import sys
import tempfile
ROOT=Path(__file__).resolve().parents[1]
CSV=ROOT/'data/unit05_belonging_study.csv'
SPEC=ROOT/'model-specs/unit05_regression.json'
EXPECTED=ROOT/'expected-results/unit05_regression.json'
def near(a,b,tol=1e-8):
    assert abs(a-b)<tol,(a,b)
def inverse(matrix):
    n=len(matrix);a=[[float(v) for v in row]+[float(i==j) for j in range(n)] for i,row in enumerate(matrix)]
    for j in range(n):
        p=max(range(j,n),key=lambda i:abs(a[i][j]))
        a[j],a[p]=a[p],a[j]
        q=a[j][j];assert abs(q)>1e-12
        a[j]=[v/q for v in a[j]]
        for i in range(n):
            if i!=j:
                q=a[i][j];a[i]=[v-q*w for v,w in zip(a[i],a[j])]
    return [row[n:] for row in a]
def mv(a,v):return [sum(x*y for x,y in zip(row,v)) for row in a]
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def cross(X):return [[sum(row[i]*row[j] for row in X) for j in range(len(X[0]))] for i in range(len(X[0]))]
def tcdf(t,df):
    c=gamma((df+1)/2)/(sqrt(df*pi)*gamma(df/2));n=2048;h=t/n
    return .5+h/3*sum((1 if j in (0,n) else 4 if j%2 else 2)*c*(1+(j*h)**2/df)**(-(df+1)/2) for j in range(n+1))
def tcrit(df):
    lo,hi=0.,5.
    for _ in range(38):
        m=(lo+hi)/2
        if tcdf(m,df)<.975:lo=m
        else:hi=m
    return (lo+hi)/2
with CSV.open(newline='', encoding='utf-8') as inp:
    records=list(csv.DictReader(inp))
assert len(records)==72 and list(records[0])==['participant_id','study_hours_per_week','long_commute','belonging_score']
assert [r['participant_id'] for r in records]==[f'B{i:03d}' for i in range(1,73)]
h=[float(r['study_hours_per_week']) for r in records];z=[int(r['long_commute']) for r in records];y=[float(r['belonging_score']) for r in records]
assert set(z)=={0,1} and sum(z)==36 and all(0<=a<=100 for a in y)
X=[[1.,a-5,float(b)] for a,b in zip(h,z)]
sp=json.loads(SPEC.read_text(encoding='utf-8'));ref=json.loads(EXPECTED.read_text(encoding='utf-8'))
assert sp['spec_id']==ref['spec_id']=='unit05-cross-sectional-regression-v1'
assert sp['primary_prior']['beta']==[0,2] and sp['bayesian_likelihood'].find('sigma=8')!=-1
assert ref['data_sha256']==hashlib.sha256(CSV.read_bytes()).hexdigest()
assert ref['n']==72 and ref['long_commute_n']==36
for a,b in zip(ref['ranges']['hours'],[min(h),max(h)]):near(a,b)
for a,b in zip(ref['ranges']['belonging'],[min(y),max(y)]):near(a,b)
with tempfile.TemporaryDirectory() as tmp:
    output=[]
    for j in range(2):
        path=Path(tmp)/f'copy{j}.csv'
        subprocess.run([sys.executable,str(ROOT/'scripts/generate_unit05.py'),'--output',str(path)],check=True,capture_output=True)
        output.append(path.read_bytes())
    assert output[0]==output[1]==CSV.read_bytes()
for name,A in [('ols_unadjusted',[row[:2] for row in X]),('ols_adjusted',X)]:
    C=inverse(cross(A));b=mv(C,[sum(row[j]*v for row,v in zip(A,y)) for j in range(len(A[0]))])
    s=sqrt(sum((v-dot(row,b))**2 for row,v in zip(A,y))/(len(y)-len(A[0])))
    se=s*sqrt(C[1][1]);q=tcrit(len(y)-len(A[0]))
    for a,c in zip(ref[name]['coefficients'],b):near(a,c)
    near(ref[name]['residual_sd'],s)
    for a,c in zip(ref[name]['slope_95_t'],[b[1]-q*se,b[1]+q*se]):near(a,c,1e-7)
N=NormalDist()
for name,sd,sigma in [('posterior_primary',2.,8.),('posterior_skeptical',.4,8.),('posterior_scale_10',2.,10.)]:
    crossprod=cross(X);P=[[crossprod[i][j]/sigma**2+(1/[225,sd*sd,100][i] if i==j else 0.) for j in range(3)] for i in range(3)]
    V=inverse(P)
    score=[sum(row[j]*v for row,v in zip(X,y))/sigma**2+[55/225,0,0][j] for j in range(3)]
    b=mv(V,score);unc=sqrt(V[1][1]);result=ref[name]
    for a,c in zip(result['mean'],b):near(a,c)
    near(result['slope_sd'],unc)
    for a,c in zip(result['slope_95'],[b[1]+N.inv_cdf(.025)*unc,b[1]+N.inv_cdf(.975)*unc]):near(a,c)
    near(result['prob_slope_above_zero'],1-N.cdf(-b[1]/unc))
    near(result['prob_slope_above_one'],1-N.cdf((1-b[1])/unc))
    near(result['two_hour_difference'],2*b[1])
    tails=sum(abs(v-dot(row,b))>16 for row,v in zip(X,y))
    expected=sum(2*N.cdf(-16/sqrt(sigma*sigma+dot(row,mv(V,row)))) for row in X)
    assert result['observed_abs_residual_gt_16']==tails
    near(result['replicated_expected_abs_residual_gt_16'],expected)
for recorded,(x0,z0) in zip(ref['prior_predictive_outside_0_100'],[(0,0),(1.5,1)]):
    scale=sqrt(225+4*x0*x0+100*z0*z0+64)
    near(recorded,N.cdf(-55/scale)+1-N.cdf(45/scale))
assert ref['computation']=='exact multivariate Normal conditional on fixed sigma; no MCMC'
assert (ROOT/'figures/unit05_slope_intervals.svg').exists()
output=subprocess.check_output([sys.executable,str(ROOT/'scripts/analyze_unit05.py')],text=True)
for term in ('n=72 long_commute=36','OLS_slope_unadjusted=3.274','OLS_slope_adjusted=2.766',
             'posterior_primary_slope=2.554','P(beta>1)=0.996','posterior_skeptical_slope=0.852',
             'P(beta>1)=0.329','prior_predictive_outside=[0.005,0.015]',
             'observed_tail=3 replicated_expected_tail=3.589','exact no MCMC'):
    assert term in output,term
print('BAYESIAN_UNIT05_INDEPENDENT_OK n=72 models=OLS2+Normal3')
