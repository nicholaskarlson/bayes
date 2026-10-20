#!/usr/bin/env python3
"""Raw-decimal, analytic Normal and independent composite-Simpson t audit."""
import csv
from decimal import Decimal
from hashlib import sha256
import json
import math
from pathlib import Path
import statistics
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/unit04_reflection_attention_paired.csv'
PLAN=ROOT/'model-specs/unit04_likelihood_plan.json'
PARAMS=ROOT/'model-specs/unit04_likelihood_parameters.csv'
BASELINE=ROOT/'expected-results/unit04_exact.json'
REF=ROOT/'expected-results/unit04_likelihood_reference.json'
ANCHORS=ROOT/'expected-results/unit04_likelihood_anchors.csv'


def check(ok,why):
    if not ok:raise AssertionError('UNIT04_LIKELIHOOD_ERROR '+why)


def near(a,b,why,tol):
    check(math.isfinite(a) and math.isfinite(b) and abs(a-b)<=tol,why+' '+repr((a,b)))


plan=json.loads(PLAN.read_text(encoding='utf-8'))
reference=json.loads(REF.read_text(encoding='utf-8'))
baseline=json.loads(BASELINE.read_text(encoding='utf-8'))
with PARAMS.open(newline='',encoding='utf-8') as h:params={r['parameter']:float(r['value']) for r in csv.DictReader(h)}
expected_params={'mu_prior_mean_points':float(plan['prior']['mu_mean_points']),
 'mu_prior_sd_points':float(plan['prior']['mu_sd_points']),
 'normal_sd_points':float(plan['matched_likelihoods'][0]['known_scale_points']),
 'student_t_scale_points':float(plan['matched_likelihoods'][1]['known_scale_points']),
 'student_degrees_of_freedom':float(plan['matched_likelihoods'][1]['degrees_of_freedom']),
 'two_point_threshold':float(plan['two_point_threshold']),
 'integration_lower':float(plan['integration_interval_points'][0]),
 'integration_upper':float(plan['integration_interval_points'][1]),
 'enlarged_lower':float(plan['enlarged_interval_points'][0]),
 'enlarged_upper':float(plan['enlarged_interval_points'][1])}
check(params==expected_params,'R/Python plan mismatch')
check(plan['matched_likelihoods'][0]['known_scale_points']==4 and
      abs(plan['matched_likelihoods'][1]['known_scale_points']**2*2-16)<1e-12,
      'not a matched-variance comparison')
for key,path in (('csv_sha256',DATA),('plan_sha256',PLAN),('baseline_exact_sha256',BASELINE)):
    check(reference['source'][key]==sha256(path.read_bytes()).hexdigest(),'source '+key)
check(reference['source']['source_base_commit']=='9edadbe5385a8012658accc9b0c00f74a4debe30','base SHA')
check(reference['source']['artifact_commit'] is None,'invented artifact SHA')
check(reference['method']=='normal_exact_and_student_t_adaptive_quadrature_not_mcmc','method label')
near(reference['baseline_normal_inverse_gamma']['p_mu_gt_2'],baseline['bayesian_exact']['primary']['p_mu_gt_two'],
     'baseline model',1e-14)
with DATA.open(newline='',encoding='utf-8') as h:rows=list(csv.DictReader(h))
check(len(rows)==64 and [r['participant_id'] for r in rows]==[f'P{i:03d}' for i in range(1,65)],'rows')
observations=[(r['participant_id'],float(Decimal(r['after_score'])-Decimal(r['before_score']))) for r in rows]
check(plan['cases']=={'all':None,'omit_P059':'P059','omit_P042':'P042'},'case plan')
with ANCHORS.open(newline='',encoding='utf-8') as h:
    anchor_rows=list(csv.DictReader(h))
check(len(anchor_rows)==24,'anchor count')
anchors={(r['case'],r['model'],r['field']):float(r['value']) for r in anchor_rows}
check(len(anchors)==24,'duplicate anchors')


def simpson_t(values):
    # Independent pure-Python log-likelihood; composite Simpson on 16,000 subintervals.
    nu=4.;s=math.sqrt(8);tau=5.;lo=-20.;hi=20.;intervals=16000;h=(hi-lo)/intervals
    const=math.lgamma((nu+1)/2)-math.lgamma(nu/2)-.5*math.log(nu*math.pi)-math.log(s)
    def log_kernel(mu):
        return -.5*(mu/tau)**2 + math.fsum(
            const-.5*(nu+1)*math.log1p(((v-mu)/s)**2/nu) for v in values)
    offset=log_kernel(2.2)
    weights=[math.exp(log_kernel(lo+i*h)-offset) for i in range(intervals+1)]
    def weighted(integrand):
        return h/3*(integrand(0)*weights[0]+integrand(intervals)*weights[-1]+4*math.fsum(integrand(i)*weights[i] for i in range(1,intervals,2))+
                    2*math.fsum(integrand(i)*weights[i] for i in range(2,intervals,2)))
    z=weighted(lambda i:1.)
    mean=weighted(lambda i:lo+i*h)/z
    cum=0.;cdf=[0.]
    for i in range(0,intervals,2):
        cum+=h/3*(weights[i]+4*weights[i+1]+weights[i+2])
        cdf.append(cum/z)
    def quantile(p):
        i=next(i for i in range(1,len(cdf)) if cdf[i]>=p)
        return lo+2*h*(i-1+(p-cdf[i-1])/(cdf[i]-cdf[i-1]))
    threshold_index=int((2-lo)/(2*h))
    return {'posterior_mean':mean,'ci95':[quantile(.025),quantile(.975)],
            'p_mu_gt_2':1-cdf[threshold_index]}


normal=statistics.NormalDist()
for case,removed in plan['cases'].items():
    vals=[v for pid,v in observations if pid!=removed]
    check(reference['cases'][case]['n']==len(vals),'case n')
    check(reference['cases'][case]['omitted_participant']==removed,'case omitted id')
    prior_precision=1/25;postvar=1/(prior_precision+len(vals)/16)
    postmean=postvar*(math.fsum(vals)/16)
    sd=math.sqrt(postvar)
    expected_normal={'posterior_mean':postmean,
      'ci95':[postmean+sd*normal.inv_cdf(.025),postmean+sd*normal.inv_cdf(.975)],
      'p_mu_gt_2':normal.cdf((postmean-2)/sd)}
    expected_t=simpson_t(vals)
    for label,computed,tolerance in (('normal_fixed',expected_normal,2e-10),('student_t4_equal_variance',expected_t,3e-4)):
        target=reference['cases'][case][label]
        for key in ('posterior_mean','p_mu_gt_2'):
            near(target[key],computed[key],case+' '+label+' '+key,tolerance)
            near(anchors[(case,label,key)],target[key],case+' anchor '+key,2e-15)
        for idx,key in enumerate(('ci95_lower','ci95_upper')):
            near(target['ci95'][idx],computed['ci95'][idx],case+' '+label+' '+key,tolerance)
            near(anchors[(case,label,key)],target['ci95'][idx],case+' anchor '+key,2e-15)
        if label=='student_t4_equal_variance':
            quad=target['quadrature']
            check(quad['bounds']==[-20,20] and quad['enlarged_bounds']==[-30,30],'quadrature bounds')
            check(quad['mass_fraction_outside_bounds'] < 1e-8,'truncation')
            near(quad['p_threshold_enlarged'],target['p_mu_gt_2'],'enlarged result',1e-8)

full=reference['cases']['all']
near(full['student_t4_equal_variance']['p_mu_gt_2']-full['normal_fixed']['p_mu_gt_2'],.1138514950677598,
     'matched likelihood difference',2e-6)
for label in ('normal_fixed','student_t4_equal_variance'):
    check(reference['cases']['omit_P059'][label]['p_mu_gt_2']<full[label]['p_mu_gt_2']<
          reference['cases']['omit_P042'][label]['p_mu_gt_2'],'omission direction '+label)
note=(ROOT/'technical-notes/unit04_deep_dive.md').read_text(encoding='utf-8')
for phrase in ('conditional variance','post hoc','0.636','0.750','0.553','figures/unit04_likelihood.svg',
               'not a Bayes factor','not MCMC'):
    check(phrase in note,'technical audit safeguard '+phrase)
with tempfile.TemporaryDirectory() as tmp:
    made=Path(tmp)/'figure.svg'
    subprocess.check_call([sys.executable,str(ROOT/'scripts/make_unit04_likelihood_figure.py'),str(made)])
    figure=ROOT/'figures/unit04_likelihood.svg'
    check(made.read_bytes()==figure.read_bytes(),'source figure drift')
    check('aria-labelledby="title desc"' in figure.read_text(encoding='utf-8'),'accessible SVG')
print('BAYESIAN_UNIT04_LIKELIHOOD_INDEPENDENT_OK cases=3 models=2 adaptive_vs_Simpson=yes normal_exact=yes')
