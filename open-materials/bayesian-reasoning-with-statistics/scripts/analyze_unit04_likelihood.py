#!/usr/bin/env python3
"""Matched fixed-scale normal/Student-t likelihood sensitivity; adaptive quadrature."""
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar
from scipy.special import gammaln
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data/unit04_reflection_attention_paired.csv'
PLAN = ROOT / 'model-specs/unit04_likelihood_plan.json'
BASELINE = ROOT / 'expected-results/unit04_exact.json'


def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def read():
    with DATA.open(newline='', encoding='utf-8') as handle: rows = list(csv.DictReader(handle))
    assert len(rows) == 64 and [r['participant_id'] for r in rows] == [f'P{i:03d}' for i in range(1, 65)]
    return [(r['participant_id'], float(r['after_score']) - float(r['before_score'])) for r in rows]


def normal_exact(values, prior_mean, prior_sd, scale, threshold):
    precision = 1/prior_sd**2 + len(values)/scale**2
    v = 1/precision
    mean = v*(prior_mean/prior_sd**2 + math.fsum(values)/scale**2)
    sd = math.sqrt(v)
    return dict(posterior_mean=mean, ci95=[float(x) for x in norm.ppf((.025, .975), loc=mean, scale=sd)],
                p_mu_gt_2=float(norm.sf(threshold, loc=mean, scale=sd)))


def robust_quadrature(values, prior_mean, prior_sd, scale, df, threshold, bounds, larger):
    prior_const = -.5*math.log(2*math.pi*prior_sd**2)
    student_const = gammaln((df+1)/2) - gammaln(df/2) - .5*math.log(df*math.pi) - math.log(scale)
    exponent = (df+1)/2
    def log_density(mu):
        return prior_const - .5*((mu-prior_mean)/prior_sd)**2 + math.fsum(
            student_const - exponent*math.log1p(((value-mu)/scale)**2/df) for value in values)
    mode = minimize_scalar(lambda x: -log_density(x), bounds=bounds, method='bounded',
                           options={'xatol': 1e-11})
    assert mode.success and bounds[0] < mode.x < bounds[1]
    peak = log_density(mode.x)
    def weight(mu): return math.exp(log_density(mu)-peak)
    def integral(fun, low, high):
        val, err = quad(fun, low, high, epsabs=1e-12, epsrel=1e-11, limit=500,
                        points=[mode.x] if low < mode.x < high else None)
        assert math.isfinite(val) and err < 1e-8
        return val
    z = integral(weight, *bounds)
    assert z > 0
    postmean = integral(lambda mu: mu*weight(mu), *bounds)/z
    probability = integral(weight, threshold, bounds[1])/z
    def cdf(x): return integral(weight, bounds[0], x)/z
    interval = [brentq(lambda x: cdf(x)-p, bounds[0], bounds[1], xtol=1e-11)
                for p in (.025,.975)]
    z_large = integral(weight, *larger)
    p_large = integral(weight, threshold, larger[1])/z_large
    outside_fraction = (z_large-z)/z_large
    assert outside_fraction >= -1e-10 and outside_fraction < 1e-8
    assert abs(p_large-probability) < 1e-8
    return dict(posterior_mean=postmean, ci95=interval, p_mu_gt_2=probability,
                quadrature=dict(bounds=bounds, enlarged_bounds=larger,
                                mass_fraction_outside_bounds=max(0.,outside_fraction),
                                p_threshold_enlarged=p_large, mode=mode.x))


def analyze():
    spec=json.loads(PLAN.read_text(encoding='utf-8'))
    assert spec['schema_version']=='unit04-likelihood-sensitivity-1'
    assert spec['data_rows']==64 and spec['two_point_threshold']==2
    assert [p['id'] for p in spec['matched_likelihoods']]==['normal_fixed','student_t4_equal_variance']
    normal_spec,student_spec=spec['matched_likelihoods']
    assert normal_spec['known_scale_points']==4
    assert student_spec['degrees_of_freedom']==4
    assert math.isclose(student_spec['known_scale_points']**2*4/(4-2),
                        normal_spec['known_scale_points']**2,rel_tol=0,abs_tol=1e-12)
    records=read(); prior=spec['prior'];m=prior['mu_mean_points'];tau=prior['mu_sd_points']
    threshold=spec['two_point_threshold'];scale=normal_spec['known_scale_points'];t_scale=student_spec['known_scale_points'];df=student_spec['degrees_of_freedom']
    all_results={}
    for name,removed in spec['cases'].items():
        values=[v for pid,v in records if pid!=removed]
        assert len(values)==(64 if removed is None else 63)
        all_results[name]={'n':len(values),'omitted_participant':removed,
            'normal_fixed':normal_exact(values,m,tau,scale,threshold),
            'student_t4_equal_variance':robust_quadrature(values,m,tau,t_scale,df,threshold,
                spec['integration_interval_points'],spec['enlarged_interval_points'])}
    baseline=json.loads(BASELINE.read_text(encoding='utf-8'))
    return {'schema_version':'unit04-likelihood-reference-1',
            'source':{'csv_sha256':digest(DATA),'plan_sha256':digest(PLAN),'baseline_exact_sha256':digest(BASELINE),
                      'source_base_commit':'9edadbe5385a8012658accc9b0c00f74a4debe30','artifact_commit':None},
            'method':'normal_exact_and_student_t_adaptive_quadrature_not_mcmc',
            'model_note':'Both sensitivity models use mu ~ Normal(0,5^2) and conditional variance 16: Normal SD 4; Student-t(4) scale sqrt(8). Baseline NIG estimates unknown variance and has a different prior, so contrasting it with either side model does not isolate the likelihood.',
            'cases':all_results,
            'baseline_normal_inverse_gamma':{'p_mu_gt_2':baseline['bayesian_exact']['primary']['p_mu_gt_two']},
            'interpretation':'Synthetic cohort, post hoc model sensitivity. Both fixed dispersion parameters are illustrative, not known in a real study. These are conditional comparisons, not model probabilities, Bayes factors, predictive validation, or causal effects.'}


if __name__ == '__main__':
    output=json.dumps(analyze(),indent=2,sort_keys=True,allow_nan=False)+'\n'
    if len(sys.argv)==1: print(output,end='')
    elif len(sys.argv)==2: Path(sys.argv[1]).write_text(output,encoding='utf-8')
    else: raise SystemExit('usage: analyze_unit04_likelihood.py [OUTPUT.json]')
