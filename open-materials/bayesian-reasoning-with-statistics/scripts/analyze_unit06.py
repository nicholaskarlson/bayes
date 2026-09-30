#!/usr/bin/env python3
"""Unit 6 exact joint Normal hierarchy, with fixed within/between group scales."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[1]
SIZES = (4, 5, 6, 8, 10, 12, 15, 20)
GROUPS = tuple(f'W{j:02d}' for j in range(1, 9))

def calculate():
    spec = json.loads((ROOT/'model-specs/unit06_partial_pooling.json').read_text())
    if spec['spec_id'] != 'unit06-workshop-partial-pooling-v1':
        raise ValueError('Wrong Unit 6 specification')
    csv_path = ROOT/spec['source']
    with csv_path.open(newline='', encoding='utf-8') as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != sum(SIZES) or list(rows[0]) != ['student_id', 'workshop_id', 'confidence_score']:
        raise ValueError('Wrong Unit 6 layout')
    scores = []
    for group, count in zip(GROUPS, SIZES):
        group_rows = [r for r in rows if r['workshop_id'] == group]
        if [r['student_id'] for r in group_rows] != [f'{group}-S{k:02d}' for k in range(1, count+1)]:
            raise ValueError('Missing, duplicate, or out-of-order student')
        scores.append(np.array([float(r['confidence_score']) for r in group_rows]))
    if not all(np.isfinite(y).all() and ((0 <= y) & (y <= 100)).all() for y in scores):
        raise ValueError('Invalid score')
    n = np.array(SIZES)
    raw = np.array([np.mean(y) for y in scores]); pooled = float(np.mean(np.concatenate(scores)))
    within_sse = float(sum(np.sum((y-y.mean())**2) for y in scores))
    sigma = 8.
    def posterior(tau):
        precision = np.zeros((9, 9)); rhs = np.zeros(9)
        precision[0, 0] = 1/100 + 8/tau**2
        precision[0, 1:] = precision[1:, 0] = -1/tau**2
        precision[np.arange(1, 9), np.arange(1, 9)] = 1/tau**2+n/sigma**2
        rhs[0] = 55/100
        rhs[1:] = n*raw/sigma**2
        covariance = np.linalg.inv(precision)
        mean = np.linalg.solve(precision, rhs)
        sds = np.sqrt(np.diag(covariance)[1:])
        intervals = norm.ppf([.025, .975], loc=mean[1:, None], scale=sds[:, None])
        difference = mean[1]-mean[8]
        difference_sd = np.sqrt(covariance[1, 1]+covariance[8, 8]-2*covariance[1, 8])
        return {'mu_mean':float(mean[0]), 'mu_sd':float(np.sqrt(covariance[0, 0])),
                'theta_mean':mean[1:].tolist(), 'theta_sd':sds.tolist(),
                'theta_95':intervals.tolist(),
                'weights':(n*tau**2/(sigma**2+n*tau**2)).tolist(),
                'A_minus_H_mean':float(difference), 'A_minus_H_sd':float(difference_sd),
                'P_A_gt_H':float(norm.sf(0, loc=difference, scale=difference_sd))}
    prior_sd = np.sqrt(100+5**2+sigma**2/SIZES[0])
    single_sd = np.sqrt(100+5**2+sigma**2)
    return {'spec_id':spec['spec_id'], 'data_sha256':hashlib.sha256(csv_path.read_bytes()).hexdigest(),
            'n':len(rows), 'sizes':dict(zip(GROUPS, SIZES)), 'pooled_mean':pooled,
            'raw_means':dict(zip(GROUPS, raw.tolist())), 'score_range':[float(min(y.min() for y in scores)),float(max(y.max() for y in scores))],
            'within_sse':within_sse, 'within_df':len(rows)-len(GROUPS),
            'within_sd':float(np.sqrt(within_sse/(len(rows)-len(GROUPS)))),
            'replicated_expected_within_sse':float((len(rows)-len(GROUPS))*sigma**2),
            'prior_predictive_W01_group_mean_95':norm.ppf([.025,.975],loc=55,scale=prior_sd).tolist(),
            'prior_predictive_score_outside_0_100':float(norm.cdf(0,loc=55,scale=single_sd)+norm.sf(100,loc=55,scale=single_sd)),
            'conditional_prior_icc':float(5**2/(5**2+sigma**2)),
            'primary_tau5':posterior(5.), 'sensitivity_tau2':posterior(2.), 'sensitivity_tau10':posterior(10.),
            'computation':'exact joint Normal integration with fixed sigma and tau; no MCMC',
            'claim_limit':'fictional group association only; no workshop effects, individual change or population inference'}

if __name__ == '__main__':
    v = calculate()
    print(f"n={v['n']} groups=8 sizes={','.join(str(k) for k in v['sizes'].values())} pooled={v['pooled_mean']:.3f} score_range=[{v['score_range'][0]:.1f},{v['score_range'][1]:.1f}]")
    p = v['primary_tau5']
    print(f"posterior_mu={p['mu_mean']:.3f} mu_95=[{p['mu_mean']-1.959963984540054*p['mu_sd']:.3f},{p['mu_mean']+1.959963984540054*p['mu_sd']:.3f}]")
    for j, (g, size) in enumerate(v['sizes'].items()):
        print(f"group_{g} n={size} raw={v['raw_means'][g]:.3f} posterior={p['theta_mean'][j]:.3f} 95=[{p['theta_95'][j][0]:.3f},{p['theta_95'][j][1]:.3f}] weight={p['weights'][j]:.3f}")
    for tau, key in ((5,'primary_tau5'),(2,'sensitivity_tau2'),(10,'sensitivity_tau10')):
        s=v[key]
        print(f"tau={tau} P_A_gt_H={s['P_A_gt_H']:.3f} A={s['theta_mean'][0]:.3f} H={s['theta_mean'][7]:.3f}")
    print(f"prior_W01_mean_95=[{v['prior_predictive_W01_group_mean_95'][0]:.3f},{v['prior_predictive_W01_group_mean_95'][1]:.3f}] outside_score={v['prior_predictive_score_outside_0_100']:.4f} prior_ICC={v['conditional_prior_icc']:.3f}")
    print(f"within_sd={v['within_sd']:.3f} observed_SSE={v['within_sse']:.3f} replicated_expected_SSE={v['replicated_expected_within_sse']:.3f} df={v['within_df']}")
    print('scope=fictional one-time grouped scores; no workshop effects or population claim; fixed scales; exact no MCMC')
