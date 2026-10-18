#!/usr/bin/env python3
"""Unit 5 primary NumPy/SciPy exact fixed-scale regression calculation."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.stats import norm, t
ROOT=Path(__file__).resolve().parents[1]

def calculate():
    spec=json.loads((ROOT/'model-specs/unit05_regression.json').read_text())
    if spec['spec_id']!='unit05-cross-sectional-regression-v1':raise ValueError('spec id')
    path=ROOT/spec['source']
    with path.open(newline='') as handle: rows=list(csv.DictReader(handle))
    if len(rows)!=72 or [r['participant_id'] for r in rows]!=[f'B{i:03d}' for i in range(1,73)]:raise ValueError('row ids')
    hours=np.array([float(r['study_hours_per_week']) for r in rows]);z=np.array([int(r['long_commute']) for r in rows]);y=np.array([float(r['belonging_score']) for r in rows])
    if not (np.isfinite(hours).all() and np.isfinite(y).all() and set(z)=={0,1} and
            np.all((hours>=0)&(hours<=168)) and np.all((y>=0)&(y<=100))):raise ValueError('invalid data')
    x=hours-5; X=np.column_stack([np.ones(len(y)),x,z]);U=X[:,:2]
    def ols(A):
        beta=np.linalg.solve(A.T@A,A.T@y)
        s=np.sqrt(np.sum((y-A@beta)**2)/(len(y)-A.shape[1]))
        cov=s*s*np.linalg.inv(A.T@A)
        return {'coefficients':beta.tolist(),'residual_sd':float(s),
                'slope_95_t':(beta[1]+t.ppf([.025,.975],len(y)-A.shape[1])*np.sqrt(cov[1,1])).tolist()}
    def posterior(beta_sd,sigma):
        prior_mean=np.array([55.,0.,0.]); prior_var=np.array([225.,beta_sd**2,100.])
        P0=np.diag(1/prior_var)
        covariance=np.linalg.inv(P0+X.T@X/sigma**2)
        mean=covariance@(P0@prior_mean+X.T@y/sigma**2)
        slope_sd=np.sqrt(covariance[1,1])
        result={'mean':mean.tolist(),'slope_sd':float(slope_sd),
                'slope_95':norm.ppf([.025,.975],loc=mean[1],scale=slope_sd).tolist(),
                'prob_slope_above_zero':float(norm.sf(0,loc=mean[1],scale=slope_sd)),
                'prob_slope_above_one':float(norm.sf(1,loc=mean[1],scale=slope_sd)),
                'two_hour_difference':float(2*mean[1])}
        residual=y-X@mean
        result['observed_abs_residual_gt_16']=int(np.sum(np.abs(residual)>16))
        predictive_sd=np.sqrt(sigma**2+np.einsum('ij,jk,ik->i',X,covariance,X))
        result['replicated_expected_abs_residual_gt_16']=float(np.sum(2*norm.sf(16,scale=predictive_sd)))
        return result
    prior_checks=[]
    for x0,z0 in [(0,0),(1.5,1)]:
        prior_sd=np.sqrt(225+4*x0*x0+100*z0*z0+64)
        prior_checks.append(float(norm.cdf(0,loc=55,scale=prior_sd)+norm.sf(100,loc=55,scale=prior_sd)))
    return {'spec_id':spec['spec_id'],'data_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'n':len(y),'long_commute_n':int(z.sum()),'ranges':{'hours':[float(hours.min()),float(hours.max())],'belonging':[float(y.min()),float(y.max())]},
            'ols_unadjusted':ols(U),'ols_adjusted':ols(X),
            'posterior_primary':posterior(2,8),'posterior_skeptical':posterior(.4,8),'posterior_scale_10':posterior(2,10),
            'prior_predictive_outside_0_100':prior_checks,
            'computation':'exact multivariate Normal conditional on fixed sigma; no MCMC',
            'claim_limit':'fictional cross-sectional association; no effect of increasing hours or population inference'}
if __name__=='__main__':
    v=calculate()
    print(f"n={v['n']} long_commute={v['long_commute_n']} hours=[{v['ranges']['hours'][0]:.1f},{v['ranges']['hours'][1]:.1f}] score=[{v['ranges']['belonging'][0]:.1f},{v['ranges']['belonging'][1]:.1f}]")
    print(f"OLS_slope_unadjusted={v['ols_unadjusted']['coefficients'][1]:.3f} 95_t=[{v['ols_unadjusted']['slope_95_t'][0]:.3f},{v['ols_unadjusted']['slope_95_t'][1]:.3f}]")
    print(f"OLS_slope_adjusted={v['ols_adjusted']['coefficients'][1]:.3f} 95_t=[{v['ols_adjusted']['slope_95_t'][0]:.3f},{v['ols_adjusted']['slope_95_t'][1]:.3f}] commute_coefficient={v['ols_adjusted']['coefficients'][2]:.3f}")
    for label in ('primary','skeptical','scale_10'):
        p=v['posterior_'+label]
        print(f"posterior_{label}_slope={p['mean'][1]:.3f} 95_credible=[{p['slope_95'][0]:.3f},{p['slope_95'][1]:.3f}] P(beta>1)={p['prob_slope_above_one']:.3f} two_hour_difference={p['two_hour_difference']:.3f}")
    p=v['posterior_primary'];print(f"prior_predictive_outside=[{v['prior_predictive_outside_0_100'][0]:.3f},{v['prior_predictive_outside_0_100'][1]:.3f}] observed_tail={p['observed_abs_residual_gt_16']} replicated_expected_tail={p['replicated_expected_abs_residual_gt_16']:.3f}")
    print('scope=fictional cross-sectional association, not an effect or population estimate; fixed-sigma working model; exact no MCMC')
