#!/usr/bin/env python3
"""Exact independent-Beta posteriors and deterministic difference integrals."""
import csv
import hashlib
from pathlib import Path
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import betaln
from scipy.stats import beta
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/unit07_help_seeking_routes.csv'
def read_counts():
    with DATA.open(newline='',encoding='utf-8') as f:rows=list(csv.DictReader(f))
    assert list(rows[0])==['volunteer_id','recruitment_route','reported_contact_7d']
    assert [(r['volunteer_id'],r['recruitment_route']) for r in rows]==[(f'D{i:03d}','desk') for i in range(1,25)]+[(f'P{i:03d}','poster') for i in range(1,13)]
    assert all(r['reported_contact_7d'] in ('0','1') for r in rows)
    return len(rows),sum(int(r['reported_contact_7d']) for r in rows[:24]),sum(int(r['reported_contact_7d']) for r in rows[24:])
def tail(a,b,c,d,cut):
    # P(theta_desk - theta_poster > cut); full support [-1,1].
    if cut<=-1:return 1.
    if cut>=1:return 0.
    lo=max(0.,cut);hi=min(1.,cut+1.)
    interior=quad(lambda x:beta.pdf(x,a,b)*beta.cdf(x-cut,c,d),lo,hi,epsabs=2e-12,points=[lo,hi])[0] if hi>lo else 0.
    above=beta.sf(1.+cut,a,b) if cut<0 else 0.
    return interior+above
def case(prior):
    a0,b0=prior; nD,yD=24,8;nP,yP=12,0
    a,b=a0+yD,b0+nD-yD;c,d=a0+yP,b0+nP-yP
    p0=tail(a,b,c,d,0.)
    result={'prior':[a0,b0],'desk_posterior':[a,b],'poster_posterior':[c,d],
            'desk_mean':a/(a+b),'poster_mean':c/(c+d),
            'desk_95':list(beta.ppf([.025,.975],a,b)),
            'poster_95':list(beta.ppf([.025,.975],c,d)),
            'difference_mean':a/(a+b)-c/(c+d),
            'difference_95':[brentq(lambda z:1-tail(a,b,c,d,z)-p,-1.,1.,xtol=2e-12) for p in (.025,.975)],
            'P_difference_gt_zero':p0,'P_difference_gt_0_10':tail(a,b,c,d,.1),
            'prior_predictive_poster_zero':__import__('math').exp(betaln(a0,b0+nP)-betaln(a0,b0)),
            'posterior_predictive_poster_zero':__import__('math').exp(betaln(c,d+nP)-betaln(c,d))}
    return result
def calculate():
    n,desk,poster=read_counts()
    assert (n,desk,poster)==(36,8,0)
    return {'spec_id':'unit07-help-seeking-routes-v1','data_sha256':hashlib.sha256(DATA.read_bytes()).hexdigest(),
            'n':n,'desk_n':24,'desk_y':desk,'poster_n':12,'poster_y':poster,
            'primary':case((1,1)),'sensitivity':case((2,2)),
            'scope':'fictional convenience volunteers; route association only; no causal effect or campus rate; exact no MCMC'}
if __name__=='__main__':
    r=calculate();print('n=36 desk=8/24 poster=0/12 raw_difference=0.333')
    for name in ('primary','sensitivity'):
        v=r[name]
        print(f"{name} prior=Beta({v['prior'][0]},{v['prior'][1]}) desk=Beta({v['desk_posterior'][0]},{v['desk_posterior'][1]}) poster=Beta({v['poster_posterior'][0]},{v['poster_posterior'][1]})")
        print(f"{name} means=[{v['desk_mean']:.3f},{v['poster_mean']:.3f}] difference={v['difference_mean']:.3f} difference_95=[{v['difference_95'][0]:.3f},{v['difference_95'][1]:.3f}] P_diff_gt_zero={v['P_difference_gt_zero']:.3f} P_diff_gt_0.10={v['P_difference_gt_0_10']:.3f}")
        print(f"{name} prior_zero={v['prior_predictive_poster_zero']:.3f} replicated_zero={v['posterior_predictive_poster_zero']:.3f}")
    print('scope='+r['scope'])
