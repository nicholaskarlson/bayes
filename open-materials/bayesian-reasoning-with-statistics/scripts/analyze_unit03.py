#!/usr/bin/env python3
"""Matched Wilson and exact beta intervals; illustrative threshold decision."""
import csv
import json
from math import sqrt
from pathlib import Path
from scipy.stats import beta, norm

ROOT=Path(__file__).resolve().parents[1]
spec=json.loads((ROOT/'model-specs/unit03_intervals_decision.json').read_text())
assert spec['spec_id']=='unit03-intervals-decision-v1'
assert spec['source']=='data/unit02_help_seeking.csv'
with (ROOT/spec['source']).open(newline='',encoding='utf-8') as stream:
    rows=list(csv.DictReader(stream))
assert [row['participant_id'] for row in rows]==[f'H{i:03d}' for i in range(1,13)]
assert all(row['help_sought_within_7_days'] in ('0','1') for row in rows)
n=len(rows);y=sum(int(row['help_sought_within_7_days']) for row in rows)
assert (n,y)==(spec['n'],spec['y'])
a,b=spec['bayesian']['prior']
assert spec['bayesian']['posterior']==[a+y,b+n-y]
p=y/n;z=norm.ppf((1+spec['frequentist']['confidence_level'])/2)
den=1+z*z/n
center=(p+z*z/(2*n))/den
half=z*sqrt(p*(1-p)/n+z*z/(4*n*n))/den
posterior=beta(a+y,b+n-y)
t=spec['illustrative_decision']['threshold_theta']
u_yes=spec['illustrative_decision']['utility_if_above']
u_no=spec['illustrative_decision']['utility_if_not_above']
u_skip=spec['illustrative_decision']['utility_no_action']
def decision(prior):
    aa,bb=prior
    fit=beta(aa+y,bb+n-y)
    mass=float(fit.sf(t))
    value=u_yes*mass+u_no*(1-mass)
    return {'prior':prior,'posterior':[aa+y,bb+n-y],
            'probability_above_0_6':mass,'expected_action_utility':value,
            'higher_expected_utility':'action' if value>u_skip else 'no_action'}
result={'spec_id':spec['spec_id'],'n':n,'y':y,'observed_proportion':p,
        'wilson_95':[center-half,center+half],
        'posterior':[a+y,b+n-y],
        'posterior_equal_tail_95':[float(v) for v in posterior.ppf([.025,.975])],
        'posterior_mean':float(posterior.mean()),
        'primary_decision':decision([a,b]),
        'sensitivity':[decision(q) for q in spec['bayesian']['sensitivity_priors']],
        'computation':'exact; no MCMC'}
if __name__=='__main__':
    print(f"n={n} y={y} observed={p:.3f} target=reported_contact_probability")
    print(f"wilson_95=[{center-half:.3f},{center+half:.3f}] procedure=repeated_sample_coverage")
    print(f"posterior=Beta({a+y},{b+n-y}) mean={posterior.mean():.3f} credible_95=[{result['posterior_equal_tail_95'][0]:.3f},{result['posterior_equal_tail_95'][1]:.3f}]")
    for item in [result['primary_decision']]+result['sensitivity']:
        aa,bb=item['prior']
        print(f"Beta({aa},{bb}) P(theta>0.6)={item['probability_above_0_6']:.3f} toy_expected_utility={item['expected_action_utility']:.3f} toy_choice={item['higher_expected_utility']}")
    print('scope=fictional conditional model and hypothetical utility; no orientation effect or real service recommendation')
