#!/usr/bin/env python3
"""Exact reader entry point: a beta-binomial count, no Monte Carlo."""
import csv
import json
from pathlib import Path
from scipy.stats import beta, betabinom

ROOT = Path(__file__).resolve().parents[1]
spec = json.loads((ROOT / 'model-specs/unit02_model.json').read_text())
assert spec['spec_id'] == 'unit02-beta-binomial-v1'
with (ROOT / 'data/unit02_help_seeking.csv').open(newline='') as stream:
    records = list(csv.DictReader(stream))
assert [r['participant_id'] for r in records] == [f'H{i:03d}' for i in range(1,13)]
assert all(r['help_sought_within_7_days'] in ('0','1') for r in records)
n = len(records)
y = sum(int(r['help_sought_within_7_days']) for r in records)
assert (n,y) == (12,8)
a,b = (spec['primary_prior'][key] for key in ('alpha','beta'))
posterior = beta(a+y,b+n-y)
result = {
    'spec_id': spec['spec_id'], 'n': n, 'y': y, 'sample_proportion': y/n,
    'prior': [a,b], 'posterior': [a+y,b+n-y],
    'prior_predictive_ge_10_of_12': float(betabinom.sf(9,n,a,b)),
    'prior_predictive_outer_0_to_2_or_10_to_12': float(betabinom.cdf(2,n,a,b)+betabinom.sf(9,n,a,b)),
    'posterior_mean': float(posterior.mean()),
    'posterior_equal_tail_95': [float(v) for v in posterior.ppf([.025,.975])],
    'posterior_probability_above_half': float(posterior.sf(.5)),
    'posterior_predictive_ge_8_of_12': float(betabinom.sf(7,n,a+y,b+n-y)),
    'sensitivity': [], 'computation': 'exact; no MCMC',
}
for prior in spec['sensitivity_priors']:
    a,b = prior['alpha'],prior['beta']
    fitted = beta(a+y,b+n-y)
    result['sensitivity'].append({'prior':[a,b], 'posterior':[a+y,b+n-y],
                                  'mean':float(fitted.mean()), 'probability_above_half':float(fitted.sf(.5))})
if __name__ == '__main__':
    print(f"n={n} y={y} observed={y/n:.3f} prior=Beta(2,2) posterior=Beta(10,6)")
    print(f"prior_P(Y>=10_of_12)={result['prior_predictive_ge_10_of_12']:.3f} prior_P(outer)={result['prior_predictive_outer_0_to_2_or_10_to_12']:.3f}")
    print(f"posterior_mean={result['posterior_mean']:.3f} equal_tail_95=[{result['posterior_equal_tail_95'][0]:.3f},{result['posterior_equal_tail_95'][1]:.3f}] P(theta>0.5)={result['posterior_probability_above_half']:.3f}")
    print(f"posterior_predictive_P(Yrep>=8_of_12)={result['posterior_predictive_ge_8_of_12']:.3f}")
    for item in result['sensitivity']:
        print(f"sensitivity_Beta({item['prior'][0]},{item['prior'][1]})_mean={item['mean']:.3f} P(theta>0.5)={item['probability_above_half']:.3f}")
    print('scope=conditional model probability; no programme effect or population representation')
