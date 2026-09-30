#!/usr/bin/env python3
"""Recalculate reader-facing Unit 8 diagnostics from the existing cases."""
import csv
import json
import math
from pathlib import Path
from statistics import NormalDist

ROOT = Path(__file__).resolve().parents[1]

def rows(name):
    with (ROOT/'data'/name).open(newline='') as handle:
        return list(csv.DictReader(handle))

def prior_zero(a, b, n):
    # Beta-binomial probability of zero: B(a,b+n)/B(a,b).
    return math.exp(math.lgamma(b+n)+math.lgamma(a+b)-math.lgamma(a+b+n)-math.lgamma(b))

def normal_tail(differences):
    # Prior mu ~ N(0,25), observations ~ N(mu,16), fixed variance.
    precision = 1/25+len(differences)/16
    mean = (sum(differences)/16)/precision
    return NormalDist(mu=mean,sigma=math.sqrt(1/precision)).cdf(2)  # complemented below

def results():
    routes = rows('unit07_help_seeking_routes.csv')
    poster = [r for r in routes if r['recruitment_route']=='poster']
    desk = [r for r in routes if r['recruitment_route']=='desk']
    assert len(poster)==12 and len(desk)==24
    assert sum(int(r['reported_contact_7d']) for r in poster)==0
    assert sum(int(r['reported_contact_7d']) for r in desk)==8
    paired = rows('unit04_reflection_attention_paired.csv')
    assert len(paired)==64 and len({r['participant_id'] for r in paired})==64
    matched = {}
    anchors = {}
    with (ROOT/'expected-results/unit04_likelihood_anchors.csv').open(newline='') as handle:
        for row in csv.DictReader(handle):
            anchors[(row['case'],row['model'],row['field'])]=float(row['value'])
    for case,excluded in [('all',None),('omit_P059','P059'),('omit_P042','P042')]:
        diffs = [float(r['after_score'])-float(r['before_score']) for r in paired if r['participant_id']!=excluded]
        p_normal = 1-normal_tail(diffs)
        p_t = anchors[(case,'student_t4_equal_variance','p_mu_gt_2')]
        assert abs(p_normal-anchors[(case,'normal_fixed','p_mu_gt_2')])<1e-12
        matched[case]={'normal_fixed':p_normal,'student_t4_equal_variance':p_t}
    old = json.loads((ROOT/'expected-results/unit04_checks.json').read_text())['leave_one_out']
    return {
        'schema':'unit08-reader-diagnostics-v1',
        'unit07':{'poster_n':len(poster),'poster_reports':0,
                  'prior_predictive_zero_beta_2_2':prior_zero(2,2,len(poster)),
                  'prior_predictive_zero_beta_1_1':prior_zero(1,1,len(poster)),
                  'p_desk_minus_poster_gt_zero':json.loads((ROOT/'expected-results/unit07_proportions.json').read_text())['primary']['P_difference_gt_zero']},
        'unit04_matched_p_mu_gt_two':matched,
        'unit04_original_nig_p_mu_gt_two':{'all':old['full_p_mu_gt_two'],
            'omit_P059':old['min']['p_mu_gt_two'],'omit_P042':old['max']['p_mu_gt_two']},
        'note':'Student-t and NIG values are sourced from existing independently audited Unit 4 anchors; Normal and beta-binomial values are recalculated from raw data. No new quadrature or MCMC.'
    }

if __name__=='__main__':
    r=results()
    print('UNIT08_ZERO poster=0/12 prior_Beta22=%.6f prior_Beta11=%.6f' % (r['unit07']['prior_predictive_zero_beta_2_2'],r['unit07']['prior_predictive_zero_beta_1_1']))
    for case,p in r['unit04_matched_p_mu_gt_two'].items():
        print('UNIT08_MATCHED %s normal=%.6f t4=%.6f' % (case,p['normal_fixed'],p['student_t4_equal_variance']))
    print('UNIT08_SIGN desk_minus_poster=%.6f poster_minus_desk=%.6f' %
          (r['unit07']['p_desk_minus_poster_gt_zero'],1-r['unit07']['p_desk_minus_poster_gt_zero']))
    print('UNIT08_SCOPE toy_threshold=0.60 design_limit=yes mcmc=no')
