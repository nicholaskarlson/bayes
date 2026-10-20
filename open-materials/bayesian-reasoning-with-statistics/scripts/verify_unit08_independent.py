#!/usr/bin/env python3
"""Independent arithmetic, source and sign audit for Unit 8."""
import csv
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from math import erfc,sqrt
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
def need(test,what):
    if not test:raise SystemExit('UNIT08_INDEPENDENT_ERROR '+what)
def close(a,b,tol=2e-12):return abs(float(a)-float(b))<tol
def main():
    want=json.loads((ROOT/'expected-results/unit08_diagnostics.json').read_text(encoding='utf-8'))
    u7=want['unit07'];paired=want['unit04_matched_p_mu_gt_two']
    with (ROOT/'data/unit07_help_seeking_routes.csv').open(newline='', encoding='utf-8') as h: routes=list(csv.DictReader(h))
    poster=[r for r in routes if r['recruitment_route']=='poster']
    need(len(poster)==12 and all(r['reported_contact_7d']=='0' for r in poster),'poster records')
    # Integer beta-binomial products, independent of lgamma implementation.
    beta22=Fraction(2*3,14*15)  # B(2,14)/B(2,2) = 6/(14*15)
    beta11=Fraction(1,13)
    need(close(u7['prior_predictive_zero_beta_2_2'],beta22),'Beta(2,2) zero')
    need(close(u7['prior_predictive_zero_beta_1_1'],beta11),'Beta(1,1) zero')
    with (ROOT/'data/unit04_reflection_attention_paired.csv').open(newline='', encoding='utf-8') as h: records=list(csv.DictReader(h))
    need(len(records)==64,'paired records')
    with (ROOT/'expected-results/unit04_likelihood_anchors.csv').open(newline='', encoding='utf-8') as h:
        anchors={(r['case'],r['model'],r['field']):float(r['value']) for r in csv.DictReader(h)}
    ref=json.loads((ROOT/'expected-results/unit04_likelihood_reference.json').read_text(encoding='utf-8'))
    for case,exclude in [('all',None),('omit_P059','P059'),('omit_P042','P042')]:
        diffs=[Decimal(r['after_score'])-Decimal(r['before_score']) for r in records if r['participant_id']!=exclude]
        n=len(diffs);v=1/(1/25+n/16)
        m=(float(sum(diffs))/16)*v
        p=0.5*erfc((2-m)/sqrt(2*v))
        need(close(p,paired[case]['normal_fixed']),'Normal from raw '+case)
        for model in ('normal_fixed','student_t4_equal_variance'):
            val=paired[case][model]
            need(close(val,anchors[(case,model,'p_mu_gt_2')]),'existing anchor '+case+model)
            need(close(val,ref['cases'][case][model]['p_mu_gt_2']),'existing reference '+case+model)
    old=json.loads((ROOT/'expected-results/unit04_checks.json').read_text(encoding='utf-8'))['leave_one_out']
    nig=want['unit04_original_nig_p_mu_gt_two']
    need(close(nig['all'],old['full_p_mu_gt_two']) and close(nig['omit_P059'],old['min']['p_mu_gt_two']) and
         close(nig['omit_P042'],old['max']['p_mu_gt_two']),'NIG trace')
    p=u7['p_desk_minus_poster_gt_zero']
    u7ref=json.loads((ROOT/'expected-results/unit07_proportions.json').read_text(encoding='utf-8'))
    need(close(p,u7ref['primary']['P_difference_gt_zero']) and close(1-p,0.012532693984307),'sign reversal')
    # A deliberate wrong-sign transcription must fail this guard.
    need(not close(p,1-p) and p>0.98 and 1-p<0.02,'wrong-sign trap')
    need(len(json.loads((ROOT/'model-specs/unit08_diagnostic_plan.json').read_text(encoding='utf-8'))['cases'])==2,'model card')
    print('BAYESIAN_UNIT08_INDEPENDENT_OK zero=exact normal=raw matched_t=prior_audited sign=guarded mcmc=no')
if __name__=='__main__':main()
