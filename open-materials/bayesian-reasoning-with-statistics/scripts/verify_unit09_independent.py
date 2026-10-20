#!/usr/bin/env python3
"""Independent gate for the Unit 9 runnable handoff."""
import csv
from fractions import Fraction
from hashlib import sha256
import json
from math import comb,exp,lgamma
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[1]
def need(ok,why):
    if not ok:raise SystemExit('UNIT09_INDEPENDENT_ERROR '+why)
def near(x,y,tol=2e-12):return abs(float(x)-float(y))<tol
def finite_tail(a,b,c,d):
    n=c+d-1
    return sum(comb(n,k)*exp(lgamma(a+k)+lgamma(b+n-k)-lgamma(a+b+n)-lgamma(a)-lgamma(b)+lgamma(a+b))
               for k in range(c,n+1))

def main():
    data=ROOT/'data/unit07_help_seeking_routes.csv'
    reference=ROOT/'expected-results/unit07_proportions.json'
    model=ROOT/'model-specs/unit07_proportions.json'
    plan=ROOT/'model-specs/unit09_handoff.json'
    with data.open(newline='', encoding='utf-8') as h:rows=list(csv.DictReader(h))
    need(len(rows)==36 and len({r['volunteer_id'] for r in rows})==36,'rows/IDs')
    desk=[r for r in rows if r['recruitment_route']=='desk'];poster=[r for r in rows if r['recruitment_route']=='poster']
    need(len(desk)==24 and len(poster)==12 and sum(int(r['reported_contact_7d']) for r in desk)==8
         and sum(int(r['reported_contact_7d']) for r in poster)==0,'raw counts')
    ref=json.loads(reference.read_text(encoding='utf-8'));spec=json.loads(model.read_text(encoding='utf-8'));card=json.loads(plan.read_text(encoding='utf-8'))
    need(ref['data_sha256']==sha256(data.read_bytes()).hexdigest(),'data digest')
    need(card['case']==ref['spec_id']==spec['spec_id'],'model identity')
    with tempfile.TemporaryDirectory(prefix='unit09-handoff-audit-') as temp:
        subprocess.run([sys.executable,str(ROOT/'scripts/build_unit09_handoff.py'),'--outdir',temp],check=True,capture_output=True)
        evidence=json.loads((Path(temp)/'unit09_case_evidence.json').read_text(encoding='utf-8'))
        report=(Path(temp)/'unit09_case_report.md').read_text(encoding='utf-8')
        sources=evidence['sources']
        for key,path in [('data_sha256',data),('unit07_reference_sha256',reference),
                         ('unit07_model_sha256',model),('unit09_plan_sha256',plan)]:
            need(sources[key]==sha256(path.read_bytes()).hexdigest(),'source '+key)
        need(evidence['counts']=={'desk':{'n':24,'reports':8},'poster':{'n':12,'reports':0}},'evidence counts')
        need(near(evidence['raw_desk_minus_poster'],Fraction(1,3)),'raw gap')
        for label,prior in [('primary',1),('sensitivity',2)]:
            x=evidence[label];r=ref[label]
            A,B,C,D=prior+8,prior+16,prior,prior+12
            need(x['prior']==[prior,prior] and x['desk_posterior']==[A,B] and
                 x['poster_posterior']==[C,D],'shapes '+label)
            need(near(x['difference_mean'],Fraction(A,A+B)-Fraction(C,C+D)),'mean '+label)
            need(near(x['P_difference_gt_zero'],finite_tail(A,B,C,D)),'finite-tail '+label)
            need(near(x['P_difference_gt_zero'],r['P_difference_gt_zero']),'reference '+label)
            for field in ('difference_95','P_difference_gt_0_10','prior_predictive_poster_zero',
                          'posterior_predictive_poster_zero'):
                need(x[field]==r[field],'traced '+field+' '+label)
        for section in ('## Question and study card','## Data and familiar description',
                        '## Model, prior and conditional result','## Checks and sensitivity',
                        '## Reportable interpretation and next evidence','## Reproduction and provenance'):
            need(section in report,'report section '+section)
        for phrase in ('8/24','0/12','0.275','[0.042, 0.490]','0.987','[-0.017, 0.456]',
                       'no campus sampling frame','not an effect','no policy decision','does not reintegrate'):
            need(phrase.lower() in report.lower(),'report '+phrase)
    print('BAYESIAN_UNIT09_INDEPENDENT_OK cases=1 report=deterministic finite_tail=yes no_new_data=yes')
if __name__=='__main__':main()
