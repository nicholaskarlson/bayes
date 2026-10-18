#!/usr/bin/env python3
"""Create the two-panel Unit 8 diagnostic comparison without external plotting packages."""
from pathlib import Path
from analyze_unit08 import results

ROOT=Path(__file__).resolve().parents[1]
v=results()
def bar(x,y,p,color,label):
    return (f'<text x="{x}" y="{y}" class="lab">{label}</text>'
            f'<rect x="{x+145}" y="{y-15}" width="{p*310:.2f}" height="20" fill="{color}"/>'
            f'<text x="{x+160+p*310:.2f}" y="{y+1}" class="val">{p:.3f}</text>')
rows=[]
for idx,(case,label) in enumerate([('all','All 64'),('omit_P059','Omit P059'),('omit_P042','Omit P042')]):
    y=267+idx*69
    z=v['unit04_matched_p_mu_gt_two'][case]
    rows.extend([bar(34,y,z['normal_fixed'],'#287b85',label+' Normal'),
                 bar(34,y+27,z['student_t4_equal_variance'],'#d68a38',label+' t(4)')])
svg=f'''<svg xmlns="http://www.w3.org/2000/svg" width="960" height="550" viewBox="0 0 960 550" role="img" aria-labelledby="title desc">
<title id="title">Unit 8 model criticism and matched sensitivity</title>
<desc id="desc">Panel one: prior chance of zero poster reports among twelve is 0.029 with Beta two two and 0.077 with Beta one one. Panel two: probabilities that the paired mean exceeds two points under matched Normal and Student t four likelihoods. For all 64, 0.636 and 0.750; without P059, 0.530 and 0.687; without P042, 0.732 and 0.806. A classroom 0.60 threshold is crossed by the Normal result when P059 is omitted. This does not choose a real action.</desc>
<style>text{{font-family:DejaVu Sans,Arial,sans-serif;fill:#173247}}.head{{font-size:22px;font-weight:bold}}.lab{{font-size:17px}}.val{{font-size:17px;font-weight:bold}}.note{{font-size:15px}}.axis{{font-size:15px}}</style>
<rect width="960" height="550" fill="#fff"/><text x="28" y="37" class="head">A rare zero, then an influential paired record</text>
<rect x="22" y="54" width="916" height="131" fill="#f2f6f6"/><text x="34" y="78" class="lab">Prior chance of 0/12 poster reports</text>
{bar(34,114,v['unit07']['prior_predictive_zero_beta_2_2'],'#287b85','Beta(2,2)')}
{bar(34,154,v['unit07']['prior_predictive_zero_beta_1_1'],'#d68a38','Beta(1,1)')}
<text x="34" y="223" class="lab">P(paired mean &gt; 2) with the same Normal(0, 5²) prior</text>
<line x1="179" y1="236" x2="179" y2="458" stroke="#526675"/><line x1="365" y1="236" x2="365" y2="458" stroke="#805140" stroke-dasharray="5 4"/>
<text x="365" y="476" class="axis">Toy 0.60 threshold</text>
{''.join(rows)}
<text x="34" y="513" class="note">Normal SD 4; Student-t(4) scale √8; conditional variance 16. Omission is a diagnostic, not deletion.</text>
<text x="34" y="536" class="note">Post hoc sensitivity and fictional data. No causal effect, model probability, or validated action threshold.</text>
</svg>'''
(ROOT/'figures/unit08_diagnostics.svg').write_text(svg+'\n',encoding='utf-8')
print('UNIT08_FIGURE_OK')
