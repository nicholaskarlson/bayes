#!/usr/bin/env python3
"""Build an accessible model-comparison SVG from the reviewed reference record."""
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
REFERENCE=ROOT/'expected-results/unit04_likelihood_reference.json'

def make_svg():
    record=json.loads(REFERENCE.read_text(encoding='utf-8'))
    cases=record['cases']
    x=lambda p:280+(p-.45)*1200
    lines=['<svg xmlns="http://www.w3.org/2000/svg" role="img" aria-labelledby="title desc" viewBox="0 0 900 480">',
      '<title id="title">Illustrative likelihood sensitivity for a two-point threshold</title>',
      '<desc id="desc">On all 64 pairs, probability the mean exceeds two is 0.636 for a Normal and 0.750 for a Student-t with four degrees of freedom and matched conditional variance. Omitting P059 yields 0.530 and 0.687; omitting P042 yields 0.732 and 0.806. Both comparison models share a Normal mean prior and conditional variance sixteen; these are not model probabilities.</desc>',
      '<rect width="900" height="480" fill="white"/>',
      '<g font-family="DejaVu Sans, sans-serif" fill="#1b2733">',
      '<text x="60" y="48" font-size="24" font-weight="bold">How does the likelihood change the answer?</text>',
      '<text x="60" y="76" font-size="15">Same mean prior and conditional variance in both models; threshold μ &gt; 2 points</text>',
      '<text x="60" y="122" font-size="17">All 64 pairs</text>',
      '<text x="60" y="230" font-size="17">Omit P059 (+10.8)</text>',
      '<text x="60" y="338" font-size="17">Omit P042 (−6.5)</text>',
      '</g>']
    for y in (158,266,374):
        lines.append(f'<line x1="280" y1="{y}" x2="788" y2="{y}" stroke="#ccd4da" stroke-width="2"/>')
    for name,y in (('all',136),('omit_P059',244),('omit_P042',352)):
        for model,dy,color in (('normal_fixed',0,'#266a96'),('student_t4_equal_variance',32,'#b65321')):
            value=cases[name][model]['p_mu_gt_2'];xx=x(value)
            lines.append(f'<circle cx="{xx:.1f}" cy="{y+dy}" r="8" fill="{color}"/>')
            lines.append(f'<text x="{xx+14:.1f}" y="{y+dy+5}" font-family="DejaVu Sans, sans-serif" font-size="15" fill="#1b2733">{value:.3f}</text>')
    lines.append('<line x1="280" y1="405" x2="788" y2="405" stroke="#697781" stroke-width="2"/>')
    for value in (.5,.6,.7,.8,.85):
        xx=x(value)
        label=f'{value:.2f}' if value==.85 else f'{value:.1f}'
        lines.append(f'<line x1="{xx:.1f}" y1="405" x2="{xx:.1f}" y2="412" stroke="#697781"/>')
        lines.append(f'<text x="{xx:.1f}" y="434" text-anchor="middle" font-family="DejaVu Sans, sans-serif" font-size="14" fill="#405260">{label}</text>')
    lines += ['<circle cx="80" cy="455" r="6" fill="#266a96"/>',
      '<text x="94" y="460" font-family="DejaVu Sans, sans-serif" font-size="14" fill="#1b2733">Normal</text>',
      '<circle cx="200" cy="455" r="6" fill="#b65321"/>',
      '<text x="214" y="460" font-family="DejaVu Sans, sans-serif" font-size="14" fill="#1b2733">Student-t(4)</text>',
      '<text x="440" y="460" font-family="DejaVu Sans, sans-serif" font-size="13" fill="#485866">Illustrative post hoc sensitivity; no causal claim.</text>',
      '</svg>']
    return '\n'.join(lines)+'\n'

if __name__=='__main__':
    if len(sys.argv)!=2:raise SystemExit('usage: make_unit04_likelihood_figure.py OUTPUT.svg')
    Path(sys.argv[1]).write_text(make_svg(),encoding='utf-8',newline='\n')
