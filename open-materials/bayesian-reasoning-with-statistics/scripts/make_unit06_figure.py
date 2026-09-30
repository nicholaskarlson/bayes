#!/usr/bin/env python3
"""Deterministic source SVG: group means, exact intervals and pooling benchmark."""
import argparse
from pathlib import Path
from analyze_unit06 import calculate

ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'figures/unit06_pooling.svg')
    a=parser.parse_args();v=calculate();p=v['primary_tau5']
    left,right=43.,70.
    x=lambda value: 355+(float(value)-left)*36
    parts=['<svg xmlns="http://www.w3.org/2000/svg" width="1500" height="870" viewBox="0 0 1500 870" role="img" aria-labelledby="title desc">',
           '<title id="title">Eight workshop means and partial pooling</title>',
           '<desc id="desc">Workshop W01 has four students and raw mean 63.625; its exact posterior mean is 60.272 with interval 53.964 to 66.579. Workshop W08 has twenty students and raw mean 59.345; posterior mean 58.856 with interval 55.525 to 62.186. Other group means lie between about 49.5 and 58.1. Raw crosses, posterior circles and 95 percent credible intervals share a score-point axis. A dashed vertical line marks pooled mean 55.796. Intervals concern latent group means under fixed scales; no workshop effects are identified.</desc>',
           '<rect width="1500" height="820" fill="#ffffff"/>',
           '<text x="70" y="64" fill="#1c3040" font-family="DejaVu Sans,Arial,sans-serif" font-size="42" font-weight="700">Group means with partial pooling</text>',
           '<text x="70" y="98" fill="#475e6b" font-family="DejaVu Sans,Arial,sans-serif" font-size="27">Fictional confidence scores, 0–100; each student measured once</text>']
    for value in range(45,71,5):
        loc=x(value)
        parts.append(f'<line x1="{loc:.1f}" y1="145" x2="{loc:.1f}" y2="700" stroke="#e0e9ec" stroke-width="2"/>')
        parts.append(f'<text x="{loc:.1f}" y="740" text-anchor="middle" font-size="27" fill="#304653" font-family="DejaVu Sans,Arial,sans-serif">{value}</text>')
    pool=x(v['pooled_mean'])
    parts.append(f'<line x1="{pool:.1f}" y1="145" x2="{pool:.1f}" y2="700" stroke="#af651e" stroke-width="3" stroke-dasharray="9 7"/>')
    parts.append(f'<text x="{pool+8:.1f}" y="158" fill="#874e16" font-size="24" font-family="DejaVu Sans,Arial,sans-serif">Pooled 55.796</text>')
    for j,(group,n) in enumerate(v['sizes'].items()):
        y=205+j*66
        raw=v['raw_means'][group];post=p['theta_mean'][j];lo,hi=p['theta_95'][j]
        parts.append(f'<line x1="68" y1="{y+33}" x2="1430" y2="{y+33}" stroke="#eef2f3"/>')
        parts.append(f'<text x="75" y="{y+7}" fill="#243c4b" font-size="30" font-family="DejaVu Sans,Arial,sans-serif">{group}  (n={n})</text>')
        parts.append(f'<line x1="{x(lo):.1f}" y1="{y}" x2="{x(hi):.1f}" y2="{y}" stroke="#087e84" stroke-width="6"/>')
        for endpoint in (lo,hi):parts.append(f'<line x1="{x(endpoint):.1f}" y1="{y-10}" x2="{x(endpoint):.1f}" y2="{y+10}" stroke="#087e84" stroke-width="3"/>')
        parts.append(f'<line x1="{x(raw)-8:.1f}" y1="{y-8}" x2="{x(raw)+8:.1f}" y2="{y+8}" stroke="#656d76" stroke-width="3"/>')
        parts.append(f'<line x1="{x(raw)-8:.1f}" y1="{y+8}" x2="{x(raw)+8:.1f}" y2="{y-8}" stroke="#656d76" stroke-width="3"/>')
        parts.append(f'<circle cx="{x(post):.1f}" cy="{y}" r="9" fill="#0a4b71" stroke="white" stroke-width="2"/>')
        parts.append(f'<text x="1470" y="{y+7}" text-anchor="end" fill="#344955" font-size="27" font-family="DejaVu Sans,Arial,sans-serif">{raw:.1f} → {post:.1f}</text>')
    parts.extend(['<text x="820" y="780" text-anchor="middle" fill="#243c4b" font-size="30" font-family="DejaVu Sans,Arial,sans-serif">Mean confidence score (points)</text>',
                  '<text x="75" y="821" fill="#4b626c" font-size="24" font-family="DejaVu Sans,Arial,sans-serif">× Raw mean     ● Posterior mean     Teal bar: 95% credible interval</text>',
                  '<text x="75" y="851" fill="#874e16" font-size="24" font-family="DejaVu Sans,Arial,sans-serif">Dashed line: student-weighted pooled raw mean</text>',
                  '</svg>'])
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text('\n'.join(parts)+'\n',encoding='utf-8')
    print('UNIT06_FIGURE_OK '+str(a.output))
if __name__=='__main__':main()
