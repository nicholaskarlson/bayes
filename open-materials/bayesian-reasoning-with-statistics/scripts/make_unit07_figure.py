#!/usr/bin/env python3
"""Generate two-panel probability-scale intervals from the exact reference."""
import argparse
from pathlib import Path
from analyze_unit07 import calculate
ROOT=Path(__file__).resolve().parents[1]
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=ROOT/'figures/unit07_proportions.svg')
    path=parser.parse_args().output
    r=calculate();p=r['primary'];s=r['sensitivity']
    x=lambda v: 310+(v+.1)*1350
    parts=['<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="620" viewBox="0 0 1200 620" role="img" aria-labelledby="title desc">',
    '<title id="title">Help-seeking proportions and route difference</title>',
    '<desc id="desc">For fictional convenience volunteers, desk has 8 reports among 24, poster zero among 12. Under independent Beta one one priors the desk posterior mean is 0.346 and poster mean 0.071; their credible intervals are shown separately. The desk minus poster probability difference has mean 0.275 and 95 percent credible interval 0.042 to 0.490. Under Beta two two priors the contrast interval is minus 0.017 to 0.456. Zero reports does not imply zero probability; intervals do not establish a route effect or a campus rate.</desc>',
    '<rect width="1200" height="620" fill="white"/>',
    '<text x="65" y="55" fill="#17354a" font-family="DejaVu Sans,Arial" font-size="31" font-weight="bold">A zero count still leaves uncertainty</text>',
    '<text x="65" y="90" fill="#455d69" font-family="DejaVu Sans,Arial" font-size="19">Fictional volunteers; exact Beta updates; intervals conditional on the stated priors</text>',
    '<text x="65" y="145" fill="#17354a" font-family="DejaVu Sans,Arial" font-size="22" font-weight="bold">Reporting probability by route · Beta(1,1)</text>',
    '<text x="65" y="374" fill="#17354a" font-family="DejaVu Sans,Arial" font-size="22" font-weight="bold">Desk minus poster probability · prior sensitivity</text>']
    # Separate scales share x mapping for display but use labeled tick spans.
    def axis(y,lo,hi,xx,step):
        parts.append(f'<line x1="{xx(lo):.1f}" y1="{y}" x2="{xx(hi):.1f}" y2="{y}" stroke="#738d99" stroke-width="2"/>')
        k=lo
        while k<hi+1e-9:
            z=xx(k);parts.append(f'<line x1="{z:.1f}" y1="{y-5}" x2="{z:.1f}" y2="{y+5}" stroke="#738d99"/>')
            parts.append(f'<text x="{z:.1f}" y="{y+29}" text-anchor="middle" fill="#334d5b" font-family="DejaVu Sans,Arial" font-size="17">{k:.1f}</text>');k+=step
    xm=lambda v:330+v*1250
    axis(294,0,.6,xm,.1)
    for label,mean,interval,y,color in [('Desk 8/24',p['desk_mean'],p['desk_95'],205,'#0a6581'),('Poster 0/12',p['poster_mean'],p['poster_95'],255,'#a86122')]:
        parts.extend([f'<text x="65" y="{y+6}" fill="#263f4d" font-family="DejaVu Sans,Arial" font-size="20">{label}</text>',
        f'<line x1="{xm(interval[0]):.1f}" y1="{y}" x2="{xm(interval[1]):.1f}" y2="{y}" stroke="{color}" stroke-width="7"/>',
        f'<circle cx="{xm(mean):.1f}" cy="{y}" r="9" fill="{color}" stroke="white" stroke-width="2"/>'])
    xd=lambda v:410+(v+.1)*1050
    axis(548,-.1,.6,xd,.1)
    zero=xd(0);parts.append(f'<line x1="{zero:.1f}" y1="410" x2="{zero:.1f}" y2="548" stroke="#9a4141" stroke-width="2" stroke-dasharray="7 5"/>')
    for label,v,y,color in [('Beta(1,1)',p,445,'#0a6581'),('Beta(2,2)',s,502,'#a86122')]:
        lo,hi=v['difference_95'];m=v['difference_mean']
        parts.extend([f'<text x="65" y="{y+6}" fill="#263f4d" font-family="DejaVu Sans,Arial" font-size="20">{label}</text>',
        f'<line x1="{xd(lo):.1f}" y1="{y}" x2="{xd(hi):.1f}" y2="{y}" stroke="{color}" stroke-width="7"/>',
        f'<circle cx="{xd(m):.1f}" cy="{y}" r="9" fill="{color}" stroke="white" stroke-width="2"/>'])
    parts.extend(['<text x="1120" y="605" text-anchor="end" fill="#405966" font-family="DejaVu Sans,Arial" font-size="17">Dot: posterior mean · bar: central 95% credible interval</text>','</svg>'])
    path.parent.mkdir(parents=True,exist_ok=True);path.write_text('\n'.join(parts)+'\n',encoding='utf-8')
    print('UNIT07_FIGURE_OK '+str(path))
if __name__=='__main__':main()
