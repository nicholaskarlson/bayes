#!/usr/bin/env python3
"""Draw three deterministic Primer figures from checked companion anchors."""
import csv
import json
from math import comb, gamma
from pathlib import Path
from xml.sax.saxutils import escape

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / 'figures'
NAVY, TEAL, ORANGE, GRAY = '#173e59', '#087d84', '#a64f12', '#7b8790'


def base(title, desc):
    return [f'<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="650" viewBox="0 0 1200 650" role="img" aria-labelledby="title desc">',
            f'<title id="title">{escape(title)}</title><desc id="desc">{escape(desc)}</desc>',
            '<rect width="1200" height="650" fill="white"/>',
            f'<text x="50" y="56" font-size="34" font-weight="bold" fill="{NAVY}" font-family="DejaVu Sans,Arial">{escape(title)}</text>']


def line(x1, y1, x2, y2, color='#bdc9ce', width=2, extra=''):
    return f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{color}" stroke-width="{width}" {extra}/>'


def label(x, y, content, size=21, color=NAVY, anchor='start'):
    return f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" fill="{color}" font-family="DejaVu Sans,Arial" text-anchor="{anchor}">{escape(content)}</text>'


def panel_ticks(parts, x0, y0, width, height, ymax):
    parts.append(line(x0, y0, x0+width, y0))
    parts.append(line(x0, y0, x0, y0-height))
    for k in range(0, 11, 2):
        x = x0 + width*k/10
        parts.append(line(x, y0, x, y0+7))
        parts.append(label(x, y0+32, f'{k/10:g}', 19, anchor='middle'))
    parts.append(label(x0+width/2, y0+68, 'θ or sample proportion', 19, anchor='middle'))
    parts.append(label(x0, y0-height-10, f'vertical scale: {ymax:g}', 18, GRAY))


def beta_pdf(t):
    if t <= 0 or t >= 1:
        return 0.0
    return gamma(16)/(gamma(10)*gamma(6))*t**9*(1-t)**5


def density_panel(parts, x0, y0, width, height, title):
    panel_ticks(parts, x0, y0, width, height, 3.5)
    parts.append(label(x0+width/2, 126, title, 25, NAVY, 'middle'))
    xhalf = x0 + width/2
    parts.append(line(xhalf, y0, xhalf, y0-height, GRAY, 2, 'stroke-dasharray="5 6"'))
    points = [(x0+width*i/250, y0-height*beta_pdf(i/250)/3.5) for i in range(251)]
    right = points[125:]
    path = f'M {xhalf:.1f} {y0:.1f} ' + ' '.join(f'L {x:.1f} {y:.1f}' for x,y in right) + f' L {x0+width:.1f} {y0:.1f} Z'
    parts.append(f'<path d="{path}" fill="#bce5df" opacity="0.8"/>')
    parts.append('<polyline fill="none" stroke="'+TEAL+'" stroke-width="5" points="'+' '.join(f'{x:.1f},{y:.1f}' for x,y in points)+'"/>')
    parts.append(label(x0+width/2, 574, 'Area above 0.5 = 0.849121', 20, TEAL, 'middle'))


def grid_density():
    with (ROOT/'expected-results/prelude_grid.csv').open(newline='') as f:
        rows = list(csv.DictReader(f))
    parts = base('Discrete grid and continuous Beta posterior',
                 'Nine discrete posterior bars place 0.183 mass exactly at theta 0.5; bars strictly above 0.5 total 0.745. The Beta(10,6) density has shaded continuous area above 0.5 of 0.849. Both horizontal axes are zero to one, but the left vertical axis is probability mass and the right is density height.')
    x0, y0, w, h = 80, 480, 430, 280
    panel_ticks(parts,x0,y0,w,h,.35)
    parts.append(label(x0+w/2,126,'Nine-value grid: probability mass',25,NAVY,'middle'))
    for r in rows:
        t,m = float(r['theta']),float(r['posterior'])
        x = x0+w*t
        parts.append(f'<rect x="{x-13:.1f}" y="{y0-h*m/.35:.1f}" width="26" height="{h*m/.35:.1f}" fill="{ORANGE if t==.5 else TEAL if t>.5 else GRAY}"/>')
    parts.append(label(x0+w/2,574,'Mass above 0.5 = 0.745038',20,TEAL,'middle'))
    density_panel(parts,680,480,430,280,'Beta(10,6): density')
    parts.append('</svg>')
    (FIG/'appendix_primer_grid_density.svg').write_text('\n'.join(parts)+'\n')


def sampling_posterior():
    parts=base('Possible sample proportions and a posterior',
               'Left, a Binomial(12, two-thirds) sampling distribution over possible sample proportions with the model theta fixed at two-thirds. Right, the Beta(10,6) posterior density over theta conditional on eight of twelve reports and a Beta(2,2) prior. The horizontal scales match but the random quantities and vertical scales differ.')
    x0,y0,w,h=80,480,430,280
    panel_ticks(parts,x0,y0,w,h,.25)
    parts.append(label(x0+w/2,126,'Sampling: fixed θ = 2/3',25,NAVY,'middle'))
    for k in range(13):
        p=comb(12,k)*(2/3)**k*(1/3)**(12-k)
        x=x0+w*k/12
        parts.append(f'<rect x="{x-8:.1f}" y="{y0-h*p/.25:.1f}" width="16" height="{h*p/.25:.1f}" fill="{ORANGE if k==8 else GRAY}"/>')
    parts.append(label(x0+w/2,574,'Possible Y/12 across repeated samples',20,GRAY,'middle'))
    density_panel(parts,680,480,430,280,'Posterior: observed Y = 8')
    parts.append('</svg>')
    (FIG/'appendix_primer_sampling_posterior.svg').write_text('\n'.join(parts)+'\n')


def shrinkage():
    raw=json.loads((ROOT/'expected-results/unit06_partial_pooling.json').read_text())
    rows=[(key,raw['raw_means'][key],raw['primary_tau5']['theta_mean'][i]) for i,key in enumerate(raw['raw_means'])]
    parts=base('Eight workshops: raw and partially pooled means',
               'Eight horizontal line segments connect gray raw workshop means to blue posterior means under Unit 6 fixed-scale tau equals five model. W01 with n four moves from 63.625 to 60.272, while W08 with n twenty moves from 59.345 to 58.856. The learned common center near 55.032 is marked by a dashed vertical line; movement is model conditional, not an intervention effect.')
    x0,w=260,800
    scale=lambda v:x0+(v-45)/25*w
    center=raw['primary_tau5']['mu_mean']
    parts.append(line(scale(center),145,scale(center),565,ORANGE,3,'stroke-dasharray="8 7"'))
    parts.append(label(scale(center)+10,136,'shared center 55.032',19,ORANGE))
    for tick in (45,50,55,60,65,70):
        x=scale(tick)
        parts.append(line(x,565,x,575))
        parts.append(label(x,601,str(tick),20,anchor='middle'))
    parts.append(line(x0,565,x0+w,565))
    sizes=[4,5,6,8,10,12,15,20]
    for i,(key,observed,post) in enumerate(rows):
        y=180+i*50
        parts.append(label(65,y+7,f'{key}  n={sizes[i]}',21))
        parts.append(line(scale(observed),y,scale(post),y,TEAL,4))
        parts.append(f'<circle cx="{scale(observed):.1f}" cy="{y}" r="7" fill="{GRAY}"/>')
        parts.append(f'<circle cx="{scale(post):.1f}" cy="{y}" r="8" fill="{NAVY}"/>')
        parts.append(label(1140,y+7,f'{observed:.1f} → {post:.1f}',19,anchor='end'))
    parts.append(label(290,634,'Gray: observed group mean',19,GRAY))
    parts.append(label(700,634,'Blue: conditional posterior mean',19,NAVY))
    parts.append('</svg>')
    (FIG/'appendix_primer_shrinkage.svg').write_text('\n'.join(parts)+'\n')


if __name__=='__main__':
    grid_density(); sampling_posterior(); shrinkage()
    print('PRIMER_FIGURES_OK files=3')
