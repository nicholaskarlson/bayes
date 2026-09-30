#!/usr/bin/env python3
"""Rebuild the exact prior/posterior density figure as accessible SVG."""

import json
import math
import sys
from pathlib import Path
from scipy.stats import t

ROOT = Path(__file__).resolve().parents[1]
record = json.loads((ROOT / "expected-results" / "unit04_exact.json").read_text(encoding="utf-8"))
primary = record["bayesian_exact"]["primary"]
prior = record["bayesian_exact"]["prior"]
destination = (Path(sys.argv[1]) if len(sys.argv)==2 else ROOT / "figures" / "unit04_prior_posterior.svg")
if len(sys.argv)>2:
    raise SystemExit("usage: make_unit04_figure.py [OUTPUT.svg]")

width, height = 820, 440
left, right, top, bottom = 72, 785, 42, 354
xmin, xmax, ymax = -12.0, 12.0, 0.84
prior_scale = math.sqrt(float(prior["b0"])/(float(prior["a0"])*float(prior["kappa0"])))
prior_dist = t(df=2*float(prior["a0"]), loc=float(prior["m0"]), scale=prior_scale)
post_dist = t(df=primary["df"], loc=primary["m_n"], scale=primary["posterior_mu_scale"])


def px(x): return left + (x-xmin)/(xmax-xmin)*(right-left)
def py(y): return bottom - y/ymax*(bottom-top)
def path(dist):
    return " ".join(("M" if i == 0 else "L") + f"{px(x):.2f},{py(float(dist.pdf(x))):.2f}"
                    for i,x in enumerate(xmin+(xmax-xmin)*j/360 for j in range(361)))


svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-labelledby="title desc">
<title id="title">Unit 4 prior and posterior for mean paired change</title>
<desc id="desc">Prior density centered at zero is broad. Exact posterior density is narrower and centered near {primary['m_n']:.2f} points. A dashed vertical line marks zero and a dotted vertical line marks an illustrative two-point meaningful-change threshold. The posterior overlaps both sides of two points.</desc>
<rect width="{width}" height="{height}" fill="white"/>
<text x="72" y="27" font-family="Liberation Sans, sans-serif" font-size="18" font-weight="bold">Mean paired change: prior and exact posterior</text>
<path d="M{left},{bottom} H{right} M{left},{bottom} V{top}" fill="none" stroke="#303030" stroke-width="1.5"/>
<path d="M{px(0):.2f},{top} V{bottom}" stroke="#555555" stroke-width="1.2" stroke-dasharray="6 5"/>
<path d="M{px(2):.2f},{top} V{bottom}" stroke="#805c00" stroke-width="1.5" stroke-dasharray="2 5"/>
<path d="{path(prior_dist)}" fill="none" stroke="#336da8" stroke-width="3"/>
<path d="{path(post_dist)}" fill="none" stroke="#bd5839" stroke-width="3"/>
<text x="{px(0)-8:.1f}" y="{bottom+19}" font-family="Liberation Sans, sans-serif" font-size="12">0</text>
<text x="{px(2)-6:.1f}" y="{bottom+19}" font-family="Liberation Sans, sans-serif" font-size="12">2</text>
<text x="280" y="{height-22}" font-family="Liberation Sans, sans-serif" font-size="14">Mean after-minus-before change (attention-score points)</text>
<text x="{left}" y="{top+20}" font-family="Liberation Sans, sans-serif" font-size="13" fill="#336da8">Prior</text>
<text x="{left+68}" y="{top+20}" font-family="Liberation Sans, sans-serif" font-size="13" fill="#bd5839">Posterior</text>
<text x="{right-128}" y="{top+20}" font-family="Liberation Sans, sans-serif" font-size="12" fill="#805c00">2-point threshold</text>
</svg>\n'''
destination.write_text(svg, encoding="utf-8")
print("BAYESIAN_UNIT04_FIGURE_OK", str(destination))
