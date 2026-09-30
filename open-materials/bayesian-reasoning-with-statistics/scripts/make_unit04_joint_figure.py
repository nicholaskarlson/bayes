#!/usr/bin/env python3
"""Deterministic accessible SVG from the committed predictive reference record."""

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT/"expected-results/unit04_joint_reference.json"


def make_svg():
    record=json.loads(REFERENCE.read_text(encoding="utf-8"))
    observed=record["observed"]
    count=lambda x: 175+690*x/32
    skew=lambda x: 175+690*(x+.7)/1.4
    prior_count=record["prior"]["tail_count_95"]
    post_count=record["posterior"]["tail_count_95"]
    prior_skew=record["prior"]["moment_skewness_95"]
    post_skew=record["posterior"]["moment_skewness_95"]
    lines=[
        '<svg xmlns="http://www.w3.org/2000/svg" role="img" aria-labelledby="title desc" viewBox="0 0 1000 540">',
        '<title id="title">Whole-cohort prior and posterior predictive checks</title>',
        '<desc id="desc">For 64 fictional paired changes, the observed count above eight points is five. The prior middle 95 percent interval is zero to thirty; posterior is one to twelve. Observed skewness is minus 0.074, inside both roughly minus 0.58 to plus 0.58 intervals. Bars describe replicated cohorts, not uncertainty about the mean parameter.</desc>',
        '<rect width="1000" height="540" fill="#ffffff"/>',
        '<g font-family="DejaVu Sans, sans-serif" fill="#1b2733">',
        '<text x="80" y="50" font-size="26" font-weight="bold">Replicated cohorts under the stated model</text>',
        '<text x="80" y="80" font-size="16">Middle 95% intervals; each simulated cohort has 64 paired differences</text>',
        '<text x="80" y="125" font-size="20" font-weight="bold">Count with |change| &gt; 8 points</text>',
        '<text x="80" y="335" font-size="20" font-weight="bold">Moment skewness of changes</text>',
        '</g>',
    ]
    for label,y in (("Prior",175),("Posterior",220),("Observed",265),("Prior",385),("Posterior",430),("Observed",475)):
        lines.append(f'<text x="80" y="{y+6}" font-family="DejaVu Sans, sans-serif" font-size="17" fill="#1b2733">{label}</text>')
    for axis_y, ticks, transform in ((282,range(0,33,8),count),(492,(-.6,-.3,0,.3,.6),skew)):
        lines.append(f'<line x1="175" y1="{axis_y}" x2="865" y2="{axis_y}" stroke="#6d7884" stroke-width="1.5"/>')
        for tick in ticks:
            x=transform(tick)
            lines.append(f'<line x1="{x:.1f}" y1="{axis_y}" x2="{x:.1f}" y2="{axis_y+7}" stroke="#6d7884"/>')
            lines.append(f'<text x="{x:.1f}" y="{axis_y+25}" text-anchor="middle" font-family="DejaVu Sans, sans-serif" font-size="14" fill="#465462">{tick:g}</text>')
    for ys, values, transform in (((170,215),(prior_count,post_count),count),
                                  ((380,425),(prior_skew,post_skew),skew)):
        for y, interval, color in zip(ys,values,("#2a6f9e","#bd5b28")):
            a,b=map(transform,interval)
            lines.append(f'<line x1="{a:.1f}" y1="{y}" x2="{b:.1f}" y2="{y}" stroke="{color}" stroke-width="14" stroke-linecap="round"/>')
            for x in (a,b):
                lines.append(f'<line x1="{x:.1f}" y1="{y-12}" x2="{x:.1f}" y2="{y+12}" stroke="{color}" stroke-width="3"/>')
    for y,x,value in ((260,count(observed["count_abs_gt_8"]),"5"),
                      (470,skew(observed["moment_skewness"]),"-0.074")):
        lines.append(f'<circle cx="{x:.1f}" cy="{y}" r="7" fill="#15252d"/>')
        lines.append(f'<text x="{x+13:.1f}" y="{y+6}" font-family="DejaVu Sans, sans-serif" font-size="16" fill="#15252d">{value}</text>')
    lines += [
        '<text x="80" y="530" font-family="DejaVu Sans, sans-serif" font-size="14" fill="#465462">Illustrative thresholds; predictive fit cannot validate design or measurement.</text>',
        '</svg>'
    ]
    return "\n".join(lines)+"\n"


if __name__ == "__main__":
    if len(sys.argv)!=2:raise SystemExit("usage: make_unit04_joint_figure.py OUTPUT.svg")
    Path(sys.argv[1]).write_text(make_svg(),encoding="utf-8")
