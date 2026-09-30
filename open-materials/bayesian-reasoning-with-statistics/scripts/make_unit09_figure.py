#!/usr/bin/env python3
"""Make a traceable claim ladder from the committed Unit 7 reference."""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
r=json.loads((ROOT/'expected-results/unit07_proportions.json').read_text())
assert (r['desk_n'],r['desk_y'],r['poster_n'],r['poster_y'])==(24,8,12,0)
p=r['primary'];s=r['sensitivity']
svg=f'''<svg xmlns="http://www.w3.org/2000/svg" width="960" height="500" viewBox="0 0 960 500" role="img" aria-labelledby="title desc">
<title id="title">Unit 9 claim ladder for the fictional route case</title>
<desc id="desc">Four vertically ordered evidence steps: the file contains eight of twenty-four desk reports and zero of twelve poster reports; an independent Beta one one model gives a posterior mean route difference of {p['difference_mean']:.3f} and 95 percent interval {p['difference_95'][0]:.3f} to {p['difference_95'][1]:.3f}; a Beta two two sensitivity interval {s['difference_95'][0]:.3f} to {s['difference_95'][1]:.3f} includes zero; neither result identifies an intervention effect or a campus rate. The report recommends design and measurement work before action.</desc>
<style>text{{font-family:DejaVu Sans,Arial,sans-serif;fill:#173247}}.title{{font-size:23px;font-weight:bold}}.head{{font-size:21px;font-weight:bold}}.body{{font-size:18px}}.num{{font-size:24px;font-weight:bold}}.small{{font-size:16px}} </style>
<rect width="960" height="500" fill="white"/><text x="30" y="39" class="title">Move from records to a bounded report</text>
<rect x="26" y="57" width="908" height="88" rx="9" fill="#e8f2f3"/><text x="44" y="88" class="head">1  Describe the file</text><text x="475" y="99" class="num">8/24 desk · 0/12 poster</text><text x="44" y="123" class="small">Fictional convenience volunteers; one seven-day self-report each.</text>
<rect x="26" y="155" width="908" height="88" rx="9" fill="#f4f3e8"/><text x="44" y="186" class="head">2  State the model</text><text x="475" y="197" class="num">Δ mean {p['difference_mean']:.3f}</text><text x="44" y="221" class="small">Independent Beta(1,1); 95% credible interval [{p['difference_95'][0]:.3f}, {p['difference_95'][1]:.3f}].</text>
<rect x="26" y="253" width="908" height="88" rx="9" fill="#f9eee5"/><text x="44" y="284" class="head">3  Challenge the model</text><text x="475" y="295" class="num">Beta(2,2) includes zero</text><text x="44" y="319" class="small">Interval [{s['difference_95'][0]:.3f}, {s['difference_95'][1]:.3f}]; selection and measurement still need an audit.</text>
<rect x="26" y="351" width="908" height="99" rx="9" fill="#edf0f4"/><text x="44" y="382" class="head">4  Set the claim ceiling</text><text x="475" y="393" class="num">No route effect inferred</text><text x="44" y="417" class="body">Define an intervention, validate the outcome and obtain a suitable comparison.</text>
<text x="31" y="481" class="small">A probability of association under a model cannot supply missing assignment or a campus sampling frame.</text>
</svg>'''
(ROOT/'figures/unit09_claim_ladder.svg').write_text(svg+'\n',encoding='utf-8')
print('UNIT09_FIGURE_OK')
