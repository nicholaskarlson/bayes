# Bayesian Reasoning with Statistics — public numerical companion

This is the public companion to Nicholas Elliott Karlson's *Bayesian Reasoning
with Statistics: Models, Evidence, and Decisions with Python and R*. It is a
new-history repository containing the book's **fictional** datasets, variable
dictionaries, provenance, model specifications, expected numerical results,
figures, technical notes and runnable code. It does not contain the book text
or the private authoring history. Record this repository's Git commit when
reporting a reproduction, or the release tag (see Versions below).

## Run the checks on Ubuntu

```bash
git clone https://github.com/nicholaskarlson/bayes.git
cd bayes
python3 -m venv .venv
.venv/bin/python -m pip install -r open-materials/bayesian-reasoning-with-statistics/scripts/requirements.txt
make verify PYTHON=.venv/bin/python
make r-audit
```

`make verify` validates the SHA-256 inventory, checks raw counts and selected
closed-form anchors, runs each primary Python analysis, and compares the four
Unit 4 exact, joint and sensitivity records with their committed reference
results. `make r-audit` invokes fourteen independent base-R routes for the
prelude, Supplement A and Units 1–9; Rscript needs no extra R package. On
Ubuntu install R with `sudo apt-get install r-base` if needed. The Actions
workflow exercises Python 3.10/3.12 and the base-R route on clean Ubuntu
checkouts, and the Python route on a clean Windows checkout.

## Run the Python checks on Windows

In PowerShell, with Python 3.10 to 3.12 from python.org:

```powershell
git clone https://github.com/nicholaskarlson/bayes.git
cd bayes
py -3.12 -m venv .venv
.venv\Scripts\python -m pip install -r open-materials/bayesian-reasoning-with-statistics/scripts/requirements.txt
.venv\Scripts\python open-materials/bayesian-reasoning-with-statistics/scripts/verify_public.py
```

The repository's `.gitattributes` keeps Git from converting line endings, so
the SHA-256 inventory matches on Windows too. For the base-R route, use
Ubuntu, WSL, or read the workflow's `r-audit` job.

To follow one example, read the data dictionary and model specification, then
run `make unit04 PYTHON=.venv/bin/python`; compare its JSON result to
`expected-results/unit04_exact.json` and inspect
`technical-notes/unit04_deep_dive.md`. To produce the Unit 9 teaching report,
run `make report PYTHON=.venv/bin/python`; its output is in ignored `review/`.
No MCMC sampler is used. The Unit 4 joint check uses seeded independent
simulation, not chains, and the matched likelihood route uses quadrature.

## Contents and scope

All CSVs represent constructed cases; their IDs and measurements do not
describe actual participants. Deterministic generators in `scripts/generate_*`
disclose how several fictional files were made. Their construction parameters
are reproducibility metadata, **not** evidence for an inferential claim. The
`expected-results/` directory contains full-precision numeric reference
values, while `model-specs/` records prior/likelihood choices. The `figures/`
and `technical-notes/` directories help locate calculations used in the book.
The resource's `SHA256SUMS` covers the shipped files under that directory.

These checks establish arithmetic agreement under stated models. They cannot
validate sampling, measurement, exchangeability, causal identification or
generalization. In particular the before/after paired case does not identify
an effect and the route comparison is not a randomized intervention.

## Versions

Releases are Git tags. Cite the tag, or the commit it points to.

- `v1.0.0`: the companion as published with the book (commit `bdfb2ac`).
- `v1.0.1`: the Python check also runs on Windows (forward-slash paths in the
  inventory check; `.gitattributes` turns off line-ending conversion); the clone
  commands use HTTPS. No data, dictionary, model specification, expected result
  or analysis changed. The only `SHA256SUMS` change is the line for
  `scripts/verify_public.py`.

## Rights and attribution

Original companion data, dictionaries, model/expected records, figures,
provenance, notes and this README are CC BY 4.0, with attribution to Nicholas
Elliott Karlson; original code is MIT. See [LICENSES.md](LICENSES.md) for scope.
The book itself has separate rights. Dependencies retain their own licenses.
