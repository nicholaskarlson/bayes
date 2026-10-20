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
results. It then runs fourteen independent Python checks,
`scripts/verify_*_independent.py`, one beside each base-R check. Each reads the
raw files itself, uses exact fractions where the book says so, and never
imports the analysis scripts. `make r-audit` invokes fourteen independent base-R routes for the
prelude, Supplement A and Units 1–9; Rscript needs no extra R package, and
any R from 3.6 on works. `make r-audit` prints the R version first, so a
report shows which R ran the checks. On Ubuntu install R with
`sudo apt-get install r-base` if needed. The Actions workflow exercises
Python 3.10/3.12 and the base-R route on clean Ubuntu checkouts (Ubuntu's
current R, and R 3.6.3 on Ubuntu 20.04), and both Python routes on a clean
Windows checkout.

## Run the Python checks on Windows

In PowerShell, with Python 3.10 to 3.12 from python.org:

```powershell
git clone https://github.com/nicholaskarlson/bayes.git
cd bayes
py -3.12 -m venv .venv
.venv\Scripts\python -m pip install -r open-materials/bayesian-reasoning-with-statistics/scripts/requirements.txt
.venv\Scripts\python open-materials/bayesian-reasoning-with-statistics/scripts/verify_public.py
Get-ChildItem open-materials/bayesian-reasoning-with-statistics/scripts/verify_*_independent.py | ForEach-Object { .venv\Scripts\python $_.FullName; if ($LASTEXITCODE) { throw "failed: $($_.Name)" } }
```

The last line runs the fourteen independent Python checks that `make verify`
runs on Ubuntu.

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

## Names in the book that are not in this repository

The book was written alongside a private authoring repository, and a few things
it names belong there:

- `make review-build`, `MANUSCRIPT_BUILD.md` and the assembled-PDF verifier
  build and check the book's printed interior from its manuscript. Here,
  `make review-build` only prints this explanation.
- The private editorial suite that Unit 4 mentions also checks the chapters'
  wording: each unit's six problems and solutions, required phrases, and the
  printed tables. Those checks need the book's text, so they stay private. The
  public `verify_*_independent.py` scripts make every numeric check that the
  private versions make.

Everything else the book names is here. The Unit 8 and Unit 9 handoffs run
`scripts/verify_unit08_independent.py` and `scripts/verify_unit09_independent.py`.
Supplement A's independent Python grid and Primer checks are
`scripts/verify_prelude_independent.py` and `scripts/verify_primer_independent.py`.
The topic crosswalk in the book's closing section is
`technical-notes/topic_crosswalk.md`.

## Versions

Releases are Git tags. Cite the tag, or the commit it points to.

- `v1.0.0`: the companion as published with the book (commit `bdfb2ac`).
- `v1.0.1`: the Python check also runs on Windows. Changes:
  - the inventory check uses forward-slash paths;
  - `.gitattributes` turns off line-ending conversion;
  - the Unit 4 quadrature centering point (`quadrature.mode`) is compared at
    1e-6 instead of 1e-8. In double precision a maximizer can locate it only to
    about 1e-7, and it is a numerical diagnostic, not a reported result. Every
    reported value is still compared at 1e-8;
  - a failed comparison now prints the observed and expected values;
  - the clone commands use HTTPS.

  No data, dictionary, model specification, expected result or analysis changed.
  The only `SHA256SUMS` change is the line for `scripts/verify_public.py`.
- `v1.0.2`: the base-R checks also run on R before 4.0. Changes:
  - five R checks now pass `stringsAsFactors=FALSE` to `read.csv()`, as the
    other R scripts already did: `verify_primer_in_R.R`,
    `verify_prelude_in_R.R`, `verify_unit04_checks_in_R.R`,
    `verify_unit04_joint_in_R.R` and `verify_unit04_likelihood_in_R.R`.
    Before R 4.0, `read.csv()` turns text columns into factors by default, so
    on an older R the participant-ID check in `verify_unit04_checks_in_R.R`
    stopped with
    `identical(raw$participant_id, sprintf("P%03d", 1:64)) is not TRUE`;
  - `make r-audit` prints the R version before the checks;
  - the Actions workflow also runs `make r-audit` on R 3.6.3 (Ubuntu 20.04).

  No data, dictionary, model specification, expected result, analysis or
  reported number changed. The only `SHA256SUMS` changes are the lines for
  those five R checks.
- `v1.1.0`: the independent Python checks the book names are public. Changes:
  - fourteen `scripts/verify_*_independent.py` checks, one beside each base-R
    check. They include `verify_unit08_independent.py` and
    `verify_unit09_independent.py` from the Unit 8 and Unit 9 handoffs, and the
    Prelude grid and Primer checks that Supplement A describes. They are the
    author's own verifiers with their chapter-text assertions removed. Their
    numeric checks are unchanged, except that every file is now read as UTF-8,
    so they also run on Windows;
  - `make verify` runs them after `verify_public.py`, and the Windows job in
    the Actions workflow runs them too;
  - `technical-notes/topic_crosswalk.md`, the topic crosswalk that the book's
    closing section names;
  - the three Unit 4 figure scripts now write LF line endings on every system,
    so the figure byte checks hold on Windows. The figures are unchanged;
  - `make review-build` explains that it is an author-side tool, and a new
    README section lists the other author-side names in the book.

  No data, dictionary, model specification, expected result, figure or
  reported number changed. `SHA256SUMS` gains fifteen lines and changes three
  (the figure scripts), so `verify_public.py` now reports `files=130`.

## Rights and attribution

Original companion data, dictionaries, model/expected records, figures,
provenance, notes and this README are CC BY 4.0, with attribution to Nicholas
Elliott Karlson; original code is MIT. See [LICENSES.md](LICENSES.md) for scope.
The book itself has separate rights. Dependencies retain their own licenses.
