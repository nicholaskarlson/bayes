RESOURCE := open-materials/bayesian-reasoning-with-statistics
PYTHON ?= python3
.PHONY: verify r-audit unit04 report review-build

# verify: the inventory and anchors, then the fourteen independent Python checks.
# Each reads the raw files directly and never imports the analysis scripts.
verify:
	$(PYTHON) $(RESOURCE)/scripts/verify_public.py
	$(PYTHON) $(RESOURCE)/scripts/verify_primer_independent.py
	$(PYTHON) $(RESOURCE)/scripts/verify_prelude_independent.py
	$(PYTHON) $(RESOURCE)/scripts/verify_unit01_independent.py
	$(PYTHON) $(RESOURCE)/scripts/verify_unit02_independent.py
	$(PYTHON) $(RESOURCE)/scripts/verify_unit03_independent.py
	$(PYTHON) $(RESOURCE)/scripts/verify_unit04_independent.py
	$(PYTHON) $(RESOURCE)/scripts/verify_unit04_checks_independent.py
	$(PYTHON) $(RESOURCE)/scripts/verify_unit04_joint_independent.py
	$(PYTHON) $(RESOURCE)/scripts/verify_unit04_likelihood_independent.py
	$(PYTHON) $(RESOURCE)/scripts/verify_unit05_independent.py
	$(PYTHON) $(RESOURCE)/scripts/verify_unit06_independent.py
	$(PYTHON) $(RESOURCE)/scripts/verify_unit07_independent.py
	$(PYTHON) $(RESOURCE)/scripts/verify_unit08_independent.py
	$(PYTHON) $(RESOURCE)/scripts/verify_unit09_independent.py

r-audit:
	@Rscript -e 'cat(R.version.string, "\n")'
	Rscript $(RESOURCE)/scripts/verify_primer_in_R.R
	Rscript $(RESOURCE)/scripts/verify_prelude_in_R.R
	Rscript $(RESOURCE)/scripts/verify_unit01_in_R.R
	Rscript $(RESOURCE)/scripts/verify_unit02_in_R.R
	Rscript $(RESOURCE)/scripts/verify_unit03_in_R.R
	Rscript $(RESOURCE)/scripts/verify_unit04_in_R.R $(RESOURCE)
	Rscript $(RESOURCE)/scripts/verify_unit04_checks_in_R.R $(RESOURCE)
	Rscript $(RESOURCE)/scripts/verify_unit04_joint_in_R.R $(RESOURCE)
	Rscript $(RESOURCE)/scripts/verify_unit04_likelihood_in_R.R $(RESOURCE)
	Rscript $(RESOURCE)/scripts/verify_unit05_in_R.R
	Rscript $(RESOURCE)/scripts/verify_unit06_in_R.R
	Rscript $(RESOURCE)/scripts/verify_unit07_in_R.R
	Rscript $(RESOURCE)/scripts/verify_unit08_in_R.R
	Rscript $(RESOURCE)/scripts/verify_unit09_in_R.R

unit04:
	$(PYTHON) $(RESOURCE)/scripts/analyze_unit04.py

report:
	$(PYTHON) $(RESOURCE)/scripts/build_unit09_handoff.py

# The book also names make review-build. It assembles the book's interior from its
# private manuscript source, which is not part of this public companion.
review-build:
	@echo "review-build assembles the book's interior from its private manuscript source."
	@echo "It is an author-side tool and is not part of this public companion; nothing here needs it."
	@echo "Run make verify and make r-audit to check the numbers."
