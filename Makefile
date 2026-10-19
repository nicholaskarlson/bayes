RESOURCE := open-materials/bayesian-reasoning-with-statistics
PYTHON ?= python3
.PHONY: verify r-audit unit04 report

verify:
	$(PYTHON) $(RESOURCE)/scripts/verify_public.py

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
