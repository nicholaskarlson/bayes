# Unit 8 diagnostic audit (optional)

This note reproduces selected previously audited Unit 4 and Unit 7 quantities for the reader chapter. Neither case has a sampler, and no diagnostic below is an MCMC convergence statistic. Files remain fictional. See `model-specs/unit08_diagnostic_plan.json`, `expected-results/unit08_diagnostics.json`, and the original Unit 4 and Unit 7 technical notes.

## Zero count and a prior predictive warning

For θ~Beta(a,b) and a prospective twelve-person count Y|θ~Binomial(12,θ), P(Y=0)=B(a,b+12)/B(a,b). For Beta(2,2), this is **1/35=0.0285714285714**; for Beta(1,1), it is **1/13=0.0769230769231**. Unit 7 observed zero in its poster route. The check asks whether a chosen pre-data model would often produce that feature; it does not assign an ordinary null-hypothesis p value to the observed data, and the two priors were not selected as posterior model probabilities. Post-data changes prompted by zero must be described as sensitivity, justified externally, and checked for other consequences. A low probability suggests an audit of recruitment, coding, independence, outcome measurement and prior information. It does not by itself identify which component failed.

Unit 7's exact independent-Beta posterior yields P(θ_D−θ_P>0)=0.987467306015693. Because the two continuous posteriors assign zero probability to equality, P(θ_P−θ_D>0)=1−0.987467306015693=0.012532693984307. Reporting 0.987 for the reversed target is a deliberately planted transcription error in the chapter exercise, not a real result. The independent audit explicitly tests both directions against the declared target.

## Paired observation and likelihood

For Unit 4, let d_i=after_i−before_i. In this *matched comparison*, μ~Normal(0,25). For the Normal likelihood, d_i|μ~Normal(μ,16), the posterior variance is v=(1/25+n/16)^−1 and the posterior mean is m=v Σd_i/16. Thus P(μ>2)=1−Φ((2−m)/√v). The R and independent Python audits recompute this from the 64 CSV rows, then after omitting P059 and P042 separately. A second likelihood uses Student-t with four degrees of freedom and scale √8, which also has conditional variance 16. Its posterior probabilities come from the existing independently checked deterministic quadrature in `unit04_likelihood_reference.json` and `unit04_likelihood_anchors.csv`; Unit 8 reads those audited results and does not implement a second integrator.

| Records | Normal fixed P(μ>2) | Student-t4 matched P(μ>2) | Original NIG P(μ>2) |
| --- | ---: | ---: | ---: |
| All 64 | 0.636412 | 0.750264 | 0.552522 |
| Omit P059 | 0.530226 | 0.687279 | 0.448079 |
| Omit P042 | 0.731709 | 0.805715 | 0.653158 |

The original Normal-inverse-gamma model estimates variance and has a different prior. Its column is an influence illustration **within that model**; differences between this column and the matched Normal column are not isolated likelihood effects. Omitting a row reveals sensitivity; it is not license to erase a valid observation. Inspect the P059 record and measurement source; retain it unless a documented data error warrants correction. The illustrative 0.60 cutoff switches for the matched Normal model when P059 is omitted; its choice is not validated by consequences. The Student-t likelihood changes the tail under a shared prior and conditional variance, but these are post hoc alternatives, not posterior odds for model truth. A robust claim should report the conflict and avoid a binary action built on this toy line.

## Computation triage and claim boundary

Exact conjugate updates, deterministic quadrature and independent numeric audits do not provide sampling-chain R-hat, effective sample size or MCSE. If a future model genuinely uses MCMC, declare chains, warmup, seed, target and tolerances in its own plan; inspect divergences, traces, mixing, effective sample size, R-hat and MCSE relative to claimed precision. A warning is a stop-and-diagnose signal, not a reason to hide a chain or increase decimal places. A converged sampler still cannot repair poor measurement, self-selection, dependence, or a causal design ceiling. These prior and posterior checks are in-sample criticism, not held-out calibration or generalization; the latter belong to the later *Predictive Reasoning* book.
