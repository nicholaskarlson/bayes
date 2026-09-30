# Unit 6 technical audit — exact Normal partial pooling

This optional note derives the reported calculation. It is **not** required to read the student chapter. Input is the unchanged fictional 80-row CSV, eight known group sizes `n=(4,5,6,8,10,12,15,20)`, group averages and within-group SSE. No sampler, group effect, or population representative estimate is asserted. Generator centers in `scripts/generate_unit06.py` are author provenance, not inferred truths.

## Model and the exact update

For group `j` and volunteer `i`, `y_ij | theta_j ~ Normal(theta_j, sigma²)` with **sigma=8 fixed**. Conditional on an unknown common center `mu`, the eight `theta_j` are independent `Normal(mu, tau²)` with **tau=5 fixed**. `mu ~ Normal(55, 10²)`. All scales are in score points and the group labels are assumed exchangeable under this prior. The model is unbounded even though the reported scale is 0–100.

For group mean `a_j` from `n_j` rows, marginalizing `theta_j` gives `a_j | mu ~ Normal(mu, tau²+sigma²/n_j)`. Put `v_j=tau²+sigma²/n_j`. Conjugacy then gives

`V_mu = [1/10² + sum_j 1/v_j]^-1` and
`M_mu = V_mu [55/10² + sum_j a_j/v_j]`.

For tau=5 this gives posterior mean **55.032** and SD about **1.986** for `mu`. Conditional on `mu` and the data, the group mean has Normal variance `c_j=[1/tau²+n_j/sigma²]^-1` and conditional mean `w_j a_j+(1-w_j)mu`, where `w_j=n_j tau²/(sigma²+n_j tau²)`. Integrating uncertainty in `mu` produces posterior group mean `w_j a_j+(1-w_j)M_mu`, variance `c_j+(1-w_j)² V_mu`, and between-group covariance `(1-w_j)(1-w_k)V_mu`. **Dropping that covariance makes the W01–W08 comparison too uncertain or otherwise wrong.** Its contrast variance is `Var(theta_1)+Var(theta_8)-2Cov(theta_1,theta_8)`. A Normal tail integral yields the comparison probabilities 0.653 at tau=5, 0.401 at tau=2 and 0.783 at tau=10. These are probabilities of **latent score means**, not intervention or predictive-choice probabilities.

The primary Python route constructs a nine-parameter precision matrix for `(mu,theta_1,...,theta_8)` and solves its joint posterior. The independent Python audit first integrates each `theta_j`, then uses the scalar expressions above. The base-R entry point uses the scalar expressions, and the independent base-R audit solves the joint matrix system. `expected-results/unit06_partial_pooling.json` holds full precision reference values, with a CSV of selected anchors suitable for base R. Both routes use the same declared data and model but different algebra; this checks implementation arithmetic, not assumptions.

## Familiar descriptions and predictive implications

The **student-weighted** raw pooled mean 55.796 is not `mu`, the group-level exchangeable center. The individual raw group means answer another descriptive question again. W01's n=4 raw 63.625 receives self weight 0.610; W08's n=20 raw 59.345 receives 0.887. Credible intervals integrate the uncertainty in the learned common center. The model fixes tau and sigma instead of estimating them, which understates uncertainty if their selection is doubtful; the three tau cases demonstrate material sensitivity without constituting a posterior for tau.

Before observing any scores, a new group of four has mean `Normal(55, 10²+tau²+sigma²/4)`, central 95% range [31.727,78.273]. An individual prior predictive score has variance `10²+tau²+sigma²`, giving P(score outside 0–100) about 0.0006. Those bounds are judgments about a hypothetical exchangeable setting, not an empirical check of an instrument.

For a replicated cohort at the same group sizes, center *each replicated workshop* on its own replicated sample mean. Conditional on any values of the latent means, the expected within-group SSE is `(80−8) sigma² = 4608`; the observed SSE is 4540.640, corresponding to within-group SD 7.941 on 72 degrees of freedom. The expectation remains 4608 after integrating the posterior. This check is exact, in-sample, and narrow. It neither checks between-workshop variation nor establishes normal tails, independent students, reliable scores, a representative sampling frame, or a workshop effect.

For further conceptual treatment of pooling, see the [Stan education multilevel tutorial](https://mc-stan.org/learn-stan/case-studies/tutorial_rstanarm.html). The chapter's fixed-scale conjugate teaching model is simpler and does **not** claim to have run that tutorial or fitted Stan. A real analysis may need estimated group variation, predictors, dependence, bounded/robust outcomes and computational diagnostics if a sampler is used.
