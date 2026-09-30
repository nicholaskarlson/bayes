# Unit 5 technical note: fixed-scale regression and matched descriptive slopes

This note supports, but is not required for, Unit 5 of the book. All 72 records are fictional and cross-sectional. Read `model-specs/unit05_regression.json`, the dictionary and provenance before running the code. `scripts/generate_unit05.py` reproduces fictional data and discloses their simulation choices; its parameter values do not answer the reader's inferential question.

Let `x_i=study_hours_per_week_i−5`, `z_i=long_commute_i`, `X_i=(1,x_i,z_i)`, and `y_i=belonging_score_i`. The common conditional-mean slope `β` is in score points per extra *reported* hour at the same `z`, among comparable modeled volunteers. All ordinary and Bayesian adjusted statements use this target. The unadjusted `y~1+x` slope is deliberately a different description; it aggregates over commute composition. It cannot be called an effect or a confounding correction that is known to suffice.

## Ordinary route

For `A=(1,x)` or `A=X`, compute `b=(A' A)^{-1}A'y`, `s²=||y−Ab||²/(n−p)` and `SE(b_j)=s sqrt([(A'A)^{-1}]_{jj})`. The 95% t interval is `b_beta ± t_.975,n−p SE(b_beta)`. Here the unadjusted slope is 3.27396531748, interval [2.13935611537, 4.40857451959]; the commute-adjusted slope is 2.76553546840, interval [1.57184175971, 3.95922917708], with fitted residual SD 7.91783302142. These are model-based repeated-sampling intervals, not posterior probabilities. The Bayesian model uses a fixed scale of 8 as a working assumption, so its uncertainty is not exactly the OLS t procedure's uncertainty.

## Exact Bayesian route

Given the fixed `σ=8`, `y|b ~ Normal(Xb, σ²I)` and independent Normal priors with means `m0=(55,0,0)` and SDs `(15,2,10)`, write `P0=diag(1/15²,1/2²,1/10²)`. Then `V=(P0+X'X/σ²)^{-1}` and `m=V(P0 m0+X'y/σ²)`. The posterior `b|y ~ Normal(m,V)` is exact under these declared assumptions. Primary `β`: mean 2.55416429294, SD 0.57698926091, central 95% credible interval [1.42328612449, 3.68504246140], `P(β>1)=0.99646549621`. A two-hour contrast at the same `z` is `2β`: posterior mean 5.10832858589 points; this is a comparison of conditional means, not a predicted person's score or a causal response to adding two hours.

A skeptical Normal(0,0.4²) slope prior changes the posterior to mean 0.85209828087, SD 0.33326381379, 95% [0.19891320622, 1.50528335553], and `P(β>1)=0.32859442022`. It is a stress test. External evidence would be needed to defend such strong prior concentration. Keeping the primary prior but changing the fixed likelihood SD to 10 gives mean slope 2.45039237443, 95% [1.07083206681, 3.82995268205] and `P(β>1)=0.98032921331`. No prior may be chosen after these rows merely to produce a desired tail probability.

## Checks and limitations

At `(hours,z)=(5,0)` the prior predictive distribution for a new score has mean 55 and variance `15²+8²`, with probability 0.00466735 outside the 0–100 score bounds. At `(6.5,1)` it has variance `15²+(1.5×2)²+10²+8²`, giving probability 0.01496400 outside the bounds. The normal model has nonzero mass outside the score's possible range. The check exposes that approximation; it does not certify the score instrument.

With the primary posterior mean `m`, define the descriptive discrepancy `T(y)=#{i: |y_i−X_i m|>16}`. The observed file has 3 such rows. For a replicated score at each observed `X_i`, integrating coefficient uncertainty gives predictive variance `σ²+X_i V X_i'`; summing its two-sided tail probabilities yields an expected count of 3.58896223. This is an in-sample fitted-model diagnostic, *not* a calibrated p value or held-out accuracy. The alternate tight prior gives a different fitted mean and observed discrepancy. The 72 rows may have unrecorded programme or tutor dependence, and both predictors and outcome are self-reported. Neither cross-language agreement nor this check fixes those design problems.

The posterior calculation is exact matrix algebra with Normal quantiles, and uses **no MCMC**. The independent standard-library Python verifier implements Gauss–Jordan inversion, NormalDist and a separately integrated t quantile; the base-R audit uses `qr.solve`, Normal and t functions. Both start from the CSV and check the frozen expected record and teaching anchors. They verify arithmetic, not prior credibility, selection or causal identification. Python: `python3 scripts/analyze_unit05.py`; R: `Rscript scripts/analyze_unit05.R` from the resource directory, or use the paths shown in Unit 5 from the repo root.
