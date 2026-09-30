#!/usr/bin/env Rscript
# Base-R path reads the raw CSV and numerical parameter table; no Python analyzer.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("usage: verify_unit04_in_R.R RESOURCE_DIR")
root <- normalizePath(args[[1]], mustWork = TRUE)
data <- read.csv(file.path(root, "data", "unit04_reflection_attention_paired.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(identical(names(data), c("participant_id", "before_score", "after_score")))
stopifnot(nrow(data) == 64L, identical(data$participant_id, sprintf("P%03d", 1:64)))
stopifnot(!anyNA(data), all(data$before_score >= 0), all(data$after_score <= 100))
params_table <- read.csv(file.path(root, "model-specs", "unit04_parameters.csv"),
                         stringsAsFactors = FALSE)
p <- setNames(params_table$value, params_table$parameter)
expected <- read.csv(file.path(root, "expected-results", "unit04_anchors.csv"),
                     stringsAsFactors = FALSE)
target <- setNames(expected$value, expected$field)
d <- data$after_score - data$before_score
n <- length(d)
sample_mean <- mean(d)
sample_sd <- sd(d)
sample_se <- sample_sd / sqrt(n)
paired <- t.test(d, mu = 0, conf.level = .95)
stopifnot(abs(as.numeric(paired$statistic) - sample_mean/sample_se) < 1e-12)
k0 <- p[["kappa0"]]; m0 <- p[["m0"]]; a0 <- p[["a0"]]; b0 <- p[["b0"]]
kn <- k0 + n
mn <- (k0*m0 + n*sample_mean)/kn
an <- a0+n/2
bn <- b0 + sum((d-sample_mean)^2)/2 + k0*n*(sample_mean-m0)^2/(2*kn)
post_scale <- sqrt(bn/(an*kn))
limits <- mn + qt(c(.025,.975), df=2*an)*post_scale
pr_above_two <- pt((p[["meaningful_change"]]-mn)/post_scale,
                    df=2*an, lower.tail=FALSE)
actual <- c(n=n, mean_difference=sample_mean, t=as.numeric(paired$statistic),
            ci95_lower=paired$conf.int[[1]], ci95_upper=paired$conf.int[[2]],
            m_n=mn, posterior_ci95_lower=limits[[1]], posterior_ci95_upper=limits[[2]],
            p_mu_gt_two=pr_above_two)
stopifnot(identical(sort(names(target)), sort(names(actual))))
for (name in names(actual)) {
  if (!is.finite(actual[[name]]) || abs(actual[[name]]-target[[name]]) > 2e-8)
    stop(sprintf("R verification differs for %s: actual %.12f expected %.12f",
                 name, actual[[name]], target[[name]]))
}
cat(sprintf("R_UNIT04_RESULT mean=%.3f t=%.3f ci95=[%.3f,%.3f] posterior_location=%.3f posterior_ci95=[%.3f,%.3f] p_mu_gt_two=%.3f\n",
            sample_mean, as.numeric(paired$statistic), paired$conf.int[[1]],
            paired$conf.int[[2]], mn, limits[[1]], limits[[2]], pr_above_two))
cat("BAYESIAN_UNIT04_R_VERIFY_OK n=64 exact_posterior=yes frequentist=yes\n")
