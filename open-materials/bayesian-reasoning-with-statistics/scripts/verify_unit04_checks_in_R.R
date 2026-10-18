#!/usr/bin/env Rscript
# Independent base-R calculation; reads raw CSV, numeric prior table and scalar anchors.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1L) stop("usage: verify_unit04_checks_in_R.R RESOURCE_DIR")
root <- normalizePath(args[[1]], mustWork=TRUE)
raw <- read.csv(file.path(root,"data","unit04_reflection_attention_paired.csv"))
stopifnot(nrow(raw)==64L, identical(raw$participant_id,sprintf("P%03d",1:64)))
prior_table <- read.csv(file.path(root,"model-specs","unit04_parameters.csv"))
p <- setNames(prior_table$value,prior_table$parameter)
anchors <- read.csv(file.path(root,"expected-results","unit04_checks_anchors.csv"))
expected <- setNames(anchors$value,anchors$field)
d <- raw$after_score-raw$before_score
n <- length(d); avg <- mean(d); variance <- var(d)
m0 <- p[["m0"]]; k0 <- p[["kappa0"]]; a0 <- p[["a0"]]; b0 <- p[["b0"]]
kn <- k0+n; mn <- (k0*m0+n*avg)/kn; an <- a0+n/2
bn <- b0 + (n-1)*variance/2 + k0*n*(avg-m0)^2/(2*kn)
predictive <- function(m,k,a,b) {
  mean_scale <- sqrt((b/a)*(1/k+1/n)); var_scale <- b/a
  c(mean_lower=m+qt(.025,2*a)*mean_scale,
    mean_upper=m+qt(.975,2*a)*mean_scale,
    mean_cdf=pt((avg-m)/mean_scale,2*a),
    var_lower=var_scale*qf(.025,n-1,2*a),
    var_upper=var_scale*qf(.975,n-1,2*a),
    var_cdf=pf(variance/var_scale,n-1,2*a))
}
pr <- predictive(m0,k0,a0,b0); po <- predictive(mn,kn,an,bn)
scale_new <- sqrt(bn/an*(1+1/kn))
prob <- pt((-8-mn)/scale_new,2*an) + pt((8-mn)/scale_new,2*an,lower.tail=FALSE)
omission <- vapply(seq_along(d), function(i) {
  x <- d[-i]; z <- length(x); xbar <- mean(x); ks <- k0+z
  ms <- (k0*m0+z*xbar)/ks; as <- a0+z/2
  bs <- b0+(z-1)*var(x)/2+k0*z*(xbar-m0)^2/(2*ks)
  pt((2-ms)/sqrt(bs/(as*ks)),2*as,lower.tail=FALSE)
},numeric(1))
actual <- c(observed_mean=avg, observed_sample_variance=variance,
  observed_tail_count=sum(abs(d)>8),
  prior_mean_lower=pr[["mean_lower"]],prior_mean_upper=pr[["mean_upper"]],
  posterior_mean_lower=po[["mean_lower"]],posterior_mean_upper=po[["mean_upper"]],
  posterior_var_lower=po[["var_lower"]],posterior_var_upper=po[["var_upper"]],
  posterior_mean_cdf=po[["mean_cdf"]],posterior_var_cdf=po[["var_cdf"]],
  individual_tail_probability=prob,expected_tail_count=n*prob,
  loo_min_id_number=which.min(omission),loo_min_probability=min(omission),
  loo_max_id_number=which.max(omission),loo_max_probability=max(omission))
stopifnot(identical(sort(names(actual)),sort(names(expected))))
for (key in names(actual)) {
  if (!is.finite(actual[[key]]) || abs(actual[[key]]-expected[[key]])>2e-8)
    stop(sprintf("R checks discrepancy %s: %.12f versus %.12f",key,actual[[key]],expected[[key]]))
}
cat(sprintf("R_UNIT04_CHECKS mean95=[%.3f,%.3f] variance95=[%.3f,%.3f] tail=%d expected_tail=%.3f loo=[%.3f,%.3f]\n",
  po[["mean_lower"]],po[["mean_upper"]],po[["var_lower"]],po[["var_upper"]],
  sum(abs(d)>8),n*prob,min(omission),max(omission)))
cat("BAYESIAN_UNIT04_CHECKS_R_OK exact_replications=mean,variance influence=64\n")
