#!/usr/bin/env Rscript
# Independent base-R direct conjugate simulation; no Python, JSON library or MCMC.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=1L) stop("usage: verify_unit04_joint_in_R.R RESOURCE_DIR")
root <- normalizePath(args[[1]],mustWork=TRUE)
raw <- read.csv(file.path(root,"data","unit04_reflection_attention_paired.csv"),stringsAsFactors=FALSE)
stopifnot(nrow(raw)==64L,identical(raw$participant_id,sprintf("P%03d",1:64)))
prior_table <- read.csv(file.path(root,"model-specs","unit04_parameters.csv"),stringsAsFactors=FALSE)
p <- setNames(prior_table$value,prior_table$parameter)
plan_table <- read.csv(file.path(root,"model-specs","unit04_joint_parameters.csv"),stringsAsFactors=FALSE)
plan <- setNames(plan_table$value,plan_table$parameter)
reference_table <- read.csv(file.path(root,"expected-results","unit04_joint_anchors.csv"),stringsAsFactors=FALSE)
reference <- setNames(reference_table$value,reference_table$field)
exact_table <- read.csv(file.path(root,"expected-results","unit04_checks_anchors.csv"),stringsAsFactors=FALSE)
exact <- setNames(exact_table$value,exact_table$field)
stopifnot(plan[["replicates_per_distribution"]]==30000,plan[["absolute_tail_cutoff_points"]]==8,
          plan[["large_change_cutoff_points"]]==15)
set.seed(as.integer(plan[["r_seed"]]))
d <- raw$after_score-raw$before_score
n <- length(d);avg <- mean(d);ss <- (n-1)*var(d)
observed_tail <- sum(abs(d)>8);observed_negative <- sum(d<0)
centered <- d-avg
observed_skew <- mean(centered^3)/mean(centered^2)^1.5
m0 <- p[["m0"]];k0 <- p[["kappa0"]];a0 <- p[["a0"]];b0 <- p[["b0"]]
kn <- k0+n;mn <- (k0*m0+n*avg)/kn;an <- a0+n/2
bn <- b0+ss/2+k0*n*(avg-m0)^2/(2*kn)
reps <- as.integer(plan[["replicates_per_distribution"]])

check_mean <- function(label,key,values,control=NULL) {
  current <- mean(values);se <- sd(values)/sqrt(length(values))
  target <- reference[[paste0(label,"_",key)]]
  target_se <- reference[[paste0(label,"_se_",key)]]
  if (is.null(target)||is.null(target_se) ||
      !is.finite(current)||abs(current-target)>6*sqrt(se^2+target_se^2)+0.05)
    stop(paste("R joint mean mismatch",label,key,current,target))
  if (!is.null(control) && abs(current-control)>6*se+0.05)
    stop(paste("R joint analytic control mismatch",label,key,current,control))
  invisible(current)
}

check_probability <- function(label,key,condition) {
  current <- mean(condition);se <- sqrt(current*(1-current)/length(condition))
  target <- reference[[paste0(label,"_",key)]]
  target_se <- reference[[paste0(label,"_",key,"_mcse")]]
  if (is.null(target)||is.null(target_se)||!is.finite(current)||
      abs(current-target)>6*sqrt(se^2+target_se^2)+0.01)
    stop(paste("R joint probability mismatch",label,key,current,target))
  current
}

simulate <- function(label,m,k,a,b) {
  # rgamma(shape=a, rate=b) is precision 1/v when v ~ IG(shape=a, scale=b).
  v <- 1/rgamma(reps,shape=a,rate=b)
  mu <- rnorm(reps,mean=m,sd=sqrt(v/k))
  y <- matrix(rnorm(reps*n),nrow=reps,ncol=n)
  y <- y*sqrt(v)+mu
  cohort_mean <- rowMeans(y)
  residual <- y-cohort_mean
  sample_variance <- rowSums(residual^2)/(n-1)
  tail_count <- rowSums(abs(y)>8)
  negative_count <- rowSums(y<0)
  skew <- rowMeans(residual^3)/rowMeans(residual^2)^1.5
  large <- rowSums(abs(y)>15)>0
  check_mean(label,"mean_of_replicated_means",cohort_mean,m)
  check_mean(label,"mean_of_replicated_sample_variances",sample_variance,b/(a-1))
  check_mean(label,"mean_tail_count",tail_count)
  check_mean(label,"mean_negative_count",negative_count)
  check_probability(label,"p_tail_count_ge_observed",tail_count>=observed_tail)
  check_probability(label,"p_abs_skew_ge_observed",abs(skew)>=abs(observed_skew))
  check_probability(label,"p_negative_count_ge_observed",negative_count>=observed_negative)
  check_probability(label,"p_any_abs_change_gt_15",large)
  c(tail_probability=mean(tail_count>=observed_tail),
    skew_probability=mean(abs(skew)>=abs(observed_skew)),
    expected_tail=mean(tail_count))
}

prior <- simulate("prior",m0,k0,a0,b0)
posterior <- simulate("posterior",mn,kn,an,bn)
if (abs(posterior[["expected_tail"]]-exact[["expected_tail_count"]])>
    6*reference[["posterior_se_mean_tail_count"]]+0.15)
  stop("posterior joint tail expectation disagrees with exact marginal expectation")
cat(sprintf("R_UNIT04_JOINT prior_tail_ge=%0.3f posterior_tail_ge=%0.3f posterior_skew_abs_ge=%0.3f\n",
            prior[["tail_probability"]],posterior[["tail_probability"]],posterior[["skew_probability"]]))
cat("BAYESIAN_UNIT04_JOINT_R_OK cohorts=30000+30000 direct_conjugate=yes mcmc=no\n")
