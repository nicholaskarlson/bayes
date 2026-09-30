#!/usr/bin/env Rscript
# Independent base-R adaptive integration under a common fixed-scale prior/likelihood comparison.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=1L) stop('usage: verify_unit04_likelihood_in_R.R RESOURCE_DIR')
root <- normalizePath(args[[1]],mustWork=TRUE)
raw <- read.csv(file.path(root,'data','unit04_reflection_attention_paired.csv'))
stopifnot(nrow(raw)==64L,identical(raw$participant_id,sprintf('P%03d',1:64)))
model <- read.csv(file.path(root,'model-specs','unit04_likelihood_parameters.csv'))
p <- setNames(model$value,model$parameter)
anchors <- read.csv(file.path(root,'expected-results','unit04_likelihood_anchors.csv'))
stopifnot(nrow(anchors)==24L,length(unique(paste(anchors$case,anchors$model,anchors$field)))==24L)
stopifnot(identical(sort(unique(anchors$case)),sort(c('all','omit_P042','omit_P059'))))
stopifnot(p[['mu_prior_mean_points']]==0,p[['mu_prior_sd_points']]==5,
          p[['normal_sd_points']]==4,abs(p[['student_t_scale_points']]^2*2-16)<1e-12,p[['student_degrees_of_freedom']]==4,
          p[['two_point_threshold']]==2,p[['integration_lower']]==-20,p[['integration_upper']]==20,
          p[['enlarged_lower']]==-30,p[['enlarged_upper']]==30)
d <- raw$after_score-raw$before_score

anchor <- function(case,model,field) {
  rows <- anchors[anchors$case==case & anchors$model==model & anchors$field==field,]
  if (nrow(rows)!=1L) stop('missing/duplicate anchor')
  rows$value[[1]]
}
check <- function(case,model,field,value,tol) {
  expected <- anchor(case,model,field)
  if (!is.finite(value)||abs(value-expected)>tol)
    stop(sprintf('R likelihood mismatch %s %s %s got %.10f expected %.10f',
                 case,model,field,value,expected))
}

normal_exact <- function(values,case) {
  n <- length(values)
  v <- 1/(1/25+n/16)
  m <- v*sum(values)/16
  sd <- sqrt(v)
  check(case,'normal_fixed','posterior_mean',m,2e-10)
  check(case,'normal_fixed','ci95_lower',qnorm(.025,m,sd),2e-10)
  check(case,'normal_fixed','ci95_upper',qnorm(.975,m,sd),2e-10)
  check(case,'normal_fixed','p_mu_gt_2',pnorm(2,m,sd,lower.tail=FALSE),2e-10)
  invisible(m)
}

student_integrate <- function(values,case) {
  df <- p[['student_degrees_of_freedom']]
  s <- p[['student_t_scale_points']]
  prior_sd <- p[['mu_prior_sd_points']]
  student_const <- lgamma((df+1)/2)-lgamma(df/2)-.5*log(df*pi)-log(s)
  log_kernel_one <- function(mu) {
    -.5*(mu/prior_sd)^2+sum(student_const-(df+1)/2*log1p(((values-mu)/s)^2/df))
  }
  mode <- optimize(function(mu) -log_kernel_one(mu),c(-10,15),tol=1e-10)$minimum
  peak <- log_kernel_one(mode)
  weight <- function(x) vapply(x,function(mu) exp(log_kernel_one(mu)-peak),numeric(1))
  integral <- function(fun,a,b) integrate(fun,a,b,rel.tol=1e-10,
                                           abs.tol=1e-12,subdivisions=500L)$value
  z <- integral(weight,-20,20)
  if (!(z>0)) stop('zero normalizer')
  m <- integral(function(x) x*weight(x),-20,20)/z
  tail <- integral(weight,2,20)/z
  quantile <- function(prob) uniroot(function(x) integral(weight,-20,x)/z-prob,
                                      interval=c(-10,15),tol=1e-8)$root
  low <- quantile(.025);high <- quantile(.975)
  z_large <- integral(weight,-30,30)
  tail_large <- integral(weight,2,30)/z_large
  if (abs(tail_large-tail)>1e-8||(z_large-z)/z_large>1e-8) stop('integration interval sensitivity')
  check(case,'student_t4_equal_variance','posterior_mean',m,3e-5)
  check(case,'student_t4_equal_variance','ci95_lower',low,3e-5)
  check(case,'student_t4_equal_variance','ci95_upper',high,3e-5)
  check(case,'student_t4_equal_variance','p_mu_gt_2',tail,3e-5)
  invisible(tail)
}

for (case in c('all','omit_P059','omit_P042')) {
  values <- switch(case,all=d,omit_P059=d[raw$participant_id!='P059'],
                   omit_P042=d[raw$participant_id!='P042'])
  stopifnot(length(values)==if (case=='all') 64L else 63L)
  normal_exact(values,case)
  student_integrate(values,case)
}
cat(sprintf('R_UNIT04_LIKELIHOOD normal_P_gt_two=%.3f student_t4_P_gt_two=%.3f\n',
            anchor('all','normal_fixed','p_mu_gt_2'),anchor('all','student_t4_equal_variance','p_mu_gt_2')))
cat('BAYESIAN_UNIT04_LIKELIHOOD_R_OK cases=3 models=2 independent_adaptive_integration=yes mcmc=no\n')
