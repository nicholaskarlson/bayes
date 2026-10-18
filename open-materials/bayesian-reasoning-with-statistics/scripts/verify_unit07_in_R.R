#!/usr/bin/env Rscript
# Independent base-R route checks full-precision anchors and a finite-sum P(Delta>0).
args <- commandArgs(trailingOnly=FALSE)
self <- sub('^--file=', '', args[grepl('^--file=',args)])
if (length(self)!=1L) stop('Run with Rscript.')
root <- normalizePath(file.path(dirname(self),'..'))
d <- read.csv(file.path(root,'data','unit07_help_seeking_routes.csv'),stringsAsFactors=FALSE)
a <- read.csv(file.path(root,'expected-results','unit07_anchors.csv'),stringsAsFactors=FALSE)
stopifnot(nrow(d)==36L,nrow(a)==2L,identical(a$case,c('primary','sensitivity')),
          identical(d$volunteer_id,c(sprintf('D%03d',1:24),sprintf('P%03d',1:12))),
          identical(d$recruitment_route,rep(c('desk','poster'),c(24,12))),
          !anyNA(d),all(d$reported_contact_7d %in% c(0,1)),
          sum(d$reported_contact_7d[1:24])==8,sum(d$reported_contact_7d[25:36])==0)
near <- function(actual,expected,tol=2e-7,label='value') {
  if(any(!is.finite(actual)) || length(actual)!=length(expected) || any(abs(actual-expected)>tol))
    stop(label,' mismatch: actual=',paste(actual,collapse=','),' expected=',paste(expected,collapse=','))
}
for (i in 1:2) {
  prior <- i
  A <- prior+8;B <- prior+24-8;C <- prior;D <- prior+12
  row <- a[i,]
  near(as.numeric(row[c('a_desk','b_desk','a_poster','b_poster')]),c(A,B,C,D),0,'shapes')
  m <- c(A/(A+B),C/(C+D))
  near(c(row$desk_mean,row$poster_mean,row$difference_mean),c(m,m[1]-m[2]),label='mean')
  n <- C+D-1L
  # Integer-Beta CDF is a finite binomial upper tail. Integrate termwise via beta ratios.
  p0 <- sum(vapply(C:n,function(k) exp(lchoose(n,k)+lbeta(A+k,B+n-k)-lbeta(A,B)),numeric(1)))
  near(p0,row$P_diff_gt_zero,label='finite-sum P(Delta>0)')
  tdiff <- function(h) integrate(function(x) dbeta(x,A,B)*pbeta(x-h,C,D),
                                max(0,h),1,subdivisions=500L,rel.tol=1e-11)$value
  near(tdiff(.1),row$P_diff_gt_0_10,label='ten-point tail')
  for (j in c(.025,.975)) {
    target <- if(j==.025) row$difference_lower else row$difference_upper
    val <- integrate(function(x) dbeta(x,A,B)*pbeta(x-target,C,D,lower.tail=FALSE),
                     0,1,subdivisions=500L,rel.tol=1e-11)$value
    near(val,j,3e-7,'difference quantile')
  }
  near(exp(lbeta(prior,prior+12)-lbeta(prior,prior)),row$prior_zero,label='prior predictive zero')
  near(exp(lbeta(C,D+12)-lbeta(C,D)),row$replicated_zero,label='posterior predictive zero')
}
out <- suppressWarnings(system2('Rscript',file.path(root,'scripts','analyze_unit07.R'),stdout=TRUE))
if(length(attr(out,'status')) && attr(out,'status')!=0L) stop('Unit 7 R analysis failed.')
out <- trimws(out);out <- out[nzchar(out)]
stopifnot(length(out)==8L,identical(out[1],'n=36 desk=8/24 poster=0/12 raw_difference=0.333'),
          identical(out[8],'scope=fictional convenience volunteers; route association only; no causal effect or campus rate; exact no MCMC'))
for(i in 1:2) {
  row <- a[i,];name <- as.character(row$case);j <- 2+(i-1)*3
  stopifnot(identical(out[j],sprintf('%s prior=Beta(%d,%d) desk=Beta(%d,%d) poster=Beta(%d,%d)',
                  name,i,i,row$a_desk,row$b_desk,row$a_poster,row$b_poster)))
  expected <- sprintf('%s means=[%.3f,%.3f] difference=%.3f difference_95=[%.3f,%.3f] P_diff_gt_zero=%.3f P_diff_gt_0.10=%.3f',
                  name,row$desk_mean,row$poster_mean,row$difference_mean,row$difference_lower,row$difference_upper,row$P_diff_gt_zero,row$P_diff_gt_0_10)
  stopifnot(identical(out[j+1],expected))
  stopifnot(identical(out[j+2],sprintf('%s prior_zero=%.3f replicated_zero=%.3f',name,row$prior_zero,row$replicated_zero)))
}
cat('BAYESIAN_UNIT07_R_VERIFY_OK n=36 finite_sum_vs_integral=yes exact=yes\n')
