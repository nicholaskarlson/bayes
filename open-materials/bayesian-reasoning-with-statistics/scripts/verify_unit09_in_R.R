#!/usr/bin/env Rscript
# Independent end-to-end route handoff check in base R; no JSON package.
args <- commandArgs(trailingOnly=FALSE)
self <- sub('^--file=', '', args[grepl('^--file=',args)])
if(length(self)!=1L) stop('Run with Rscript.')
root <- normalizePath(file.path(dirname(self),'..'))
d <- read.csv(file.path(root,'data','unit07_help_seeking_routes.csv'),stringsAsFactors=FALSE)
a <- read.csv(file.path(root,'expected-results','unit07_anchors.csv'),stringsAsFactors=FALSE)
stopifnot(nrow(d)==36L,length(unique(d$volunteer_id))==36L,
          sum(d$reported_contact_7d[d$recruitment_route=='desk'])==8L,
          sum(d$reported_contact_7d[d$recruitment_route=='poster'])==0L,
          sum(d$recruitment_route=='desk')==24L,sum(d$recruitment_route=='poster')==12L)
for(i in 1:2) {
  q <- i; A <- q+8; B <- q+16; C <- q; D <- q+12
  row <- a[i,]
  ppositive <- integrate(function(x) dbeta(x,A,B)*pbeta(x,C,D),0,1,
                         subdivisions=500L,rel.tol=1e-11)$value
  mean_difference <- A/(A+B)-C/(C+D)
  prior_zero <- exp(lbeta(q,q+12)-lbeta(q,q))
  repeat_zero <- exp(lbeta(C,D+12)-lbeta(C,D))
  stopifnot(abs(ppositive-row$P_diff_gt_zero)<2e-8,
            abs(mean_difference-row$difference_mean)<1e-12,
            abs(prior_zero-row$prior_zero)<1e-12,
            abs(repeat_zero-row$replicated_zero)<1e-12)
  cat(sprintf('R_UNIT09_HANDOFF %s difference=%.3f P_positive=%.3f prior_zero=%.3f repeat_zero=%.3f\n',
              row$case,mean_difference,ppositive,prior_zero,repeat_zero))
}
cat('BAYESIAN_UNIT09_R_VERIFY_OK raw=yes prior=yes posterior=yes numeric_scope_only=yes no_sampler=yes\n')
