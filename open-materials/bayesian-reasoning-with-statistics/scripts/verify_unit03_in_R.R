#!/usr/bin/env Rscript
# Independent base-R audit, from raw CSV and analytic beta/Wilson functions.
args <- commandArgs(trailingOnly=FALSE)
self <- sub('^--file=', '', args[grepl('^--file=',args)])
if (length(self)!=1L) stop('Run with Rscript.')
root <- normalizePath(file.path(dirname(self),'..'))
d <- read.csv(file.path(root,'data','unit02_help_seeking.csv'),stringsAsFactors=FALSE)
stopifnot(identical(names(d),c('participant_id','help_sought_within_7_days')),
 identical(d$participant_id,sprintf('H%03d',1:12)),nrow(d)==12L,!anyNA(d),
 all(d$help_sought_within_7_days %in% c(0L,1L)),sum(d$help_sought_within_7_days)==8L)
near <- function(x,y,tol=1e-9) stopifnot(is.finite(x),abs(x-y)<tol)
p <- sum(d$help_sought_within_7_days)/nrow(d);n <- nrow(d);z <- qnorm(.975)
center <- (p+z*z/(2*n))/(1+z*z/n)
spread <- z*sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n)
near(center-spread,.39062208887279953)
near(center+spread,.8618799089087867)
stopifnot(max(abs(qbeta(c(.025,.975),10,6)-c(.38380373254115413,.8366356761404865)))<1e-8)
for (case in list(c(10,6,.5967844495851521),c(9,5,.646958153728),c(16,12,.38726771570698787))) {
 mass <- pbeta(.6,case[1],case[2],lower.tail=FALSE)
 near(mass,case[3]);near(2*mass-1,2*case[3]-1)
}
out <- suppressWarnings(system2('Rscript',file.path(root,'scripts','analyze_unit03.R'),stdout=TRUE))
if (length(attr(out,'status')) && attr(out,'status')!=0L) stop('Unit 3 R route failed.')
needed <- c('n=12 y=8 observed=0.667 target=reported_contact_probability',
 'wilson_95=[0.391,0.862] procedure=repeated_sample_coverage',
 'posterior=Beta(10,6) mean=0.625 credible_95=[0.384,0.837]',
 'Beta(2,2) P(theta>0.6)=0.597 toy_expected_utility=0.194 toy_choice=action',
 'Beta(1,1) P(theta>0.6)=0.647 toy_expected_utility=0.294 toy_choice=action',
 'Beta(8,8) P(theta>0.6)=0.387 toy_expected_utility=-0.225 toy_choice=no_action',
 'scope=fictional conditional model and hypothetical utility; no orientation effect or real service recommendation')
stopifnot(setequal(out,needed))
cat('BAYESIAN_UNIT03_R_VERIFY_OK n=12 y=8 wilson=yes posterior=yes decision_sensitivity=yes\n')
