#!/usr/bin/env Rscript
# Base-R entry point for the same exact beta-binomial model.
args <- commandArgs(trailingOnly=FALSE)
self <- sub('^--file=', '', args[grepl('^--file=',args)])
if (length(self)!=1L) stop('Run with Rscript.')
resource <- normalizePath(file.path(dirname(self),'..'))
spec <- readLines(file.path(resource,'model-specs','unit02_model.json'),warn=FALSE)
if (!any(grepl('"spec_id": "unit02-beta-binomial-v1"',spec,fixed=TRUE))) stop('Wrong model spec.')
d <- read.csv(file.path(resource,'data','unit02_help_seeking.csv'),stringsAsFactors=FALSE)
if (nrow(d)!=12L || !identical(d$participant_id,sprintf('H%03d',1:12)) ||
    anyNA(d) || !all(d$help_sought_within_7_days %in% c(0L,1L))) stop('Invalid data.')
n <- nrow(d); y <- sum(d$help_sought_within_7_days)
if (y!=8L) stop('Unexpected count.')
bb <- function(k,n,a,b) choose(n,k)*exp(lbeta(a+k,b+n-k)-lbeta(a,b))
cat(sprintf('n=%d y=%d observed=%.3f prior=Beta(2,2) posterior=Beta(10,6)\n',n,y,y/n))
cat(sprintf('prior_P(Y>=10_of_12)=%.3f prior_P(outer)=%.3f\n',
            sum(bb(10:12,n,2,2)),sum(bb(c(0:2,10:12),n,2,2))))
interval <- qbeta(c(.025,.975),10,6)
cat(sprintf('posterior_mean=%.3f equal_tail_95=[%.3f,%.3f] P(theta>0.5)=%.3f\n',
            10/16,interval[1],interval[2],pbeta(.5,10,6,lower.tail=FALSE)))
cat(sprintf('posterior_predictive_P(Yrep>=8_of_12)=%.3f\n',sum(bb(8:12,n,10,6))))
for (prior in list(c(1,1),c(8,8))) {
 a <- prior[1]+y; b <- prior[2]+n-y
 cat(sprintf('sensitivity_Beta(%d,%d)_mean=%.3f P(theta>0.5)=%.3f\n',
             prior[1],prior[2],a/(a+b),pbeta(.5,a,b,lower.tail=FALSE)))
}
cat('scope=conditional model probability; no programme effect or population representation\n')
