#!/usr/bin/env Rscript
# Base-R entry point, same binary CSV and explicit Unit 3 comparison record.
args <- commandArgs(trailingOnly=FALSE)
self <- sub('^--file=', '', args[grepl('^--file=',args)])
if (length(self)!=1L) stop('Run with Rscript.')
root <- normalizePath(file.path(dirname(self),'..'))
spec <- readLines(file.path(root,'model-specs','unit03_intervals_decision.json'),warn=FALSE)
if (!any(grepl('"spec_id": "unit03-intervals-decision-v1"',spec,fixed=TRUE))) stop('Wrong model spec.')
d <- read.csv(file.path(root,'data','unit02_help_seeking.csv'),stringsAsFactors=FALSE)
if (nrow(d)!=12L || !identical(d$participant_id,sprintf('H%03d',1:12)) ||
    anyNA(d) || !all(d$help_sought_within_7_days %in% c(0L,1L))) stop('Invalid source rows.')
n <- nrow(d); y <- sum(d$help_sought_within_7_days)
if (y!=8L) stop('Unexpected count.')
p <- y/n; z <- qnorm(.975)
center <- (p+z^2/(2*n))/(1+z^2/n)
half <- z*sqrt(p*(1-p)/n+z^2/(4*n^2))/(1+z^2/n)
cat(sprintf('n=%d y=%d observed=%.3f target=reported_contact_probability\n',n,y,p))
cat(sprintf('wilson_95=[%.3f,%.3f] procedure=repeated_sample_coverage\n',center-half,center+half))
ci <- qbeta(c(.025,.975),10,6)
cat(sprintf('posterior=Beta(10,6) mean=%.3f credible_95=[%.3f,%.3f]\n',10/16,ci[1],ci[2]))
for (prior in list(c(2,2),c(1,1),c(8,8))) {
 mass <- pbeta(.6,prior[1]+y,prior[2]+n-y,lower.tail=FALSE)
 utility <- 2*mass-1
 choice <- if (utility>0) 'action' else 'no_action'
 cat(sprintf('Beta(%d,%d) P(theta>0.6)=%.3f toy_expected_utility=%.3f toy_choice=%s\n',
             prior[1],prior[2],mass,utility,choice))
}
cat('scope=fictional conditional model and hypothetical utility; no orientation effect or real service recommendation\n')
