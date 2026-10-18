#!/usr/bin/env Rscript
# Base R exact beta updates and one-dimensional quadrature; no packages.
args <- commandArgs(trailingOnly=FALSE)
self <- sub('^--file=', '', args[grepl('^--file=',args)])
if (length(self)!=1L) stop('Run with Rscript.')
root <- normalizePath(file.path(dirname(self),'..'))
d <- read.csv(file.path(root,'data','unit07_help_seeking_routes.csv'),stringsAsFactors=FALSE)
ids <- c(sprintf('D%03d',1:24),sprintf('P%03d',1:12))
stopifnot(identical(names(d),c('volunteer_id','recruitment_route','reported_contact_7d')),
          identical(d$volunteer_id,ids),identical(d$recruitment_route,rep(c('desk','poster'),c(24,12))),
          !anyNA(d),all(d$reported_contact_7d %in% c(0,1)),
          sum(d$reported_contact_7d[1:24])==8,sum(d$reported_contact_7d[25:36])==0)
tail_diff <- function(a,b,c,d,h) {
  if (h<=-1) return(1)
  if (h>=1) return(0)
  lo <- max(0,h);hi <- min(1,1+h)
  center <- if (hi>lo) integrate(function(x) dbeta(x,a,b)*pbeta(x-h,c,d),lo,hi,
                                 rel.tol=1e-11,subdivisions=500L)$value else 0
  center + if (h<0) pbeta(1+h,a,b,lower.tail=FALSE) else 0
}
fit <- function(alpha,beta0) {
  a <- alpha+8;b <- beta0+16;c <- alpha;d <- beta0+12
  quant <- vapply(c(.025,.975),function(p) uniroot(function(h) 1-tail_diff(a,b,c,d,h)-p,
                                               c(-1,1),tol=1e-10)$root,numeric(1))
  list(shapes=c(a,b,c,d),means=c(a/(a+b),c/(c+d)),difference=a/(a+b)-c/(c+d),
       interval=quant,p0=tail_diff(a,b,c,d,0),p10=tail_diff(a,b,c,d,.1),
       prior_zero=exp(lbeta(alpha,beta0+12)-lbeta(alpha,beta0)),
       replicated_zero=exp(lbeta(c,d+12)-lbeta(c,d)))
}
cat('n=36 desk=8/24 poster=0/12 raw_difference=0.333\n')
for (name in c('primary','sensitivity')) {
  prior <- if(name=='primary') c(1,1) else c(2,2)
  v <- fit(prior[1],prior[2])
  cat(sprintf('%s prior=Beta(%d,%d) desk=Beta(%d,%d) poster=Beta(%d,%d)\n',
    name,prior[1],prior[2],v$shapes[1],v$shapes[2],v$shapes[3],v$shapes[4]))
  cat(sprintf('%s means=[%.3f,%.3f] difference=%.3f difference_95=[%.3f,%.3f] P_diff_gt_zero=%.3f P_diff_gt_0.10=%.3f\n',
    name,v$means[1],v$means[2],v$difference,v$interval[1],v$interval[2],v$p0,v$p10))
  cat(sprintf('%s prior_zero=%.3f replicated_zero=%.3f\n',name,v$prior_zero,v$replicated_zero))
}
cat('scope=fictional convenience volunteers; route association only; no causal effect or campus rate; exact no MCMC\n')
