#!/usr/bin/env Rscript
# Base-R scalar conjugate calculation; fixed sigma and tau, no R packages.
args <- commandArgs(trailingOnly=FALSE)
self <- sub('^--file=', '', args[grepl('^--file=',args)])
if (length(self)!=1L) stop('Run with Rscript.')
root <- normalizePath(file.path(dirname(self),'..'))
spec <- readLines(file.path(root,'model-specs','unit06_partial_pooling.json'),warn=FALSE)
if (!any(grepl('unit06-workshop-partial-pooling-v1',spec,fixed=TRUE))) stop('Wrong model spec.')
d <- read.csv(file.path(root,'data','unit06_workshop_confidence.csv'),stringsAsFactors=FALSE)
groups <- sprintf('W%02d',1:8)
sizes <- c(4L,5L,6L,8L,10L,12L,15L,20L)
expected_ids <- unlist(lapply(seq_along(groups),function(j) sprintf('%s-S%02d',groups[j],seq_len(sizes[j]))),use.names=FALSE)
if (!identical(names(d),c('student_id','workshop_id','confidence_score')) ||
    nrow(d)!=80L || !identical(d$student_id,expected_ids) || anyNA(d) ||
    any(!is.finite(d$confidence_score)) || any(d$confidence_score<0|d$confidence_score>100) ||
    !identical(d$workshop_id,rep(groups,sizes))) stop('Invalid Unit 6 rows.')
y <- split(d$confidence_score,factor(d$workshop_id,levels=groups))
n <- unname(vapply(y,length,integer(1)))
raw <- unname(vapply(y,mean,numeric(1)))
sse <- sum(vapply(y,function(x) sum((x-mean(x))^2),numeric(1)))
sigma <- 8
fit <- function(tau) {
  vbar <- tau^2+sigma^2/n
  vmu <- 1/(1/100+sum(1/vbar))
  mmu <- vmu*(55/100+sum(raw/vbar))
  weight <- tau^2/vbar
  m <- weight*raw+(1-weight)*mmu
  vtheta <- 1/(1/tau^2+n/sigma^2)+(1-weight)^2*vmu
  diff <- m[1]-m[8]
  diff_var <- vtheta[1]+vtheta[8]-2*(1-weight[1])*(1-weight[8])*vmu
  list(mu=mmu,mu_sd=sqrt(vmu),mean=m,sd=sqrt(vtheta),weight=weight,
       p=pnorm(0,mean=diff,sd=sqrt(diff_var),lower.tail=FALSE))
}
p <- fit(5);tight <- fit(2);wide <- fit(10)
if (!all(is.finite(unlist(p))) || !all(is.finite(unlist(tight))) ||
    !all(is.finite(unlist(wide)))) stop('Non-finite posterior.')
prior_sd <- sqrt(100+25+64/4)
single_sd <- sqrt(100+25+64)
outside <- pnorm(0,mean=55,sd=single_sd)+pnorm(100,mean=55,sd=single_sd,lower.tail=FALSE)
cat(sprintf('n=%d groups=8 sizes=%s pooled=%.3f score_range=[%.1f,%.1f]\n',nrow(d),paste(sizes,collapse=','),mean(d$confidence_score),min(d$confidence_score),max(d$confidence_score)))
cat(sprintf('posterior_mu=%.3f mu_95=[%.3f,%.3f]\n',p$mu,qnorm(.025,p$mu,p$mu_sd),qnorm(.975,p$mu,p$mu_sd)))
for (j in seq_along(groups))
  cat(sprintf('group_%s n=%d raw=%.3f posterior=%.3f 95=[%.3f,%.3f] weight=%.3f\n',
      groups[j],n[j],raw[j],p$mean[j],qnorm(.025,p$mean[j],p$sd[j]),qnorm(.975,p$mean[j],p$sd[j]),p$weight[j]))
for (tau in c(5,2,10)) {
  s <- switch(as.character(tau),'5'=p,'2'=tight,'10'=wide)
  cat(sprintf('tau=%d P_A_gt_H=%.3f A=%.3f H=%.3f\n',tau,s$p,s$mean[1],s$mean[8]))
}
cat(sprintf('prior_W01_mean_95=[%.3f,%.3f] outside_score=%.4f prior_ICC=%.3f\n',qnorm(.025,55,prior_sd),qnorm(.975,55,prior_sd),outside,25/(25+64)))
cat(sprintf('within_sd=%.3f observed_SSE=%.3f replicated_expected_SSE=%.3f df=%d\n',sqrt(sse/72),sse,72*sigma^2,72))
cat('scope=fictional one-time grouped scores; no workshop effects or population claim; fixed scales; exact no MCMC\n')
