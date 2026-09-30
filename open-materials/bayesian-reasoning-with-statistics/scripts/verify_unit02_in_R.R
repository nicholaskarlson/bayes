#!/usr/bin/env Rscript
# Base-R independent numerical audit: read raw rows and recompute anchors.
args <- commandArgs(trailingOnly=FALSE)
self <- sub('^--file=', '', args[grepl('^--file=',args)])
if (length(self)!=1L) stop('Run with Rscript.')
resource <- normalizePath(file.path(dirname(self),'..'))
d <- read.csv(file.path(resource,'data','unit02_help_seeking.csv'),stringsAsFactors=FALSE)
stopifnot(identical(names(d),c('participant_id','help_sought_within_7_days')),
          identical(d$participant_id,sprintf('H%03d',1:12)),
          nrow(d)==12L, !anyNA(d), all(d$help_sought_within_7_days %in% c(0L,1L)),
          sum(d$help_sought_within_7_days)==8L)
# Exact beta-binomial probabilities from a different route than the student script.
beta_ratio <- function(k,a,b) choose(12,k)*beta(a+k,b+12-k)/beta(a,b)
near <- function(x,y) stopifnot(is.finite(x),abs(x-y)<1e-9)
near(sum(vapply(10:12,beta_ratio,0.0,2,2)),2/13)
near(sum(vapply(c(0:2,10:12),beta_ratio,0.0,2,2)),4/13)
near(10/16,.625)
near(pbeta(.5,10,6,lower.tail=FALSE),.84912109375)
near(sum(vapply(8:12,beta_ratio,0.0,10,6)),.5192982456140351)
near(pbeta(.5,9,5,lower.tail=FALSE),.8665771484375)
near(pbeta(.5,16,12,lower.tail=FALSE),.7789658308029175)
stopifnot(max(abs(qbeta(c(.025,.975),10,6)-c(.38380373,.83663568)))<1e-7)
out <- suppressWarnings(system2('Rscript',file.path(resource,'scripts','analyze_unit02.R'),stdout=TRUE))
if (length(attr(out,'status')) && attr(out,'status')!=0L) stop('Student script failed.')
needed <- c('n=12 y=8 observed=0.667 prior=Beta(2,2) posterior=Beta(10,6)',
 'prior_P(Y>=10_of_12)=0.154 prior_P(outer)=0.308',
 'posterior_mean=0.625 equal_tail_95=[0.384,0.837] P(theta>0.5)=0.849',
 'posterior_predictive_P(Yrep>=8_of_12)=0.519',
 'sensitivity_Beta(1,1)_mean=0.643 P(theta>0.5)=0.867',
 'sensitivity_Beta(8,8)_mean=0.571 P(theta>0.5)=0.779',
 'scope=conditional model probability; no programme effect or population representation')
stopifnot(setequal(out,needed))
cat('BAYESIAN_UNIT02_R_VERIFY_OK n=12 y=8 exact=yes\n')
