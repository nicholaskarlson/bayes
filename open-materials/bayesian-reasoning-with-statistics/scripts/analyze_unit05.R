#!/usr/bin/env Rscript
# Base-R route from the unchanged fictional CSV; no R package dependencies.
args <- commandArgs(trailingOnly=FALSE)
self <- sub('^--file=', '', args[grepl('^--file=',args)])
if (length(self)!=1L) stop('Run with Rscript.')
root <- normalizePath(file.path(dirname(self),'..'))
spec <- readLines(file.path(root,'model-specs','unit05_regression.json'),warn=FALSE)
if (!any(grepl('unit05-cross-sectional-regression-v1',spec,fixed=TRUE))) stop('Wrong model spec.')
d <- read.csv(file.path(root,'data','unit05_belonging_study.csv'),stringsAsFactors=FALSE)
if (!identical(names(d),c('participant_id','study_hours_per_week','long_commute','belonging_score')) ||
    nrow(d)!=72L || !identical(d$participant_id,sprintf('B%03d',1:72)) ||
    anyNA(d) || any(!d$long_commute %in% c(0L,1L))) stop('Invalid Unit 5 rows.')
y <- d$belonging_score;h <- d$study_hours_per_week;z <- d$long_commute
if (any(h<0|h>168) || any(y<0|y>100)) stop('Out of range.')
x <- h-5;X <- cbind(1,x,z)
fit <- function(M) {
  b <- as.vector(solve(crossprod(M),crossprod(M,y)))
  sd <- sqrt(sum((y-M%*%b)^2)/(length(y)-ncol(M)))
  se <- sd*sqrt(diag(solve(crossprod(M))))
  # diag() can carry predictor names into se. A named c() would then create
  # lower.x/upper.x, so indexing it by "lower"/"upper" would produce NA.
  list(slope=b[2],
       lower=b[2]+qt(.025,length(y)-ncol(M))*se[2],
       upper=b[2]+qt(.975,length(y)-ncol(M))*se[2],
       commute=if(ncol(M)==3L)b[3] else NA_real_)
}
u <- fit(X[,1:2]);a <- fit(X)
if (!all(is.finite(unlist(u[c('slope','lower','upper')]))) ||
    !all(is.finite(unlist(a)))) stop('Non-finite ordinary-regression estimate or interval.')
posterior <- function(prior_slope_sd, sigma) {
  P0 <- diag(1/c(225,prior_slope_sd^2,100))
  V <- solve(P0+crossprod(X)/sigma^2)
  b <- as.vector(V%*%(P0%*%c(55,0,0)+crossprod(X,y)/sigma^2))
  s <- sqrt(V[2,2]);ci <- qnorm(c(.025,.975),b[2],s)
  residuals <- y-as.vector(X%*%b)
  predicted_sd <- sqrt(sigma^2+rowSums((X%*%V)*X))
  c(slope=b[2],lower=ci[1],upper=ci[2],prob=pnorm(1,b[2],s,lower.tail=FALSE),
    two_hours=2*b[2],observed=sum(abs(residuals)>16),expected=sum(2*pnorm(-16/predicted_sd)))
}
p <- posterior(2,8);skeptical <- posterior(.4,8);scale_10 <- posterior(2,10)
prior_outside <- function(xx,zz) {
  sd <- sqrt(225+4*xx^2+100*zz^2+64)
  pnorm(0,55,sd)+pnorm(100,55,sd,lower.tail=FALSE)
}
cat(sprintf('n=%d long_commute=%d hours=[%.1f,%.1f] score=[%.1f,%.1f]\n',nrow(d),sum(z),min(h),max(h),min(y),max(y)))
cat(sprintf('OLS_slope_unadjusted=%.3f 95_t=[%.3f,%.3f]\n',u$slope,u$lower,u$upper))
cat(sprintf('OLS_slope_adjusted=%.3f 95_t=[%.3f,%.3f] commute_coefficient=%.3f\n',a$slope,a$lower,a$upper,a$commute))
for (name in c('primary','skeptical','scale_10')) {
  v <- switch(name,primary=p,skeptical=skeptical,scale_10=scale_10)
  cat(sprintf('posterior_%s_slope=%.3f 95_credible=[%.3f,%.3f] P(beta>1)=%.3f two_hour_difference=%.3f\n',
              name,v['slope'],v['lower'],v['upper'],v['prob'],v['two_hours']))
}
cat(sprintf('prior_predictive_outside=[%.3f,%.3f] observed_tail=%d replicated_expected_tail=%.3f\n',
            prior_outside(0,0),prior_outside(1.5,1),p['observed'],p['expected']))
cat('scope=fictional cross-sectional association, not an effect or population estimate; fixed-sigma working model; exact no MCMC\n')
