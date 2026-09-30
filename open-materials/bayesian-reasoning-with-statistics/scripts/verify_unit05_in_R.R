#!/usr/bin/env Rscript
# Independent base-R matrix audit, not sourcing the main R entry point.
args <- commandArgs(trailingOnly=FALSE)
self <- sub('^--file=', '', args[grepl('^--file=',args)])
if (length(self)!=1L) stop('Run with Rscript.')
root <- normalizePath(file.path(dirname(self),'..'))
d <- read.csv(file.path(root,'data','unit05_belonging_study.csv'),stringsAsFactors=FALSE)
stopifnot(identical(names(d),c('participant_id','study_hours_per_week','long_commute','belonging_score')),
          nrow(d)==72L,identical(d$participant_id,sprintf('B%03d',1:72)),!anyNA(d),
          sum(d$long_commute)==36,all(d$long_commute %in% c(0L,1L)),
          min(d$study_hours_per_week)==1.8,max(d$study_hours_per_week)==8.5,
          min(d$belonging_score)==25.9,max(d$belonging_score)==76.7)
near <- function(actual,expected,tol=1e-7) stopifnot(is.finite(actual),abs(actual-expected)<tol)
y <- d$belonging_score; M <- cbind(1,d$study_hours_per_week-5,d$long_commute)
ordinary <- function(A) {
  b <- qr.solve(A,y);r <- y-as.vector(A%*%b)
  s <- sqrt(sum(r*r)/(length(y)-ncol(A)))
  v <- s*s*solve(crossprod(A))
  c(b, b[2]+qt(c(.025,.975),length(y)-ncol(A))*sqrt(v[2,2]))
}
u <- ordinary(M[,1:2]);a <- ordinary(M)
near(u[2],3.27396531818725);near(u[3],2.13935611971153);near(u[4],4.40857451666297)
near(a[2],2.7655354702025);near(a[3],-4.50645169723869)
near(a[4],1.57184175772652);near(a[5],3.95922918267847)
model <- function(slope_sd,sigma) {
  prior <- c(55,0,0);V0 <- diag(c(225,slope_sd^2,100))
  cov <- solve(solve(V0)+crossprod(M)/sigma^2)
  mean <- as.vector(cov%*%(solve(V0,prior)+crossprod(M,y)/sigma^2))
  sd <- sqrt(cov[2,2])
  tail <- sum(abs(y-as.vector(M%*%mean))>16)
  vpred <- sigma^2+rowSums((M%*%cov)*M)
  list(mean=mean,sd=sd,ci=qnorm(c(.025,.975),mean[2],sd),prob=1-pnorm(1,mean[2],sd),
       tail=tail,rep=sum(2*pnorm(-16/sqrt(vpred))))
}
p <- model(2,8);tight <- model(.4,8);large <- model(2,10)
near(p$mean[2],2.55416429294408);near(p$ci[1],1.42328612449297)
near(p$ci[2],3.68504246139518);near(p$prob,.996465496205975)
stopifnot(p$tail==3);near(p$rep,3.58896223097279)
near(tight$mean[2],.852098280873906);near(tight$prob,.328594420216753)
near(large$mean[2],2.45039237199425);near(large$prob,.980329213309264)
pred <- function(x,z) {v <- 225+4*x*x+100*z*z+64;pnorm(0,55,sqrt(v))+1-pnorm(100,55,sqrt(v))}
near(pred(0,0),.00466734882508108);near(pred(1.5,1),.0149640010240850)
out <- suppressWarnings(system2('Rscript',file.path(root,'scripts','analyze_unit05.R'),stdout=TRUE))
if (length(attr(out,'status')) && attr(out,'status')!=0L) stop('Unit 5 R entry failed.')
out <- trimws(out)
out <- out[nzchar(out)]
if (length(out)!=8L) stop('Unit 5 R route must print exactly eight labeled lines:\n',paste(out,collapse='\n'))
prefix <- c('n=','OLS_slope_unadjusted=','OLS_slope_adjusted=',
            'posterior_primary_slope=','posterior_skeptical_slope=',
            'posterior_scale_10_slope=','prior_predictive_outside=','scope=')
for (i in seq_along(prefix)) {
  if (!startsWith(out[i],prefix[i])) stop('Unit 5 R output label at line ',i,': ',out[i])
}
if (!all(grepl(' 95_t=[',out[2:3],fixed=TRUE)) ||
    !all(grepl(' 95_credible=[',out[4:6],fixed=TRUE)) ||
    !all(grepl(' P(beta>1)=',out[4:6],fixed=TRUE)) ||
    !all(grepl(' two_hour_difference=',out[4:6],fixed=TRUE)) ||
    !grepl(' observed_tail=',out[7],fixed=TRUE) ||
    !grepl(' replicated_expected_tail=',out[7],fixed=TRUE))
  stop('Unit 5 R output fields changed:\n',paste(out,collapse='\n'))
if (!identical(out[[8L]],
  'scope=fictional cross-sectional association, not an effect or population estimate; fixed-sigma working model; exact no MCMC'))
  stop('Unit 5 R claim boundary changed: ',out[8])
# Compare printed three-decimal numbers against the independent raw-data audit.
# The 95, beta>1, and scale_10 tokens are also checked in their labeled rows.
expected <- list(
  c(72,36,1.8,8.5,25.9,76.7),
  c(u[2],95,u[3],u[4]),
  c(a[2],95,a[4],a[5],a[3]),
  c(p$mean[2],95,p$ci,1,p$prob,2*p$mean[2]),
  c(tight$mean[2],95,tight$ci,1,tight$prob,2*tight$mean[2]),
  c(10,large$mean[2],95,large$ci,1,large$prob,2*large$mean[2]),
  c(pred(0,0),pred(1.5,1),p$tail,p$rep))
for (i in seq_along(expected)) {
  printed <- as.numeric(regmatches(out[i],gregexpr('-?[0-9]+(\\.[0-9]+)?',out[i],perl=TRUE))[[1]])
  if (length(printed)!=length(expected[[i]]) ||
      any(!is.finite(printed)) || any(abs(printed-expected[[i]])>0.00051))
    stop('Unit 5 R route differs from the independent audit on line ',i,
         ': ',out[i],'\nExpected numeric fields: ',paste(expected[[i]],collapse=', '))
}
cat('BAYESIAN_UNIT05_R_VERIFY_OK n=72 OLS=yes posterior=yes checks=yes\n')
