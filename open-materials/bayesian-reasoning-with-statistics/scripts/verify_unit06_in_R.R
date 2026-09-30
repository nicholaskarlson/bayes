#!/usr/bin/env Rscript
# Independent matrix-precision audit, not sourcing Unit 6's scalar R route.
args <- commandArgs(trailingOnly=FALSE)
self <- sub('^--file=', '', args[grepl('^--file=',args)])
if (length(self)!=1L) stop('Run with Rscript.')
root <- normalizePath(file.path(dirname(self),'..'))
d <- read.csv(file.path(root,'data','unit06_workshop_confidence.csv'),stringsAsFactors=FALSE)
groups <- sprintf('W%02d',1:8);n <- c(4L,5L,6L,8L,10L,12L,15L,20L)
ids <- unlist(lapply(seq_along(n),function(j) sprintf('%s-S%02d',groups[j],seq_len(n[j]))),use.names=FALSE)
stopifnot(identical(names(d),c('student_id','workshop_id','confidence_score')),
          nrow(d)==80L,identical(d$student_id,ids),identical(d$workshop_id,rep(groups,n)),
          !anyNA(d),all(is.finite(d$confidence_score)),all(d$confidence_score>=0),all(d$confidence_score<=100))
y <- split(d$confidence_score,factor(d$workshop_id,levels=groups))
raw <- unname(vapply(y,mean,numeric(1)))
sse <- sum(vapply(y,function(a) sum((a-mean(a))^2),numeric(1)))
near <- function(a,b,tol=1e-7) stopifnot(all(is.finite(a)),all(abs(a-b)<tol))
fit <- function(tau) {
  P <- matrix(0,9,9)
  P[1,1] <- 1/100+8/tau^2
  P[1,2:9] <- -1/tau^2
  P[2:9,1] <- -1/tau^2
  for (j in seq_along(n)) P[j+1,j+1] <- 1/tau^2+n[j]/64
  rhs <- c(55/100,n*raw/64)
  V <- solve(P)
  b <- as.numeric(solve(P,rhs))
  vs <- diag(V)[2:9]
  delta <- b[2]-b[9]
  ds <- sqrt(V[2,2]+V[9,9]-2*V[2,9])
  list(mu=b[1],mu_sd=sqrt(V[1,1]),m=b[2:9],sd=sqrt(vs),
       lo=qnorm(.025,b[2:9],sqrt(vs)),hi=qnorm(.975,b[2:9],sqrt(vs)),
       w=n*tau^2/(64+n*tau^2),p=pnorm(0,mean=delta,sd=ds,lower.tail=FALSE))
}
anchors <- read.csv(file.path(root,'expected-results','unit06_anchors.csv'),stringsAsFactors=FALSE)
stopifnot(nrow(anchors)==24L,identical(unique(anchors$case),c('primary_tau5','sensitivity_tau2','sensitivity_tau10')))
results <- list(primary_tau5=fit(5),sensitivity_tau2=fit(2),sensitivity_tau10=fit(10))
for (case in names(results)) {
  a <- anchors[anchors$case==case,]
  v <- results[[case]]
  stopifnot(identical(a$group,groups),identical(a$n,n))
  for (pair in list(cbind(a$raw,raw),cbind(a$posterior,v$m),cbind(a$lower,v$lo),
                    cbind(a$upper,v$hi),cbind(a$weight,v$w),cbind(a$mu,rep(v$mu,8)),
                    cbind(a$p_A_gt_H,rep(v$p,8)))) near(pair[,1],pair[,2])
}
p <- results$primary_tau5
near(p$mu,55.032,0.001);near(p$m[1],60.272,0.001);near(p$m[8],58.856,0.001)
# The displayed 4540.640 is rounded to three decimals; the CSV's exact
# decimal within-group sum of squares is 4540.64025.
near(sse,4540.64025);near(results$sensitivity_tau2$p,.401,0.001)
near(results$sensitivity_tau10$p,.783,0.001)

out <- suppressWarnings(system2('Rscript',file.path(root,'scripts','analyze_unit06.R'),stdout=TRUE))
if (length(attr(out,'status')) && attr(out,'status')!=0L) stop('Unit 6 R route failed.')
out <- trimws(out);out <- out[nzchar(out)]
if (length(out)!=16L) stop('Unit 6 R route must print 16 labeled lines:\n',paste(out,collapse='\n'))
prefix <- c('n=','posterior_mu=',paste0('group_',groups,' n='),
            'tau=5 ','tau=2 ','tau=10 ','prior_W01_mean_95=','within_sd=','scope=')
if (!all(startsWith(out,prefix))) stop('Unit 6 R route labels differ:\n',paste(out,collapse='\n'))
if (!identical(out[16],'scope=fictional one-time grouped scores; no workshop effects or population claim; fixed scales; exact no MCMC'))
  stop('Unit 6 claim boundary changed.')
num <- function(line) as.numeric(regmatches(line,gregexpr('-?[0-9]+(\\.[0-9]+)?',line,perl=TRUE))[[1]])
expected <- list(c(80,8,n,mean(d$confidence_score),min(d$confidence_score),max(d$confidence_score)),
                 c(p$mu,95,qnorm(.025,p$mu,p$mu_sd),qnorm(.975,p$mu,p$mu_sd)))
for (j in seq_along(groups)) expected[[length(expected)+1L]] <- c(j,n[j],raw[j],p$m[j],95,p$lo[j],p$hi[j],p$w[j])
for (tau in c(5,2,10)) {
  v <- results[[paste0(if(tau==5) 'primary_tau' else 'sensitivity_tau',tau)]]
  expected[[length(expected)+1L]] <- c(tau,v$p,v$m[1],v$m[8])
}
prior_sd <- sqrt(100+25+64/4);score_sd <- sqrt(100+25+64)
outside <- pnorm(0,55,score_sd)+pnorm(100,55,score_sd,lower.tail=FALSE)
expected[[14]] <- c(1,95,qnorm(.025,55,prior_sd),qnorm(.975,55,prior_sd),outside,25/89)
expected[[15]] <- c(sqrt(sse/72),sse,72*64,72)
for (i in seq_along(expected)) {
  printed <- num(out[i]);gold <- expected[[i]]
  if (length(printed)!=length(gold) || any(!is.finite(printed)) ||
      any(abs(printed-gold)>.00051))
    stop('Unit 6 R route differs from independent matrix audit on line ',i,': ',out[i],
         '\nExpected: ',paste(gold,collapse=', '))
}
# The prior outside-score field prints four decimals and has a tighter bound.
near(num(out[14])[5],outside,.000051)
cat('BAYESIAN_UNIT06_R_VERIFY_OK n=80 groups=8 matrix_vs_scalar=yes exact=yes\n')
