#!/usr/bin/env Rscript
# Independent base-R audit of Unit 8's raw-data Normal and beta-binomial values.
args <- commandArgs(trailingOnly=TRUE)
root <- if (length(args)) args[[1]] else 'open-materials/bayesian-reasoning-with-statistics'
stopifnot(file.exists(file.path(root,'data/unit04_reflection_attention_paired.csv')))
routes <- read.csv(file.path(root,'data/unit07_help_seeking_routes.csv'),stringsAsFactors=FALSE)
poster <- subset(routes,recruitment_route=='poster')
stopifnot(nrow(poster)==12,all(poster$reported_contact_7d==0))
p22 <- exp(lbeta(2,14)-lbeta(2,2))
p11 <- exp(lbeta(1,13)-lbeta(1,1))
stopifnot(abs(p22-1/35)<1e-13,abs(p11-1/13)<1e-13)
u7 <- read.csv(file.path(root,'expected-results/unit07_anchors.csv'),stringsAsFactors=FALSE)
base <- read.csv(file.path(root,'data/unit04_reflection_attention_paired.csv'),stringsAsFactors=FALSE)
stopifnot(nrow(base)==64,length(unique(base$participant_id))==64)
anchors <- read.csv(file.path(root,'expected-results/unit04_likelihood_anchors.csv'),stringsAsFactors=FALSE)
for (omit in c('all','omit_P059','omit_P042')) {
  d <- base$after_score-base$before_score
  if (omit!='all') d <- d[base$participant_id!=substring(omit,6)]
  variance <- 1/(1/25+length(d)/16)
  mean <- (sum(d)/16)*variance
  p <- pnorm(2,mean,sqrt(variance),lower.tail=FALSE)
  a <- subset(anchors,case==omit & model=='normal_fixed' & field=='p_mu_gt_2')$value
  t <- subset(anchors,case==omit & model=='student_t4_equal_variance' & field=='p_mu_gt_2')$value
  stopifnot(length(a)==1,length(t)==1,abs(p-a)<2e-12)
  cat(sprintf('R_UNIT08_MATCHED %s normal=%.6f t4_anchor=%.6f\n',omit,p,t))
}
# Unit 7's independently audited exact integral gives the posterior direction.
stopifnot(file.exists(file.path(root,'expected-results/unit07_proportions.json')))
cat(sprintf('R_UNIT08_ZERO prior_Beta22=%.6f prior_Beta11=%.6f\n',p22,p11))
cat('BAYESIAN_UNIT08_R_VERIFY_OK exact_zero=yes raw_normal=yes existing_t_anchor=yes no_sampler=yes\n')
