args <- commandArgs(trailingOnly=TRUE)
root <- if (length(args) > 0) args[[1]] else 'open-materials/bayesian-reasoning-with-statistics'
reference <- read.csv(file.path(root, 'expected-results/prelude_grid.csv'), stringsAsFactors = FALSE)
stopifnot(identical(names(reference), c('theta', 'prior', 'relative_fit', 'product', 'posterior')),
          nrow(reference) == 9L, all(abs(reference$theta - (1:9)/10) < 1e-12))
theta <- reference$theta
prior <- theta * (1 - theta)
prior <- prior / sum(prior)
# dbinom includes the same combination factor at every proposed theta; use
# relative likelihoods to theta = 0.7 so the rendered table stays readable.
fit <- dbinom(8, size=12, prob=theta) / dbinom(8, size=12, prob=0.7)
product <- prior * fit
posterior <- product / sum(product)
for (name in c('prior', 'relative_fit', 'product', 'posterior')) {
  calculated <- switch(name, prior=prior, relative_fit=fit, product=product,
                       posterior=posterior)
  stopifnot(all(abs(reference[[name]] - calculated) <= 0.0000005))
}
stopifnot(abs(sum(product) - 0.432173) <= 0.0000005,
          abs(sum(posterior[theta > 0.5]) - 0.745038) <= 0.0000005,
          abs((1 - pbeta(0.5, 10, 6)) - 0.84912109375) < 1e-12,
          abs(sum(posterior[theta > 0.5]) - (1 - pbeta(0.5, 10, 6))) > 0.10)
cat('PRELUDE_GRID_R_VERIFY_OK candidates=9 grid_gt_half=0.745038 continuous_gt_half=0.849121\n')
