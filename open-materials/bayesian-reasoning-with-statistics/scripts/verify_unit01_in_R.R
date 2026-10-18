#!/usr/bin/env Rscript
# Independent base-R audit of the Unit 1 entry point and printed anchors.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
if (length(file_arg) != 1L) stop("Run with Rscript.")
resource <- normalizePath(file.path(dirname(file_arg), ".."))
target <- file.path(resource, "scripts", "inspect_unit01.R")
output <- suppressWarnings(system2("Rscript", target, stdout = TRUE))
if (length(attr(output, "status")) && attr(output, "status") != 0L) stop("Unit 1 R route failed.")
data <- read.csv(file.path(resource, "data", "unit04_reflection_attention_paired.csv"),
                 stringsAsFactors = FALSE)
# Recompute the declared mean from integer tenths to avoid depending on the route.
to_tenths <- function(x) as.integer(round(10 * x))
delta <- to_tenths(data$after_score) - to_tenths(data$before_score)
stopifnot(nrow(data) == 64L, sum(delta) == 1405L, sum(delta < 0) == 20L,
          sum(delta > 0) == 44L, delta[1] == -2L)
needed <- c("pairs=64", "mean_after_minus_before=2.1953125",
            "negative_differences=20", "positive_differences=44",
            "P001_after_minus_before=-0.2",
            "scope=description_only; no causal effect or Bayesian model estimated")
stopifnot(setequal(output, needed))
cat("BAYESIAN_UNIT01_R_VERIFY_OK pairs=64 mean=2.1953125 negative=20\n")
