#!/usr/bin/env Rscript
# Base-R descriptive route for Unit 1; reads the original paired CSV.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
if (length(file_arg) != 1L) stop("Run this script with Rscript.")
resource <- normalizePath(file.path(dirname(file_arg), ".."))
data <- read.csv(file.path(resource, "data", "unit04_reflection_attention_paired.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(identical(names(data), c("participant_id", "before_score", "after_score")))
stopifnot(nrow(data) == 64L, identical(data$participant_id, sprintf("P%03d", 1:64)))
stopifnot(all(is.finite(data$before_score)), all(is.finite(data$after_score)))
stopifnot(all(data$before_score >= 0 & data$before_score <= 100))
stopifnot(all(data$after_score >= 0 & data$after_score <= 100))
differences <- data$after_score - data$before_score
cat(sprintf("pairs=%d\n", nrow(data)))
cat(sprintf("mean_after_minus_before=%.7f\n", mean(differences)))
cat(sprintf("negative_differences=%d\n", sum(differences < 0)))
cat(sprintf("positive_differences=%d\n", sum(differences > 0)))
cat(sprintf("P001_after_minus_before=%.1f\n", differences[1]))
cat("scope=description_only; no causal effect or Bayesian model estimated\n")
