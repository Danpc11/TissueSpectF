#!/usr/bin/env Rscript
# calibrate_null.R -- family-wise error of the permutation nulls under noise
# that is autocorrelated but NOT periodic.
#
# WHY
# ---
# The permutation null of maxT assumes the observed values are exchangeable
# across positions, i.e. white noise. Expression along a chromosome is not
# white: neighbouring genes are co-expressed and expression tracks GC content,
# gene density and replication timing over long stretches. A signal with that
# kind of autocorrelation concentrates power at low frequency, and a test that
# compares it against white-noise permutations will call the lowest frequencies
# significant whether or not any periodicity is there.
#
# selfcheck.R injects a sinusoid into white noise and recovers it. That is the
# easy direction. This script is the other direction: no periodicity at all,
# and the question is how often each scheme finds one anyway. A scheme whose
# rate at nominal alpha is far above alpha is not testing periodicity.
#
# Usage:
#   Rscript scripts/calibrate_null.R [--phi 0.3,0.6] [--reps 40] [--B 200]
#                                    [--N 1200] [--coverage 0.7] [--alpha 0.05]
#
# Prints a table: for each phi, the fraction of replicates with at least one
# frequency below alpha under `full`, under each block scheme, and under `all`.

args <- commandArgs(trailingOnly = TRUE)
flag <- function(name, default) {
  i <- match(name, args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}
phis     <- as.numeric(strsplit(flag("--phi", "0,0.3,0.6"), ",")[[1]])
reps     <- as.integer(flag("--reps", "40"))
B        <- as.integer(flag("--B", "200"))
N        <- as.integer(flag("--N", "1200"))
coverage <- as.numeric(flag("--coverage", "0.7"))
alpha    <- as.numeric(flag("--alpha", "0.05"))
block_sizes <- c(10L, 20L, 50L)

root <- normalizePath(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(), value = TRUE)[1])), ".."))
for (f in c("utils_io.R", "grid.R", "multitaper.R", "maxt.R")) {
  source(file.path(root, "R", f))
}
if (!exists("tsf_spectrum")) {
  tsf_spectrum <- function(y, terms, estimator = NULL) gls_spectrum(y, terms)
}

set.seed(1)
t <- sort(sample.int(N, round(coverage * N)))
terms <- gls_prepare(t, N)

cat(sprintf("N = %d, coverage = %.2f, B = %d, %d replicate(s), alpha = %.3f\n",
            N, coverage, B, reps, alpha))
cat(sprintf("%-8s %-8s %-10s %-10s %-10s %-8s %-8s\n", "phi", "full",
            "block10", "block20", "block50", "all", "k_top"))

for (phi in phis) {
  res <- t(vapply(seq_len(reps), function(i) {
    x <- if (phi == 0) stats::rnorm(N) else
      as.numeric(stats::arima.sim(list(ar = phi), n = N))
    y <- x[t]
    r <- suppressWarnings(permutation_gls_test(
      y, terms, B = B, seed = 1000L + i, block_sizes = block_sizes,
      primary_scheme = "all"))
    hit <- function(col) if (col %in% names(r) && any(is.finite(r[[col]])))
      any(r[[col]] <= alpha, na.rm = TRUE) else NA
    c(full    = hit("p_empirical_maxT_full"),
      block10 = hit("p_empirical_maxT_block10"),
      block20 = hit("p_empirical_maxT_block20"),
      block50 = hit("p_empirical_maxT_block50"),
      all     = hit("p_empirical_maxT_all"),
      k_top   = r$k[which.min(r$p_empirical_maxT_full)])
  }, numeric(6)))
  cat(sprintf("%-8.2f %-8.3f %-10.3f %-10.3f %-10.3f %-8.3f %-8.0f\n", phi,
              mean(res[, "full"], na.rm = TRUE),
              mean(res[, "block10"], na.rm = TRUE),
              mean(res[, "block20"], na.rm = TRUE),
              mean(res[, "block50"], na.rm = TRUE),
              mean(res[, "all"], na.rm = TRUE),
              stats::median(res[, "k_top"])))
}

cat("\nRead the `full` column against alpha. A value far above it means the\n",
    "scheme reacts to autocorrelation, not to periodicity; `all` is the\n",
    "default for that reason. `k_top` is the frequency index of the most\n",
    "significant false peak: low k = long period, which is where a red-noise\n",
    "artefact lands and where chromosomal structure is claimed.\n", sep = "")
