# stability.R -- which peaks are stable within a condition.
#
# A peak (chr, N, k) is stable in a condition when maxT called it significant in
# at least `stable_frac` of that condition's samples, and every sample produced
# a result for it. This is the definition the whole downstream analysis rests
# on, so it lives in one function with one threshold, read from project.R.

stable_peaks_maxt <- function(maxt_individual, expected_samples,
                              alpha = 0.05, stable_frac = 0.9) {
  if (is.null(maxt_individual) || !nrow(maxt_individual)) return(NULL)
  n_expected <- length(unique(expected_samples))

  d <- maxt_individual
  d$is_sig <- d$p_empirical_maxT <= alpha
  # Both definitions are always computed and written, whichever is primary, so
  # the sensitivity of a result to the choice of null is visible without a rerun.
  d$is_sig_full <- if ("p_empirical_maxT_full" %in% colnames(d))
    d$p_empirical_maxT_full <= alpha else d$is_sig
  d$is_sig_all <- if ("p_empirical_maxT_all" %in% colnames(d))
    !is.na(d$p_empirical_maxT_all) & d$p_empirical_maxT_all <= alpha else NA
  by <- list(chr = as.character(d$chr), N = as.integer(d$N), k = as.integer(d$k))

  n_res <- stats::aggregate(d$sample, by, function(x) length(unique(x)))
  n_sig <- stats::aggregate(d$sample[d$is_sig],
                            lapply(by, function(v) v[d$is_sig]),
                            function(x) length(unique(x)))
  pw_mean <- stats::aggregate(d$power, by, mean)
  pw_med  <- stats::aggregate(d$power, by, stats::median)

  out <- n_res
  names(out)[names(out) == "x"] <- "n_samples_with_result"
  out <- merge(out, stats::setNames(pw_mean, c("chr", "N", "k", "mean_power")),
               by = c("chr", "N", "k"))
  out <- merge(out, stats::setNames(pw_med, c("chr", "N", "k", "median_power")),
               by = c("chr", "N", "k"))
  out <- merge(out, stats::setNames(n_sig, c("chr", "N", "k", "n_samples_significant")),
               by = c("chr", "N", "k"), all.x = TRUE)
  out$n_samples_significant[is.na(out$n_samples_significant)] <- 0L

  count_scheme <- function(flag, label) {
    sub <- d[flag %in% TRUE, , drop = FALSE]
    if (!nrow(sub)) { out[[label]] <<- 0L; return(invisible(NULL)) }
    agg <- stats::aggregate(sub$sample,
                            list(chr = as.character(sub$chr), N = as.integer(sub$N),
                                 k = as.integer(sub$k)),
                            function(x) length(unique(x)))
    v <- agg$x[match(paste(out$chr, out$N, out$k),
                     paste(agg$chr, agg$N, agg$k))]
    v[is.na(v)] <- 0L
    out[[label]] <<- as.integer(v)
    invisible(NULL)
  }
  count_scheme(d$is_sig_full, "n_samples_significant_full")
  count_scheme(d$is_sig_all,  "n_samples_significant_all")

  out$n_samples_expected <- n_expected
  out$freq <- out$k / out$N
  out$period <- out$N / out$k
  out$pct_samples_significant <- 100 * out$n_samples_significant / out$n_samples_expected
  complete <- out$n_samples_with_result == n_expected
  out$is_stable      <- complete & out$n_samples_significant >= stable_frac * n_expected
  out$is_stable_full <- complete & out$n_samples_significant_full >= stable_frac * n_expected
  out$is_stable_all  <- complete & out$n_samples_significant_all  >= stable_frac * n_expected
  out[order(-out$pct_samples_significant, out$chr, out$k), ]
}

#' Condition-level peak table for one branch (average or median).
#'
#' Keys are the maxT-stable peaks; amplitude, phase and the Fisher p-value come
#' from the condition's summary spectrum. Peaks that maxT calls stable but that
#' the Fisher test did not flag are kept -- that union is what the original
#' pipeline built with union_fisher_maxt(), and dropping it would silently
#' shrink the peak set.
condition_peak_table <- function(stability, spectra, branch) {
  if (is.null(stability) || !nrow(stability)) return(NULL)
  stable <- stability[stability$is_stable %in% TRUE, c("chr", "N", "k",
                                                       "n_samples_expected",
                                                       "n_samples_significant",
                                                       "pct_samples_significant",
                                                       "mean_power", "median_power")]
  if (!nrow(stable)) return(NULL)
  spec <- branch_spectrum(spectra, branch)
  amp_col <- if (branch == "average") "amplitude_mean" else "amplitude_median"
  cols <- intersect(c("chr", "N", "k", "freq", "period", "phase", "power",
                      "power_normalised", "window_power", "coverage",
                      amp_col, "p_value"), colnames(spec))
  spec <- spec[, cols, drop = FALSE]
  colnames(spec)[colnames(spec) == "p_value"] <- "p_value_fisher"

  out <- merge(stable, spec, by = c("chr", "N", "k"), all.x = TRUE)
  if ("p_value_fisher" %in% colnames(out)) {
    out$p_fdr_fisher <- stats::p.adjust(out$p_value_fisher, method = "BH")
  }
  out$branch <- branch
  # How much of this frequency the sampling pattern alone can explain. A peak
  # whose window_rank is near 1 sits exactly where the gaps are most periodic
  # and must not be read as biological without further evidence.
  if ("window_power" %in% colnames(out)) {
    out$window_rank <- rank(-out$window_power)
  }
  n_missing <- sum(is.na(out[[amp_col]]))
  if (n_missing) {
    tsf_warn(n_missing, " stable peak(s) absent from the ", branch,
             " summary spectrum; amplitude is NA for those rows")
  }
  out[order(-out$pct_samples_significant, out$chr, out$k), ]
}
