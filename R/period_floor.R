# period_floor.R -- the smallest period worth testing, technical and biological.
#
# TWO FLOORS, ONE MAXIMUM
# -----------------------
# Technical. On a chromosome with coverage rho the observed genes sit, on
# average, 1/rho grid positions apart. A period shorter than 2/rho is below the
# Nyquist limit of the sampling that actually happened: nothing at that period
# can be distinguished from an alias of something slower, no matter how strong
# it looks. `--min-period auto` computes 2/rho per row and adds a margin, since
# coverage is an average and the local gaps around any given peak may be worse.
#
# Biological. Independently of resolution, a period of two or three genes is not
# a claim anyone would make about chromosome organisation. `--min-period-
# biological` sets that floor in genes, flat across chromosomes.
#
# The applied floor is the larger of the two: a frequency has to clear both.
#
# WHY THIS FILE EXISTS
# --------------------
# This logic was defined inside scripts/build_final_condition_spectra.R and
# nowhere else, which meant the library build filtered frequencies that the
# consensus stage had already tested. That ordering is the problem: the
# multiplicity correction in consensus ran over every frequency including the
# unresolvable ones, so the family was inflated by frequencies that were never
# candidates. Filtering afterwards cleans the output but cannot give back the
# statistical power the inflated family cost.
#
# One definition, applied before the null. Both callers use this function.

#' Restrict a frequency table to periods above the technical and biological floors.
#'
#' @param long data frame with `period`, and `coverage` when min_period is "auto"
#' @param label what to call this table in the log
#' @param min_period "off", "auto", or a number of genes
#' @param period_margin added to (or multiplied into) the technical floor
#' @param margin_mode "add" or "mult"
#' @param min_period_biological flat floor in genes
#' @param quiet suppress the log line
#' @return `long` with rows below the floor removed
apply_period_floor <- function(long,
                               label = "table",
                               min_period = "off",
                               period_margin = 2,
                               margin_mode = "add",
                               min_period_biological = 0,
                               quiet = FALSE) {
  if (is.null(long) || !nrow(long)) return(long)

  min_period <- tolower(as.character(min_period))
  margin_mode <- tolower(as.character(margin_mode))
  if (!margin_mode %in% c("add", "mult")) {
    tsf_abort("margin_mode must be 'add' or 'mult', got '", margin_mode, "'")
  }
  period_margin <- as.numeric(period_margin)
  min_period_biological <- as.numeric(min_period_biological)

  # "off" still honours a biological floor if one was asked for: the two are
  # independent claims, and silently dropping the biological one because the
  # technical one is off would be surprising.
  if (identical(min_period, "off") && !(min_period_biological > 0)) {
    return(long)
  }

  cov <- if ("coverage" %in% names(long)) {
    suppressWarnings(as.numeric(long$coverage))
  } else {
    rep(NA_real_, nrow(long))
  }
  cov[!is.finite(cov) | cov <= 0] <- NA_real_

  if (identical(min_period, "auto") && all(is.na(cov))) {
    tsf_abort("--min-period auto needs per-chromosome coverage, and this table ",
              "has none. Give --min-period a number of genes instead.")
  }

  technical <- if (identical(min_period, "auto")) {
    base <- 2 / cov
    if (identical(margin_mode, "mult")) base * (1 + period_margin)
    else base + period_margin
  } else if (identical(min_period, "off")) {
    rep(NA_real_, nrow(long))
  } else {
    v <- suppressWarnings(as.numeric(min_period))
    if (!is.finite(v) || v < 0) {
      tsf_abort("min_period must be 'off', 'auto' or a non-negative number, ",
                "got '", min_period, "'")
    }
    rep(v, nrow(long))
  }

  floor_period <- pmax(technical, min_period_biological, na.rm = TRUE)

  # A row with no period is kept: it is not evidence that the period is short.
  keep <- !is.finite(long$period) | long$period >= floor_period
  n_drop <- sum(!keep)

  if (n_drop && !quiet) {
    tsf_log(sprintf(
      "  period floor for %s: %d of %d frequencies below the floor (%.1f%%), %d remain",
      label, n_drop, nrow(long), 100 * n_drop / nrow(long), nrow(long) - n_drop))
  }
  long[keep, , drop = FALSE]
}

#' Smallest BH q-value any frequency could reach, given the family size.
#'
#' At rank 1 the BH-adjusted p is p_min * m, and a permutation p cannot fall
#' below 1/(B+1). So the floor is m/(B+1) -- and a q-value is capped at 1, which
#' the reported number must respect: printing "q is 10" is arithmetic, not a
#' q-value, and invites the reader to think a threshold was missed by a factor
#' rather than being unreachable outright.
reachable_bh_q <- function(m, n_draws) {
  if (!is.finite(m) || !is.finite(n_draws) || n_draws < 1) return(1)
  min(1, m / (n_draws + 1))
}

#' Draws needed for the pointwise BH route to be able to reach `target`.
draws_for_bh <- function(m, target = 0.05) {
  if (!is.finite(m) || !is.finite(target) || target <= 0) return(NA_integer_)
  as.integer(ceiling(m / target)) - 1L
}
