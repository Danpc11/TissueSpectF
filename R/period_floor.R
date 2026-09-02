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
    tsf_abort("--min-period auto needs a coverage column, and the table for ",
              label, " has none or has no finite value in it.\n",
              "  The technical floor is 2/coverage, so coverage per frequency ",
              "is what makes it computable.\n",
              "  Either pass --min-period <genes> as a flat floor, or check ",
              "that the per-sample spectra carry `coverage` -- the consensus ",
              "stage joins it from there, and an older results tree written ",
              "before that column existed will not have it. Re-running the ",
              "`spectra` stage regenerates it.")
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

#' Conservative rank-1 BH diagnostic: m/(B+1), capped at 1.
#'
#' NOT A LOWER BOUND. This was documented as "the smallest q any frequency
#' could reach", and that is wrong. BH takes q_i = min over j >= i of
#' m * p_j / j, so the achievable minimum depends on how many hypotheses tie at
#' the permutation floor. If all m of them sit at p = 1/(B+1), then at rank m
#' the adjusted value is m * p / m = 1/(B+1) -- a factor of m below what this
#' function returns. With k hypotheses tied at the floor it is m/(k(B+1)).
#'
#' So this is the value at rank 1 with no ties: the worst case, useful as a
#' warning that the pointwise route is hostile at this family size, and useless
#' as grounds for saying nothing can pass. Reporting it as a bound understates
#' the method's power by up to a factor of m, which for m = 10,030 is not a
#' rounding difference -- it could be the difference between reporting nothing
#' and reporting a real set of components.
#'
#' Use it to decide whether to look at the pointwise route. Do not use it to
#' decide that the pointwise route is closed.
bh_rank1_diagnostic <- function(m, n_draws) {
  if (!is.finite(m) || !is.finite(n_draws) || n_draws < 1) return(1)
  min(1, m / (n_draws + 1))
}

#' The actual achievable minimum BH q, given how many p-values tie at the floor.
#'
#' `n_tied` is how many of the m hypotheses reach the permutation floor
#' 1/(n_draws+1). At rank n_tied the BH value is m/(n_tied * (n_draws+1)), and
#' that is the honest answer to "could anything pass".
bh_achievable_q <- function(m, n_draws, n_tied = 1L) {
  if (!is.finite(m) || !is.finite(n_draws) || n_draws < 1) return(1)
  n_tied <- max(1L, min(as.integer(n_tied), as.integer(m)))
  min(1, m / (n_tied * (n_draws + 1)))
}

#' Draws needed for the rank-1 (no-ties) case to reach `target`.
#'
#' The worst case, so an upper bound on what is needed -- with ties at the
#' floor, fewer draws suffice.
draws_for_bh <- function(m, target = 0.05) {
  if (!is.finite(m) || !is.finite(target) || target <= 0) return(NA_integer_)
  as.integer(ceiling(m / target)) - 1L
}

#' Deprecated name for bh_rank1_diagnostic(), kept so nothing breaks silently.
reachable_bh_q <- function(m, n_draws) bh_rank1_diagnostic(m, n_draws)
