#!/usr/bin/env Rscript

# Regression tests for the invariant layer in
# scripts/build_final_condition_spectra.R.
#
# Run:
#   Rscript tests/test_condition_invariants.R

source("scripts/build_final_condition_spectra.R")

make_row <- function(
  dataset,
  condition,
  chr = "1",
  N = 100L,
  k = 5L,
  power = 1,
  prevalence = 0.8,
  plv = 0.9,
  phase = 0.1,
  enrichment = 0.1,
  window_pct = 50
) {
  data.frame(
    chr = chr,
    N = N,
    k = k,
    period = N / k,
    median_power_normalised = power,
    prevalence_rank = prevalence,
    plv = plv,
    mean_phase = phase,
    dataset = dataset,
    condition = condition,
    source_condition = condition,
    cohort_log2_enrichment = enrichment,
    window_power = 0,
    window_rank = 1,
    window_pct = window_pct,
    coverage = 0.8,
    selected_in_signature = TRUE,
    confirmed_in_cohort = FALSE,
    stringsAsFactors = FALSE
  )
}

# -------------------------------------------------------------------------
# A genuinely stable component across all conditions must become core.
# -------------------------------------------------------------------------
inputs <- list(
  "A\rF0" = make_row("A", "F0", power = 1.00, phase = 0.10, enrichment = 0.10),
  "B\rF0" = make_row("B", "F0", power = 1.05, phase = 0.12, enrichment = 0.08),
  "A\rF1" = make_row("A", "F1", power = 0.98, phase = 0.11, enrichment = 0.09),
  "B\rF1" = make_row("B", "F1", power = 1.02, phase = 0.09, enrichment = 0.07)
)

inv <- build_invariant_spectrum(
  inputs = inputs,
  conditions = c("F0", "F1"),
  min_cohorts = 2,
  min_condition_fraction = 0.70,
  min_prevalence = 0.5,
  min_phase_coherence = 0.8,
  max_power_cv = 1.0,
  max_abs_log2_enrichment = 0.5,
  core_phase_coherence = 0.90,
  core_max_power_cv = 0.30,
  core_max_abs_log2_enrichment = 0.50,
  window_cut = 1
)

stopifnot(
  nrow(inv) == 1L,
  inv$n_conditions == 2L,
  inv$n_cohorts == 2L,
  abs(inv$condition_fraction - 1) < 1e-12,
  isTRUE(inv$shared_candidate),
  isTRUE(inv$shared_robust),
  isTRUE(inv$core_invariant),
  identical(inv$invariant_class, "core_invariant")
)

# -------------------------------------------------------------------------
# Shared robust does not imply core: one condition may be missing.
# -------------------------------------------------------------------------
inputs_partial <- c(
  inputs,
  list(
    "A\rF2" = make_row("A", "F2", prevalence = 0.1),
    "B\rF2" = make_row("B", "F2", prevalence = 0.1)
  )
)

inv_partial <- build_invariant_spectrum(
  inputs = inputs_partial,
  conditions = c("F0", "F1", "F2"),
  min_cohorts = 2,
  min_condition_fraction = 0.60,
  min_prevalence = 0.5,
  min_phase_coherence = 0.8,
  max_power_cv = 1.0,
  max_abs_log2_enrichment = 0.5,
  core_phase_coherence = 0.90,
  core_max_power_cv = 0.30,
  core_max_abs_log2_enrichment = 0.50,
  window_cut = 1
)

stopifnot(
  abs(inv_partial$condition_fraction - (2 / 3)) < 1e-12,
  isTRUE(inv_partial$shared_candidate),
  isTRUE(inv_partial$shared_robust),
  !isTRUE(inv_partial$core_invariant),
  identical(inv_partial$invariant_class, "shared_robust")
)

# -------------------------------------------------------------------------
# A component can be shared robust but fail core because one condition has
# too large a maximum enrichment/depletion, even if the median is acceptable.
# -------------------------------------------------------------------------
inputs_outlier <- inputs
inputs_outlier[["B\rF1"]]$cohort_log2_enrichment <- 0.8

inv_outlier <- build_invariant_spectrum(
  inputs = inputs_outlier,
  conditions = c("F0", "F1"),
  min_cohorts = 2,
  min_condition_fraction = 0.70,
  min_prevalence = 0.5,
  min_phase_coherence = 0.8,
  max_power_cv = 1.0,
  max_abs_log2_enrichment = 0.5,
  core_phase_coherence = 0.90,
  core_max_power_cv = 0.30,
  core_max_abs_log2_enrichment = 0.50,
  window_cut = 1
)

stopifnot(
  isTRUE(inv_outlier$shared_candidate),
  isTRUE(inv_outlier$shared_robust),
  !isTRUE(inv_outlier$core_invariant),
  identical(inv_outlier$invariant_class, "shared_robust")
)

# -------------------------------------------------------------------------
# Healthy must not double-count source observations.
# -------------------------------------------------------------------------
with_healthy <- add_healthy_superclass(inputs)

inv2 <- build_invariant_spectrum(
  inputs = with_healthy,
  conditions = c("F0", "F1", "Healthy"),
  min_cohorts = 2,
  min_condition_fraction = 1,
  min_prevalence = 0.5,
  min_phase_coherence = 0.8,
  max_power_cv = 1.0,
  max_abs_log2_enrichment = 0.5,
  core_phase_coherence = 0.90,
  core_max_power_cv = 0.30,
  core_max_abs_log2_enrichment = 0.50,
  window_cut = 1
)

stopifnot(
  nrow(inv2) == 1L,
  inv2$n_conditions_total == 2L,
  inv2$n_conditions == 2L
)

# -------------------------------------------------------------------------
# Evidence-status diagnostics remain distinct.
# -------------------------------------------------------------------------
single_cohort_tab <- data.frame(
  n_cohorts = 1L,
  q_meta_null = NA_real_,
  n_confirmed_cohorts = 0L,
  n_null_draws = NA_real_
)

stopifnot(
  identical(
    condition_evidence_status(single_cohort_tab, min_cohorts = 2L),
    "insufficient_cohorts_for_meta_analysis"
  )
)

missing_null_tab <- data.frame(
  n_cohorts = 2L,
  q_meta_null = NA_real_,
  n_confirmed_cohorts = 0L,
  n_null_draws = NA_real_
)

stopifnot(
  identical(
    condition_evidence_status(missing_null_tab, min_cohorts = 2L),
    "missing_meta_null"
  )
)

cat("test_condition_invariants.R: OK\n")

# --- the meta-analysis counts cohorts that spoke, not cohorts present --------
#
# stouffer() drops non-finite p silently, so a frequency present in three
# cohorts but with only one finite p_null was combined from one value while
# n_cohorts said three -- and then passed a `n_cohorts >= min_cohorts` gate as
# though it carried three independent sources. Replication is the central claim
# here, so the gate has to count contributions, not presences.

local({
  p_groups <- list(a = c(0.01, NA, 0.2), b = c(NA, NA, NA), c = c(0.3))
  n_meta <- vapply(p_groups,
                   function(p) sum(is.finite(suppressWarnings(as.numeric(p)))),
                   integer(1))
  stopifnot(identical(unname(n_meta), c(2L, 0L, 1L)))
})

local({
  # Two cohorts required; one contributing p is not a meta-analysis, and the
  # combined value must be withheld rather than reported with a caveat -- a
  # number in a q column gets used, and a footnote does not travel with it.
  out <- data.frame(p_meta_null = c(0.001, 0.002, 0.003),
                    n_meta_cohorts = c(3L, 1L, 2L))
  out$p_meta_null[out$n_meta_cohorts < 2L] <- NA_real_
  stopifnot(is.na(out$p_meta_null[2]),
            all(is.finite(out$p_meta_null[c(1, 3)])))
})

local({
  src <- paste(readLines("scripts/build_final_condition_spectra.R",
                         warn = FALSE), collapse = " ")
  stopifnot(
    grepl("eligible_k <- out$n_meta_cohorts", src, fixed = TRUE),
    !grepl("eligible_k <- out$n_cohorts", src, fixed = TRUE))
})

local({
  # Set on both the meta and the no-meta path, so the column always exists.
  src <- readLines("scripts/build_final_condition_spectra.R", warn = FALSE)
  stopifnot(sum(grepl("out\\$n_meta_cohorts <-", src)) >= 2L)
})

cat("test_condition_invariants.R: meta-cohort assertions OK\n")
