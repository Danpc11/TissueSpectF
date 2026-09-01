#!/usr/bin/env Rscript

# Lightweight regression tests for the invariant layer in
# scripts/build_final_condition_spectra.R.
#
# Run:
#   Rscript tests/test_condition_invariants.R
#
# The constructor is sourced without executing main() because its bottom guard
# uses sys.nframe() == 0L.

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

inputs <- list(
  "A\rF0" = make_row(
    "A", "F0",
    power = 1.00,
    phase = 0.10
  ),
  "B\rF0" = make_row(
    "B", "F0",
    power = 1.05,
    phase = 0.12
  ),
  "A\rF1" = make_row(
    "A", "F1",
    power = 0.98,
    phase = 0.11
  ),
  "B\rF1" = make_row(
    "B", "F1",
    power = 1.02,
    phase = 0.09
  )
)

inv <- build_invariant_spectrum(
  inputs = inputs,
  conditions = c("F0", "F1"),
  min_cohorts = 2,
  min_condition_fraction = 1,
  min_prevalence = 0.5,
  min_phase_coherence = 0.8,
  max_power_cv = 0.2,
  max_abs_log2_enrichment = 0.5,
  window_cut = 1
)

stopifnot(
  nrow(inv) == 1L,
  inv$n_conditions == 2L,
  inv$n_cohorts == 2L,
  abs(inv$condition_fraction - 1) < 1e-12,
  identical(
    inv$invariant_class,
    "robust"
  )
)

# Healthy must not double-count source observations.
with_healthy <- add_healthy_superclass(inputs)

inv2 <- build_invariant_spectrum(
  inputs = with_healthy,
  conditions = c("F0", "F1", "Healthy"),
  min_cohorts = 2,
  min_condition_fraction = 1,
  min_prevalence = 0.5,
  min_phase_coherence = 0.8,
  max_power_cv = 0.2,
  max_abs_log2_enrichment = 0.5,
  window_cut = 1
)

stopifnot(
  nrow(inv2) == 1L,
  inv2$n_conditions_total == 2L,
  inv2$n_conditions == 2L
)

# A single-cohort condition must be reported as insufficient for meta-analysis,
# not as merely missing a permutation null.
single_cohort_tab <- data.frame(
  n_cohorts = 1L,
  q_meta_null = NA_real_,
  n_confirmed_cohorts = 0L,
  n_null_draws = NA_real_
)

stopifnot(
  identical(
    condition_evidence_status(
      single_cohort_tab,
      min_cohorts = 2L
    ),
    "insufficient_cohorts_for_meta_analysis"
  )
)

# Missing null is a separate status when enough cohorts exist.
missing_null_tab <- data.frame(
  n_cohorts = 2L,
  q_meta_null = NA_real_,
  n_confirmed_cohorts = 0L,
  n_null_draws = NA_real_
)

stopifnot(
  identical(
    condition_evidence_status(
      missing_null_tab,
      min_cohorts = 2L
    ),
    "missing_meta_null"
  )
)

cat("test_condition_invariants.R: OK\n")
