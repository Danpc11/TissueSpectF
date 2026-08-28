#!/usr/bin/env Rscript
# Stable peaks per condition, and the condition-level peak tables per branch.
#
#   Rscript scripts/04_stability.R
#
# Fast: reads maxT and spectra from disk. Re-run this after changing
# maxt$alpha or maxt$stable_frac in config/project.R -- no need to redo maxT.

suppressPackageStartupMessages(source("R/paths.R"))
source("R/utils_io.R"); tsf_source_pipeline()

cli <- tsf_args()
project <- load_project_config("config/project.R")

for (id in tsf_dataset_ids(cli)) {
  inp <- tsf_stage_inputs(project, id, need = c("maxt", "spectra"))
  tsf_log("=== ", id, ": stability (alpha = ", project$maxt$alpha,
          ", stable_frac = ", project$maxt$stable_frac, ") ===")

  for (cond in tsf_conditions(inp$conditions, cli)) {
    m <- inp$maxt[[cond]]
    if (is.null(m)) { tsf_warn("  ", cond, ": no maxT file, run 03 first"); next }
    expected <- read_tsv_tsf(file.path(inp$paths$maxt,
                                       sprintf("expected_samples_%s.tsv", cond)),
                             required = FALSE)
    expected_samples <- if (!is.null(expected)) expected$sample_id else unique(m$sample)

    st <- stable_peaks_maxt(m, expected_samples,
                            alpha = project$maxt$alpha,
                            stable_frac = project$maxt$stable_frac)
    write_tsv_tsf(st, p_stability(inp$paths, cond))
    tsf_log("  ", cond, ": ", sum(st$is_stable), " stable of ", nrow(st), " peaks")

    for (branch in c("average", "median")) {
      pk <- condition_peak_table(st, inp$spectra[[cond]], branch)
      if (is.null(pk)) next
      write_tsv_tsf(pk, p_peaks(inp$paths, branch, cond))
    }
  }
}
tsf_log("Done. Next: scripts/05_peaks_genes.R")
