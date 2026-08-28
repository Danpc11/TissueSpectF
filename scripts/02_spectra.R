#!/usr/bin/env Rscript
# FFT of chromosome-ordered expression, per sample and per condition.
#
#   Rscript scripts/02_spectra.R                 # all datasets, all conditions
#   Rscript scripts/02_spectra.R GSE135251 --cond=F2
#
# Cheap to re-run (seconds to minutes). Reads only the common format.

suppressPackageStartupMessages(source("R/paths.R"))
source("R/utils_io.R"); tsf_source_pipeline()

cli <- tsf_args()
project <- load_project_config("config/project.R")

for (id in tsf_dataset_ids(cli)) {
  inp <- tsf_stage_inputs(project, id)
  tsf_log("=== ", id, ": spectra over ", length(inp$chrom_idx), " chromosomes ===")
  for (cond in tsf_conditions(inp$conditions, cli)) {
    spec <- compute_condition_spectra(inp$dataset, cond, inp$chrom_idx)
    if (is.null(spec)) next
    is_summary <- spec$sample %in% c("avg_signal", "median_signal")
    write_tsv_tsf(spec[is_summary, ], p_spectra_condition(inp$paths, cond))
    write_tsv_tsf(spec[!is_summary, ], p_spectra_samples(inp$paths, cond))
    tsf_log("  ", cond, ": ", sum(is_summary), " summary rows, ",
            sum(!is_summary), " sample rows")
  }
}
tsf_log("Done. Next: scripts/03_maxt.R")
