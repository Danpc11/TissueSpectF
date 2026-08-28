#!/usr/bin/env Rscript
# maxT permutation test, per sample and per chromosome. THIS IS THE SLOW STAGE
# (B x n_samples x n_chromosomes spectra). Parallel over chromosomes.
#
#   Rscript scripts/03_maxt.R GSE135251
#   Rscript scripts/03_maxt.R GSE135251 --cond=F3        # one condition
#   Rscript scripts/03_maxt.R GSE135251 --force          # recompute existing
#
# Without --force, a condition whose maxt_individual_<cond>.tsv already exists is
# skipped and says so. Reuse is opt-out here and nowhere else, because this is
# the only stage where recomputation costs hours.

suppressPackageStartupMessages(source("R/paths.R"))
source("R/utils_io.R"); tsf_source_pipeline()

cli <- tsf_args()
project <- load_project_config("config/project.R")

for (id in tsf_dataset_ids(cli)) {
  inp <- tsf_stage_inputs(project, id)
  n_cores <- maxt_cores(inp$chrom_idx)
  tsf_log("=== ", id, ": maxT (B = ", project$maxt$B, ", ", n_cores, " cores) ===")

  for (cond in tsf_conditions(inp$conditions, cli)) {
    out_path <- p_maxt(inp$paths, cond)
    if (!cli$force && file.exists(out_path) && file.size(out_path) > 0) {
      tsf_log("  ", cond, ": reusing existing maxT (use --force to recompute)")
      next
    }
    t0 <- Sys.time()
    res <- maxt_condition(inp$dataset, cond, inp$chrom_idx, project$maxt, n_cores)
    if (is.null(res)) { tsf_warn("  ", cond, ": no maxT result"); next }
    write_tsv_tsf(res, out_path)
    write_tsv_tsf(data.frame(condition = cond,
                             sample_id = attr(res, "expected_samples")),
                  file.path(inp$paths$maxt, sprintf("expected_samples_%s.tsv", cond)))
    tsf_log("  ", cond, ": ", nrow(res), " rows in ",
            round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
  }
}
tsf_log("Done. Next: scripts/04_stability.R")
