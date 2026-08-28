#!/usr/bin/env Rscript
# Gene-level reconstruction for every stable peak.
#
#   Rscript scripts/05_peaks_genes.R                        # everything
#   Rscript scripts/05_peaks_genes.R GSE162694 --cond=F2    # one condition
#   Rscript scripts/05_peaks_genes.R --branch=median        # one branch
#
# This replaces regenerar_picos_genes_condicion.R. That script existed because
# the reconstruction could only be rebuilt from objects left in a finished R
# session, so a change of threshold meant re-running the whole pipeline or
# keeping the session alive. Here every input is a file: after changing
# stable_frac, run 04 and then this, and nothing else is recomputed.

suppressPackageStartupMessages(source("R/paths.R"))
source("R/utils_io.R"); tsf_source_pipeline()

cli <- tsf_args()
project <- load_project_config("config/project.R")
branches <- if (is.null(cli$branch)) c("average", "median") else cli$branch

for (id in tsf_dataset_ids(cli)) {
  inp <- tsf_stage_inputs(project, id)
  tsf_log("=== ", id, ": peak-gene reconstruction ===")

  for (cond in tsf_conditions(inp$conditions, cli)) {
    for (branch in branches) {
      peaks <- read_tsv_tsf(p_peaks(inp$paths, branch, cond), required = FALSE)
      if (is.null(peaks) || !nrow(peaks)) {
        tsf_warn("  ", cond, "/", branch, ": no peak table, run 04 first")
        next
      }
      out_dir <- p_peak_genes_dir(inp$paths, branch, cond)
      unlink(list.files(out_dir, pattern = "^pico_chr.*\\.tsv$", full.names = TRUE))
      n <- write_peak_gene_tables(peaks, inp$dataset$genes, inp$chrom_idx,
                                  out_dir, branch, cond)
      tsf_log("  ", cond, "/", branch, ": ", n, " peak file(s)")
    }
  }
}
tsf_log("Done. Next: scripts/06_compare_datasets.R")
