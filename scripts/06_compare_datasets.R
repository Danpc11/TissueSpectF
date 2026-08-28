#!/usr/bin/env Rscript
# Constant signature, transitions, and the cross-dataset comparison.
#
#   Rscript scripts/06_compare_datasets.R
#   Rscript scripts/06_compare_datasets.R --branch=average
#
# Runs over comparable_conditions(): conditions present in every dataset. A
# condition missing on one side is reported, not silently turned into an empty
# intersection. Nothing is pooled across datasets -- the point is replication.

suppressPackageStartupMessages(source("R/paths.R"))
source("R/utils_io.R"); tsf_source_pipeline()

cli <- tsf_args()
project <- load_project_config("config/project.R")
branches <- if (is.null(cli$branch)) c("average", "median") else cli$branch
dataset_ids <- tsf_dataset_ids(cli)
if (length(dataset_ids) < 2) tsf_abort("Need at least two datasets to compare")

compare_dir <- file.path(project$results_dir, "comparison")
ensure_dir(compare_dir)

loaded <- lapply(dataset_ids, function(id) {
  inp <- tsf_stage_inputs(project, id, need = c("stability", "maxt"))
  inp$peaks <- lapply(c(average = "average", median = "median"), function(br)
    stats::setNames(lapply(inp$conditions, function(c)
      read_tsv_tsf(p_peaks(inp$paths, br, c), required = FALSE)), inp$conditions))
  inp
})
names(loaded) <- dataset_ids

present <- lapply(loaded, function(x) x$conditions)
common <- comparable_conditions(present)

for (branch in branches) {
  tsf_log("=== branch: ", branch, " ===")

  # ---- per dataset ---------------------------------------------------------
  signature_by_ds <- list(); transitions_by_ds <- list(); conditions_by_ds <- list()

  for (id in dataset_ids) {
    x <- loaded[[id]]
    sig <- constant_signature(x$stability, common)
    if (!is.null(sig)) {
      sig$branch <- branch
      write_tsv_tsf(sig, file.path(x$paths$base, "comparison",
                                   sprintf("constant_signature_%s.tsv", branch)))
      tsf_log(id, ": constant signature = ", nrow(sig), " peak(s) over ",
              paste(common, collapse = ", "))
    } else {
      tsf_log(id, ": constant signature is empty")
    }
    signature_by_ds[[id]] <- sig

    cond_tbl <- do.call(rbind, lapply(common, function(c) {
      p <- x$peaks[[branch]][[c]]
      if (is.null(p) || !nrow(p)) return(NULL)
      p$condition <- c
      p
    }))
    conditions_by_ds[[id]] <- cond_tbl

    trans <- list()
    for (i in seq_len(length(common) - 1L)) {
      t <- transition_table(common[i], common[i + 1L], x$peaks[[branch]],
                            x$stability, x$maxt, branch)
      if (is.null(t)) {
        tsf_log(id, ": ", common[i], " -> ", common[i + 1L], ": no shared stable peak")
        next
      }
      tsf_log(id, ": ", common[i], " -> ", common[i + 1L], ": ", nrow(t),
              " shared peak(s), ", sum(t$p_power_fdr <= 0.05, na.rm = TRUE),
              " with a significant power change")
      trans[[length(trans) + 1]] <- t
    }
    trans <- if (length(trans)) do.call(rbind, trans) else NULL
    if (!is.null(trans)) {
      write_tsv_tsf(trans, file.path(x$paths$base, "comparison",
                                     sprintf("transitions_%s.tsv", branch)))
    }
    transitions_by_ds[[id]] <- trans
  }

  # ---- across datasets -----------------------------------------------------
  sig_cross <- cross_datasets(signature_by_ds, c("chr", "N", "k"))
  if (!is.null(sig_cross)) {
    write_tsv_tsf(sig_cross, file.path(compare_dir,
                                       sprintf("constant_signature_shared_%s.tsv", branch)))
  }

  cond_cross <- cross_datasets(conditions_by_ds, c("chr", "N", "k", "condition"))
  if (!is.null(cond_cross)) {
    write_tsv_tsf(cond_cross, file.path(compare_dir,
                                        sprintf("conditions_shared_%s.tsv", branch)))
  }

  trans_cross <- cross_datasets(transitions_by_ds, c("chr", "N", "k", "transition"))
  if (!is.null(trans_cross)) {
    trans_cross <- add_replication_flags(trans_cross, dataset_ids)
    write_tsv_tsf(trans_cross, file.path(compare_dir,
                                         sprintf("transitions_shared_%s.tsv", branch)))
  }
}

tsf_log("Comparison written to ", compare_dir)
