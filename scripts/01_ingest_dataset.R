#!/usr/bin/env Rscript
# Usage:
#   Rscript scripts/01_ingest_dataset.R GSE135251
#   Rscript scripts/01_ingest_dataset.R            # all configured datasets
#
# Run from the repository root.

suppressPackageStartupMessages({
  source("R/utils_io.R")
  source("R/config.R")
  source("R/labels.R")
  source("R/ingest.R")
})

args <- commandArgs(trailingOnly = TRUE)
project <- load_project_config("config/project.R")

dataset_ids <- if (length(args)) {
  args
} else {
  sub("\\.R$", "", list.files("config/datasets", pattern = "\\.R$"))
}

tsf_log("Datasets: ", paste(dataset_ids, collapse = ", "))

results <- lapply(dataset_ids, function(id) {
  tryCatch(ingest_dataset(id, project),
           error = function(e) {
             tsf_warn("FAILED ", id, ": ", conditionMessage(e))
             NULL
           })
})
names(results) <- dataset_ids

ok <- !vapply(results, is.null, logical(1))
present <- lapply(results[ok], function(r) r$audit$present)
if (length(present) > 1) invisible(comparable_conditions(present))

if (!all(ok)) {
  tsf_abort("Ingest failed for: ", paste(dataset_ids[!ok], collapse = ", "))
}
tsf_log("Done.")
