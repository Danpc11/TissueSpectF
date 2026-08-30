#!/usr/bin/env Rscript
# Check that every file the configs expect is present, before running ingest.
#
#   Rscript scripts/00_check_inputs.R
#
# Reports what is missing and, for each missing file, the closest names actually
# found in geo_dir -- GEO downloads often differ by a suffix (_raw_counts_GRCh38,
# _norm_counts_TPM, a date stamp), and that mismatch is the most common reason
# ingest fails on a new machine.

# When invoked as `./tsf check`, the CLI has already loaded the config and
# applied every --geo-dir / --results-dir override; reloading it here would
# quietly discard them and check the wrong directory.
if (!exists("project", inherits = TRUE) || !is.list(get0("project"))) {
  suppressPackageStartupMessages({
    source("R/utils_io.R"); source("R/config.R"); source("R/labels.R")
  })
  project <- load_project_config("config/project.R")
}
tsf_log("geo_dir: ", project$geo_dir)

if (!dir.exists(project$geo_dir)) {
  tsf_abort("geo_dir does not exist or is not readable: ", project$geo_dir)
}

available <- list.files(project$geo_dir)
tsf_log(length(available), " file(s) found in geo_dir")

closest <- function(target, pool, n = 3L) {
  if (!length(pool)) return(character(0))
  d <- utils::adist(tolower(target), tolower(pool))[1, ]
  pool[order(d)][seq_len(min(n, length(pool)))]
}

report <- function(label, filename) {
  path <- file.path(project$geo_dir, filename)
  if (file.exists(path)) {
    size_mb <- round(file.size(path) / 1024^2, 1)
    tsf_log("  OK      ", label, ": ", filename, " (", size_mb, " MB)")
    return(TRUE)
  }
  tsf_warn("  MISSING ", label, ": ", filename)
  cand <- closest(filename, available)
  if (length(cand)) {
    tsf_warn("          closest names present: ", paste(cand, collapse = ", "))
  }
  FALSE
}

ok <- report("annotation", project$annotation_file)

dataset_ids <- sub("\\.R$", "", list.files("config/datasets", pattern = "\\.R$"))
for (id in dataset_ids) {
  cfg <- load_dataset_config(id)
  tsf_log("dataset ", cfg$id, " (has_control_cohort = ", cfg$has_control_cohort, ")")
  ok <- report("counts", cfg$counts_file) && ok
  ok <- report("series matrix", cfg$series_matrix) && ok
}

# Writability of the output trees, checked now rather than after an hour of work.
for (d in c(project$interim_dir, project$results_dir)) {
  ensure_dir(d)
  probe <- file.path(d, ".tsf_write_probe")
  can_write <- tryCatch({ file.create(probe); file.remove(probe); TRUE },
                        warning = function(w) FALSE, error = function(e) FALSE)
  if (can_write) tsf_log("  OK      writable: ", d)
  else { tsf_warn("  NOT WRITABLE: ", d); ok <- FALSE }
}

if (!ok) {
  tsf_abort("Inputs incomplete. Fix the file names in config/datasets/<GSE>.R ",
            "(counts_file, series_matrix) or config/project.R (annotation_file, paths).")
}
tsf_log("All inputs present. Run: Rscript scripts/01_ingest_dataset.R")
