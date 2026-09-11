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

# --help imprime la cabecera del propio archivo y sale con 0. Antes `--help` se
# tomaba como el VALOR de la bandera anterior --"Missing value for --help",
# "No such file: --help"-- o el script corria con la configuracion vacia. Un
# script que no sabe explicarse es un script que nadie usa bien.
if (any(commandArgs(TRUE) %in% c("-h", "--help"))) {
  self <- sub("^--file=", "",
              grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  if (!is.na(self) && file.exists(self)) {
    hdr <- readLines(self, warn = FALSE)
    hdr <- hdr[!grepl("^#!", hdr)]                    # fuera el shebang
    stop_at <- which(!grepl("^#", hdr) & nzchar(hdr))[1]
    if (is.na(stop_at)) stop_at <- length(hdr) + 1L
    hdr <- hdr[seq_len(stop_at - 1L)]
    cat(paste(sub("^#[ ]?", "", hdr[grepl("^#", hdr)]), collapse = "\n"), "\n")
  }
  quit(save = "no", status = 0)
}

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

#' Human-readable size. A series matrix is tens of kilobytes, and rounding it
#' to "0 MB" reads like a zero-byte file, which is alarming for no reason.
human_size <- function(bytes) {
  if (bytes >= 1024^2) return(paste0(round(bytes / 1024^2, 1), " MB"))
  if (bytes >= 1024) return(paste0(round(bytes / 1024), " KB"))
  paste0(bytes, " B")
}

report <- function(label, filename) {
  path <- file.path(project$geo_dir, filename)
  if (file.exists(path)) {
    tsf_log("  OK      ", label, ": ", filename, " (", human_size(file.size(path)), ")")
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

# Only the datasets asked for (positional, via ./tsf check <ids>); every config
# in the directory when none is named. A stale config for a cohort this run
# does not use must not fail the check.
dataset_ids <- if (exists("opt", inherits = TRUE) && length(get0("opt")$datasets))
  get0("opt")$datasets else sub("\\.R$", "", list.files("config/datasets", pattern = "\\.R$"))
for (id in dataset_ids) {
  cfg <- load_dataset_config(id)
  tsf_log("dataset ", cfg$id, " (has_control_cohort = ", cfg$has_control_cohort, ")")
  ok <- report("counts", cfg$counts_file) && ok
  # GEO configs declare series_matrix; recount3 / matrix configs declare
  # metadata_file. One of the two has to be there, and file.exists(NULL) is an
  # R error, not a MISSING.
  pheno <- cfg$series_matrix %||% cfg$metadata_file
  if (is.null(pheno)) {
    tsf_warn("  MISSING phenotype: config declares neither series_matrix nor metadata_file"); ok <- FALSE
  } else {
    ok <- report(if (is.null(cfg$series_matrix)) "metadata" else "series matrix", pheno) && ok
  }
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
  tsf_abort("Inputs incomplete. Fix the file names in config/datasets/<id>.R ",
            "(counts_file, series_matrix or metadata_file) or config/project.R ",
            "(annotation_file, paths). To check only some datasets: ./tsf check <id> ...")
}
tsf_log("All inputs present. Next: ./tsf ingest")
