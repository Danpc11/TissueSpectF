#!/usr/bin/env Rscript
# validate_recount3_configs.R -- check (and optionally migrate) config/datasets/R3_*.R
#
# Usage:
#   Rscript scripts/validate_recount3_configs.R --datasets R3_LIVER,R3_SRP217231 [--vocabulary liver_fibrosis]
#   Rscript scripts/validate_recount3_configs.R --datasets ... --migrate   # fix header fields in place
#   Rscript scripts/validate_recount3_configs.R --all                      # every R3_*.R in the repo
# Only the configs a run uses are validated (--datasets); --all is for auditing
# the repository, so an experimental config of another tissue never stops a run.
#
# fetch_recount3_project() never overwrites an existing config, so a file
# written by an earlier version keeps whatever that version emitted. Configs
# from before the current generator can have:
#   - series_matrix = <pheno> and phenotype_format = "tsv" instead of
#     source = "matrix" + metadata_file = <pheno>   -> ingest reads the pheno as a GEO matrix
#   - vocabulary = "case_control" while the GTEx rule assigns
#     Control_external_study                         -> every GTEx sample dropped at labelling
#   - count_id_type other than ENSEMBL              -> no gene matches
#
# --migrate rewrites ONLY those header fields, as text, and leaves the
# condition_rules block (which may carry hand-written SRA rules) untouched. A
# backup <file>.bak is written first. Nothing is ever deleted.
if (any(commandArgs(TRUE) %in% c("-h", "--help"))) {
  cat(paste(sub("^# ?", "", grep("^#", readLines(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))[-1], value = TRUE)), collapse = "\n"), "\n")
  quit(save = "no")
}
suppressWarnings({ source("R/utils_io.R"); tsf_load_all("R") })
args <- commandArgs(trailingOnly = TRUE)
flag <- function(n, d = NULL) {
  h <- grep(paste0("^", n, "="), args, value = TRUE)
  if (length(h)) return(sub(paste0("^", n, "="), "", h[1]))
  i <- match(n, args); if (!is.na(i) && length(args) > i) args[i + 1] else d
}
migrate <- "--migrate" %in% args
want_vocab <- flag("--vocabulary", NULL)
ds <- flag("--datasets", NULL)
if (is.null(ds) && !"--all" %in% args) {
  tsf_abort("Pass --datasets R3_LIVER,R3_SRP... (the configs THIS run uses) or --all. ",
            "An experimental config for another tissue must not stop a liver run.")
}
files <- if (!is.null(ds)) {
  f <- file.path("config", "datasets", paste0(trimws(strsplit(ds, ",")[[1]]), ".R"))
  miss <- f[!file.exists(f)]
  if (length(miss)) tsf_abort("config(s) not found: ", paste(miss, collapse = ", "))
  f
} else Sys.glob("config/datasets/R3_*.R")
if (!length(files)) { tsf_log("no recount3 configs to validate"); quit(save = "no") }

vocab_levels <- function(id) {
  f <- file.path("config", "vocabularies", paste0(id, ".R"))
  if (!file.exists(f)) return(NULL)
  v <- tryCatch(source(f, local = TRUE)$value, error = function(e) NULL)
  v$levels
}

bad_any <- FALSE
for (f in files) {
  txt <- readLines(f, warn = FALSE)
  cfg <- tryCatch(source(f, local = TRUE)$value, error = function(e) NULL)
  probs <- character(0); fixes <- list()
  if (is.null(cfg)) { tsf_warn(f, ": does not parse"); bad_any <- TRUE; next }

  if (!identical(cfg$source, "matrix")) {
    probs <- c(probs, "source != \"matrix\"")
    fixes$source <- TRUE
  }
  if (is.null(cfg$metadata_file)) {
    probs <- c(probs, "no metadata_file")
    fixes$metadata_file <- cfg$series_matrix %||% NA_character_
  }
  if (!is.null(cfg$phenotype_format)) { probs <- c(probs, "obsolete phenotype_format"); fixes$drop_pf <- TRUE }
  if (!identical(cfg$count_id_type, "ENSEMBL")) probs <- c(probs, "count_id_type != ENSEMBL (recount3 ids are Ensembl)")
  # the vocabulary must contain every label the rules assign
  assigned <- unlist(lapply(cfg$condition_rules, function(r) r$assign))
  lv <- vocab_levels(cfg$vocabulary)
  if (!is.null(assigned) && !is.null(lv) && !all(assigned %in% lv)) {
    probs <- c(probs, sprintf("rules assign %s but vocabulary '%s' lacks it",
                              paste(setdiff(assigned, lv), collapse = ","), cfg$vocabulary))
    if (!is.null(want_vocab)) fixes$vocabulary <- want_vocab
  }
  if (!is.null(want_vocab) && !identical(cfg$vocabulary, want_vocab)) {
    probs <- c(probs, sprintf("vocabulary '%s' != requested '%s'", cfg$vocabulary, want_vocab))
    fixes$vocabulary <- want_vocab
  }
  if (is.null(cfg$condition_rules) || !length(cfg$condition_rules)) {
    probs <- c(probs, "condition_rules empty: every sample will be unlabelled (SRA configs need hand-written rules)")
  }

  if (!length(probs)) { tsf_log("OK   ", f); next }
  bad_any <- TRUE
  tsf_warn("INVALID ", f, ": ", paste(probs, collapse = "; "))
  if (!migrate || !length(fixes)) next

  file.copy(f, paste0(f, ".bak"), overwrite = TRUE)
  if (isTRUE(fixes$source) && !any(grepl("^\\s*source\\s*=", txt))) {
    i <- grep("^\\s*counts_file\\s*=", txt)[1]
    txt <- append(txt, '  source        = "matrix",', after = i)
  }
  if (!is.null(fixes$metadata_file)) {
    i <- grep("^\\s*series_matrix\\s*=", txt)
    if (length(i)) txt[i[1]] <- sub("series_matrix\\s*=", "metadata_file =", txt[i[1]])
    else txt <- append(txt, sprintf('  metadata_file = "%s",', fixes$metadata_file),
                       after = grep("^\\s*source\\s*=", txt)[1])
  }
  if (isTRUE(fixes$drop_pf)) txt <- txt[!grepl("^\\s*phenotype_format\\s*=", txt)]
  if (!is.null(fixes$vocabulary)) {
    i <- grep("^\\s*vocabulary\\s*=", txt)[1]
    txt[i] <- sprintf('  vocabulary    = "%s",', fixes$vocabulary)
  }
  writeLines(txt, f)
  chk <- tryCatch(source(f, local = TRUE)$value, error = function(e) NULL)
  if (is.null(chk)) { file.copy(paste0(f, ".bak"), f, overwrite = TRUE); tsf_abort(f, ": migration broke the file; restored") }
  tsf_log("MIGRATED ", f, " (backup ", basename(f), ".bak)")
}
if (bad_any && !migrate) {
  tsf_log("Run with --migrate to fix header fields; condition_rules are never touched.")
  quit(status = 1)
}
if (bad_any && migrate) {
  tsf_log("Re-run without --migrate to confirm; anything still INVALID needs a hand edit.")
  quit(status = 1)
}
