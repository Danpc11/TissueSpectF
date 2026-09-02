#!/usr/bin/env Rscript
# clean_results.R -- delete the contents of results_dir, and nothing else.
#
# This used to be a backslash-continued one-liner in the Makefile. That does not
# work: make hands the continuation to the shell, the shell leaves it alone
# because it sits inside single quotes, and R then receives a stray backslash
# and refuses to parse the expression. A destructive operation guarded by code
# that never ran is worse than no guard, so it lives in a file that can be read
# and tested.
#
#   Rscript scripts/clean_results.R [--results-dir <dir>] [--dry-run]
#
# The guard is structural rather than a length heuristic: it refuses the
# filesystem root, the home directory and the repository itself. A results_dir
# that resolves to one of those means the configuration is wrong, and clearing
# it would take the code or the user's home with it.

suppressWarnings({
  source("R/utils_io.R")
  source("R/config.R")
})

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args

flag_value <- function(name) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (length(hit)) return(sub(paste0("^", name, "="), "", hit[1]))
  i <- match(name, args)
  if (!is.na(i) && length(args) > i) args[i + 1] else NULL
}

target <- flag_value("--results-dir") %||% flag_value("--results_dir")
if (is.null(target)) {
  target <- load_project_config("config/project.R")$results_dir
}

if (is.null(target) || !nzchar(target)) {
  tsf_abort("results_dir is unset: nothing to clean, and nothing safe to guess.")
}

absolute <- normalizePath(target, mustWork = FALSE)
forbidden <- unique(normalizePath(c("/", "~", "."), mustWork = FALSE))
if (absolute %in% forbidden) {
  tsf_abort("Refusing to clean ", absolute, ": that is the filesystem root, ",
            "the home directory or the repository itself. Point results_dir ",
            "somewhere of its own.")
}

if (!dir.exists(absolute)) {
  tsf_log("Nothing to clean: ", absolute, " does not exist")
  quit(save = "no", status = 0L)
}

entries <- list.files(absolute, all.files = TRUE, no.. = TRUE, full.names = TRUE)
if (!length(entries)) {
  tsf_log("Nothing to clean: ", absolute, " is already empty")
  quit(save = "no", status = 0L)
}

if (dry_run) {
  tsf_log("Would remove ", length(entries), " entr(y/ies) from ", absolute, ":")
  for (e in entries) tsf_log("  ", basename(e))
  quit(save = "no", status = 0L)
}

tsf_log("Removing ", length(entries), " entr(y/ies) from ", absolute)
unlink(entries, recursive = TRUE, force = TRUE)

still <- list.files(absolute, all.files = TRUE, no.. = TRUE)
if (length(still)) {
  tsf_abort(length(still), " entr(y/ies) survived removal in ", absolute,
            "; check permissions rather than assuming the tree is clean.")
}
tsf_log("Clean: ", absolute)
