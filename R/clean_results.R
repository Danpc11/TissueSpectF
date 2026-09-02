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
#   Rscript scripts/clean_results.R [--results-dir <dir>] [--dry-run] [--force]
#
# Deleting a non-empty tree asks first. config/project.R points results_dir at
# the active run tree, so the default target is real work -- a consensus stage
# is tens of minutes and is not reproducible from what is left behind. An
# interactive session gets a prompt; a script or a CI job has no one to answer,
# so there --force is required and the default is to refuse.
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
forced <- "--force" %in% args

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

# Size is reported before the question, because "3 entries" and "3 entries,
# 44 GB" deserve different answers and the count alone hides which one this is.
bytes <- suppressWarnings(sum(file.size(list.files(
  absolute, recursive = TRUE, all.files = TRUE, full.names = TRUE)),
  na.rm = TRUE))
size <- if (is.finite(bytes) && bytes > 0) {
  sprintf(" (%.1f GB)", bytes / 1024^3)
} else ""

if (!forced) {
  cat(sprintf("\nAbout to delete %d entr%s%s from:\n  %s\n\n",
              length(entries), if (length(entries) == 1L) "y" else "ies",
              size, absolute))
  cat("This is the tree config/project.R points at.\n")
  cat("Type the directory name to confirm: ")

  # Read from the terminal rather than testing for one. isatty() does not
  # distinguish reliably under Rscript, and guessing wrong in either direction
  # is bad: a false negative refuses a legitimate `make clean`, a false
  # positive deletes a tree in a CI job with nobody to answer. Reading and
  # requiring a specific answer handles all three cases without a guess --
  # a closed or empty stdin yields nothing and is treated as a refusal.
  answer <- tryCatch(readLines(con = "stdin", n = 1L, warn = FALSE),
                     error = function(e) character(0))
  answer <- if (length(answer)) trimws(answer[1]) else ""

  if (!identical(answer, basename(absolute))) {
    cat("\n")
    tsf_log("Not confirmed: expected \"", basename(absolute), "\", got \"",
            answer, "\". Nothing was deleted.")
    tsf_log("Use --dry-run to list what would go, or --force to skip this prompt.")
    quit(save = "no", status = 1L)
  }
}

tsf_log("Removing ", length(entries), " entr(y/ies) from ", absolute)
unlink(entries, recursive = TRUE, force = TRUE)

still <- list.files(absolute, all.files = TRUE, no.. = TRUE)
if (length(still)) {
  tsf_abort(length(still), " entr(y/ies) survived removal in ", absolute,
            "; check permissions rather than assuming the tree is clean.")
}
tsf_log("Clean: ", absolute)
