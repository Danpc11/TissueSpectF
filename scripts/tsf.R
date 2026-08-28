#!/usr/bin/env Rscript
# tsf -- one entry point for the whole pipeline.
#
#   ./tsf run                        ingest -> spectra -> maxt -> stability
#                                    -> peaks -> compare
#   ./tsf run --from=stability       resume from a stage
#   ./tsf maxt GSE135251 --cond=F3   one stage, one dataset, one condition
#   ./tsf status                     what exists on disk
#   ./tsf selfcheck                  end-to-end correctness on synthetic data
#
# Run from the repository root.

if (!dir.exists("R") || !dir.exists("config")) {
  stop("Run tsf from the repository root (the directory holding R/ and config/).",
       call. = FALSE)
}
suppressPackageStartupMessages({
  for (f in c("utils_io", "config", "labels", "ingest", "paths", "grid", "spectrum",
              "maxt", "stability", "peaks_genes", "compare", "stages")) {
    source(file.path("R", paste0(f, ".R")))
  }
})

USAGE <- "
tsf -- TissueSpectF pipeline

Usage:
  ./tsf <command> [<dataset>...] [options]

Commands:
  check       verify the GEO inputs are where the configs expect them
  ingest      GEO -> common format (labels, expression, genes)
  spectra     FFT per chromosome, per sample and per condition
  maxt        permutation test for peak significance   [the slow stage]
  stability   stable peaks per condition + peak tables
  peaks       gene-level reconstruction of every stable peak
  compare     constant signature, transitions, cross-dataset replication
  window      spectral window: what the gap pattern alone can produce
  run         all of the above, in order
  status      what each dataset has on disk
  selfcheck   run the full pipeline on synthetic data with a known peak

Options:
  --from=<stage>    start at this stage (run only)
  --to=<stage>      stop after this stage (run only)
  --cond=F2,F3      restrict to these conditions
  --branch=median   restrict to one branch (default: average and median)
  --force           recompute maxT even if output exists
  --dry-run         print the plan without doing anything
  --log=<file>      also append all output to this file
  --help            this message

Stages, in order: ingest, spectra, maxt, stability, peaks, compare
Paths come from config/project.R (override with TSF_GEO_DIR, TSF_INTERIM_DIR,
TSF_RESULTS_DIR).
"

parse_cli <- function(args) {
  flag <- function(name) {
    hit <- grep(paste0("^--", name, "="), args, value = TRUE)
    if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else NULL
  }
  positional <- args[!grepl("^--", args)]
  branch <- flag("branch")
  list(
    command  = if (length(positional)) positional[1] else "help",
    datasets = if (length(positional) > 1) positional[-1] else character(0),
    from     = flag("from"),
    to       = flag("to"),
    cond     = flag("cond"),
    branch   = branch,
    branches = if (is.null(branch)) c("average", "median") else branch,
    force    = any(args == "--force"),
    dry_run  = any(args == "--dry-run"),
    log      = flag("log"),
    help     = any(args %in% c("--help", "-h"))
  )
}

opt <- parse_cli(commandArgs(trailingOnly = TRUE))
if (opt$help || opt$command %in% c("help", "--help")) { cat(USAGE); quit(status = 0L) }

known <- c(stage_names, "check", "run", "status", "selfcheck", "window")
if (!opt$command %in% known) {
  cat(USAGE)
  tsf_abort("Unknown command: ", opt$command)
}
if (!is.null(opt$branch) && !opt$branch %in% c("average", "median")) {
  tsf_abort("--branch must be average or median")
}

if (!is.null(opt$log)) {
  ensure_dir(dirname(opt$log))
  con <- file(opt$log, open = "at")
  sink(con, split = TRUE); sink(con, type = "message")
  on.exit({ sink(type = "message"); sink(); close(con) }, add = TRUE)
}

project <- load_project_config("config/project.R")

# --- commands that are not stages -------------------------------------------
if (opt$command == "check") {
  source("scripts/00_check_inputs.R")
  quit(status = 0L)
}

if (opt$command == "status") {
  st <- pipeline_status(project, opt)
  tsf_log("results_dir: ", project$results_dir)
  print(st, row.names = FALSE)
  quit(status = 0L)
}

if (opt$command == "selfcheck") {
  source("scripts/selfcheck.R")
  quit(status = run_selfcheck())
}

# --- stage selection ---------------------------------------------------------
plan <- if (opt$command == "run") {
  from <- if (is.null(opt$from)) stage_names[1] else opt$from
  to   <- if (is.null(opt$to)) stage_names[length(stage_names)] else opt$to
  for (s in c(from, to)) if (!s %in% stage_names) {
    tsf_abort("Unknown stage '", s, "'. Stages: ", paste(stage_names, collapse = ", "))
  }
  i <- match(from, stage_names); j <- match(to, stage_names)
  if (i > j) tsf_abort("--from comes after --to")
  stage_names[i:j]
} else {
  opt$command
}

datasets <- stage_datasets(opt)
tsf_log("datasets: ", paste(datasets, collapse = ", "))
tsf_log("stages:   ", paste(plan, collapse = " -> "))
if (!is.null(opt$cond))   tsf_log("conditions: ", opt$cond)
if (!identical(opt$branches, c("average", "median"))) tsf_log("branch: ", opt$branch)

if (opt$dry_run) {
  tsf_log("dry run: nothing was executed")
  quit(status = 0L)
}

# --- run ---------------------------------------------------------------------
started <- Sys.time()
summary_rows <- list()

for (stage in plan) {
  tsf_log(strrep("-", 62))
  tsf_log("stage: ", stage)
  t0 <- Sys.time()
  result <- tryCatch(stage_functions[[stage]](project, opt),
                     error = function(e) structure(conditionMessage(e), class = "tsf_error"))
  mins <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)

  if (inherits(result, "tsf_error")) {
    summary_rows[[stage]] <- data.frame(stage = stage, minutes = mins,
                                        result = paste("FAILED:", result),
                                        stringsAsFactors = FALSE)
    tsf_warn("stage '", stage, "' failed: ", result)
    print(do.call(rbind, summary_rows), row.names = FALSE)
    tsf_abort("Pipeline stopped at '", stage, "'.")
  }
  summary_rows[[stage]] <- data.frame(stage = stage, minutes = mins,
                                      result = as.character(result),
                                      stringsAsFactors = FALSE)
  tsf_log("stage '", stage, "' done in ", mins, " min: ", result)
}

tsf_log(strrep("=", 62))
print(do.call(rbind, summary_rows), row.names = FALSE)
tsf_log("total: ", round(as.numeric(difftime(Sys.time(), started, units = "mins")), 2),
        " min | results in ", project$results_dir)
