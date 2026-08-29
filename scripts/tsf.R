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
  source(file.path("R", "utils_io.R"))
  tsf_load_all("R")
})

USAGE <- "
tsf -- TissueSpectF pipeline

Usage:
  ./tsf <command> [<dataset>...] [options]

Commands:
  fetch       download the GEO inputs the configs declare
  check       verify the GEO inputs are where the configs expect them
  ingest      GEO -> common format (labels, expression, genes)
  spectra     FFT per chromosome, per sample and per condition
  maxt        per-sample permutation test              [the slow stage]
  condition   condition-level test on the summary signal  [~1/n the cost]
  consensus   characteristic spectrum of each condition (power, prevalence, PLV)
  clean       CLEAN decomposition: components by EBIC, no threshold
  stability   stable peaks per condition + peak tables
  peaks       gene-level reconstruction of every stable peak
  compare     constant signature, transitions, cross-dataset replication
  window      spectral window: what the gap pattern alone can produce
  reference   build the fingerprint library + out-of-cohort validation
  match       identify one expression profile against the reference
  app         open the local desktop app in a browser (needs shiny)
  run         all of the above, in order
  status      what each dataset has on disk
  selfcheck   run the full pipeline on synthetic data with a known peak

Options take `--key value` or `--key=value`; names are case-insensitive and
`-` and `_` are interchangeable (--results-dir = --RESULTS_DIR).

Scope:
  --from <stage>        start at this stage            (run only)
  --to <stage>          stop after this stage          (run only)
  --cond F2,F3          restrict to these conditions
  --branch median       one branch (default: both)
  --force               recompute even if output exists
  --dry-run             print the plan, do nothing

Paths:
  --geo-dir <dir>       raw GEO downloads
  --interim-dir <dir>   the common format
  --results-dir <dir>   everything downstream

Parameters (override config/project.R):
  --gene-universe <re>  biotypes on the grid, e.g. '^protein-coding$'
  --maxt-b <n>          permutations, per-sample maxT
  --condition-b <n>     permutations, condition-level test
  --stable-frac <f>     consistency threshold
  --primary-scheme <s>  full | all
  --criterion <s>       condition | consistency
  --ebic-gamma <f>      CLEAN selection penalty
  --k-max <n>           frequencies per chromosome in a fingerprint
  --target <s>          condition | tissue   (reference only)

Matching:
  --query <file>        counts TSV to identify
  --reference <file>    reference .rds to match against
  --input-unit <u>      counts (default) | cpm | tpm | logged

Other:
  --log <file>          append all output to this file
  --help                this message

Stages, in order: ingest, spectra, maxt, condition, clean, stability, peaks,
compare. Precedence: command line > environment (TSF_*) > config/project.R.
"

# Options accept both `--key value` and `--key=value`, and key names are
# matched case-insensitively with `-` and `_` interchangeable, so --results-dir,
# --results_dir and --RESULTS_DIR are the same flag. Every path and tuning
# parameter that used to require an environment variable is a flag; the
# environment variables still work, and the precedence is
#
#     command line  >  environment  >  config/project.R
#
# so a flag never has to fight a stale export.

OPTION_ALIASES <- c(
  from = "from", to = "to", cond = "cond", condition = "cond",
  branch = "branch", log = "log", query = "query", reference = "reference",
  geodir = "geo_dir", geo_dir = "geo_dir",
  interimdir = "interim_dir", interim_dir = "interim_dir",
  resultsdir = "results_dir", results_dir = "results_dir",
  geneuniverse = "gene_universe", gene_universe = "gene_universe",
  maxtb = "maxt_b", maxt_b = "maxt_b", b = "maxt_b",
  conditionb = "condition_b", condition_b = "condition_b",
  stablefrac = "stable_frac", stable_frac = "stable_frac",
  primaryscheme = "primary_scheme", primary_scheme = "primary_scheme",
  criterion = "criterion", ebicgamma = "ebic_gamma", ebic_gamma = "ebic_gamma",
  kmax = "k_max", k_max = "k_max", target = "target",
  cores = "cores",
  inputunit = "input_unit", input_unit = "input_unit", unit = "input_unit"
)

FLAG_ALIASES <- c(force = "force", dryrun = "dry_run", dry_run = "dry_run",
                  help = "help", h = "help", persample = "per_sample",
                  per_sample = "per_sample")

normalise_key <- function(k) tolower(gsub("-", "_", k))

parse_cli <- function(args) {
  opts <- list(); positional <- character(0); i <- 1L
  while (i <= length(args)) {
    a <- args[i]
    if (!startsWith(a, "-")) { positional <- c(positional, a); i <- i + 1L; next }

    key <- normalise_key(sub("^--?", "", sub("=.*$", "", a)))
    inline <- if (grepl("=", a, fixed = TRUE)) sub("^[^=]*=", "", a) else NULL

    if (key %in% names(FLAG_ALIASES) && is.null(inline)) {
      opts[[FLAG_ALIASES[[key]]]] <- TRUE; i <- i + 1L; next
    }
    if (!key %in% names(OPTION_ALIASES)) {
      cat(USAGE); tsf_abort("Unknown option: ", a)
    }
    value <- inline
    if (is.null(value)) {
      if (i + 1L > length(args) || startsWith(args[i + 1L], "--")) {
        cat(USAGE); tsf_abort("Option --", key, " needs a value")
      }
      value <- args[i + 1L]; i <- i + 1L
    }
    opts[[OPTION_ALIASES[[key]]]] <- value
    i <- i + 1L
  }

  # R's `$` does partial matching on lists, so opt$cond would silently resolve to
  # opt$condition_b when no exact "cond" element exists. Every option name is
  # therefore pre-created as an explicit NULL entry, which makes `$` exact.
  all_names <- unique(c(unname(OPTION_ALIASES), unname(FLAG_ALIASES)))
  full <- stats::setNames(vector("list", length(all_names)), all_names)
  full[names(opts)] <- opts
  opts <- full

  branch <- opts$branch
  c(list(
    command  = if (length(positional)) positional[1] else "help",
    datasets = if (length(positional) > 1) positional[-1] else character(0),
    branches = if (is.null(branch)) c("average", "median") else branch,
    force    = isTRUE(opts$force),
    dry_run  = isTRUE(opts$dry_run),
    help     = isTRUE(opts$help)
  ), opts[setdiff(names(opts), c("force", "dry_run", "help"))])
}

#' Apply command-line overrides to the loaded project config.
apply_cli_overrides <- function(project, opt) {
  set <- function(path, value, cast = identity) {
    if (is.null(value)) return(invisible(NULL))
    v <- cast(value)
    if (length(path) == 1L) project[[path]] <<- v
    else project[[path[1]]][[path[2]]] <<- v
    tsf_log("override: ", paste(path, collapse = "$"), " = ", value)
  }
  set("geo_dir", opt$geo_dir)
  set("interim_dir", opt$interim_dir)
  set("results_dir", opt$results_dir)
  set("gene_universe", opt$gene_universe)
  set("stability_criterion", opt$criterion)
  set(c("maxt", "B"), opt$maxt_b, as.integer)
  set(c("maxt", "condition_B"), opt$condition_b, as.integer)
  set(c("maxt", "stable_frac"), opt$stable_frac, as.numeric)
  set(c("maxt", "primary_scheme"), opt$primary_scheme)
  set(c("clean", "ebic_gamma"), opt$ebic_gamma, as.numeric)
  set(c("clean", "per_sample"), opt$per_sample, isTRUE)
  set(c("fingerprint", "k_max"), opt$k_max, as.integer)
  set(c("fingerprint", "target"), opt$target)
  project
}

opt <- parse_cli(commandArgs(trailingOnly = TRUE))
if (opt$help || opt$command %in% c("help", "--help")) { cat(USAGE); quit(status = 0L) }

known <- c(stage_names, "check", "run", "status", "selfcheck", "window",
           "fetch", "reference", "match", "app")
if (!opt$command %in% known) {
  cat(USAGE)
  tsf_abort("Unknown command: ", opt$command)
}
if (!all(opt$branches %in% c("average", "median"))) {
  tsf_abort("--branch must be average or median")
}

if (!is.null(opt$log)) {
  ensure_dir(dirname(opt$log))
  con <- file(opt$log, open = "at")
  sink(con, split = TRUE); sink(con, type = "message")
  on.exit({ sink(type = "message"); sink(); close(con) }, add = TRUE)
}

project <- apply_cli_overrides(load_project_config("config/project.R"), opt)

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

if (opt$command == "app") {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    tsf_abort("The app needs shiny. Install it once with:\n",
              "  Rscript -e 'install.packages(\"shiny\")'\n",
              "Everything else in TissueSpectF works without it.")
  }
  # Hand the resolved path to the app, so --results-dir / --reference on the
  # command line reach it. The app must not re-read config/project.R: doing so
  # made it look on the cluster default while the CLI pointed elsewhere.
  ref_path <- opt$reference %||%
    file.path(project$results_dir, "reference", "reference.rds")
  Sys.setenv(TSF_APP_REFERENCE = normalizePath(ref_path, mustWork = FALSE))
  if (!file.exists(ref_path)) {
    tsf_warn("No reference at ", ref_path,
             " -- the app will say so. Build one with ./tsf reference")
  } else {
    tsf_log("Reference: ", ref_path)
  }
  tsf_log("Starting the app. It runs locally; nothing leaves this machine.")
  shiny::runApp("app", launch.browser = TRUE)
  quit(status = 0L)
}

if (opt$command == "match") {
  source("scripts/match_query.R")
  quit(status = run_match(project, opt))
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
if (!identical(opt$branches, c("average", "median")))
  tsf_log("branch: ", paste(opt$branches, collapse = ", "))

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
