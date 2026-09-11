#!/usr/bin/env Rscript
# config_digest.R -- md5 of the EFFECTIVE project configuration.
#
# Usage:  Rscript scripts/config_digest.R [--show]
#
# Loads config/project.R exactly as the pipeline does (so every TSF_* override
# from the environment is applied), drops the fields that are locations rather
# than method (geo_dir, interim_dir, results_dir), serialises the rest
# deterministically and prints its md5. run_differential.sh puts this in every
# artefact's input digest, so a change to TSF_ESTIMATOR, TSF_PRIMARY_SCHEME,
# TSF_MT_NW, TSF_BIN_AGGREGATE, ... invalidates what was built with the old
# value. Grepping project.R could not see any of those.
suppressWarnings({ source("R/utils_io.R"); tsf_load_all("R") })
project <- load_project_config(Sys.getenv("TSF_CONFIG", "config/project.R"))
drop <- c("geo_dir", "interim_dir", "results_dir", "log_file", "n_workers", "cores")
eff <- project[setdiff(names(project), drop)]
# functions and environments have no stable text form: keep their names only
eff <- rapply(eff, function(x) if (is.function(x) || is.environment(x)) "<fn>" else x, how = "replace")
eff <- eff[order(names(eff))]
txt <- paste(deparse(eff, control = c("keepNA", "keepInteger", "niceNames", "showAttributes")), collapse = "\n")
if ("--show" %in% commandArgs(TRUE)) cat(txt, "\n")
tf <- tempfile(); writeLines(txt, tf)
cat(unname(tools::md5sum(tf)), "\n")
