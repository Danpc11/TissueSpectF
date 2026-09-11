#!/usr/bin/env Rscript
# apply_reference_profile.R -- switch datasets to the differential signal.
#
# Usage:
#   Rscript scripts/apply_reference_profile.R --datasets A,B,C --profile <tsv>
#   Rscript scripts/apply_reference_profile.R --datasets A,B,C --restore
#
# Rewrites <interim>/<dataset>/expression.tsv as (asinh TPM - profile median);
# keeps expression_raw.tsv; idempotent. Run AFTER ingest and BEFORE spectra.
# Every dataset that enters the same library must be corrected with the same
# profile -- including the reference dataset itself, whose deviations are the
# "healthy" class.
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
# --geo-dir / --interim-dir / --results-dir override the environment, same as ./tsf
for (pair in list(c("--geo-dir", "TSF_GEO_DIR"), c("--interim-dir", "TSF_INTERIM_DIR"),
                  c("--results-dir", "TSF_RESULTS_DIR"), c("--config", "TSF_CONFIG"))) {
  v <- flag(pair[1], NULL); if (!is.null(v)) do.call(Sys.setenv, stats::setNames(list(v), pair[2]))
}
project <- load_project_config(Sys.getenv("TSF_CONFIG", "config/project.R"))
datasets <- trimws(strsplit(flag("--datasets", ""), ",")[[1]])
if (!length(datasets)) tsf_abort("Pasa --datasets A,B,C")
if ("--restore" %in% args) {
  for (d in datasets) { restore_raw_expression(d, project); tsf_log(d, ": restored") }
  quit(save = "no")
}
profile <- flag("--profile", "")
if (!nzchar(profile)) tsf_abort("Pasa --profile <reference_profile.tsv> o --restore")
ref <- read_reference_profile(profile)
profile <- normalizePath(profile)
for (d in datasets) apply_reference_profile(d, project, ref, profile_name = basename(profile),
                                            profile_path = profile)
tsf_log("Now run ./tsf run --from spectra ... then ./tsf reference. The reference stage ",
        "stores this profile (md5 ", substr(unname(tools::md5sum(profile)), 1, 8),
        ") inside reference.rds and `match` corrects queries from there.")
