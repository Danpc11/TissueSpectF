#!/usr/bin/env Rscript
# build_tissue_reference.R -- median asinh(TPM) per gene over reference samples.
#
# Usage:
#   Rscript scripts/build_tissue_reference.R --datasets R3_LIVER \
#       --out $TSF_INTERIM_DIR/reference_profile_liver.tsv [--tissue liver] [--condition ...]
#
# The output is the profile apply_reference_profile.R subtracts from every
# cohort, and the one TSF_REFERENCE_PROFILE must name when a query is matched.
if (any(commandArgs(TRUE) %in% c("-h", "--help"))) {
  cat(paste(sub("^# ?", "", grep("^#", readLines(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))[-1], value = TRUE)), collapse = "\n"), "\n")
  quit(save = "no")
}
suppressWarnings({ source("R/utils_io.R"); tsf_load_all("R") })
project <- load_project_config(Sys.getenv("TSF_CONFIG", "config/project.R"))
args <- commandArgs(trailingOnly = TRUE)
flag <- function(n, d = NULL) {
  h <- grep(paste0("^", n, "="), args, value = TRUE)
  if (length(h)) return(sub(paste0("^", n, "="), "", h[1]))
  i <- match(n, args); if (!is.na(i) && length(args) > i) args[i + 1] else d
}
datasets <- trimws(strsplit(flag("--datasets", ""), ",")[[1]])
if (!length(datasets)) tsf_abort("Pasa --datasets R3_LIVER")
out <- flag("--out", file.path(project$interim_dir, "reference_profile.tsv"))
cond <- flag("--condition", NULL)
tissue <- flag("--tissue", NULL)

ref <- build_reference_profile(datasets, project, condition = cond, tissue = tissue)
write_reference_profile(ref, out, project)
tsf_log("reference profile: ", nrow(ref), " positions from ", ref$n_ref_samples[1],
        " sample(s) of ", paste(datasets, collapse = ","), " -> ", out,
        " (+ ", basename(sub("\\.tsv$", "_manifest.tsv", out)), ", md5 ",
        substr(unname(tools::md5sum(out)), 1, 8), ")")
