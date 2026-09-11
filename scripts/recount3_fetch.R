#!/usr/bin/env Rscript
# recount3_fetch.R -- download recount3 projects into the ingest format.
#
# Usage:
#   Rscript scripts/recount3_fetch.R --projects LIVER[,SRP217231,...] \
#       [--tissue liver] [--vocabulary case_control] [--max-samples N]
#
# LIVER (and any other GTEx tissue id: KIDNEY_CORTEX, LUNG, ...) is found under
# recount3/gtex; SRPxxxxxx under recount3/sra. Each project yields
#   $TSF_GEO_DIR/R3_<P>_reads.tsv.gz, R3_<P>_pheno.tsv, config/datasets/R3_<P>.R
# then `./tsf ingest R3_<P>` runs the ordinary path.
#
# GEO series are NOT accepted here on purpose: recount3 is indexed by SRP.
# scripts/resolve_recount3.R turns GSE ids into SRPs and tells you which ones
# recount3 actually holds.
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
projects <- strsplit(flag("--projects", ""), ",")[[1]]
if (!length(projects)) tsf_abort("Pasa --projects LIVER[,SRP...]")
geo_dir <- flag("--geo-dir", Sys.getenv("TSF_GEO_DIR", ""))
if (!nzchar(geo_dir)) tsf_abort("Pasa --geo-dir o define TSF_GEO_DIR")
tissue <- flag("--tissue", NA_character_)
vocab <- flag("--vocabulary", "case_control")
maxn <- flag("--max-samples", NULL); if (!is.null(maxn)) maxn <- as.integer(maxn)

for (p in trimws(projects)) {
  r <- fetch_recount3_project(p, geo_dir = geo_dir, tissue = tissue,
                              vocabulary = vocab, max_samples = maxn)
  tsf_log("ready: ./tsf ingest ", r$id)
}
