#!/usr/bin/env Rscript
# resolve_recount3.R -- GSE -> SRP via ENA, then "is it in recount3?".
#
# Usage:
#   Rscript scripts/resolve_recount3.R --gse GSE135251,GSE130970,GSE162694,GSE276114,GSE142530
#
# Writes config/recount3_sources.tsv: gse, srp, in_recount3, gene_sums_url.
# A cohort with in_recount3 = FALSE stays on the GEO path (./tsf fetch/ingest);
# the write-up must then say its quantification differs from the reference's.
#
# ENA's portal API accepts GEO series accessions and returns the SRA study
# (secondary_study_accession). recount3 holds SRA studies made public before
# its 2019-2020 snapshot; later releases (GSE276114) are not there, and a
# series deposited early but released late may be missing too.
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
gses <- trimws(strsplit(flag("--gse", ""), ",")[[1]])
if (!length(gses)) tsf_abort("Pasa --gse GSE1,GSE2")
out <- flag("--out", "config/recount3_sources.tsv")

ena_srp <- function(gse) {
  u <- sprintf(paste0("https://www.ebi.ac.uk/ena/portal/api/search?result=study",
                      "&query=geo_accession%%3D%%22%s%%22&fields=secondary_study_accession",
                      "&format=tsv"), gse)
  txt <- tryCatch(readLines(url(u), warn = FALSE), error = function(e) character(0))
  if (length(txt) < 2L) return(NA_character_)
  srp <- strsplit(txt[2], "\t")[[1]]
  srp[grepl("^[SED]RP", srp)][1]
}

rows <- lapply(gses, function(g) {
  srp <- ena_srp(g)
  found <- if (is.na(srp)) NA_character_ else recount3_locate(srp)
  url <- if (!is.na(found)) recount3_urls(srp, found)$gene_sums else NA_character_
  tsf_log(g, " -> ", if (is.na(srp)) "no SRP in ENA" else srp,
          " | recount3: ", if (is.na(found)) "NOT FOUND" else found)
  data.frame(gse = g, srp = srp, in_recount3 = !is.na(found),
             gene_sums_url = url, stringsAsFactors = FALSE)
})
tab <- do.call(rbind, rows)
ensure_dir(dirname(out)); write_tsv_tsf(tab, out)
tsf_log("wrote ", out)
hits <- tab$srp[tab$in_recount3 %in% TRUE]
if (length(hits)) tsf_log("fetch with: Rscript scripts/recount3_fetch.R --projects ",
                          paste(hits, collapse = ","), " --tissue liver --vocabulary liver_fibrosis")
