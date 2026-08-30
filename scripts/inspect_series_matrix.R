#!/usr/bin/env Rscript
# Print what a GEO series matrix actually contains, before writing rules for it.
#
#   Rscript scripts/inspect_series_matrix.R data/GSE276114_series_matrix.txt.gz
#   Rscript scripts/inspect_series_matrix.R <file> "etiology" "fibrosis"
#
# The two optional arguments cross-tabulate two fields against each other, which
# is the one view that settles whether a class is what its name suggests. It is
# how GSE162694's 31 "normal liver histology" samples turned out to have NAS = 0
# while its 35 stage-0 samples had NAS 1-5 -- two groups the pipeline had been
# merging.
#
# Write no config until this output is in front of you.

source("R/utils_io.R")
source("R/config.R")
source("R/labels.R")
source("R/grid.R")
source("R/ingest.R")

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) tsf_abort("Usage: inspect_series_matrix.R <series_matrix> [fieldA] [fieldB]")
path <- args[1]
if (!file.exists(path)) tsf_abort("No such file: ", path)

pheno <- parse_series_matrix_pheno(path)
tsf_log(basename(path), ": ", nrow(pheno), " samples, ", ncol(pheno), " fields")

cat("\n--- fields ---\n")
for (nm in colnames(pheno)) {
  v <- clean_pheno_value(pheno[[nm]])
  u <- sort(unique(v[!is.na(v)]))
  n_na <- sum(is.na(v))
  if (length(u) <= 15) {
    tab <- table(v, useNA = "no")
    cat(sprintf("%-34s %s%s\n", nm,
                paste(sprintf("%s=%d", names(tab), as.integer(tab)), collapse = "  "),
                if (n_na) sprintf("  [NA=%d]", n_na) else ""))
  } else {
    cat(sprintf("%-34s (%d distinct) e.g. %s\n", nm, length(u),
                paste(utils::head(u, 4), collapse = ", ")))
  }
}

if (length(args) >= 3) {
  a <- find_pheno_column(pheno, args[2]); b <- find_pheno_column(pheno, args[3])
  if (is.null(a) || is.null(b)) {
    tsf_abort("Could not find both fields: '", args[2], "' and '", args[3], "'")
  }
  cat("\n--- ", a, " x ", b, " ---\n", sep = "")
  print(table(clean_pheno_value(pheno[[a]]), clean_pheno_value(pheno[[b]]),
              useNA = "ifany"))
}

cat("\nWrite the config from this, not from the paper's description of it.\n")
