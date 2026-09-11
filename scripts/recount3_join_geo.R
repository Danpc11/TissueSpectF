#!/usr/bin/env Rscript
# recount3_join_geo.R -- give a recount3 SRA cohort its GEO labels.
#
# Usage:
#   Rscript scripts/recount3_join_geo.R --dataset R3_SRP197353 --geo GSE130970 --geo-dir data
#
# recount3's SRA metadata does not always carry the GEO sample characteristics
# (for SRP197353 `sra.sample_attributes` is empty), so the exploded pheno has no
# "fibrosis stage" column to write rules on. The labels exist, in the GEO
# series matrix the GEO config already uses, keyed by GSM; and each GSM names
# its SRA experiment (SRX) in !Sample_relation. This joins the two on SRX
# (falling back to the sample title), adds every series-matrix column to the
# recount3 pheno, and rewrites config/datasets/R3_<SRP>.R so that its
# condition_rules, covariate_columns, exclude_samples, tissue and vocabulary
# are the GEO config's, verbatim. Same samples, same labels, same rules --
# only the quantification differs, which is the point of the recount3 path.
#
# The GEO config's rules stay the single source of truth: this script copies
# them, it does not fork them. Re-run it after editing the GEO config.
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
ds  <- flag("--dataset", ""); geo <- flag("--geo", "")
geo_dir <- flag("--geo-dir", Sys.getenv("TSF_GEO_DIR", ""))
if (!nzchar(ds) || !nzchar(geo) || !nzchar(geo_dir)) tsf_abort("Pasa --dataset R3_SRP... --geo GSExxxxx --geo-dir <dir>")

cfg_r3_path  <- file.path("config", "datasets", paste0(ds, ".R"))
cfg_geo_path <- file.path("config", "datasets", paste0(geo, ".R"))
for (f in c(cfg_r3_path, cfg_geo_path)) if (!file.exists(f)) tsf_abort("missing ", f)
cfg_r3  <- source(cfg_r3_path, local = TRUE)$value
cfg_geo <- source(cfg_geo_path, local = TRUE)$value

pheno_path <- file.path(geo_dir, cfg_r3$metadata_file %||% paste0(ds, "_pheno.tsv"))
pheno <- utils::read.delim(pheno_path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE,
                           na.strings = c("", "NA"), quote = "")
sm_path <- file.path(geo_dir, cfg_geo$series_matrix)
if (!file.exists(sm_path)) tsf_abort("series matrix not found: ", sm_path, " (run ./tsf fetch ", geo, ")")
sm <- read_series_pheno(sm_path)
tsf_log(geo, ": series matrix with ", nrow(sm), " sample(s), ", ncol(sm), " column(s)")

# SRX per GSM: any column whose values contain an SRX accession
srx_of <- function(df) {
  out <- rep(NA_character_, nrow(df))
  for (col in names(df)) {
    has <- grepl("SRX[0-9]+", as.character(df[[col]]))
    if (any(has)) out[has & is.na(out)] <- regmatches(as.character(df[[col]])[has & is.na(out)],
                                                      regexpr("SRX[0-9]+", as.character(df[[col]])[has & is.na(out)]))
  }
  out
}
sm$srx <- srx_of(sm)
gsm_col <- grep("^geo_accession$", colnames(sm), value = TRUE)[1]
if (is.na(gsm_col)) gsm_col <- colnames(sm)[grepl("^GSM", as.character(sm[[1]]))][1] %||% NA
title_col <- grep("^title$", colnames(sm), value = TRUE)[1]

# join: SRX first, title second
idx <- match(pheno$experiment, sm$srx)
via <- rep("srx", nrow(pheno)); via[is.na(idx)] <- NA
if (any(is.na(idx)) && !is.na(title_col) && "sample_title" %in% colnames(pheno)) {
  alt <- match(trimws(pheno$sample_title), trimws(as.character(sm[[title_col]])))
  fill <- is.na(idx) & !is.na(alt)
  idx[fill] <- alt[fill]; via[fill] <- "title"
}
n_ok <- sum(!is.na(idx))
tsf_log(ds, ": ", n_ok, "/", nrow(pheno), " recount3 sample(s) matched to a GSM (",
        sum(via == "srx", na.rm = TRUE), " by SRX, ", sum(via == "title", na.rm = TRUE), " by title)")
if (n_ok < 0.9 * nrow(pheno)) {
  tsf_abort("fewer than 90% matched. Inspect: head data/", basename(pheno_path),
            " and the !Sample_relation lines of ", basename(sm_path))
}
if (any(duplicated(idx[!is.na(idx)]))) tsf_abort("two recount3 runs map to the same GSM; a GSM with several runs must be summed before this step")

# add every series-matrix column; keep recount3's own first
sm_cols <- setdiff(colnames(sm), "srx")
add <- sm[idx, sm_cols, drop = FALSE]
colnames(add) <- ifelse(colnames(add) %in% colnames(pheno), paste0("geo_", colnames(add)), colnames(add))
out <- cbind(pheno, add, matched_via = via, stringsAsFactors = FALSE)
if (!is.na(gsm_col)) out$donor <- ifelse(is.na(idx), out$donor, as.character(sm[[gsm_col]])[idx])
out$condition <- NULL   # the rules assign it
out <- out[!is.na(idx), , drop = FALSE]
utils::write.table(out, pheno_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
tsf_log("wrote ", pheno_path, " (", nrow(out), " sample(s), ", ncol(out), " column(s))")

# rewrite the R3 config: header from recount3, rules & friends from the GEO config
geo_txt <- readLines(cfg_geo_path, warn = FALSE)
grab <- function(field) {
  i <- grep(paste0("^\\s*", field, "\\s*="), geo_txt)
  if (!length(i)) return(NULL)
  i <- i[1]
  # a field ends where the next top-level field starts (two-space indent, name =)
  nxt <- grep("^  [a-z_]+\\s*=", geo_txt); nxt <- nxt[nxt > i]
  j <- if (length(nxt)) nxt[1] - 1L else length(geo_txt) - 1L
  block <- geo_txt[i:j]
  block <- sub("\\s*$", "", block); block <- block[nzchar(trimws(block)) | seq_along(block) < length(block)]
  # make sure the block ends with a comma so it can sit mid-list
  last <- max(which(nzchar(trimws(block))))
  if (!grepl(",\\s*$", block[last])) block[last] <- paste0(block[last], ",")
  block
}
copied <- unlist(lapply(c("covariate_columns", "condition_rules", "exclude_samples"), grab))
tissue <- cfg_geo$tissue %||% cfg_r3$tissue
vocab  <- cfg_geo$vocabulary %||% cfg_r3$vocabulary
invisible(file.copy(cfg_r3_path, paste0(cfg_r3_path, ".bak"), overwrite = TRUE))
txt <- c(
  sprintf("# %s: recount3 quantification of %s (%s). Header from R/recount3.R;", ds, geo, cfg_r3$description %||% ""),
  sprintf("# labels, covariates, exclusions and rules copied VERBATIM from %s by", basename(cfg_geo_path)),
  "# scripts/recount3_join_geo.R, which joined the series matrix on SRX. The GEO",
  "# config is the source of truth for the labels; re-run the join after editing it.",
  "list(",
  sprintf('  id            = "%s",', ds),
  sprintf('  tissue        = "%s",', tissue),
  sprintf('  vocabulary    = "%s",', vocab),
  sprintf('  description   = "recount3 sra/%s = %s (%s)",', sub("^R3_", "", ds), geo, gsub('"', "'", cfg_geo$description %||% "")),
  sprintf('  counts_file   = "%s",', cfg_r3$counts_file),
  '  source        = "matrix",',
  sprintf('  metadata_file = "%s",', basename(pheno_path)),
  '  sample_id_column = "sample_id",',
  '  counts_spec   = list(sep = "\\t", id_column = "Name"),',
  '  count_id_type = "ENSEMBL",',
  '  expression_unit = "counts",',
  sprintf('  has_control_cohort = %s,', if (isTRUE(cfg_geo$has_control_cohort)) "TRUE" else "FALSE"),
  '  donor_column  = "donor",',
  sprintf('  geo_series    = "%s",', geo),
  copied,
  sprintf('  notes = "Same samples and labels as %s; only the quantification (recount3/Monorail) differs."', geo),
  ")"
)
writeLines(txt, cfg_r3_path)
chk <- tryCatch(source(cfg_r3_path, local = TRUE)$value, error = function(e) NULL)
if (is.null(chk)) { file.copy(paste0(cfg_r3_path, ".bak"), cfg_r3_path, overwrite = TRUE)
  tsf_abort("generated config does not parse; restored the backup. Copy the rules by hand from ", cfg_geo_path) }
tsf_log("wrote ", cfg_r3_path, " with ", length(chk$condition_rules), " condition rule(s) from ", geo)
tsf_log("next: Rscript scripts/validate_recount3_configs.R --datasets ", ds, " --vocabulary ", vocab,
        "  then re-run run_differential.sh")
