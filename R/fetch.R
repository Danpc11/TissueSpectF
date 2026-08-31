# fetch.R -- download the GEO inputs a dataset config declares.
#
# Exists so a run can start from nothing (a fresh Colab VM, a new machine)
# without anyone hand-copying files. Downloads only what the configs name, skips
# what is already present, and never deletes.

#' Candidate URLs for a supplementary file.
#'
#' GEO names RNA-seq count tables inconsistently across series (raw_counts_GRCh38,
#' norm_counts_TPM, or the bare accession), so the config names the file and
#' this only builds the location.
geo_series_url <- function(gse, file) {
  stub <- paste0(substr(gse, 1, nchar(gse) - 3), "nnn")
  if (grepl("series_matrix", file)) {
    sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s", stub, gse, file)
  } else {
    sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/%s", stub, gse, file)
  }
}

download_if_missing <- function(url, dest, force = FALSE) {
  if (!force && file.exists(dest) && file.size(dest) > 0) {
    sz <- file.size(dest)
    tsf_log("  have: ", basename(dest), " (",
            if (sz >= 1024^2) paste0(round(sz / 1024^2, 1), " MB")
            else if (sz >= 1024) paste0(round(sz / 1024), " KB")
            else paste0(sz, " B"), ")")
    return(invisible(TRUE))
  }
  ensure_dir(dirname(dest))
  tsf_log("  get:  ", basename(dest))
  ok <- tryCatch({
    utils::download.file(url, dest, mode = "wb", quiet = TRUE)
    file.exists(dest) && file.size(dest) > 0
  }, error = function(e) { tsf_warn("    failed: ", conditionMessage(e)); FALSE })
  if (!ok && file.exists(dest)) unlink(dest)
  invisible(ok)
}

#' Download counts, series matrices and the annotation for every configured dataset.
stage_fetch <- function(project, opt) {
  ensure_dir(project$geo_dir)
  ok <- TRUE
  for (id in stage_datasets(opt)) {
    cfg <- load_dataset_config(id)
    tsf_log(cfg$id, ":")
    ok <- download_if_missing(geo_series_url(cfg$id, cfg$series_matrix),
                              file.path(project$geo_dir, cfg$series_matrix),
                              opt$force) && ok
    ok <- download_if_missing(geo_series_url(cfg$id, cfg$counts_file),
                              file.path(project$geo_dir, cfg$counts_file),
                              opt$force) && ok
  }
  annot <- file.path(project$geo_dir, project$annotation_file)
  if (!file.exists(annot)) {
    tsf_log("annotation:")
    # The NCBI-generated annotation ships with the RNA-seq counts of any GEO
    # series that has them; take it from the first configured dataset.
    first <- load_dataset_config(stage_datasets(opt)[1])
    ok <- download_if_missing(geo_series_url(first$id, project$annotation_file),
                              annot, opt$force) && ok
  } else {
    tsf_log("annotation: have ", project$annotation_file)
  }
  if (!ok) tsf_warn("Some downloads failed; run ./tsf check to see what is missing")
  sprintf("%d dataset(s) fetched", length(stage_datasets(opt)))
}
