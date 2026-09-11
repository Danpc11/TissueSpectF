# reference_profile.R -- the differential spectrum.
#
# THE PROBLEM IT SOLVES
# ---------------------
# The spectrum of asinh(TPM) along a chromosome is dominated by structure every
# human tissue shares: gene density, GC isochores, housekeeping clusters,
# replication domains. That is the "invariant layer" the pipeline finds, and it
# is most of the power. A condition -- fibrosis stage, tumour, anything -- is a
# small modulation on top of that carrier. A classifier fed the raw spectrum
# learns the carrier and the cohort's spectral window, and very little else.
#
# THE FIX
# -------
# Subtract, position by position, a robust reference for the tissue:
#     d_j = asinh(TPM_j)  -  median over reference samples of asinh(TPM_j)
# and run every downstream stage on d. The carrier cancels. What is left is
# "how this sample deviates from a healthy liver along the genome", and its
# spectrum is the spectrum of the deviation. The maxT null becomes far more
# reasonable too: the residual has much less autocorrelation than the raw
# signal, because the autocorrelated part IS the carrier.
#
# WHERE IT PLUGS IN
# -----------------
# The contract between ingest and everything downstream is
# <interim>/<dataset>/expression.tsv. apply_reference_profile() rewrites that
# file (keeping expression_raw.tsv) so spectra / condition / consensus /
# reference / fingerprint run unchanged. A query sample has to go through the
# same subtraction: the profile is STORED in reference.rds by stage_reference
# (with its md5) and fingerprint_query() applies it from there, after the
# query has been put on the library's positions (bins or genes). A raw query
# against a differential library, or the reverse, is refused.
#
# The reference and the cohorts MUST come from the same quantification
# pipeline (recount3 for all of them, see R/recount3.R). Subtracting a GTEx
# median computed with one aligner from a GEO cohort counted with another
# subtracts pipeline as much as biology, and the residual spectrum would be a
# spectrum of the pipeline difference.

#' Build a tissue reference profile from one or more ingested datasets.
#'
#' @param dataset_ids  interim datasets whose samples are the reference
#'   (e.g. "R3_LIVER"). All samples are used unless `condition` is given.
#' @param condition  optional condition level to restrict to.
#' @return data.frame gene_id, ref_median, ref_mad, n_ref
build_reference_profile <- function(dataset_ids, project, condition = NULL,
                                    min_samples = 20L, tissue = NULL) {
  mats <- lapply(dataset_ids, function(id) {
    d <- load_dataset(id, project)
    m <- d$expression
    rownames(m) <- sub("\\..*$", "", rownames(m))   # unversioned Ensembl ids
    m <- m[!duplicated(rownames(m)), , drop = FALSE]
    if (!is.null(condition)) {
      keep <- d$samples$sample_id[as.character(d$samples$condition) %in% condition]
      m <- m[, colnames(m) %in% keep, drop = FALSE]
    }
    m
  })
  genes <- Reduce(intersect, lapply(mats, rownames))
  if (!length(genes)) tsf_abort("reference datasets share no gene ids")
  m <- do.call(cbind, lapply(mats, function(x) x[genes, , drop = FALSE]))
  if (ncol(m) < min_samples) {
    tsf_abort("only ", ncol(m), " reference sample(s); a median over fewer than ",
              min_samples, " is not a tissue reference, it is another cohort")
  }
  ref <- data.frame(
    gene_id    = genes,
    ref_median = apply(m, 1, stats::median, na.rm = TRUE),
    ref_mad    = apply(m, 1, stats::mad, na.rm = TRUE),
    n_ref      = rowSums(is.finite(m)),
    # Provenance as COLUMNS, so it survives the TSV round trip. R attributes
    # do not, and reference_profile_applied.tsv used to record n_ref_samples
    # = NA for that reason.
    n_ref_samples   = ncol(m),
    source_datasets = paste(dataset_ids, collapse = ","),
    tissue          = tissue %||% NA_character_,
    stringsAsFactors = FALSE)
  attr(ref, "n_samples") <- ncol(m)
  attr(ref, "datasets") <- dataset_ids
  ref
}

#' Write the profile plus a sidecar manifest with its digest and provenance.
write_reference_profile <- function(ref, path, project = NULL) {
  ensure_dir(dirname(path))
  write_tsv_tsf(ref, path)
  man <- data.frame(
    key = c("profile", "md5", "n_positions", "n_ref_samples", "source_datasets", "tissue",
            "grid_axis", "bin_size", "annotation_file", "created"),
    value = c(basename(path), unname(tools::md5sum(path)), nrow(ref),
              ref$n_ref_samples[1] %||% NA, ref$source_datasets[1] %||% NA, ref$tissue[1] %||% NA,
              project$grid_axis %||% NA, project$bin_size %||% NA, project$annotation_file %||% NA,
              format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    stringsAsFactors = FALSE)
  write_tsv_tsf(man, sub("\\.tsv$", "_manifest.tsv", path))
  invisible(path)
}

read_reference_profile <- function(path) {
  ref <- read_tsv_tsf(path)
  need <- c("gene_id", "ref_median")
  if (!all(need %in% colnames(ref))) {
    tsf_abort("reference profile lacks columns ", paste(setdiff(need, colnames(ref)), collapse = ", "))
  }
  ref$gene_id <- sub("\\..*$", "", as.character(ref$gene_id))
  ref[!duplicated(ref$gene_id), , drop = FALSE]
}

#' Subtract the profile from a vector or matrix on the asinh(TPM) scale.
#'
#' Genes absent from the profile become NA: the pipeline treats NA as
#' unmeasured (gls_observed drops the position), never as zero deviation.
subtract_reference_profile <- function(x, ref) {
  ids <- if (is.matrix(x)) rownames(x) else names(x)
  if (is.null(ids)) tsf_abort("subtract_reference_profile: x has no gene ids")
  ids <- sub("\\..*$", "", ids)
  r <- ref$ref_median[match(ids, ref$gene_id)]
  out <- if (is.matrix(x)) x - r else x - r     # r recycles down rows / along the vector
  attr(out, "reference_profile") <- TRUE
  attr(out, "n_unmatched") <- sum(is.na(r))
  out
}

#' Rewrite <interim>/<dataset>/expression.tsv as the deviation from the profile.
#'
#' Idempotent: the original is kept as expression_raw.tsv and re-read on each
#' call, so applying a different profile later does not stack subtractions.
#' Genes not in the profile are DROPPED from expression.tsv and genes.tsv (both
#' files must describe the same positions), and the drop is logged.
apply_reference_profile <- function(dataset_id, project, ref, profile_name = "reference",
                                    profile_path = NA_character_) {
  dir <- file.path(project$interim_dir, dataset_id)
  raw_path <- file.path(dir, "expression_raw.tsv")
  expr_path <- file.path(dir, "expression.tsv")
  genes_path <- file.path(dir, "genes.tsv")
  genes_raw_path <- file.path(dir, "genes_raw.tsv")
  applied_path <- file.path(dir, "reference_profile_applied.tsv")
  # WHICH expression.tsv is on disk is decided by content, not by mtime: the
  # marker records the md5 of the corrected file it wrote and of the ingest
  # output it started from. If the current file is the corrected one, the raw
  # copy is still valid; anything else means ingest rewrote it (--force, a new
  # mask, a new annotation) and the raw copy is refreshed from it. Timestamps
  # on cluster filesystems are not reliable enough to carry this decision.
  cur_md5 <- unname(tools::md5sum(expr_path))
  prev <- if (file.exists(applied_path)) read_reference_profile_applied(dir, check = FALSE) else NULL
  still_corrected <- !is.null(prev) && identical(prev$corrected_md5, cur_md5) &&
    file.exists(raw_path) && identical(prev$raw_md5, unname(tools::md5sum(raw_path)))
  if (!still_corrected) {
    file.copy(expr_path, raw_path, overwrite = TRUE)
    file.copy(genes_path, genes_raw_path, overwrite = TRUE)
    if (!is.null(prev)) tsf_log(dataset_id, ": expression.tsv is not the file the last ",
                                "correction wrote; taking it as the new ingest output")
  }
  raw_md5 <- unname(tools::md5sum(raw_path))
  genes_raw_md5 <- unname(tools::md5sum(genes_raw_path))
  ingest_manifest <- file.path(dir, "manifest.tsv")
  ingest_md5 <- if (file.exists(ingest_manifest)) unname(tools::md5sum(ingest_manifest)) else NA_character_
  expr <- read_tsv_tsf(raw_path)
  genes <- read_tsv_tsf(genes_raw_path)
  mat <- as.matrix(expr[, setdiff(colnames(expr), "gene_id"), drop = FALSE])
  rownames(mat) <- sub("\\..*$", "", expr$gene_id)

  d <- subtract_reference_profile(mat, ref)
  in_ref <- rownames(mat) %in% ref$gene_id
  if (!any(in_ref)) {
    tsf_abort(dataset_id, ": no gene of expression.tsv is in the reference profile. ",
              "Id systems differ (the profile must be Ensembl ids, unversioned).")
  }
  d <- d[in_ref, , drop = FALSE]
  genes_keep <- genes[sub("\\..*$", "", genes$gene_id) %in% rownames(d), , drop = FALSE]

  write_tsv_tsf(data.frame(gene_id = rownames(d), d, check.names = FALSE), expr_path)
  write_tsv_tsf(genes_keep, genes_path)
  digest <- if (!is.na(profile_path) && file.exists(profile_path))
    unname(tools::md5sum(profile_path)) else NA_character_
  write_tsv_tsf(data.frame(
    key = c("profile", "profile_path", "profile_digest", "n_ref_samples", "source_datasets",
            "tissue", "n_genes_before", "n_genes_after", "applied",
            "raw_md5", "genes_raw_md5", "corrected_md5", "genes_md5", "ingest_manifest_md5"),
    value = c(profile_name, profile_path, digest,
              ref$n_ref_samples[1] %||% attr(ref, "n_samples") %||% NA,
              ref$source_datasets[1] %||% NA, ref$tissue[1] %||% NA,
              nrow(mat), nrow(d), format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              raw_md5, genes_raw_md5, unname(tools::md5sum(expr_path)),
              unname(tools::md5sum(genes_path)), ingest_md5),
    stringsAsFactors = FALSE), applied_path)
  tsf_log(dataset_id, ": expression.tsv is now the deviation from '", profile_name,
          "' (", nrow(d), "/", nrow(mat), " genes kept; ", sum(!in_ref),
          " not in the profile, dropped)")
  invisible(list(n_before = nrow(mat), n_after = nrow(d)))
}

#' The marker apply_reference_profile() leaves, as a list, or NULL if the
#' dataset is raw: no marker, or expression.tsv / genes.tsv are not the files
#' the marker says it wrote (content check by md5, not mtime). With
#' `check = FALSE` the marker is returned as recorded, for callers that decide
#' themselves.
read_reference_profile_applied <- function(dir, check = TRUE) {
  p <- file.path(dir, "reference_profile_applied.tsv")
  if (!file.exists(p)) return(NULL)
  t <- read_tsv_tsf(p)
  out <- as.list(stats::setNames(as.character(t$value), t$key))
  out$n_ref_samples <- suppressWarnings(as.integer(out$n_ref_samples))
  if (!check) return(out)
  cur <- unname(tools::md5sum(file.path(dir, "expression.tsv")))
  cur_g <- unname(tools::md5sum(file.path(dir, "genes.tsv")))
  if (!identical(out$corrected_md5, cur) || !identical(out$genes_md5, cur_g)) {
    tsf_warn(basename(dir), ": expression.tsv/genes.tsv are not the files the reference profile ",
             "was applied to (md5 differs); treating the dataset as RAW. Re-run ",
             "apply_reference_profile.R.")
    return(NULL)
  }
  out
}

#' Restore expression.tsv / genes.tsv to the ingest output.
restore_raw_expression <- function(dataset_id, project) {
  dir <- file.path(project$interim_dir, dataset_id)
  for (f in c("expression", "genes")) {
    raw <- file.path(dir, paste0(f, "_raw.tsv"))
    if (file.exists(raw)) file.rename(raw, file.path(dir, paste0(f, ".tsv")))
  }
  unlink(file.path(dir, "reference_profile_applied.tsv"))
  invisible(TRUE)
}
