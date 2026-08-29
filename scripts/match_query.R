# match_query.R -- identify one expression profile against a reference.
#
#   ./tsf match --query counts.tsv
#   ./tsf match --query counts.tsv --reference /path/to/reference.rds
#
# The reference is self-contained: it carries its own grid, identifiers,
# frequency ceiling and feature representation. Nothing here loads a previously
# ingested dataset, so a reference.rds is portable to a machine that has never
# seen the cohorts it was built from.

run_match <- function(project, opt) {
  ref_path <- opt$reference %||%
    file.path(project$results_dir, "reference", "reference.rds")
  if (!file.exists(ref_path)) {
    tsf_abort("No reference at ", ref_path, ". Build one first: ./tsf reference")
  }
  ref <- readRDS(ref_path)
  if (is.null(ref$grid)) {
    tsf_abort("This reference predates self-contained references (no grid). ",
              "Rebuild it with ./tsf reference.")
  }

  if (is.null(opt$query)) tsf_abort("Give a query: --query <counts.tsv>")
  if (!file.exists(opt$query)) tsf_abort("No such file: ", opt$query)

  q <- read_tsv_tsf(opt$query)
  id_col <- colnames(q)[1]
  value_cols <- setdiff(colnames(q), id_col)
  if (!length(value_cols)) tsf_abort("Query has no value column")
  ids <- sub("\\..*$", "", as.character(q[[id_col]]))

  tsf_log("Query: ", basename(opt$query), " (", nrow(q), " rows, ",
          length(value_cols), " sample column(s))")
  tsf_log("Reference: ", basename(ref_path), " | grid ", nrow(ref$grid),
          " genes | k_max ", ref$params$k_max, " | ", ref$params$features)

  cat("\n", strrep("-", 72), "\n", sep = "")
  cat("REFERENCE: ", reference_status(ref), "\n", sep = "")
  calib <- ref$validation$calibration
  if (!is.null(calib) && !is.na(calib$separability_auc)) {
    cat(sprintf("Rejection rule: similarity below the per-class %.0fth percentile
                 of correct held-out matches is reported UNKNOWN. Correct vs
                 incorrect matches separate with AUC %.2f.\n",
                100 * calib$quantile, calib$separability_auc))
  }
  cat(strrep("-", 72), "\n\n", sep = "")

  for (col in value_cols) {
    fq <- fingerprint_query(q[[col]], ids, ref)
    if (is.null(fq)) {
      tsf_warn(col, ": no usable positions on the reference grid")
      next
    }
    if (fq$coverage < 0.2) {
      cat(col, ": only ", round(100 * fq$coverage, 1),
          "% of the reference grid is present. Too little to score.\n", sep = "")
      next
    }
    proj <- project_to_reference(fq$vector, ref)
    res <- match_query(ref$model, proj$vector)
    if (is.null(res)) { tsf_warn(col, ": could not be scored"); next }
    res <- apply_rejection(res, calib)

    cat(col, "\n")
    if (identical(res$decision, "UNKNOWN")) {
      cat("  UNKNOWN -- outside the domain of this reference.\n")
      cat(sprintf("  Closest class was %s at similarity %.3f, below the %.3f
                   threshold calibrated for it.\n",
                  res$best, res$similarity, res$threshold %||% NA_real_))
    }
    for (i in seq_len(min(5, nrow(res$scores)))) {
      cat(sprintf("  %-14s similarity %.3f\n",
                  res$scores$class[i], res$scores$similarity[i]))
    }
    cat(sprintf("  margin %.3f | p(shuffled) %.3f | grid coverage %.1f%% (%s ids) | %d/%d features\n",
                res$margin, res$p_shuffle, 100 * fq$coverage, fq$id_type,
                proj$n_shared, proj$n_features))

    if (is.null(ref$validation)) {
      cat("  -> Uncalibrated reference: a suggestion, not an identification.\n")
    } else if (ref$validation$accuracy - ref$validation$baseline <= 0.02) {
      cat("  -> This reference does not beat its baseline out of cohort.\n")
    } else if (!is.na(res$margin) && res$margin < 0.02 &&
               !identical(res$decision, "UNKNOWN")) {
      cat("  -> Top two classes within 0.02: the call is not separable.\n")
    }
    cat("\n")
  }
  0L
}
