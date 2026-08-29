# match_query.R -- identify one expression profile against the reference.
#
#   ./tsf match --query=path/to/counts.tsv
#   ./tsf match --query=path/to/counts.tsv --reference=results/reference/reference.rds
#
# The query is a TSV with a gene id column and one column of counts (or several,
# scored one at a time). Gene ids are matched to the reference grid, so the
# query needs no preprocessing beyond being counts of the same annotation.

run_match <- function(project, opt) {
  ref_path <- opt$reference %||%
    file.path(project$results_dir, "reference", "reference.rds")
  if (!file.exists(ref_path)) {
    tsf_abort("No reference at ", ref_path, ". Build one first: ./tsf reference")
  }
  ref <- readRDS(ref_path)

  if (is.null(opt$query)) tsf_abort("Give a query: --query=<counts.tsv>")
  if (!file.exists(opt$query)) tsf_abort("No such file: ", opt$query)

  # The grid must be the one the reference was built on, so it is rebuilt from
  # the same annotation rather than from the query.
  dataset_id <- stage_datasets(opt)[1]
  inp <- tsf_stage_inputs(project, dataset_id)
  genes <- inp$dataset$genes

  q <- read_tsv_tsf(opt$query)
  id_col <- colnames(q)[1]
  value_cols <- setdiff(colnames(q), id_col)
  if (!length(value_cols)) tsf_abort("Query has no value column")
  tsf_log("Query: ", basename(opt$query), " (", nrow(q), " rows, ",
          length(value_cols), " sample column(s))")

  ids <- sub("\\..*$", "", as.character(q[[id_col]]))
  hit <- match(genes$gene_id, ids)
  if (sum(!is.na(hit)) < 0.2 * nrow(genes)) {
    hit <- match(genes$entrez_id %||% rep(NA, nrow(genes)), ids)
  }
  coverage <- mean(!is.na(hit))
  tsf_log("Grid coverage of the query: ", round(100 * coverage, 1), "%")
  if (coverage < 0.2) {
    tsf_abort("Fewer than 20% of grid genes found in the query. Check that the ",
              "identifiers are Ensembl or Entrez gene ids matching the annotation.")
  }

  terms_cache <- fingerprint_terms(inp$chrom_idx)
  status <- reference_status(ref)

  cat("\n", strrep("-", 70), "\n", sep = "")
  cat("REFERENCE: ", status, "\n", sep = "")
  cat(strrep("-", 70), "\n\n", sep = "")

  for (col in value_cols) {
    counts <- suppressWarnings(as.numeric(q[[col]]))[hit]
    counts[!is.finite(counts)] <- 0
    y <- asinh(counts / max(sum(counts, na.rm = TRUE), 1) * 1e6)

    fp <- fingerprint_vector(y, inp$chrom_idx, terms_cache,
                             k_max = project$fingerprint$k_max %||% 64L,
                             features = project$fingerprint$features %||% "amplitude")
    if (is.null(fp)) { tsf_warn(col, ": no fingerprint"); next }
    # Normalise exactly as the reference samples were, then project onto the
    # reference feature space. Features the query lacks stay at 0 (the mean of
    # a normalised vector), so a partially covered query degrades rather than
    # failing.
    fp <- (fp - mean(fp)) / max(stats::sd(fp), .Machine$double.eps)
    vv <- stats::setNames(rep(0, length(ref$feature_space)), ref$feature_space)
    shared <- intersect(names(fp), ref$feature_space)
    vv[shared] <- fp[shared]
    if (length(shared) < 0.5 * length(ref$model$features)) {
      tsf_warn(col, ": only ", length(shared), " of ",
               length(ref$feature_space), " reference features present")
    }

    res <- match_query(ref$model, vv)
    if (is.null(res)) { tsf_warn(col, ": could not be scored"); next }

    cat(col, "\n")
    top <- utils::head(res$scores, 5)
    for (i in seq_len(nrow(top))) {
      cat(sprintf("  %-12s similarity %.3f\n", top$class[i], top$similarity[i]))
    }
    cat(sprintf("  margin over runner-up: %.3f   |   p(shuffled query): %.3f\n",
                res$margin, res$p_shuffle))

    if (is.null(ref$validation)) {
      cat("  -> The reference is uncalibrated. This is a suggestion, not an ",
          "identification.\n", sep = "")
    } else if (ref$validation$accuracy - ref$validation$baseline <= 0.02) {
      cat("  -> The reference does not beat guessing out of cohort. ",
          "Do not read this as an identification.\n", sep = "")
    } else if (is.na(res$margin) || res$margin < 0.02) {
      cat("  -> Top two classes are within 0.02: the call is not separable.\n")
    }
    cat("\n")
  }
  0L
}
