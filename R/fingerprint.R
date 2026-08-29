# fingerprint.R -- a sample's spectrum reduced to a comparable feature vector.
#
# This is the representation the matcher works on. It is deliberately NOT the
# list of significant peaks: a fingerprint does not need any single component to
# be significant, and thresholding by q would discard exactly the information
# that makes classes separable. Selection is left to the classifier.
#
# Feature space: one entry per (chromosome, frequency index k). Because the axis
# is the annotation grid, N per chromosome is fixed by the annotation, so k means
# the same thing in every dataset and the space is shared without alignment.
#
# Phase is kept, as cos/sin of the phase weighted by amplitude. Amplitude alone
# is invariant to translation along the chromosome -- it says a 200-gene
# structure exists but not where its crest falls. Phase anchors the structure to
# genomic coordinates, and for tissue identity it is plausibly the more
# discriminative half. Whether it helps is an empirical question the validation
# answers; `features` selects which representation to build.

FINGERPRINT_FEATURES <- c("amplitude", "amplitude_phase")

#' Feature vector for one signal on one chromosome set.
#'
#' @param k_max highest frequency index kept per chromosome. Low frequencies
#'   carry chromosome-scale trends and are the most reproducible; very high ones
#'   are the noisiest and the most exposed to the sampling window.
fingerprint_vector <- function(y, chrom_idx, terms_cache, k_max = 64L,
                               features = "amplitude") {
  pieces <- list()
  for (chr_now in names(chrom_idx)) {
    ci <- chrom_idx[[chr_now]]
    terms <- terms_cache[[chr_now]]
    sp <- gls_spectrum(y[ci$rows], terms)
    keep <- sp$k <= k_max
    sp <- sp[keep, , drop = FALSE]
    if (!nrow(sp)) next

    # log1p on amplitude: spectra are heavy-tailed, and a single dominant
    # component would otherwise drive every distance.
    amp <- log1p(sp$amplitude)
    nm <- paste0("chr", chr_now, "_k", sp$k)
    if (identical(features, "amplitude")) {
      v <- stats::setNames(amp, nm)
    } else {
      v <- stats::setNames(c(amp * cos(sp$phase), amp * sin(sp$phase)),
                           c(paste0(nm, "_c"), paste0(nm, "_s")))
    }
    pieces[[chr_now]] <- v
  }
  if (!length(pieces)) return(NULL)
  unlist(pieces, use.names = TRUE)
}

#' Cache the window terms once per chromosome (they depend only on positions).
fingerprint_terms <- function(chrom_idx) {
  lapply(chrom_idx, function(ci) gls_prepare(ci$t, ci$N))
}

#' Fingerprints for every sample of a dataset.
#'
#' @return list(matrix = samples x features, labels = data.frame)
fingerprint_dataset <- function(dataset, chrom_idx, k_max = 64L,
                                features = "amplitude") {
  terms_cache <- fingerprint_terms(chrom_idx)
  samples <- dataset$samples$sample_id[dataset$samples$sample_id %in%
                                         colnames(dataset$expression)]
  if (!length(samples)) return(NULL)

  vecs <- lapply(samples, function(s) {
    fingerprint_vector(dataset$expression[, s], chrom_idx, terms_cache,
                       k_max = k_max, features = features)
  })
  keep <- !vapply(vecs, is.null, logical(1))
  vecs <- vecs[keep]; samples <- samples[keep]
  if (!length(vecs)) return(NULL)

  common <- Reduce(intersect, lapply(vecs, names))
  mat <- do.call(rbind, lapply(vecs, function(v) v[common]))
  rownames(mat) <- samples
  colnames(mat) <- common

  lab <- dataset$samples[match(samples, dataset$samples$sample_id), , drop = FALSE]
  list(matrix = mat,
       labels = data.frame(sample_id = samples,
                           dataset_id = dataset$id,
                           tissue = if ("tissue" %in% colnames(lab))
                             as.character(lab$tissue) else NA_character_,
                           condition = as.character(lab$condition),
                           stringsAsFactors = FALSE))
}

#' Per-sample normalisation.
#'
#' Centre and scale each fingerprint so that library depth, and any global shift
#' in expression level, cannot drive a match. What survives is the SHAPE of the
#' spectrum across frequencies, which is the thing claimed to be characteristic.
normalise_fingerprints <- function(mat) {
  mat[!is.finite(mat)] <- 0
  t(apply(mat, 1, function(v) {
    s <- stats::sd(v)
    if (!is.finite(s) || s == 0) return(v - mean(v))
    (v - mean(v)) / s
  }))
}


# --- query side ---------------------------------------------------------------
#
# A query file covers a different set of genes from the reference cohorts, and
# the genes it lacks must be treated exactly as unmeasured genes are everywhere
# else in this pipeline: absent from the observed positions, NEVER zero.
#
# Filling a missing gene with zero would put a deterministic value at a fixed
# grid position, and with the 20% coverage floor up to four fifths of the grid
# could be zeros. The resulting spectrum would be dominated by the pattern of
# absence rather than by expression -- the precise failure the grid design
# exists to prevent. An explicit zero count inside the file is a measurement and
# is kept; a gene not in the file is not.

#' Observed positions for a query, on the reference grid.
#'
#' @param grid the reference's own grid (gene_id, entrez_id, chr, grid_index, grid_N)
#' @param present_ids identifiers actually present in the query file
#' @return list(chrom_idx, rows_gene_id, coverage, id_type)
query_grid_index <- function(grid, present_ids, min_observed = 8L,
                             min_coverage = 0.05) {
  ids <- unique(as.character(present_ids))
  by_ensembl <- sum(grid$gene_id %in% ids)
  by_entrez <- if ("entrez_id" %in% colnames(grid))
    sum(as.character(grid$entrez_id) %in% ids) else 0L
  id_type <- if (by_entrez > by_ensembl) "entrez" else "ensembl"
  key <- if (id_type == "entrez") as.character(grid$entrez_id) else grid$gene_id

  observed <- key %in% ids
  g <- grid[observed, , drop = FALSE]
  if (!nrow(g)) return(NULL)

  chrom_idx <- list()
  for (chr_now in unique(g$chr)) {
    sel <- which(g$chr == chr_now)
    sel <- sel[order(g$grid_index[sel])]
    N <- g$grid_N[sel[1]]
    cov <- length(sel) / N
    if (length(sel) < min_observed || cov < min_coverage) next
    chrom_idx[[as.character(chr_now)]] <-
      list(rows = sel, t = as.integer(g$grid_index[sel]),
           N = as.integer(N), coverage = cov)
  }
  if (!length(chrom_idx)) return(NULL)
  list(chrom_idx = chrom_idx, genes = g, key = key[observed],
       coverage = nrow(g) / nrow(grid), id_type = id_type)
}

#' Collapse duplicate identifiers explicitly.
#'
#' A query file often carries several rows per gene (transcript-level rows,
#' duplicated symbols). match() would keep the first silently, which loses most
#' of a gene's signal and does so invisibly. Counts are summed; anything already
#' normalised cannot be summed, so duplicates there are an error the caller has
#' to resolve.
collapse_duplicate_ids <- function(values, ids, unit = "counts") {
  v <- suppressWarnings(as.numeric(values))
  ids <- as.character(ids)
  keep <- !is.na(ids) & nzchar(ids)
  v <- v[keep]; ids <- ids[keep]

  if (any(v < 0, na.rm = TRUE)) {
    tsf_abort(sum(v < 0, na.rm = TRUE), " negative value(s) in the query. ",
              "Counts and abundances cannot be negative; if this is already ",
              "log-transformed, declare it with --input-unit logged.")
  }
  dup <- duplicated(ids) | duplicated(ids, fromLast = TRUE)
  if (!any(dup)) return(list(values = v, ids = ids, collapsed = 0L))

  n_dup_ids <- length(unique(ids[dup]))
  if (!unit %in% c("counts")) {
    tsf_abort(n_dup_ids, " identifier(s) appear more than once. Rows can only ",
              "be summed for counts; with unit '", unit, "' the duplicates must ",
              "be resolved before matching.")
  }
  agg <- stats::aggregate(v, list(id = ids), sum, na.rm = TRUE)
  tsf_log("Collapsed ", sum(dup), " row(s) over ", n_dup_ids,
          " duplicated identifier(s) by summing counts")
  list(values = agg$x, ids = agg$id, collapsed = sum(dup))
}

#' Values on the observed positions, on the scale the reference was built on.
query_signal <- function(v, unit = "counts") {
  switch(unit,
    counts = asinh(v / max(sum(v, na.rm = TRUE), 1) * 1e6),
    cpm = ,
    tpm = asinh(v),
    logged = v,
    tsf_abort("Unknown input unit '", unit,
              "'. Use counts, cpm, tpm or logged."))
}

#' Fingerprint of one query column, built only on the genes it actually contains.
#'
#' Returns the fingerprint UNNORMALISED. Normalisation has to happen after the
#' intersection with the reference feature space, not before: a query missing a
#' chromosome would otherwise have its mean and standard deviation computed over
#' a different set of features from the one the centroids were trained on, and
#' the two would no longer be on the same scale.
fingerprint_query <- function(values, ids, ref, unit = "counts") {
  cl <- collapse_duplicate_ids(values, ids, unit)
  qi <- query_grid_index(ref$grid, cl$ids)
  if (is.null(qi)) return(NULL)

  v <- cl$values[match(qi$key, cl$ids)]
  # Non-finite entries are unusable measurements, not zeros: drop those
  # positions from the observed set rather than imputing them.
  if (!all(is.finite(v))) {
    qi <- query_grid_index(ref$grid, qi$key[is.finite(v)])
    if (is.null(qi)) return(NULL)
    v <- cl$values[match(qi$key, cl$ids)]
  }
  y <- query_signal(v, unit)

  terms <- fingerprint_terms(qi$chrom_idx)
  fp <- fingerprint_vector(y, qi$chrom_idx, terms,
                           k_max = ref$params$k_max,
                           features = ref$params$features)
  if (is.null(fp)) return(NULL)
  list(vector = fp, coverage = qi$coverage, id_type = qi$id_type,
       n_chromosomes = length(qi$chrom_idx), collapsed = cl$collapsed,
       unit = unit)
}

#' Which reference features a query can actually contribute.
#'
#' No zero-filling. A frequency the query never observed is ABSENT, not average:
#' filling it with 0 (the mean of a z-scored vector) would let the missing part
#' of the query pull every similarity toward the centroid mean, and would let
#' the centroid keep its full norm while the query contributes over a subset.
#' Similarity is computed over the shared features only.
project_to_reference <- function(fp, ref) {
  shared <- intersect(names(fp), ref$feature_space)
  model_shared <- intersect(shared, ref$model$features)
  list(available = shared,
       model_available = model_shared,
       n_shared = length(shared),
       n_features = length(ref$feature_space),
       feature_coverage = if (length(ref$model$features))
         length(model_shared) / length(ref$model$features) else 0,
       vector = fp[shared])
}
