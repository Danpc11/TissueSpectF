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

FINGERPRINT_FEATURES <- c("amplitude", "amplitude_phase",
                          "period_bins", "period_bins_genomic",
                          "band_ratios", "expression_baseline")

# Common period grid, in genes, shared by every chromosome.
#
# WHY THIS EXISTS
# ---------------
# "amplitude" indexes features by (chromosome, k). k is cycles per chromosome,
# so the same k means different physical scales on different chromosomes: with
# N = 2066 on chr1, k = 64 is a period of 32 genes; with N = 759 on chr21 it is
# 12. chr1_k64 and chr21_k64 are columns the model treats as parallel while
# they describe different phenomena, and 24 chromosomes' worth of that is most
# of the feature space.
#
# Indexing by period instead makes bin_p21 mean the same thing everywhere, so
# the features are comparable across chromosomes and -- with
# "period_bins_genomic" -- can be averaged into a single genome-wide spectrum.
# That is the characteristic spectrum in the sense the project is after: one
# curve over period, not a concatenation of 24 curves over incompatible axes.
#
# Log spacing because period resolution is multiplicative: the gap between 20
# and 21 genes is not the gap between 400 and 401.
#
# The range stops at 10 genes because below that the technical floor of most
# chromosomes (2/coverage + margin) puts the frequency past the resolution the
# sampling supports, and at 500 because a period that long fits only a few
# times into the shorter chromosomes.
FINGERPRINT_PERIOD_BREAKS <- exp(seq(log(10), log(500), length.out = 41L))

#' Feature vector for one signal on one chromosome set.
#'
#' @param k_max highest frequency index kept per chromosome. Low frequencies
#'   carry chromosome-scale trends and are the most reproducible; very high ones
#'   are the noisiest and the most exposed to the sampling window.
fingerprint_vector <- function(y, chrom_idx, terms_cache, k_max = 64L,
                               features = "amplitude") {
  # THE CONTROL, not a fingerprint: the gene expression itself, positioned on
  # the same grid and carried through the same leave-one-cohort-out validation,
  # the same feature selection and the same thresholds.
  #
  # Without it the accuracy of a spectral fingerprint is uninterpretable. If raw
  # expression classifies these cohorts better than the spectrum does, the
  # Fourier transform is discarding information and the work has to say so; if
  # it classifies worse, the spectrum is a genuine compression. A number with no
  # control is a number nobody can act on, and this pipeline had no control.
  #
  # It is a representation rather than a separate script so that every stage
  # downstream is bit-for-bit the same, and the comparison isolates the
  # representation instead of confounding it with the harness.
  if (identical(features, "expression_baseline")) {
    out <- list()
    for (chr_now in names(chrom_idx)) {
      ci <- chrom_idx[[chr_now]]
      v <- as.numeric(y[ci$rows])
      names(v) <- paste0("chr", chr_now, "_g", seq_along(v))
      out[[chr_now]] <- v
    }
    if (!length(out)) return(NULL)
    return(unlist(out, use.names = TRUE))
  }

  pieces <- list()
  for (chr_now in names(chrom_idx)) {
    ci <- chrom_idx[[chr_now]]
    # A query, or a masked validation sample, is missing genes the reference
    # cohorts had. Those positions leave the observed set; they are never
    # zero-filled. N is unchanged, so the feature names chr<chr>_k<k> still
    # index the same frequencies and stay comparable to the reference.
    fit <- gls_observed(y[ci$rows], ci$t, ci$N, terms = terms_cache[[chr_now]])
    if (is.null(fit)) next
    sp <- gls_spectrum(fit$y, fit$terms)

    # k_max caps cycles-per-chromosome, which is the right selection for the
    # (chromosome, k) representations and the wrong one for the period-binned
    # ones: on chr1 with N = 2066, k <= 64 means period >= 32 genes, so every
    # bin below 32 would be empty and the short-period half of the common grid
    # would exist only on the short chromosomes. The period range is itself the
    # selection there, so k_max does not apply.
    period_indexed <- features %in% c("period_bins", "period_bins_genomic")
    if (!period_indexed) {
      sp <- sp[sp$k <= k_max, , drop = FALSE]
    }
    if (!nrow(sp)) next

    # log1p on amplitude: spectra are heavy-tailed, and a single dominant
    # component would otherwise drive every distance.
    amp <- log1p(sp$amplitude)
    nm <- paste0("chr", chr_now, "_k", sp$k)
    v <- if (identical(features, "amplitude")) {
      stats::setNames(amp, nm)
    } else if (identical(features, "amplitude_phase")) {
      stats::setNames(c(amp * cos(sp$phase), amp * sin(sp$phase)),
                      c(paste0(nm, "_c"), paste0(nm, "_s")))
    } else {
      # Period-binned: average log-amplitude within each common period bin.
      # Averaging rather than interpolating, because a bin on a short
      # chromosome may contain no frequency at all and interpolation would
      # invent a value there. An empty bin stays NA and is handled below.
      # band_ratios averages across chromosomes exactly as period_bins_genomic
      # does, so it needs the unprefixed bin names too. Prefixed names would
      # make rbind align positionally while labelling every column with the
      # first chromosome's name -- right by accident, wrong on the label.
      period_bin_vector(sp, chr_now,
                        genomic = features %in% c("period_bins_genomic",
                                                  "band_ratios"),
                        transform = if (identical(features, "band_ratios"))
                          "log" else "log1p")
    }
    if (is.null(v) || !length(v)) next
    pieces[[chr_now]] <- v
  }
  if (!length(pieces)) return(NULL)

  if (identical(features, "band_ratios")) {
    # Shazam's idea, translated: what identifies a spectrum is the RELATION
    # between its landmarks, not their magnitudes. A ratio of two period bands
    # within one sample is dimensionless -- library size, sequencing depth and
    # platform scale cancel -- so it should survive the cohort effect that a
    # raw amplitude carries.
    m <- do.call(rbind, pieces)
    genomic <- colMeans(m, na.rm = TRUE)
    ok <- which(is.finite(genomic))
    if (length(ok) < 2L) return(NULL)
    pairs <- utils::combn(ok, 2L)
    v <- genomic[pairs[1, ]] - genomic[pairs[2, ]]   # log amplitudes: difference IS the ratio
    names(v) <- paste0(names(genomic)[pairs[1, ]], "/", names(genomic)[pairs[2, ]])
    return(v)
  }

  if (identical(features, "period_bins_genomic")) {
    # One spectrum for the whole genome: average each period bin across the
    # chromosomes that have a frequency in it. 40 numbers instead of ~1500.
    m <- do.call(rbind, pieces)
    return(colMeans(m, na.rm = TRUE))
  }

  out <- unlist(pieces, use.names = TRUE)
  # A bin with no frequency on that chromosome is NA, not zero: zero would
  # claim the spectrum has no power there, which is a measurement the data did
  # not make. Downstream distance code drops non-finite features pairwise.
  out
}

#' Average log-amplitude within each common period bin.
#'
#' Returns one value per bin. Bins with no frequency on this chromosome are NA
#' -- typical on short chromosomes at long periods, where the chromosome simply
#' does not contain that many genes.
#' @param transform "log1p" tames the heavy tail for distance-based features;
#'   "log" is required whenever a DIFFERENCE of two bins has to be a RATIO of
#'   amplitudes, because log1p(5a) - log1p(5b) is not log1p(a) - log1p(b) and
#'   the scale invariance band_ratios exists for would not hold.
period_bin_vector <- function(sp, chr_now, genomic = FALSE,
                              breaks = FINGERPRINT_PERIOD_BREAKS,
                              transform = c("log1p", "log")) {
  transform <- match.arg(transform)
  ok <- is.finite(sp$period) & sp$period >= min(breaks) & sp$period <= max(breaks)
  n_bins <- length(breaks) - 1L
  labels <- sprintf("p%05.1f", sqrt(breaks[-1] * breaks[-length(breaks)]))
  out <- stats::setNames(rep(NA_real_, n_bins), labels)
  if (!any(ok)) return(if (genomic) out else stats::setNames(out, paste0("chr", chr_now, "_", labels)))

  idx <- cut(sp$period[ok], breaks = breaks, include.lowest = TRUE, labels = FALSE)
  amp <- if (identical(transform, "log")) {
    # Floored at a fixed constant, not at a fraction of the sample's own
    # amplitudes: a data-dependent floor would scale with the signal and
    # reintroduce exactly the dependence this transform removes.
    log(pmax(sp$amplitude[ok], 1e-12))
  } else {
    log1p(sp$amplitude[ok])
  }
  agg <- tapply(amp, factor(idx, levels = seq_len(n_bins)), mean, na.rm = TRUE)
  out[] <- as.numeric(agg)
  if (genomic) out else stats::setNames(out, paste0("chr", chr_now, "_", labels))
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
  if (!any(dup)) return(list(values = v, ids = ids, collapsed = 0L, all_na = 0L))

  n_dup_ids <- length(unique(ids[dup]))
  if (!unit %in% c("counts")) {
    tsf_abort(n_dup_ids, " identifier(s) appear more than once. Rows can only ",
              "be summed for counts; with unit '", unit, "' the duplicates must ",
              "be resolved before matching.")
  }
  # sum(na.rm = TRUE) over an all-NA group returns 0, which would turn "this
  # gene was not measured" into "this gene was measured at zero" -- exactly the
  # substitution the grid design exists to prevent. Such a gene stays NA and is
  # therefore dropped from the observed positions later.
  agg <- stats::aggregate(v, list(id = ids), function(x)
    if (all(!is.finite(x))) NA_real_ else sum(x, na.rm = TRUE))
  n_all_na <- sum(is.na(agg$x))
  tsf_log("Collapsed ", sum(dup), " row(s) over ", n_dup_ids,
          " duplicated identifier(s) by summing counts",
          if (n_all_na) paste0("; ", n_all_na,
            " gene(s) had no finite value and stay unmeasured") else "")
  list(values = agg$x, ids = agg$id, collapsed = sum(dup), all_na = n_all_na)
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
