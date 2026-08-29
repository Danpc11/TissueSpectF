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
