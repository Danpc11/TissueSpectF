# reference.R -- build the reference library, validate it, and match a query.
#
# THE HONESTY CONSTRAINT
# ----------------------
# A matcher will always return a best match. Whether that match means anything
# is a separate question, and the only answer that counts here is out-of-cohort:
# train on one dataset, predict another, with no tuning against the second.
#
# Internal cross-validation cannot substitute. With thousands of features and
# hundreds of samples, a classifier reaches a high internal accuracy by learning
# batch, library protocol and sequencing depth. Those do not transfer between
# studies; a real spectral signature does. A reference built from a single
# dataset therefore carries `validation = NULL`, and everything downstream --
# the CLI and the app -- must say so rather than showing a confident label.

#' Nearest-centroid model on normalised fingerprints.
#'
#' Deliberately simple. With tens to hundreds of samples per class and thousands
#' of features, a flexible model would fit cohort idiosyncrasy; the centroid
#' plus a feature filter has few enough degrees of freedom that an out-of-cohort
#' result means something.
fit_centroids <- function(mat, labels, n_features = 500L) {
  # Feature filter computed on TRAINING data only: the frequencies whose
  # between-class variance is largest relative to within-class.
  cls <- factor(labels)
  keep_var <- apply(mat, 2, function(v) {
    if (!is.finite(stats::sd(v)) || stats::sd(v) == 0) return(0)
    between <- stats::var(tapply(v, cls, mean))
    within <- mean(tapply(v, cls, stats::var), na.rm = TRUE)
    if (!is.finite(between) || !is.finite(within) || within <= 0) return(0)
    between / within
  })
  sel <- names(sort(keep_var, decreasing = TRUE))[seq_len(min(n_features, ncol(mat)))]
  sel <- sel[!is.na(sel)]

  centroids <- t(vapply(levels(cls), function(l) {
    colMeans(mat[cls == l, sel, drop = FALSE])
  }, numeric(length(sel))))
  rownames(centroids) <- levels(cls)
  list(features = sel, centroids = centroids, classes = levels(cls),
       n_train = as.integer(table(cls)))
}

#' Cosine similarity of a query to every centroid, over shared features only.
#'
#' `available` names the features the query actually observed. Both sides are
#' restricted to those and re-standardised over that same subset before the
#' cosine, so query and centroid are on one scale and neither is padded with
#' zeros. Padding would be doubly wrong: an unobserved frequency is not the
#' mean, and a padded query contributes over a subset while the centroid keeps
#' its norm over everything, which shifts similarities between classes by
#' different amounts.
score_query <- function(model, query_vec, available = NULL) {
  sel <- model$features
  if (!is.null(available)) sel <- intersect(sel, available)
  sel <- intersect(sel, names(query_vec))
  if (length(sel) < 3) return(NULL)

  zs <- function(x) {
    s <- stats::sd(x)
    if (!is.finite(s) || s == 0) return(x - mean(x))
    (x - mean(x)) / s
  }
  v <- zs(query_vec[sel])
  if (all(!is.finite(v)) || all(v == 0)) return(NULL)

  cos_sim <- function(a, b) {
    d <- sqrt(sum(a^2)) * sqrt(sum(b^2))
    if (!is.finite(d) || d == 0) return(NA_real_)
    sum(a * b) / d
  }
  s <- apply(model$centroids[, sel, drop = FALSE], 1, function(c) cos_sim(v, zs(c)))
  ord <- order(-s)
  data.frame(class = names(s)[ord], similarity = unname(s[ord]),
             stringsAsFactors = FALSE)
}

#' Match a query, with a calibrated sense of how much the top score means.
#'
#' `margin` is the gap between the best and second-best similarity. `p_shuffle`
#' compares the best similarity against the distribution obtained by randomly
#' permuting the query's own feature values -- a query with no spectral shape at
#' all should not beat that. It is a sanity floor, not a class-membership test.
match_query <- function(model, query_vec, available = NULL,
                        n_shuffle = 200L, seed = 1L) {
  s <- score_query(model, query_vec, available)
  if (is.null(s)) return(NULL)
  sel <- intersect(intersect(model$features, names(query_vec)),
                   available %||% names(query_vec))
  v <- query_vec[sel]
  set.seed(seed)
  null_best <- vapply(seq_len(n_shuffle), function(i) {
    sc <- score_query(model, stats::setNames(sample(v), sel), sel)
    if (is.null(sc)) NA_real_ else max(sc$similarity, na.rm = TRUE)
  }, numeric(1))
  null_best <- null_best[is.finite(null_best)]
  if (!length(null_best)) null_best <- -Inf
  best <- s$similarity[1]
  list(scores = s,
       best = s$class[1],
       similarity = best,
       margin = if (nrow(s) > 1) best - s$similarity[2] else NA_real_,
       p_shuffle = (1 + sum(null_best >= best)) / (length(null_best) + 1),
       n_features_used = length(sel))
}

#' Leave-one-dataset-out validation. The only number worth reporting.
validate_across_datasets <- function(fps, target = c("condition", "tissue"),
                                     n_features = 500L) {
  target <- match.arg(target)
  ids <- unique(unlist(lapply(fps, function(f) f$labels$dataset_id)))
  if (length(ids) < 2) {
    tsf_warn("Only one dataset: out-of-cohort validation is not possible. ",
             "The reference will carry no validation and any match must be ",
             "reported as uncalibrated.")
    return(NULL)
  }

  common <- Reduce(intersect, lapply(fps, function(f) colnames(f$matrix)))
  mat <- do.call(rbind, lapply(fps, function(f) f$matrix[, common, drop = FALSE]))
  lab <- do.call(rbind, lapply(fps, function(f) f$labels))
  mat <- normalise_fingerprints(mat)

  rows <- list()
  excluded <- list(samples = 0L, classes = 0L, labels = character(0))
  for (held in ids) {
    tr <- lab$dataset_id != held
    te <- !tr
    y_tr <- lab[[target]][tr]
    shared <- intersect(unique(y_tr), unique(lab[[target]][te]))
    if (length(shared) < 2) {
      tsf_warn("Holding out ", held, ": fewer than two shared classes, skipped")
      next
    }
    keep_tr <- tr & lab[[target]] %in% shared
    keep_te <- te & lab[[target]] %in% shared
    dropped <- lab[[target]][te & !(lab[[target]] %in% shared)]
    excluded$samples <- excluded$samples + length(dropped)
    excluded$labels <- unique(c(excluded$labels, dropped))
    excluded$classes <- length(excluded$labels)
    model <- fit_centroids(mat[keep_tr, , drop = FALSE], lab[[target]][keep_tr],
                           n_features)
    scored <- lapply(which(keep_te), function(i) score_query(model, mat[i, ], NULL))
    pred <- vapply(scored, function(sc) if (is.null(sc)) NA_character_ else sc$class[1],
                   character(1))
    sim <- vapply(scored, function(sc) if (is.null(sc)) NA_real_ else sc$similarity[1],
                  numeric(1))
    margin <- vapply(scored, function(sc)
      if (is.null(sc) || nrow(sc) < 2) NA_real_ else sc$similarity[1] - sc$similarity[2],
      numeric(1))
    truth <- lab[[target]][keep_te]

    # Baseline: the majority class OF THE TRAINING FOLD, evaluated on the
    # held-out cohort. Taking the majority of the held-out truth would let the
    # baseline see the test labels, which is exactly what the held-out design
    # exists to prevent, and it flatters or punishes the model at random
    # depending on how the class mix differs between cohorts.
    majority_tr <- names(sort(table(lab[[target]][keep_tr]), decreasing = TRUE))[1]
    rows[[held]] <- data.frame(held_out = held, sample_id = lab$sample_id[keep_te],
                               truth = truth, predicted = pred,
                               similarity = sim, margin = margin,
                               baseline_prediction = majority_tr,
                               stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(NULL)
  pred_all <- do.call(rbind, rows)

  acc <- mean(pred_all$truth == pred_all$predicted, na.rm = TRUE)
  baseline <- mean(pred_all$truth == pred_all$baseline_prediction, na.rm = TRUE)
  cm <- table(truth = pred_all$truth, predicted = pred_all$predicted)

  tsf_log("Out-of-cohort accuracy (", target, "): ", round(100 * acc, 1),
          "% vs ", round(100 * baseline, 1), "% for the training-fold ",
          "majority class")
  if (excluded$samples > 0L) {
    tsf_log("Excluded from validation: ", excluded$samples, " sample(s) in ",
            excluded$classes, " class(es) not shared across cohorts (",
            paste(excluded$labels, collapse = ", "), ")")
  }

  calib <- calibrate_rejection(pred_all)
  list(target = target, predictions = pred_all, accuracy = acc,
       baseline = baseline, confusion = cm, excluded = excluded,
       calibration = calib, n_datasets = length(ids), datasets = ids)
}

#' Calibrate a rejection rule from the out-of-cohort predictions.
#'
#' A nearest-centroid matcher always returns a class. Nothing so far tells you
#' whether the query belongs to ANY class of the reference: a brain sample can
#' score \"liver F3\" with a comfortable margin simply because that centroid is
#' the least distant one. The fixed 0.02 margin and the shuffled-query test do
#' not address this -- the first is arbitrary, the second only asks whether the
#' query has any spectral shape at all.
#'
#' What can be calibrated from held-out data is the similarity a CORRECT match
#' typically reaches. Queries below that are reported as UNKNOWN. Per class,
#' because classes differ in how tight their centroids are.
#'
#' Honest limitation: this is calibrated on in-domain samples only. It bounds
#' how often a true member is wrongly rejected (that is the quantile), but it
#' cannot bound how often an out-of-domain sample is wrongly accepted, because
#' no out-of-domain sample was ever seen. Genuine open-set specificity needs
#' negatives in the validation -- for a tissue reference, other tissues.
calibrate_rejection <- function(pred_all, quantile_correct = 0.05) {
  ok <- !is.na(pred_all$similarity) & !is.na(pred_all$predicted)
  if (!any(ok)) return(NULL)
  correct <- ok & pred_all$truth == pred_all$predicted
  wrong <- ok & pred_all$truth != pred_all$predicted

  per_class <- do.call(rbind, lapply(sort(unique(pred_all$predicted[correct])), function(cl) {
    s <- pred_all$similarity[correct & pred_all$predicted == cl]
    if (length(s) < 3) return(NULL)
    data.frame(class = cl, n_correct = length(s),
               threshold = unname(stats::quantile(s, quantile_correct)),
               median_correct = stats::median(s), stringsAsFactors = FALSE)
  }))

  global <- if (any(correct))
    unname(stats::quantile(pred_all$similarity[correct], quantile_correct)) else NA_real_

  # How separable correct from incorrect matches are on similarity alone.
  auc <- if (any(correct) && any(wrong)) {
    a <- pred_all$similarity[correct]; b <- pred_all$similarity[wrong]
    mean(outer(a, b, ">") + 0.5 * outer(a, b, "=="))
  } else NA_real_

  list(per_class = per_class, global_threshold = global,
       quantile = quantile_correct,
       expected_false_rejection = quantile_correct,
       separability_auc = auc,
       median_correct = if (any(correct)) stats::median(pred_all$similarity[correct]) else NA_real_,
       median_incorrect = if (any(wrong)) stats::median(pred_all$similarity[wrong]) else NA_real_)
}

#' Apply the calibrated rule to one match.
apply_rejection <- function(res, calibration) {
  if (is.null(calibration)) {
    res$decision <- "UNCALIBRATED"
    return(res)
  }
  thr <- calibration$global_threshold
  pc <- calibration$per_class
  if (!is.null(pc) && res$best %in% pc$class) thr <- pc$threshold[pc$class == res$best]
  res$threshold <- thr
  res$decision <- if (is.na(thr) || res$similarity >= thr) res$best else "UNKNOWN"
  res
}

#' Build the reference: fingerprints, validation, and a model fitted on all data.
build_reference <- function(fps, target = "condition", n_features = 500L,
                            grid = NULL, params = list()) {
  common <- Reduce(intersect, lapply(fps, function(f) colnames(f$matrix)))
  mat <- normalise_fingerprints(
    do.call(rbind, lapply(fps, function(f) f$matrix[, common, drop = FALSE])))
  lab <- do.call(rbind, lapply(fps, function(f) f$labels))

  validation <- validate_across_datasets(fps, target = target, n_features = n_features)
  model <- fit_centroids(mat, lab[[target]], n_features)

  # SELF-CONTAINED. A reference must carry everything needed to fingerprint a
  # query: the grid with its identifiers, the frequency ceiling, the feature
  # representation, and the normalisation. Reconstructing that from whatever
  # dataset happens to be configured locally would silently allow a query to be
  # scored on a different grid from the one the reference was built on.
  if (is.null(grid)) tsf_abort("build_reference needs the grid it was built on")
  list(model = model, labels = lab, target = target,
       feature_space = common, validation = validation,
       grid = grid[, intersect(c("gene_id", "entrez_id", "chr", "start",
                                 "grid_index", "grid_N"), colnames(grid))],
       params = utils::modifyList(
         list(k_max = 64L, features = "amplitude", n_features = n_features,
              expression_unit = "asinh(CPM)", normalisation = "per-sample z-score",
              annotation = NA_character_, gene_universe = NA_character_,
              chrom_levels = NA_character_), params),
       built = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
       version = TSF_VERSION)
}

#' One-line honest summary of how much a match from this reference is worth.
reference_status <- function(ref) {
  v <- ref$validation
  if (is.null(v)) {
    return(paste0("UNCALIBRATED: built from a single cohort, so no out-of-cohort ",
                  "validation exists. Matches are suggestions, not identifications."))
  }
  lift <- v$accuracy - v$baseline
  verdict <- if (lift <= 0.02) {
    "NOT BETTER THAN GUESSING: matches from this reference carry no information."
  } else if (v$accuracy < 0.5) {
    "WEAK: better than chance but wrong more often than right."
  } else {
    "USABLE: correct on a majority of held-out samples."
  }
  sprintf("%s Out-of-cohort accuracy %.1f%% vs %.1f%% baseline (%d cohorts, target: %s).",
          verdict, 100 * v$accuracy, 100 * v$baseline, v$n_datasets, v$target)
}
