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

#' Cosine similarity of a query to every centroid.
score_query <- function(model, query_vec) {
  v <- query_vec[model$features]
  v[!is.finite(v)] <- 0
  if (all(v == 0)) return(NULL)
  cos_sim <- function(a, b) {
    d <- sqrt(sum(a^2)) * sqrt(sum(b^2))
    if (!is.finite(d) || d == 0) return(NA_real_)
    sum(a * b) / d
  }
  s <- apply(model$centroids, 1, function(c) cos_sim(v, c))
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
match_query <- function(model, query_vec, n_shuffle = 200L, seed = 1L) {
  s <- score_query(model, query_vec)
  if (is.null(s)) return(NULL)
  v <- query_vec[model$features]; v[!is.finite(v)] <- 0
  set.seed(seed)
  null_best <- vapply(seq_len(n_shuffle), function(i) {
    max(score_query(model, stats::setNames(sample(v), model$features))$similarity,
        na.rm = TRUE)
  }, numeric(1))
  best <- s$similarity[1]
  list(scores = s,
       best = s$class[1],
       similarity = best,
       margin = if (nrow(s) > 1) best - s$similarity[2] else NA_real_,
       p_shuffle = (1 + sum(null_best >= best)) / (n_shuffle + 1))
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
    model <- fit_centroids(mat[keep_tr, , drop = FALSE], lab[[target]][keep_tr],
                           n_features)
    pred <- vapply(which(keep_te), function(i) {
      sc <- score_query(model, mat[i, ])
      if (is.null(sc)) NA_character_ else sc$class[1]
    }, character(1))
    truth <- lab[[target]][keep_te]
    rows[[held]] <- data.frame(held_out = held, sample_id = lab$sample_id[keep_te],
                               truth = truth, predicted = pred,
                               stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(NULL)
  pred_all <- do.call(rbind, rows)

  acc <- mean(pred_all$truth == pred_all$predicted, na.rm = TRUE)
  # Baseline: always guess the commonest class of the test set.
  baseline <- max(table(pred_all$truth)) / nrow(pred_all)
  cm <- table(truth = pred_all$truth, predicted = pred_all$predicted)

  tsf_log("Out-of-cohort accuracy (", target, "): ", round(100 * acc, 1),
          "% vs ", round(100 * baseline, 1), "% for always guessing the ",
          "commonest class")
  list(target = target, predictions = pred_all, accuracy = acc,
       baseline = baseline, confusion = cm,
       n_datasets = length(ids), datasets = ids)
}

#' Build the reference: fingerprints, validation, and a model fitted on all data.
build_reference <- function(fps, target = "condition", n_features = 500L) {
  common <- Reduce(intersect, lapply(fps, function(f) colnames(f$matrix)))
  mat <- normalise_fingerprints(
    do.call(rbind, lapply(fps, function(f) f$matrix[, common, drop = FALSE])))
  lab <- do.call(rbind, lapply(fps, function(f) f$labels))

  validation <- validate_across_datasets(fps, target = target, n_features = n_features)
  model <- fit_centroids(mat, lab[[target]], n_features)

  list(model = model, labels = lab, target = target,
       feature_space = common, validation = validation,
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
