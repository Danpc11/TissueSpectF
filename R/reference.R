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
                                     n_features = 500L, n_masks = 10L,
                                     datasets = NULL, ref_params = list(),
                                     max_queries_per_mask = 25L,
                                     threshold_policy = "pooled") {
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
  band_pred <- if (is.null(datasets)) {
    tsf_warn("No expression matrices supplied: gene-level coverage cannot be ",
             "simulated, so no per-band threshold will exist and partial ",
             "queries will be reported UNCALIBRATED_COVERAGE.")
    NULL
  } else {
    calibrate_coverage_bands(datasets, fps, lab, ids, target, n_features,
                             ref_params = ref_params, n_masks = n_masks,
                             max_queries_per_mask = max_queries_per_mask)
  }
  calib$bands <- summarise_bands(band_pred, policy = threshold_policy)
  if (!is.null(calib$bands)) {
    tsf_log("Accuracy by query coverage band:")
    for (i in seq_len(nrow(calib$bands))) {
      tsf_log("  ", calib$bands$band[i], ": ",
              round(100 * calib$bands$accuracy[i], 1), "% (n = ",
              calib$bands$n[i], " over ", calib$bands$n_masks[i], " masks, sd ",
              round(calib$bands$accuracy_sd_between_masks[i], 3),
              ", threshold ",
              if (is.na(calib$bands$threshold_applied[i])) "none" else
                round(calib$bands$threshold_applied[i], 3),
              " [", calib$bands$policy[i], "], rejects ",
              round(100 * calib$bands$expected_rejection_of_members[i], 1),
              "% of members)")
    }
  }
  list(target = target, predictions = pred_all, accuracy = acc,
       band_predictions = band_pred,
       baseline = baseline, confusion = cm, excluded = excluded,
       calibration = calib, n_datasets = length(ids), datasets = ids)
}

#' Realistic GENE loss on the reference grid.
#'
#' The thing a real query loses is genes, not frequencies, and the two are not
#' interchangeable. Masking spectral features (chr1_k6, chr1_k7, ...) models a
#' panel that reports some frequencies and not others, which no experiment does.
#' Dropping half the genes of a chromosome, by contrast, leaves the GLS able to
#' estimate almost every frequency -- feature coverage stays near 100% while
#' gene coverage is 50% -- but every estimate is noisier and the spectral window
#' changes. Calibrating on feature masks therefore measures the wrong quantity
#' and, because it is the easier one, measures it optimistically.
#'
#' Loss is simulated on the grid and the fingerprint is recomputed from the
#' surviving genes.
#'
#'   random      scattered genes, a shallower library
#'   block       contiguous runs of grid positions, a capture kit
#'   chromosome  whole chromosomes, a targeted panel
mask_grid_genes <- function(chrom_idx, target_coverage,
                            mode = c("random", "block", "chromosome"), seed = 1L) {
  mode <- match.arg(mode)
  set.seed(seed)
  chrs <- names(chrom_idx)

  if (mode == "chromosome") {
    order_chr <- sample(chrs)
    total <- sum(vapply(chrom_idx, function(ci) length(ci$t), integer(1)))
    keep_chr <- character(0); kept <- 0L
    for (c in order_chr) {
      if (kept >= target_coverage * total) break
      keep_chr <- c(keep_chr, c); kept <- kept + length(chrom_idx[[c]]$t)
    }
    out <- chrom_idx[keep_chr]
    return(out)
  }

  out <- lapply(chrom_idx, function(ci) {
    n <- length(ci$t)
    m <- max(1L, floor(target_coverage * n))
    sel <- if (mode == "block") {
      start <- sample.int(max(1L, n - m + 1L), 1L)
      start:min(n, start + m - 1L)
    } else {
      sort(sample.int(n, m))
    }
    list(rows = ci$rows[sel], t = ci$t[sel], N = ci$N,
         coverage = length(sel) / ci$N)
  })
  # A chromosome left with too few genes is dropped, exactly as ingest would.
  out[vapply(out, function(ci) length(ci$t) >= 8L, logical(1))]
}

#' Recompute a sample's fingerprint from a masked grid.
fingerprint_masked <- function(y, masked_idx, k_max, features) {
  if (!length(masked_idx)) return(NULL)
  terms <- fingerprint_terms(masked_idx)
  fingerprint_vector(y, masked_idx, terms, k_max = k_max, features = features)
}

COVERAGE_BANDS <- list(
  list(name = "90-100%", lower = 0.90, upper = 1.01),
  list(name = "75-90%",  lower = 0.75, upper = 0.90),
  list(name = "50-75%",  lower = 0.50, upper = 0.75),
  list(name = "<50%",    lower = 0.00, upper = 0.50)
)

coverage_band <- function(coverage) {
  for (b in COVERAGE_BANDS) {
    if (!is.na(coverage) && coverage >= b$lower && coverage < b$upper) return(b$name)
  }
  "<50%"
}

#' Score held-out samples again at reduced GENE coverage, band by band.
#'
#' For each fold the centroids are fitted on the full training fingerprints;
#' each held-out sample is then re-fingerprinted from a masked grid and scored.
#' Bands are keyed on GENE coverage, because that is the quantity a query can
#' report about itself before anything is computed.
#'
#' Cost is one GLS fingerprint per (sample, level, mode, mask), so
#' `max_queries_per_mask` caps how many held-out samples each mask scores.
calibrate_coverage_bands <- function(datasets, fps, lab, ids, target, n_features,
                                     ref_params,
                                     coverage_levels = c(0.95, 0.8, 0.6, 0.4),
                                     modes = c("random", "block", "chromosome"),
                                     n_masks = 10L, max_queries_per_mask = 25L,
                                     seed = 7L) {
  common <- Reduce(intersect, lapply(fps, function(f) colnames(f$matrix)))
  full_mat <- normalise_fingerprints(
    do.call(rbind, lapply(fps, function(f) f$matrix[, common, drop = FALSE])))

  rows <- list()
  for (held in ids) {
    tr <- lab$dataset_id != held; te <- !tr
    shared <- intersect(unique(lab[[target]][tr]), unique(lab[[target]][te]))
    if (length(shared) < 2) next
    keep_tr <- tr & lab[[target]] %in% shared
    model <- fit_centroids(full_mat[keep_tr, , drop = FALSE],
                           lab[[target]][keep_tr], n_features)

    ds <- datasets[[held]]
    if (is.null(ds)) next
    te_idx <- which(te & lab[[target]] %in% shared)
    te_idx <- te_idx[lab$sample_id[te_idx] %in% colnames(ds$dataset$expression)]
    if (!length(te_idx)) next

    for (lev in coverage_levels) {
      for (md in modes) {
        for (m in seq_len(n_masks)) {
          mask_seed <- seed + 1000L * m + 17L * match(md, modes) +
            as.integer(round(100 * lev)) + 7L * match(held, ids)
          masked <- mask_grid_genes(ds$chrom_idx, lev, md, mask_seed)
          if (!length(masked)) next
          gene_cov <- sum(vapply(masked, function(ci) length(ci$t), integer(1))) /
            sum(vapply(ds$chrom_idx, function(ci) length(ci$t), integer(1)))

          set.seed(mask_seed)
          pick <- if (length(te_idx) > max_queries_per_mask)
            sample(te_idx, max_queries_per_mask) else te_idx

          for (i in pick) {
            y <- ds$dataset$expression[, lab$sample_id[i]]
            fp <- fingerprint_masked(y, masked, ref_params$k_max,
                                     ref_params$features)
            if (is.null(fp)) next
            avail <- intersect(names(fp), model$features)
            if (length(avail) < 3) next
            sc <- score_query(model, fp, avail)
            if (is.null(sc)) next
            rows[[length(rows) + 1]] <- data.frame(
              held_out = held, mask = m, mode = md,
              gene_coverage = gene_cov,
              feature_coverage = length(avail) / length(model$features),
              band = coverage_band(gene_cov),
              truth = lab[[target]][i], predicted = sc$class[1],
              similarity = sc$similarity[1],
              margin = if (nrow(sc) > 1) sc$similarity[1] - sc$similarity[2] else NA_real_,
              stringsAsFactors = FALSE)
          }
        }
      }
    }
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

#' Per-band accuracy and rejection threshold.
summarise_bands <- function(band_pred, quantile_correct = 0.05,
                            policy = "pooled") {
  if (is.null(band_pred)) return(NULL)
  bands <- vapply(COVERAGE_BANDS, function(b) b$name, character(1))
  out <- do.call(rbind, lapply(bands, function(b) {
    d <- band_pred[band_pred$band == b, , drop = FALSE]
    if (!nrow(d)) return(NULL)
    correct <- d$truth == d$predicted
    # Spread BETWEEN masks, not just between samples: two panels with the same
    # coverage can behave differently depending on which regions they drop, and
    # a threshold that ignores that is calibrated to one accidental panel.
    per_mask <- if ("mask" %in% colnames(d)) {
      split(seq_len(nrow(d)), paste(d$held_out, d$mode, d$mask))
    } else list(seq_len(nrow(d)))
    acc_by_mask <- vapply(per_mask, function(i)
      mean(d$truth[i] == d$predicted[i], na.rm = TRUE), numeric(1))
    thr_by_mask <- vapply(per_mask, function(i) {
      cc <- d$truth[i] == d$predicted[i]
      if (sum(cc, na.rm = TRUE) < 3) return(NA_real_)
      unname(stats::quantile(d$similarity[i][cc], quantile_correct, na.rm = TRUE))
    }, numeric(1))

    data.frame(
      band = b, n = nrow(d), n_masks = length(per_mask),
      accuracy = mean(correct, na.rm = TRUE),
      accuracy_sd_between_masks = stats::sd(acc_by_mask, na.rm = TRUE),
      accuracy_lower = unname(stats::quantile(acc_by_mask, 0.05, na.rm = TRUE)),
      median_similarity_correct = if (any(correct))
        stats::median(d$similarity[correct]) else NA_real_,
      # The conservative end of the per-mask thresholds, so a query whose
      # missing regions happen to be the awkward ones is not judged against a
      # threshold calibrated on a luckier mask.
      threshold = if (sum(correct) >= 5)
        unname(stats::quantile(d$similarity[correct], quantile_correct)) else NA_real_,
      threshold_sd_between_masks = stats::sd(thr_by_mask, na.rm = TRUE),
      threshold_conservative = if (all(is.na(thr_by_mask))) NA_real_ else
        unname(stats::quantile(thr_by_mask, 0.90, na.rm = TRUE)),
      unknown_rate_at_threshold = NA_real_,
      unknown_rate_conservative = NA_real_,
      classify = !identical(b, "<50%"),
      stringsAsFactors = FALSE)
  }))
  if (is.null(out)) return(NULL)
  # What fraction of held-out members that band's own threshold would reject.
  rate_at <- function(i, col) {
    d <- band_pred[band_pred$band == out$band[i], , drop = FALSE]
    thr <- out[[col]][i]
    if (is.na(thr) || !nrow(d)) return(NA_real_)
    mean(d$similarity < thr, na.rm = TRUE)
  }
  out$unknown_rate_at_threshold <- vapply(seq_len(nrow(out)), rate_at,
                                          numeric(1), col = "threshold")
  out$unknown_rate_conservative <- vapply(seq_len(nrow(out)), rate_at,
                                          numeric(1), col = "threshold_conservative")
  # Which of the two is actually applied is a declared policy, not an
  # accident. "conservative" uses the 90th percentile of the per-mask
  # thresholds, so a query whose missing regions happen to be the awkward ones
  # is not judged against a threshold calibrated on a luckier mask; it rejects
  # more true members, and the rate at which it does so is recorded here.
  out$threshold_applied <- switch(policy,
    conservative = ifelse(is.na(out$threshold_conservative), out$threshold,
                          out$threshold_conservative),
    pooled = out$threshold,
    tsf_abort("threshold_policy must be 'conservative' or 'pooled'"))
  out$expected_rejection_of_members <- switch(policy,
    conservative = ifelse(is.na(out$threshold_conservative),
                          out$unknown_rate_at_threshold,
                          out$unknown_rate_conservative),
    pooled = out$unknown_rate_at_threshold)
  out$policy <- policy
  out
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
#' Apply the calibrated rule, using the threshold for THIS query's coverage.
#'
#' A full-coverage threshold applied to a partial query is simply the wrong
#' number: fewer shared features shift the whole similarity distribution. When a
#' band threshold exists it wins; the per-class and global thresholds are the
#' fallback for full coverage only.
apply_rejection <- function(res, calibration, coverage = NA_real_) {
  if (is.null(calibration)) {
    res$decision <- "UNCALIBRATED"
    return(res)
  }
  band <- if (is.na(coverage)) NA_character_ else coverage_band(coverage)
  res$coverage_band <- band

  if (!is.na(band) && identical(band, "<50%")) {
    res$decision <- "LOW_COVERAGE"
    res$threshold <- NA_real_
    return(res)
  }

  thr <- NA_real_; source <- "none"
  if (!is.na(band) && !is.null(calibration$bands)) {
    b <- calibration$bands[calibration$bands$band == band, , drop = FALSE]
    col <- if ("threshold_applied" %in% colnames(b)) "threshold_applied" else "threshold"
    if (nrow(b) && !is.na(b[[col]][1])) {
      thr <- b[[col]][1]
      source <- paste0("band:", b$policy[1] %||% "pooled")
    }
  }
  if (is.na(thr)) {
    pc <- calibration$per_class
    if (!is.null(pc) && res$best %in% pc$class) {
      thr <- pc$threshold[pc$class == res$best]; source <- "class"
    } else {
      thr <- calibration$global_threshold; source <- "global"
    }
    if (!is.na(band) && !identical(band, "90-100%") && !startsWith(source, "band")) {
      # Refuse to silently reuse a full-coverage threshold on a partial query.
      res$decision <- "UNCALIBRATED_COVERAGE"
      res$threshold <- NA_real_
      res$threshold_source <- source
      return(res)
    }
  }
  res$threshold <- thr
  res$threshold_source <- source
  res$decision <- if (is.na(thr) || res$similarity >= thr) res$best else "UNKNOWN"
  res
}

#' Build the reference: fingerprints, validation, and a model fitted on all data.
build_reference <- function(fps, target = "condition", n_features = 500L,
                            grid = NULL, params = list(), n_masks = 10L,
                            datasets = NULL, max_queries_per_mask = 25L,
                            threshold_policy = "pooled") {
  common <- Reduce(intersect, lapply(fps, function(f) colnames(f$matrix)))
  mat <- normalise_fingerprints(
    do.call(rbind, lapply(fps, function(f) f$matrix[, common, drop = FALSE])))
  lab <- do.call(rbind, lapply(fps, function(f) f$labels))

  validation <- validate_across_datasets(fps, target = target,
                                         n_features = n_features,
                                         n_masks = n_masks, datasets = datasets,
                                         ref_params = params,
                                         max_queries_per_mask = max_queries_per_mask,
                                         threshold_policy = threshold_policy)
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
              expression_unit = NA_character_, normalisation = "z-score over shared features",
              annotation = NA_character_, gene_universe = NA_character_,
              chrom_levels = NA_character_), params),
       built = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
       version = TSF_VERSION)
}

#' Is a query on the same expression scale as the reference?
#'
#' A reference built on TPM and one built on CPM are different objects: length
#' normalisation changes the relative height of every gene, so it changes the
#' spectrum. Mixing them would silently compare spectra of different quantities.
#' The reference records the unit it was built on and a query must arrive on a
#' scale that reduces to it.
assert_unit_compatible <- function(ref, query_unit) {
  ref_unit <- ref$params$expression_unit
  if (is.null(ref_unit) || is.na(ref_unit)) {
    tsf_warn("This reference does not record the expression unit it was built ",
             "on. Rebuild it; matching a TPM query against a CPM reference is ",
             "silent nonsense.")
    return(invisible(TRUE))
  }
  base <- sub("^asinh\\((.*)\\)$", "\\1", ref_unit)   # "asinh(TPM)" -> "TPM"
  produced <- switch(query_unit,
                     counts = "CPM",       # counts here are depth-normalised only
                     cpm = "CPM", tpm = "TPM", logged = "PRE-TRANSFORMED",
                     toupper(query_unit))
  if (identical(produced, "PRE-TRANSFORMED")) {
    tsf_warn("The query is declared already transformed, so its scale cannot be ",
             "checked against the reference's ", ref_unit, ". You are asserting ",
             "they match.")
    return(invisible(TRUE))
  }
  if (!identical(toupper(base), produced)) {
    tsf_abort("Unit mismatch: the reference was built on ", ref_unit,
              " but this query yields ", produced, ". Length normalisation ",
              "changes the spectrum, so the two are not comparable. Build a ",
              "reference on ", produced, ", or supply gene lengths so the query ",
              "can be expressed as ", base, ".")
  }
  invisible(TRUE)
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
