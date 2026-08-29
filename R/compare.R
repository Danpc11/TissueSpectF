# compare.R -- constant signature, per-condition peaks, transitions, and the
# cross-dataset crossing.
#
# Two things changed from the original comparison script:
#
# 1. Transitions now actually measure a change. The old table named
#    *_cambian_magnitud.tsv contained peaks stable in BOTH conditions with no
#    amplitude change computed anywhere. Here transition_table() reports the
#    amplitude delta, its log2 ratio, and a per-sample Welch test on peak power,
#    and the file is named for what it holds.
#
# 2. The constant signature reports which conditions were actually available.
#    An empty intersection because one condition is missing is a different
#    finding from an empty intersection because no peak is shared, and the two
#    used to look identical.

wilson_ci <- function(n_sig, n_exp, conf = 0.95) {
  if (is.na(n_exp) || n_exp == 0 || is.na(n_sig)) return(c(NA_real_, NA_real_))
  z <- stats::qnorm(1 - (1 - conf) / 2)
  p_hat <- n_sig / n_exp
  denom <- 1 + z^2 / n_exp
  center <- (p_hat + z^2 / (2 * n_exp)) / denom
  margin <- (z / denom) * sqrt(p_hat * (1 - p_hat) / n_exp + z^2 / (4 * n_exp^2))
  c(max(0, center - margin) * 100, min(1, center + margin) * 100)
}

add_wilson <- function(df, n_sig_col, n_exp_col, prefix = "") {
  ci <- t(mapply(wilson_ci, df[[n_sig_col]], df[[n_exp_col]]))
  df[[paste0(prefix, "ci95_lower_pct")]] <- ci[, 1]
  df[[paste0(prefix, "ci95_upper_pct")]] <- ci[, 2]
  df
}

peak_key <- function(df) paste(df$chr, df$N, df$k, sep = "|")

#' Peaks stable in every available condition of one dataset.
constant_signature <- function(stability_by_cond, conditions) {
  available <- conditions[vapply(conditions, function(c)
    !is.null(stability_by_cond[[c]]) && nrow(stability_by_cond[[c]]) > 0, logical(1))]
  missing <- setdiff(conditions, available)
  if (length(missing)) {
    tsf_warn("No stable-peak table for: ", paste(missing, collapse = ", "),
             " -- the signature is over ", length(available), " condition(s), not ",
             length(conditions))
  }
  if (!length(available)) return(NULL)

  stable <- lapply(available, function(c) {
    s <- stability_by_cond[[c]]
    s[s$is_stable %in% TRUE, , drop = FALSE]
  })
  keys <- Reduce(intersect, lapply(stable, peak_key))
  if (!length(keys)) {
    tsf_log("No peak is stable in all of: ", paste(available, collapse = ", "))
    return(NULL)
  }

  base <- stable[[1]][peak_key(stable[[1]]) %in% keys, c("chr", "N", "k", "freq", "period")]
  base$n_conditions <- length(available)
  base$conditions <- paste(available, collapse = ",")
  pull <- function(col) {
    m <- vapply(stable, function(s) s[[col]][match(peak_key(base), peak_key(s))],
                numeric(nrow(base)))
    if (is.null(dim(m))) m <- matrix(m, nrow = nrow(base))
    rowSums(m)
  }
  n_sig <- pull("n_samples_significant")
  n_exp <- pull("n_samples_expected")
  base$n_samples_significant_total <- n_sig
  base$n_samples_expected_total <- n_exp
  base$pct_samples_significant <- 100 * n_sig / n_exp
  base <- add_wilson(base, "n_samples_significant_total", "n_samples_expected_total")
  base[order(-base$ci95_lower_pct), ]
}

#' Peaks stable in both conditions of a transition, WITH the change measured.
transition_table <- function(cond_from, cond_to, peaks_by_cond, stability_by_cond,
                             maxt_by_cond, branch) {
  a <- peaks_by_cond[[cond_from]]
  b <- peaks_by_cond[[cond_to]]
  if (is.null(a) || is.null(b) || !nrow(a) || !nrow(b)) return(NULL)
  keys <- intersect(peak_key(a), peak_key(b))
  if (!length(keys)) return(NULL)

  amp_col <- if (branch == "average") "amplitude_mean" else "amplitude_median"
  a <- a[match(keys, peak_key(a)), ]
  b <- b[match(keys, peak_key(b)), ]

  out <- data.frame(
    chr = a$chr, N = a$N, k = a$k, freq = a$freq, period = a$period,
    transition = paste0(cond_from, "_to_", cond_to),
    branch = branch,
    amplitude_from = a[[amp_col]], amplitude_to = b[[amp_col]],
    amplitude_delta = b[[amp_col]] - a[[amp_col]],
    log2_amplitude_ratio = log2(b[[amp_col]] / a[[amp_col]]),
    pct_significant_from = a$pct_samples_significant,
    pct_significant_to = b$pct_samples_significant,
    n_samples_significant_from = a$n_samples_significant,
    n_samples_expected_from = a$n_samples_expected,
    n_samples_significant_to = b$n_samples_significant,
    n_samples_expected_to = b$n_samples_expected,
    stringsAsFactors = FALSE)

  # Per-sample power in each condition gives a test with real degrees of
  # freedom, instead of comparing two single summary amplitudes.
  test <- lapply(seq_len(nrow(out)), function(i) {
    pw <- function(cond) {
      m <- maxt_by_cond[[cond]]
      if (is.null(m)) return(numeric(0))
      m$power[m$chr == out$chr[i] & m$N == out$N[i] & m$k == out$k[i]]
    }
    x <- pw(cond_from); y <- pw(cond_to)
    if (length(x) < 3 || length(y) < 3) return(c(NA_real_, NA_real_, NA_real_))
    tt <- suppressWarnings(stats::t.test(y, x))
    c(mean(x), mean(y), unname(tt$p.value))
  })
  test <- do.call(rbind, test)
  if (all(is.na(test[, 3]))) {
    tsf_warn("No per-sample maxT for ", cond_from, "/", cond_to,
             ": the power change is reported without a test, and `replicated` ",
             "cannot be established. Run the maxt stage for that claim.")
  }
  out$mean_power_from <- test[, 1]
  out$mean_power_to <- test[, 2]
  out$p_power_welch <- test[, 3]
  out$p_power_fdr <- stats::p.adjust(out$p_power_welch, method = "BH")

  out <- add_wilson(out, "n_samples_significant_from", "n_samples_expected_from", "from_")
  out <- add_wilson(out, "n_samples_significant_to", "n_samples_expected_to", "to_")
  out[order(out$p_power_welch, -abs(out$log2_amplitude_ratio)), ]
}

#' Cross two datasets on the peaks they share, keeping them separate.
#'
#' No pooling of proportions: the two datasets are reported side by side because
#' the question is replication, not a bigger n. A combined percentage would be
#' dominated by whichever cohort is larger and would hide heterogeneity.
cross_datasets <- function(tables_by_dataset, key_cols) {
  nms <- names(tables_by_dataset)
  ok <- vapply(tables_by_dataset, function(t) !is.null(t) && nrow(t) > 0, logical(1))
  if (sum(ok) < 2) {
    tsf_warn("Nothing to cross: fewer than two non-empty tables (",
             paste(nms[!ok], collapse = ", "), " empty)")
    return(NULL)
  }
  if (any(!ok)) tsf_warn("Excluded from the crossing (empty): ",
                         paste(nms[!ok], collapse = ", "))
  nms <- nms[ok]

  suffixed <- lapply(nms, function(nm) {
    t <- tables_by_dataset[[nm]]
    kc <- intersect(key_cols, colnames(t))
    colnames(t)[!colnames(t) %in% kc] <-
      paste0(colnames(t)[!colnames(t) %in% kc], "_", nm)
    t
  })
  names(suffixed) <- nms
  common_keys <- Reduce(intersect, lapply(suffixed, function(t)
    intersect(key_cols, colnames(t))))
  if (!length(common_keys)) {
    tsf_warn("No key column is shared by every dataset")
    return(NULL)
  }
  out <- Reduce(function(a, b) merge(a, b, by = common_keys), suffixed)
  sizes <- vapply(suffixed, nrow, integer(1))
  tsf_log("Shared peaks across ", length(nms), " dataset(s): ", nrow(out),
          " (of ", paste(sizes, collapse = ", "), ")")
  out
}

#' Do the two datasets agree on the direction of change across a transition?
add_replication_flags <- function(crossed, nms) {
  present <- nms[paste0("log2_amplitude_ratio_", nms) %in% colnames(crossed)]
  if (length(present) < 2) return(crossed)

  deltas <- vapply(present, function(nm)
    crossed[[paste0("log2_amplitude_ratio_", nm)]], numeric(nrow(crossed)))
  if (is.null(dim(deltas))) deltas <- matrix(deltas, nrow = nrow(crossed))
  finite_all <- apply(is.finite(deltas), 1, all)
  crossed$same_direction <- finite_all &
    apply(sign(deltas), 1, function(v) length(unique(v)) == 1L)

  pcols <- paste0("p_power_fdr_", present)
  if (all(pcols %in% colnames(crossed))) {
    ps <- vapply(pcols, function(cl) crossed[[cl]], numeric(nrow(crossed)))
    if (is.null(dim(ps))) ps <- matrix(ps, nrow = nrow(crossed))
    crossed$replicated <- crossed$same_direction &
      apply(ps, 1, function(v) all(!is.na(v) & v <= 0.05))
    tsf_log("Replicated in all ", length(present),
            " dataset(s) (same direction, FDR <= 0.05): ",
            sum(crossed$replicated, na.rm = TRUE))
    crossed <- crossed[order(-crossed$same_direction, crossed[[pcols[1]]]), ]
  }
  crossed
}
