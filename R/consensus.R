# consensus.R -- the characteristic spectrum of a condition, built from the
# per-sample spectra rather than from the spectrum of the mean profile.
#
# WHY NOT THE SPECTRUM OF THE MEAN
# --------------------------------
# Averaging the profiles and transforming once is not the same as transforming
# each sample and summarising. The transform is linear, so the spectrum of the
# mean equals the mean of the complex coefficients -- a VECTOR mean. Components
# present in every sample at the same frequency but with scattered phases cancel
# in that sum and disappear from the condition spectrum, however reproducible
# they are. Conversely one extreme sample can carry a peak that no other sample
# has.
#
# So a condition is summarised here on three axes that the mean profile cannot
# separate:
#
#   how strong   median power across samples (robust to one outlier sample)
#   how common   prevalence: in what fraction of samples the frequency stands
#                out within that sample's own spectrum
#   how aligned  PLV = |mean(exp(i*phase))|, the phase-locking value: 1 when
#                every sample puts the crest in the same place, ~0 when phases
#                are scattered
#
# A frequency that is strong, common AND phase-locked is a candidate signature
# of the condition. Strong but not locked means each sample has structure at
# that scale in a different place -- real, but not a shared signature, and
# invisible in the spectrum of the mean.
#
# WHAT THE SCORE IS AND IS NOT
# ----------------------------
#   consensus_score = median_power_normalised * prevalence * PLV
#
# A product, so a component must satisfy all three: any factor near zero sends
# the score to zero. Normalised power (each frequency's share of its own
# sample-chromosome spectrum) rather than raw power, because raw power differs
# by orders of magnitude between chromosomes and the score would otherwise rank
# chromosomes instead of components.
#
# The score is a RANKING statistic, not a test. It has no null distribution and
# no error rate. Bootstrap intervals say how stable it is under resampling of
# samples; the Rayleigh p-value says whether the phase alignment alone is
# unlikely under uniform phases. Neither makes the score a significance claim --
# use condition_test.R for that.

#' Phase-locking value and its Rayleigh p-value.
#'
#' Under uniformly random phases E[PLV] is about sqrt(pi)/(2*sqrt(n)), NOT zero:
#' with 8 samples a PLV of 0.3 is unremarkable. The Rayleigh test is what turns
#' a PLV into a statement, and it is reported alongside so nobody reads a raw
#' PLV as evidence of alignment.
phase_locking <- function(phase) {
  p <- phase[is.finite(phase)]
  n <- length(p)
  if (n < 2L) return(c(plv = NA_real_, rayleigh_p = NA_real_, n = n))
  plv <- abs(mean(exp(1i * p)))
  # Rayleigh: R = n*plv^2; the standard small-sample correction of Zar (1999).
  z <- n * plv^2
  p_val <- exp(-z) * (1 + (2 * z - z^2) / (4 * n) -
                        (24 * z - 132 * z^2 + 76 * z^3 - 9 * z^4) / (288 * n^2))
  c(plv = plv, rayleigh_p = min(max(p_val, 0), 1), n = n)
}

#' Circular standard deviation, the natural heterogeneity measure for phase.
circular_sd <- function(phase) {
  p <- phase[is.finite(phase)]
  if (length(p) < 2L) return(NA_real_)
  r <- abs(mean(exp(1i * p)))
  if (r <= 0) return(Inf)
  sqrt(-2 * log(r))
}

#' Prevalence: in what fraction of samples this frequency stands out.
#'
#' "Stands out" is defined within each sample against its own spectrum for that
#' chromosome (above the given quantile of that sample's normalised power), so
#' prevalence does not depend on library depth or on how strong the sample is
#' overall. When per-sample maxT results are available, significance is the
#' better definition and is used instead.
prevalence_from_rank <- function(power_norm, sample_id, chr, quantile_cut = 0.95) {
  # The threshold is per sample AND per chromosome. Pooling chromosomes would
  # make them compete: chromosomes differ in length, coverage and total power,
  # so a short chromosome whose spectrum is flatter would never clear a
  # threshold set mostly by a long one, and prevalence would encode chromosome
  # identity instead of how much a frequency stands out.
  grp <- paste(sample_id, chr, sep = "\r")
  by_group <- split(seq_along(power_norm), grp)
  thr <- vapply(by_group, function(i)
    stats::quantile(power_norm[i], quantile_cut, na.rm = TRUE), numeric(1))
  power_norm > thr[grp]
}

#' Consensus spectrum for one condition, from its per-sample spectra.
#'
#' @param spectra_samples the spectra stage output for this condition
#' @param maxt optional per-sample maxT table; when present, prevalence is the
#'   fraction of samples in which the frequency is significant
consensus_spectrum <- function(spectra_samples, maxt = NULL, n_boot = 500L,
                               alpha = 0.05, seed = 42L, quantile_cut = 0.95) {
  d <- spectra_samples
  needed <- c("chr", "N", "k", "sample", "power", "amplitude", "phase")
  missing <- setdiff(needed, colnames(d))
  if (length(missing)) tsf_abort("consensus needs columns: ",
                                 paste(missing, collapse = ", "))
  if ("power_normalised" %in% colnames(d)) {
    d$pnorm <- d$power_normalised
  } else {
    tot <- stats::ave(d$power, paste(d$sample, d$chr), FUN = function(x) sum(x, na.rm = TRUE))
    d$pnorm <- d$power / pmax(tot, .Machine$double.eps)
  }

  d$stands_out <- if (!is.null(maxt) && "p_empirical_maxT" %in% colnames(maxt)) {
    key_d <- paste(d$chr, d$N, d$k, d$sample)
    key_m <- paste(maxt$chr, maxt$N, maxt$k, maxt$sample)
    sig <- maxt$p_empirical_maxT <= alpha
    out <- sig[match(key_d, key_m)]
    out[is.na(out)] <- FALSE
    out
  } else {
    prevalence_from_rank(d$pnorm, d$sample, d$chr, quantile_cut)
  }

  key <- paste(d$chr, d$N, d$k, sep = "|")
  groups <- split(seq_len(nrow(d)), key)
  n_samples_total <- length(unique(d$sample))
  set.seed(seed)

  rows <- lapply(names(groups), function(g) {
    i <- groups[[g]]
    ok <- is.finite(d$power[i]) & is.finite(d$phase[i])
    i <- i[ok]
    n_valid <- length(unique(d$sample[i]))
    if (n_valid < 2L) return(NULL)

    pl <- phase_locking(d$phase[i])
    med_p <- stats::median(d$pnorm[i], na.rm = TRUE)
    prev <- mean(d$stands_out[i], na.rm = TRUE)
    score <- med_p * prev * pl[["plv"]]

    # Bootstrap over SAMPLES (skipped when n_boot = 0, as in the null): the unit of replication is the sample, not the
    # frequency, so resampling anything else would understate the uncertainty.
    samples_here <- unique(d$sample[i])
    boot <- if (n_boot < 1L) matrix(NA_real_, nrow = 3, ncol = 1) else
      vapply(seq_len(n_boot), function(b) {
      pick <- sample(samples_here, replace = TRUE)
      idx <- unlist(lapply(pick, function(s) i[d$sample[i] == s]), use.names = FALSE)
      if (!length(idx)) return(c(NA_real_, NA_real_, NA_real_))
      mp <- stats::median(d$pnorm[idx], na.rm = TRUE)
      pv <- abs(mean(exp(1i * d$phase[idx])))
      c(mp, pv, mp * mean(d$stands_out[idx], na.rm = TRUE) * pv)
    }, numeric(3))

    qs <- function(v) stats::quantile(v, c(0.025, 0.975), na.rm = TRUE)
    ci_p <- qs(boot[1, ]); ci_v <- qs(boot[2, ]); ci_s <- qs(boot[3, ])

    parts <- strsplit(g, "|", fixed = TRUE)[[1]]
    data.frame(
      chr = parts[1], N = as.integer(parts[2]), k = as.integer(parts[3]),
      freq = as.integer(parts[3]) / as.integer(parts[2]),
      period = as.integer(parts[2]) / as.integer(parts[3]),
      n_samples_valid = n_valid, n_samples_total = n_samples_total,
      median_power = stats::median(d$power[i], na.rm = TRUE),
      median_power_normalised = med_p,
      median_amplitude = stats::median(d$amplitude[i], na.rm = TRUE),
      prevalence = prev,
      plv = pl[["plv"]], plv_rayleigh_p = pl[["rayleigh_p"]],
      mean_phase = Arg(mean(exp(1i * d$phase[i]))),
      power_heterogeneity = stats::mad(d$pnorm[i], na.rm = TRUE) /
        max(med_p, .Machine$double.eps),
      phase_heterogeneity = circular_sd(d$phase[i]),
      consensus_score = score,
      consensus_score_ci_lower = unname(ci_s[1]),
      consensus_score_ci_upper = unname(ci_s[2]),
      median_power_ci_lower = unname(ci_p[1]),
      median_power_ci_upper = unname(ci_p[2]),
      plv_ci_lower = unname(ci_v[1]), plv_ci_upper = unname(ci_v[2]),
      stringsAsFactors = FALSE)
  })

  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  out$plv_rayleigh_q <- stats::p.adjust(out$plv_rayleigh_p, method = "BH")
  out[order(-out$consensus_score), ]
}

#' Null distribution of the consensus score, by permuting condition labels.
#'
#' `consensus_score_ci_lower > 0` is nearly automatic: the score is a product of
#' non-negative quantities, so any signal at all clears zero. It says the score
#' is stable under resampling, not that it is larger than what an arbitrary
#' group of samples of the same size would produce.
#'
#' The null here is the right one for the question: draw n samples at random
#' from the whole dataset, ignoring condition, and compute the consensus score.
#' A component of a real condition has to beat that. Prevalence and phase
#' locking both survive in the null when they reflect tissue-wide structure
#' rather than the condition, which is exactly the confound worth removing.
null_consensus_scores <- function(spectra_all, n_samples, n_null = 100L,
                                  seed = 42L, quantile_cut = 0.95,
                                  q = 0.95) {
  samples <- unique(spectra_all$sample)
  if (length(samples) <= n_samples || n_null < 10L) return(NA_real_)
  set.seed(seed)
  best <- vapply(seq_len(n_null), function(b) {
    pick <- sample(samples, n_samples)
    sub <- spectra_all[spectra_all$sample %in% pick, , drop = FALSE]
    cs <- tryCatch(consensus_spectrum(sub, n_boot = 0L, seed = seed + b,
                                      quantile_cut = quantile_cut),
                   error = function(e) NULL)
    if (is.null(cs) || !nrow(cs)) return(NA_real_)
    max(cs$consensus_score, na.rm = TRUE)
  }, numeric(1))
  best <- best[is.finite(best)]
  if (!length(best)) return(NA_real_)
  unname(stats::quantile(best, q))
}

#' The characteristic signature of a condition.
#'
#' Selection is by the lower bootstrap bound rather than the point estimate, so
#' a component ranks on what survives resampling of the samples, and by phase
#' alignment that is unlikely under uniform phases.
consensus_signature <- function(cs, max_components = 50L, min_prevalence = 0.5,
                                plv_q = 0.05, null_score = NA_real_) {
  if (is.null(cs) || !nrow(cs)) return(NULL)

  # The Rayleigh p-value cannot fall below exp(-n) (attained at PLV = 1), so
  # after BH over n_freq frequencies the smallest reachable q is
  # exp(-n) * n_freq. With few samples that exceeds the threshold no matter how
  # perfect the alignment, and the signature comes back empty for a reason that
  # has nothing to do with the data. Same shape of problem as the permutation
  # floor in condition_test.R -- say so rather than returning an empty table.
  n_min <- min(cs$n_samples_valid, na.rm = TRUE)
  reachable <- exp(-n_min) * nrow(cs)
  if (reachable > plv_q) {
    tsf_warn("With ", n_min, " sample(s) and ", nrow(cs), " frequencies the ",
             "smallest reachable phase-alignment q is ", signif(reachable, 2),
             " > ", plv_q, ": perfect alignment could not pass. About ",
             ceiling(log(nrow(cs) / plv_q)), " samples are needed for this ",
             "condition. Reporting by prevalence and score only.")
    hit <- cs[cs$prevalence >= min_prevalence, , drop = FALSE]
    if (!nrow(hit)) return(NULL)
    hit$phase_alignment_testable <- FALSE
    hit$null_score_q95 <- null_score
    hit$signature_class <- "exploratory"
    hit <- hit[order(-hit$consensus_score_ci_lower), ]
    return(utils::head(hit, max_components))
  }
  keep <- cs$prevalence >= min_prevalence &
    !is.na(cs$plv_rayleigh_q) & cs$plv_rayleigh_q <= plv_q
  hit <- cs[keep, , drop = FALSE]
  if (!nrow(hit)) return(NULL)
  hit$phase_alignment_testable <- TRUE
  hit$null_score_q95 <- null_score
  # "confirmed" requires the bootstrap lower bound to clear the label-permuted
  # null, not merely zero. Clearing zero is nearly automatic for a product of
  # non-negative quantities; clearing the null means the component belongs to
  # this condition rather than to any group of samples of the same size.
  # Without a null, the strongest available statement is exploratory.
  beats_null <- if (is.na(null_score)) rep(FALSE, nrow(hit)) else
    hit$consensus_score_ci_lower > null_score
  hit$signature_class <- ifelse(beats_null & hit$plv_rayleigh_q <= plv_q,
                                "confirmed", "exploratory")
  if (is.na(null_score)) {
    tsf_warn("No permutation null was computed, so no component can be ",
             "confirmed: clearing zero is not evidence. Set consensus$n_null.")
  }
  hit <- hit[order(-hit$consensus_score_ci_lower), ]
  utils::head(hit, max_components)
}

#' Feature names of a signature, in the fingerprint's naming scheme.
signature_features <- function(sig, features = "amplitude") {
  if (is.null(sig) || !nrow(sig)) return(character(0))
  nm <- paste0("chr", sig$chr, "_k", sig$k)
  if (identical(features, "amplitude")) nm else c(paste0(nm, "_c"), paste0(nm, "_s"))
}
