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
                               alpha = 0.05, seed = 42L, quantile_cut = 0.95,
                               n_cores = 1L) {
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

  # TWO prevalences, always both, never one standing in for the other.
  #
  # The rank definition ("above this sample's own 95th percentile for this
  # chromosome") needs nothing but the spectra, so it can be computed for the
  # observed data and for any permuted draw alike. The maxT definition
  # ("significant in that sample") is stronger evidence but exists only where
  # per-sample maxT was run, which the permutation null cannot assume.
  #
  # Mixing them would compare an observed score built on one statistic against a
  # null built on another: not a permutation test, whatever the direction of the
  # bias. So the permuted comparison uses consensus_score_rank on both sides,
  # and consensus_score_maxt is reported next to it as confirmatory evidence.
  d$stands_out_rank <- prevalence_from_rank(d$pnorm, d$sample, d$chr, quantile_cut)
  has_maxt <- !is.null(maxt) && "p_empirical_maxT" %in% colnames(maxt)
  d$stands_out_maxt <- if (has_maxt) {
    key_d <- paste(d$chr, d$N, d$k, d$sample)
    key_m <- paste(maxt$chr, maxt$N, maxt$k, maxt$sample)
    sig <- maxt$p_empirical_maxT <= alpha
    out <- sig[match(key_d, key_m)]
    out[is.na(out)] <- FALSE
    out
  } else rep(NA, nrow(d))

  key <- paste(d$chr, d$N, d$k, sep = "|")
  groups <- split(seq_len(nrow(d)), key)
  n_samples_total <- length(unique(d$sample))
  group_names <- names(groups)
  n_cores <- max(1L, min(as.integer(n_cores), length(group_names)))

  # Parallelise the expensive observed bootstrap over frequency groups. Each
  # group receives its own deterministic seed, so results do not change with
  # core count or scheduling order. The permutation null calls this function
  # with n_boot = 0 and n_cores = 1 to avoid nested parallelism.
  rows <- parallel::mclapply(seq_along(group_names), function(group_index) {
    g <- group_names[[group_index]]
    set.seed(as.integer((as.double(seed) + group_index) %% .Machine$integer.max))
    i <- groups[[g]]
    ok <- is.finite(d$power[i]) & is.finite(d$phase[i])
    i <- i[ok]
    n_valid <- length(unique(d$sample[i]))
    if (n_valid < 2L) return(NULL)

    pl <- phase_locking(d$phase[i])
    med_p <- stats::median(d$pnorm[i], na.rm = TRUE)
    prev_rank <- mean(d$stands_out_rank[i], na.rm = TRUE)
    prev_maxt <- if (has_maxt) mean(d$stands_out_maxt[i], na.rm = TRUE) else NA_real_
    score_rank <- med_p * prev_rank * pl[["plv"]]
    score_maxt <- if (has_maxt) med_p * prev_maxt * pl[["plv"]] else NA_real_

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
      c(mp, pv, mp * mean(d$stands_out_rank[idx], na.rm = TRUE) * pv)
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
      prevalence = if (has_maxt) prev_maxt else prev_rank,
      prevalence_rank = prev_rank, prevalence_maxt = prev_maxt,
      plv = pl[["plv"]], plv_rayleigh_p = pl[["rayleigh_p"]],
      mean_phase = Arg(mean(exp(1i * d$phase[i]))),
      power_heterogeneity = stats::mad(d$pnorm[i], na.rm = TRUE) /
        max(med_p, .Machine$double.eps),
      phase_heterogeneity = circular_sd(d$phase[i]),
      consensus_score = if (has_maxt) score_maxt else score_rank,
      consensus_score_rank = score_rank,
      consensus_score_maxt = score_maxt,
      consensus_score_ci_lower = unname(ci_s[1]),
      consensus_score_ci_upper = unname(ci_s[2]),
      median_power_ci_lower = unname(ci_p[1]),
      median_power_ci_upper = unname(ci_p[2]),
      plv_ci_lower = unname(ci_v[1]), plv_ci_upper = unname(ci_v[2]),
      stringsAsFactors = FALSE)
  }, mc.cores = n_cores, mc.preschedule = TRUE, mc.set.seed = FALSE)

  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  out$plv_rayleigh_q <- stats::p.adjust(out$plv_rayleigh_p, method = "BH")
  out[order(-out$consensus_score_rank), ]
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
#' @param blocks optional named vector mapping sample id -> block. When given,
#'   draws are made of whole blocks rather than of individual samples.
#'
#' Blocking matters whenever samples are not independent: several biopsies from
#' one subject, longitudinal measurements, technical batches, tumour-normal
#' pairs, multiple regions of one organ. Drawing samples freely from such a
#' dataset builds a null in which a subject's own correlated samples rarely
#' land together, while the observed condition may consist largely of them. The
#' null then looks more variable than the data and the test is anti-conservative
#' in exactly the situation where independence fails.
#' Prepare frequency x sample matrices once for the permutation null.
#'
#' The original null rebuilt a data frame, split ~10,000 frequency groups and
#' called consensus_spectrum() for every draw. None of that structure changes
#' between draws. This representation pays the reshape cost once and lets every
#' permutation reduce three shared matrices: normalised power, rank prevalence
#' and unit phase vectors.
#' @param maxt optional per-sample maxT table. When supplied, a SECOND standing
#'   matrix is built from per-sample significance, so the null can produce a
#'   maxT-based score as well as a rank-based one.
#'
#'   Why it can: the null draws random SUBSETS OF REAL SAMPLES. Significance was
#'   already decided per (frequency, sample) when maxT ran, so the prevalence of
#'   a draw is just the fraction of its samples that were significant -- no maxT
#'   is recomputed, and the cost is one extra logical matrix.
#'
#'   Why it matters: without it, a peak is SELECTED by maxT prevalence and then
#'   TESTED against a null built from rank prevalence. Selecting on one statistic
#'   and testing on another is not a permutation test of the thing selected, and
#'   it is the reason `consensus_score_maxt` could only ever be reported as
#'   "confirmatory" rather than as the quantity carrying a p-value.
prepare_null_consensus_matrices <- function(spectra_all, quantile_cut = 0.95,
                                            maxt = NULL, alpha = 0.05) {
  needed <- c("chr", "N", "k", "sample", "power", "phase")
  missing <- setdiff(needed, colnames(spectra_all))
  if (length(missing)) tsf_abort("null consensus needs columns: ",
                                 paste(missing, collapse = ", "))
  d <- spectra_all
  if ("power_normalised" %in% colnames(d)) {
    pnorm <- d$power_normalised
  } else {
    tot <- stats::ave(d$power, paste(d$sample, d$chr),
                      FUN = function(x) sum(x, na.rm = TRUE))
    pnorm <- d$power / pmax(tot, .Machine$double.eps)
  }
  stands <- prevalence_from_rank(pnorm, d$sample, d$chr, quantile_cut)

  stands_maxt <- NULL
  if (!is.null(maxt) && "p_empirical_maxT" %in% colnames(maxt)) {
    km <- paste(maxt$chr, maxt$N, maxt$k, maxt$sample)
    sig <- maxt$p_empirical_maxT <= alpha
    stands_maxt <- sig[match(paste(d$chr, d$N, d$k, d$sample), km)]
    # A frequency/sample pair maxT never scored is not significant. NA would
    # propagate into the prevalence and silently shrink the denominator.
    stands_maxt[is.na(stands_maxt)] <- FALSE
  }

  key <- paste(d$chr, d$N, d$k, sep = "|")
  keys <- sort(unique(key))             # same key order as split()
  samples <- unique(as.character(d$sample))
  ri <- match(key, keys); ci <- match(as.character(d$sample), samples)
  cell <- paste(ri, ci, sep = "|")
  if (anyDuplicated(cell)) {
    tsf_abort("The null matrix has duplicated frequency/sample rows; expected ",
              "one value per chr/N/k/sample")
  }

  dims <- c(length(keys), length(samples))
  pn <- matrix(NA_real_, nrow = dims[1], ncol = dims[2],
               dimnames = list(keys, samples))
  so <- matrix(NA_real_, nrow = dims[1], ncol = dims[2],
               dimnames = list(keys, samples))
  ph <- matrix(NA_complex_, nrow = dims[1], ncol = dims[2],
               dimnames = list(keys, samples))
  valid <- is.finite(d$power) & is.finite(d$phase)
  pos <- cbind(ri[valid], ci[valid])
  pn[pos] <- pnorm[valid]
  so[pos] <- as.numeric(stands[valid])
  ph[pos] <- exp(1i * d$phase[valid])

  sm <- NULL
  if (!is.null(stands_maxt)) {
    sm <- matrix(NA_real_, nrow = dims[1], ncol = dims[2],
                 dimnames = list(keys, samples))
    sm[pos] <- as.numeric(stands_maxt[valid])
  }

  list(pnorm = pn, stands = so, stands_maxt = sm, phase = ph,
       keys = keys, samples = samples)
}

row_medians_tsf <- function(x) {
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::rowMedians(x, na.rm = TRUE)
  } else {
    apply(x, 1L, stats::median, na.rm = TRUE)
  }
}

#' Consensus scores for one selected set of matrix columns.
null_matrix_draw <- function(prepared, picked_samples) {
  j <- match(picked_samples, prepared$samples)
  j <- j[!is.na(j)]
  if (length(j) < 2L) return(NULL)
  pn <- prepared$pnorm[, j, drop = FALSE]
  so <- prepared$stands[, j, drop = FALSE]
  ph <- prepared$phase[, j, drop = FALSE]
  n_valid <- rowSums(!is.na(ph))
  med <- row_medians_tsf(pn)
  prev <- rowMeans(so, na.rm = TRUE)
  plv <- Mod(rowMeans(ph, na.rm = TRUE))
  score <- med * prev * plv

  # The same score built on maxT prevalence, when maxT was available. This is
  # what lets a peak be tested with the statistic it was selected by.
  score_maxt <- NULL
  if (!is.null(prepared$stands_maxt)) {
    prev_m <- rowMeans(prepared$stands_maxt[, j, drop = FALSE], na.rm = TRUE)
    score_maxt <- med * prev_m * plv
    score_maxt[n_valid < 2L | !is.finite(score_maxt)] <- NA_real_
    names(score_maxt) <- prepared$keys
  }
  score[n_valid < 2L | !is.finite(score)] <- NA_real_
  plv[n_valid < 2L | !is.finite(plv)] <- NA_real_
  names(score) <- prepared$keys
  names(plv) <- prepared$keys
  # PLV under the null is returned, not discarded. The Rayleigh test that used
  # to stand alone assumes phases are iid uniform across samples; these samples
  # share one grid and one tissue, so their phases agree for structural reasons
  # and 98.9% of frequencies cleared plv_rayleigh_q <= 0.05 -- a gate that
  # passes everything is not a gate. Each null draw takes a random subset of the
  # same pool, so it carries that shared structure too, and calibrating the
  # observed PLV against it asks the question the Rayleigh test could not: is
  # this frequency more phase-coherent than the same tissue on the same grid
  # produces by itself?
  list(score = score, plv = plv, score_maxt = score_maxt)
}

#' Permutation null for the consensus score.
#'
#' RETAINED_KEYS: which frequencies form the tested family.
#'
#' The period floor removes frequencies from `cs` before any p-value is
#' computed, but the null was still built over every frequency, so each
#' permutation's global maximum was the maximum over a family larger than the
#' one being tested. p_null_fwer was therefore the family-wise error rate of a
#' family that included frequencies which were never candidates. The direction
#' is conservative -- no false positives -- but it costs real power, and the
#' claim that the filter is applied "before the null" was not true of the FWER
#' route.
#'
#' The fix is deliberately narrow. The statistic is still computed on the full
#' prepared matrices, because `prevalence_from_rank()` defines standing as the
#' top (1 - quantile_cut) within each (sample, chromosome): preparing the null
#' on a filtered pool would recompute that threshold over a smaller set, and
#' the observed `cs` came from the unfiltered per-sample spectra. Observed and
#' null would stop being comparable, which is a worse error than the one being
#' fixed.
#'
#' So: compute over everything, then restrict WHICH frequencies compete. Both
#' the pointwise null and the per-draw maximum are taken over `retained_keys`.
null_consensus_distribution <- function(spectra_all, n_samples, n_null = 50L,
                                        seed = 42L, quantile_cut = 0.95,
                                        q_global = 0.95, blocks = NULL,
                                        n_cores = 1L,
                                        engine = c("matrix", "reference"),
                                        prepared = NULL,
                                        retained_keys = NULL) {
  engine <- match.arg(engine)
  samples <- unique(spectra_all$sample)
  if (length(samples) <= n_samples || n_null < 10L) return(NULL)
  # When the condition is most of the dataset, every null draw overlaps the
  # observed set heavily and the null is nearly the observed statistic: the
  # expected overlap of a random draw of n from M with a fixed set of n is
  # n^2/M. With 42 of 55 samples that is 32 of 42, so p_null_fwer cannot get
  # small however strong the component. Say so: that is "not reachable", which
  # is a different finding from "not present".
  frac <- n_samples / length(samples)
  if (frac > 0.5) {
    tsf_warn("Consensus null: the condition holds ", n_samples, " of ",
             length(samples), " samples (", round(100 * frac), "%). A random ",
             "draw shares about ", round(n_samples * frac), " of them with the ",
             "observed set, so the null is close to the observed score and ",
             "p_null / p_null_fwer are conservative to the point of being ",
             "uninformative. Treat an unconfirmed component here as NOT ",
             "REACHABLE, not as absent.")
  }
  set.seed(seed)

  draw <- if (is.null(blocks)) {
    function() sample(samples, n_samples)
  } else {
    blk <- blocks[samples]
    blk[is.na(blk)] <- paste0("_singleton_", which(is.na(blk)))
    by_block <- split(samples, blk)
    function() {
      picked <- character(0)
      order_b <- sample(names(by_block))
      for (b in order_b) {
        if (length(picked) >= n_samples) break
        picked <- c(picked, by_block[[b]])
      }
      utils::head(picked, n_samples)
    }
  }

  # Generate every draw in the parent process. This makes the null exactly
  # reproducible for a fixed seed regardless of worker count.
  picks <- lapply(seq_len(n_null), function(b) draw())
  n_cores <- max(1L, min(as.integer(n_cores), n_null))

  prepared <- if (engine == "matrix" && is.null(prepared))
    prepare_null_consensus_matrices(spectra_all, quantile_cut) else prepared
  rows_by_sample <- if (engine == "reference")
    split(seq_len(nrow(spectra_all)), spectra_all$sample) else NULL

  draws <- parallel::mclapply(seq_len(n_null), function(b) {
    pick <- picks[[b]]
    plv_draw <- NULL
    maxt_draw <- NULL
    values <- if (engine == "matrix") {
      d <- tryCatch(null_matrix_draw(prepared, pick), error = function(e) NULL)
      if (is.null(d)) NULL else {
        plv_draw <- d$plv; maxt_draw <- d$score_maxt; d$score
      }
    } else {
      idx <- unlist(rows_by_sample[pick], use.names = FALSE)
      sub <- spectra_all[idx, , drop = FALSE]
      cs <- tryCatch(consensus_spectrum(sub, n_boot = 0L, seed = seed + b,
                                        quantile_cut = quantile_cut,
                                        n_cores = 1L),
                     error = function(e) NULL)
      if (is.null(cs) || !nrow(cs)) NULL else
        stats::setNames(cs$consensus_score_rank,
                        paste(cs$chr, cs$N, cs$k, sep = "|"))
    }
    values <- values[is.finite(values)]
    # Restrict to the tested family before the maximum is taken. Applying it
    # here rather than to the input keeps the statistic identical to the
    # observed one and makes both the pointwise and family-wise nulls describe
    # the same family as `cs`.
    if (!is.null(retained_keys)) {
      values <- values[names(values) %in% retained_keys]
      if (!is.null(plv_draw)) {
        plv_draw <- plv_draw[names(plv_draw) %in% retained_keys]
      }
      if (!is.null(maxt_draw)) {
        maxt_draw <- maxt_draw[names(maxt_draw) %in% retained_keys]
      }
    }
    if (!length(values)) return(NULL)
    # plv_null: the phase coherence this tissue, on this grid, produces from a
    # random subset of the same pool. It is the reference the Rayleigh test
    # should have been.
    list(values = values, best = max(values),
         plv = if (is.null(plv_draw)) NULL else plv_draw[is.finite(plv_draw)],
         maxt = if (is.null(maxt_draw)) NULL else maxt_draw[is.finite(maxt_draw)],
         best_maxt = if (is.null(maxt_draw) || !any(is.finite(maxt_draw))) NA_real_
                     else max(maxt_draw, na.rm = TRUE))
  }, mc.cores = n_cores, mc.preschedule = TRUE, mc.set.seed = FALSE)

  draws <- draws[!vapply(draws, is.null, logical(1))]
  if (!length(draws)) return(NULL)
  per_key <- lapply(draws, `[[`, "values")
  best <- vapply(draws, `[[`, numeric(1), "best")
  if (!length(per_key)) return(NULL)

  keys <- Reduce(union, lapply(per_key, names))
  mat <- vapply(per_key, function(v) v[keys], numeric(length(keys)))
  if (is.null(dim(mat))) mat <- matrix(mat, nrow = length(keys))
  rownames(mat) <- keys

  # Per-key PLV null, assembled the same way as the score null.
  plv_draws <- lapply(draws, `[[`, "plv")
  plv_draws <- plv_draws[!vapply(plv_draws, is.null, logical(1))]
  plv_mat <- NULL
  if (length(plv_draws)) {
    plv_keys <- Reduce(union, lapply(plv_draws, names))
    plv_mat <- vapply(plv_draws, function(v) v[plv_keys],
                      numeric(length(plv_keys)))
    if (is.null(dim(plv_mat))) plv_mat <- matrix(plv_mat, nrow = length(plv_keys))
    rownames(plv_mat) <- plv_keys
  }

  # The maxT-based null, assembled the same way. Present only when maxT reached
  # the consensus stage; NULL otherwise, and the p-values fall back to the
  # rank route.
  maxt_draws <- lapply(draws, `[[`, "maxt")
  maxt_draws <- maxt_draws[!vapply(maxt_draws, is.null, logical(1))]
  maxt_mat <- NULL
  best_maxt <- NULL
  if (length(maxt_draws)) {
    mk <- Reduce(union, lapply(maxt_draws, names))
    maxt_mat <- vapply(maxt_draws, function(v) v[mk], numeric(length(mk)))
    if (is.null(dim(maxt_mat))) maxt_mat <- matrix(maxt_mat, nrow = length(mk))
    rownames(maxt_mat) <- mk
    best_maxt <- vapply(draws, function(d) d$best_maxt %||% NA_real_, numeric(1))
    best_maxt <- best_maxt[is.finite(best_maxt)]
  }

  list(per_key = mat, global = unname(stats::quantile(best, q_global)),
       global_max_draws = best, per_key_plv = plv_mat,
       per_key_maxt = maxt_mat, global_max_draws_maxt = best_maxt,
       n_null = ncol(mat), n_samples = n_samples)
}

#' Empirical p-value per component against its own null, then BH.
#'
#' The global maximum null controls the error rate over the whole signature and
#' is very strict -- in practice nothing is confirmed under it, which is correct
#' for "is there any component at all" and useless for "which components".
#' Keeping the per-(chr, k) null as well gives a p-value per component; BH over
#' those is the intermediate that lets a signature be localised. Both are
#' reported, and which one a claim rests on is a stated choice.
#' PLV calibrated against the permutation null rather than against Rayleigh.
#'
#' plv_rayleigh_p assumes phases iid uniform across samples. These samples share
#' one reference grid and one tissue, so their phases agree for structural
#' reasons: on a real run the PLV distribution over 1995 frequencies had median
#' 0.971 and q25 0.914, and 98.9% of frequencies cleared plv_rayleigh_q <= 0.05.
#' A gate that passes 98.9% of candidates is not a gate, and a Rayleigh p of
#' 1e-12 attached to a median frequency is not evidence of anything.
#'
#' Each null draw is a random subset of the same pool, so it reproduces the
#' shared grid and the shared tissue. Calibrating against it asks the question
#' the Rayleigh test cannot: is this frequency more phase-coherent than this
#' tissue on this grid produces on its own?
#'
#' plv_rayleigh_p is kept, unchanged, because removing a column silently changes
#' what old results mean. It should not be used as a selection gate.
plv_null_pvalues <- function(cs, null_dist) {
  cs$p_plv_null <- NA_real_
  cs$q_plv_null <- NA_real_
  cs$plv_null_median <- NA_real_
  if (is.null(null_dist) || is.null(null_dist$per_key_plv)) return(cs)

  key <- paste(cs$chr, cs$N, cs$k, sep = "|")
  m <- null_dist$per_key_plv
  idx <- match(key, rownames(m))

  for (i in which(!is.na(idx))) {
    draws <- m[idx[i], ]
    draws <- draws[is.finite(draws)]
    if (length(draws) < 10L) next
    obs <- cs$plv[i]
    if (!is.finite(obs)) next
    cs$p_plv_null[i] <- (1 + sum(draws >= obs)) / (length(draws) + 1)
    cs$plv_null_median[i] <- stats::median(draws)
  }
  ok <- is.finite(cs$p_plv_null)
  if (any(ok)) cs$q_plv_null[ok] <- stats::p.adjust(cs$p_plv_null[ok], "BH")
  cs
}

#' p-values for the maxT-based consensus score, against the maxT-based null.
#'
#' Selecting a peak by maxT prevalence and testing it against a rank-based null
#' compares an observed score built on one statistic with a null built on
#' another. It is conservative or anticonservative depending on how the two
#' statistics differ on that frequency, and either way it is not a permutation
#' test of the thing that was selected. With both nulls available the peak can
#' finally be tested with the statistic it was chosen by.
#'
#' The rank-based columns are kept unchanged beside these, so a result computed
#' before this existed still means what it said.
maxt_null_pvalues <- function(cs, null_dist) {
  cs$p_null_maxt <- NA_real_
  cs$q_null_maxt <- NA_real_
  cs$p_null_fwer_maxt <- NA_real_
  if (is.null(null_dist) || is.null(null_dist$per_key_maxt)) return(cs)
  if (!"consensus_score_maxt" %in% names(cs)) return(cs)

  m <- null_dist$per_key_maxt
  idx <- match(paste(cs$chr, cs$N, cs$k, sep = "|"), rownames(m))
  for (i in which(!is.na(idx))) {
    dr <- m[idx[i], ]; dr <- dr[is.finite(dr)]
    obs <- cs$consensus_score_maxt[i]
    if (!length(dr) || !is.finite(obs)) next
    cs$p_null_maxt[i] <- (1 + sum(dr >= obs)) / (length(dr) + 1)
  }
  ok <- is.finite(cs$p_null_maxt)
  if (any(ok)) cs$q_null_maxt[ok] <- stats::p.adjust(cs$p_null_maxt[ok], "BH")

  # Family-wise: against the per-draw maximum, exactly as the rank route does.
  gm <- null_dist$global_max_draws_maxt
  if (!is.null(gm) && length(gm)) {
    obs <- cs$consensus_score_maxt
    cs$p_null_fwer_maxt <- vapply(obs, function(o)
      if (!is.finite(o)) NA_real_ else (1 + sum(gm >= o)) / (length(gm) + 1),
      numeric(1))
  }
  cs
}

null_component_pvalues <- function(cs, null_dist, null_q = 0.05) {
  cs <- plv_null_pvalues(cs, null_dist)
  cs <- maxt_null_pvalues(cs, null_dist)
  if (is.null(null_dist)) {
    cs$p_null_fwer <- NA_real_
    cs$p_null <- NA_real_; cs$q_null <- NA_real_
    cs$beats_global_null <- NA
    return(cs)
  }
  key <- paste(cs$chr, cs$N, cs$k, sep = "|")
  n_b <- null_dist$n_null
  cs$p_null <- vapply(seq_len(nrow(cs)), function(i) {
    row <- null_dist$per_key[key[i], ]
    row <- row[is.finite(row)]
    if (!length(row)) return(NA_real_)
    (1 + sum(row >= cs$consensus_score_rank[i])) / (length(row) + 1)
  }, numeric(1))
  cs$q_null <- stats::p.adjust(cs$p_null, method = "BH")

  # PRIMARY: the maxT-style p-value against the distribution of the null's
  # global maximum. Comparing an observed score against the largest score any
  # random draw produced anywhere controls the family-wise error rate across all
  # frequencies by construction, so nothing is adjusted afterwards and the floor
  # is 1/(n_null+1) -- reachable with 50 draws.
  #
  # The pointwise p-value above cannot be used this way. Its own floor is also
  # 1/(n_null+1), but BH across ~n_f frequencies pushes the smallest reachable q
  # to n_f/(n_null+1): with 297 frequencies and 50 draws that is 5.8, so no
  # component could ever be confirmed however strong. Confirming on it would
  # need n_null >= n_f/null_q, which is thousands of draws for one chromosome
  # and far more for a genome. It is kept as a secondary, localising statistic,
  # with the reachability warning below.
  gmax <- null_dist$global_max_draws
  cs$p_null_fwer <- if (is.null(gmax) || !length(gmax)) NA_real_ else
    vapply(cs$consensus_score_rank, function(x)
      (1 + sum(gmax >= x)) / (length(gmax) + 1), numeric(1))

  # A conservative diagnostic, not a bound. m/(n_b+1) is the BH value at rank 1
  # with no ties; with k p-values tied at the permutation floor the achievable
  # value is m/(k*(n_b+1)), which for a large k is orders of magnitude smaller.
  # An earlier version reported this as "the smallest reachable q" and concluded
  # that nothing could pass, which could understate the method's power by a
  # factor of m. So: report both, and say which is which.
  rank1 <- bh_rank1_diagnostic(nrow(cs), n_b)
  n_at_floor <- sum(cs$p_null <= 1 / (n_b + 1), na.rm = TRUE)
  achievable <- bh_achievable_q(nrow(cs), n_b, max(1L, n_at_floor))
  if (rank1 > null_q) {
    tsf_log("Pointwise null over ", nrow(cs), " frequencies at ", n_b, " draws: ",
            "rank-1 BH diagnostic ", signif(rank1, 3), " (conservative, assumes ",
            "no ties). ", n_at_floor, " p-value(s) sit at the permutation floor, ",
            "so the achievable BH q is ", signif(achievable, 3),
            if (achievable <= null_q)
              paste0(" -- at or below ", null_q, ", the pointwise route is usable.")
            else
              paste0(" -- still above ", null_q, ". Family-wise p_null_fwer is ",
                     "reported alongside; ", draws_for_bh(nrow(cs), null_q),
                     " draws would clear the rank-1 case, or reduce the family ",
                     "with --min-period."))
  }

  cs$beats_global_null <- cs$consensus_score_ci_lower > null_dist$global
  cs$null_global_q95 <- null_dist$global
  cs$n_null <- n_b
  cs
}

#' The characteristic signature of a condition.
#'
#' Selection is by the lower bootstrap bound rather than the point estimate, so
#' a component ranks on what survives resampling of the samples, and by phase
#' alignment that is unlikely under uniform phases.
consensus_signature <- function(cs, max_components = 50L, min_prevalence = 0.5,
                                plv_q = 0.05, null_q = 0.05) {
  if (is.null(cs) || !nrow(cs)) return(NULL)

  # The Rayleigh p-value cannot fall below exp(-n) (attained at PLV = 1), so
  # after BH over n_freq frequencies the smallest reachable q is
  # exp(-n) * n_freq. With few samples that exceeds the threshold no matter how
  # perfect the alignment, and the signature comes back empty for a reason that
  # has nothing to do with the data. Same shape of problem as the permutation
  # floor in condition_test.R -- say so rather than returning an empty table.
  has_null <- "p_null_fwer" %in% colnames(cs) && any(is.finite(cs$p_null_fwer))
  n_min <- min(cs$n_samples_valid, na.rm = TRUE)
  reachable <- exp(-n_min) * nrow(cs)
  # The reachability argument is about the RAYLEIGH q. When the permutation-
  # calibrated p_plv_null exists it is the gate (see below) and this branch
  # would wrongly bypass it, labelling everything exploratory by prevalence.
  rayleigh_is_gate <- !("p_plv_null" %in% colnames(cs) && any(is.finite(cs$p_plv_null)))
  if (rayleigh_is_gate && reachable > plv_q) {
    tsf_warn("With ", n_min, " sample(s) and ", nrow(cs), " frequencies the ",
             "smallest reachable phase-alignment q is ", signif(reachable, 2),
             " > ", plv_q, ": perfect alignment could not pass. About ",
             ceiling(log(nrow(cs) / plv_q)), " samples are needed for this ",
             "condition. Reporting by prevalence and score only.")
    hit <- cs[cs$prevalence >= min_prevalence, , drop = FALSE]
    if (!nrow(hit)) return(NULL)
    hit$phase_alignment_testable <- FALSE
    hit$signature_class <- "exploratory"
    hit <- hit[order(-hit$consensus_score_ci_lower), ]
    return(utils::head(hit, max_components))
  }
  # PHASE GATE: the permutation-calibrated PLV when it exists, Rayleigh only as
  # the fallback. Rayleigh assumes phases iid uniform across samples, which is
  # false on a shared grid and shared tissue -- on a real cohort it admitted
  # 98.9% of frequencies (THEORY.md 5.6). p_plv_null asks the question that
  # matters: more phase-coherent than a random group of the same size from the
  # same tissue.
  #
  # The UNADJUSTED p_plv_null is used, deliberately. This is a filter, like
  # `prevalence >= min_prevalence` next to it, not the inferential claim: the
  # claim is p_null_fwer, which is family-wise by construction and whose score
  # already contains the PLV as a factor. BH-adjusting a pointwise permutation
  # p over n_f frequencies puts its floor at n_f/(B+1), unreachable at any
  # realistic B (5.5b) -- so q_plv_null is reported but cannot gate.
  # `phase_gate` records which statistic decided each row.
  has_plv_null <- "p_plv_null" %in% colnames(cs) && any(is.finite(cs$p_plv_null))
  phase_ok <- if (has_plv_null) {
    !is.na(cs$p_plv_null) & cs$p_plv_null <= plv_q
  } else {
    !is.na(cs$plv_rayleigh_q) & cs$plv_rayleigh_q <= plv_q
  }
  if (!has_plv_null) {
    tsf_warn("No permutation-calibrated PLV (q_plv_null) available; the phase ",
             "gate falls back to Rayleigh, which over-admits on a shared grid.")
  }
  keep <- cs$prevalence >= min_prevalence & phase_ok
  hit <- cs[keep, , drop = FALSE]
  if (!nrow(hit)) return(NULL)
  hit$phase_alignment_testable <- TRUE
  hit$phase_gate <- if (has_plv_null) "plv_null" else "rayleigh"
  # "confirmed" needs the component to beat its OWN permuted null (BH-adjusted
  # across components), not merely to clear zero, which a product of
  # non-negative quantities does automatically. The global-maximum null is
  # kept as the strict flag: it controls error over the whole signature and
  # confirms only components that dominate every frequency of every random
  # draw. Without any null, the strongest available statement is exploratory.
  # Family-wise by construction: no further adjustment, and the floor is
  # reachable with the default number of draws.
  beats_null <- if (has_null)
    !is.na(hit$p_null_fwer) & hit$p_null_fwer <= null_q else rep(FALSE, nrow(hit))
  hit$signature_class <- ifelse(beats_null, "confirmed", "exploratory")
  if (!has_null) {
    tsf_warn("No permutation null was computed, so no component can be ",
             "confirmed: clearing zero is not evidence. Set consensus$n_null.")
  }
  hit <- hit[order(-hit$consensus_score_ci_lower), ]
  utils::head(hit, max_components)
}

#' Classify a condition's components by how common they are WITHIN that
#' condition alone, using nothing but per-sample prevalence.
#'
#' A different question from consensus_signature(): that one asks whether a
#' component is a CONFIRMED signature (prevalent AND phase-locked AND beats a
#' permutation null against the whole dataset). This one asks only "in what
#' fraction of this condition's own samples does the frequency stand out",
#' using `prevalence` -- already computed condition-blind, per sample, per
#' chromosome, against that sample's own spectrum alone
#' (see prevalence_from_rank()) -- with no requirement on phase or power.
#'
#' `prevalence` is deliberately used here rather than `prevalence_rank`: it
#' already takes the best evidence available per sample -- `prevalence_maxt`
#' (the fraction of samples where the peak was significant under the maxT
#' permutation test, when the maxt stage was run) when it exists, and falls
#' back to `prevalence_rank` (a same-sample top-quantile heuristic, no
#' significance test) only when it does not. So the tiers below are backed by
#' real per-sample significance whenever maxT ran, without any change here.
#'
#' Three nested tiers, from `thresholds` (default 0.60/0.80/0.90). Nested by
#' construction: prevalence >= 0.90 implies >= 0.80 implies >= 0.60, so a
#' component is labelled by the HIGHEST threshold it clears.
#'
#' @param cs the consensus_spectrum() table for ONE condition (needs the
#'   `prevalence` column it already produces)
#' @param thresholds increasing prevalence cutoffs (default 0.60, 0.80, 0.90)
#' @return `cs` with one new logical column per threshold
#'   (`condition_invariant_60`, `_80`, `_90`) and `condition_invariant_class`,
#'   the highest tier label reached ("60", "80", "90") or "none".
classify_condition_invariants <- function(cs, thresholds = c(0.60, 0.80, 0.90)) {
  if (is.null(cs) || !nrow(cs)) return(cs)
  if (!"prevalence" %in% colnames(cs)) {
    tsf_abort("classify_condition_invariants needs a prevalence column ",
              "(from consensus_spectrum())")
  }
  thresholds <- sort(as.numeric(thresholds))
  labels <- as.character(round(thresholds * 100))
  tier_cols <- paste0("condition_invariant_", labels)
  for (i in seq_along(thresholds)) {
    cs[[tier_cols[i]]] <-
      is.finite(cs$prevalence) & cs$prevalence >= thresholds[i]
  }
  best <- rep("none", nrow(cs))
  for (i in seq_along(thresholds)) {
    best[cs[[tier_cols[i]]] %in% TRUE] <- labels[i]
  }
  cs$condition_invariant_class <- best
  cs
}

#' Per-sample "stands out" flags for a WHOLE dataset, computed ONCE on the
#' pooled per-sample spectra -- every condition's samples together -- before
#' any split by condition.
#'
#' prevalence_from_rank() (and the maxT flag below) were already
#' condition-blind in their MATH: each sample's own threshold depends only on
#' that sample's own spectrum, never on which other samples happen to be in
#' the table. But calling them once per condition-filtered subset, the way
#' stage_consensus()'s per-condition loop calls consensus_spectrum(), still
#' means the code reads a condition label before anything is decided. This
#' function makes that impossible: it takes the pool assembled BEFORE the
#' per-condition loop begins, and condition is not a parameter here -- the
#' pool does not carry one that this function looks at. A caller groups the
#' result by condition afterwards, in condition_invariants_from_pool().
#'
#' @param pool per-sample spectra across every condition of a dataset (chr,
#'   N, k, sample, power[, power_normalised], phase)
#' @param maxt optional POOLED per-sample maxT table (every condition's,
#'   already row-bound -- the same table stage_consensus() builds for the
#'   null, reused here)
#' @param quantile_cut passed to prevalence_from_rank()
#' @param alpha the maxT significance level (only used when `maxt` is given)
#' @return one row per (chr, N, k, sample): `stands_out` (maxT-based when
#'   `maxt` is supplied and covers that row, rank-based otherwise), plus
#'   `stands_out_rank` and `stands_out_maxt` kept separately for audit, and
#'   `pnorm`/`phase` carried through so a caller can compute median power and
#'   phase-locking WITHIN a group (e.g. a condition) without re-reading the
#'   spectra -- see condition_invariants_from_pool().
pooled_stands_out <- function(pool, maxt = NULL, quantile_cut = 0.95, alpha = 0.05) {
  needed <- c("chr", "N", "k", "sample", "power", "phase")
  missing <- setdiff(needed, colnames(pool))
  if (length(missing)) tsf_abort("pooled_stands_out needs columns: ",
                                 paste(missing, collapse = ", "))
  d <- pool
  if ("power_normalised" %in% colnames(d)) {
    d$pnorm <- d$power_normalised
  } else {
    tot <- stats::ave(d$power, paste(d$sample, d$chr), FUN = function(x) sum(x, na.rm = TRUE))
    d$pnorm <- d$power / pmax(tot, .Machine$double.eps)
  }
  d$stands_out_rank <- prevalence_from_rank(d$pnorm, d$sample, d$chr, quantile_cut)
  has_maxt <- !is.null(maxt) && "p_empirical_maxT" %in% colnames(maxt)
  d$stands_out_maxt <- if (has_maxt) {
    key_d <- paste(d$chr, d$N, d$k, d$sample)
    key_m <- paste(maxt$chr, maxt$N, maxt$k, maxt$sample)
    sig <- maxt$p_empirical_maxT <= alpha
    out <- sig[match(key_d, key_m)]
    out[is.na(out)] <- FALSE
    out
  } else rep(NA, nrow(d))
  d$stands_out <- if (has_maxt) d$stands_out_maxt else d$stands_out_rank
  d[, c("chr", "N", "k", "sample", "pnorm", "phase",
        "stands_out", "stands_out_rank", "stands_out_maxt")]
}

#' Condition-own invariant components, from per-sample "stands out" flags
#' that were computed on the whole pool before any condition was consulted
#' (see pooled_stands_out()). Condition enters ONLY here, as a filter on
#' flags that already exist -- never inside pooled_stands_out().
#'
#' A component is not just "present or absent": `prevalence` alone would call
#' two components the same just because they share (chr, N, k), even if their
#' phase and power look nothing alike across this condition's samples. So
#' alongside prevalence, this reports -- WITHIN the condition's own samples,
#' no other data, no permutation, no bootstrap -- the same two axes
#' consensus_spectrum() reports for the whole-dataset question:
#'   median_power_normalised  median of pnorm across the condition's samples
#'   plv                      phase_locking()'s |mean(exp(i*phase))|, 1 when
#'                            every sample puts the crest in the same place
#' These two are computed here, deterministically, from this call alone.
#' Confidence intervals (bootstrap) and a null p-value (permutation against
#' the whole pool) are NOT computed here -- see condition_invariant_bootstrap()
#' and stage_consensus()'s use of null_consensus_distribution(), which attach
#' them afterwards as DESCRIPTIVE columns. None of it gates the tier: no
#' threshold on power, phase, the CI or the null p-value is applied here,
#' because none was specified for this tier system -- inventing one would be
#' a filter nobody asked for. `prevalence` alone decides the 60/80/90% tier;
#' everything else rides along as evidence to judge each component by.
#'
#' @param pooled the table pooled_stands_out() returns
#' @param condition_samples sample ids belonging to ONE condition
#' @param thresholds forwarded to classify_condition_invariants()
#' @return one row per (chr, N, k) reached by this condition's samples, with
#'   `prevalence`, `n_samples_condition`, `median_power_normalised`, `plv`,
#'   `consensus_score_rank` (their product, so null_component_pvalues() can
#'   be reused unchanged), and the tier columns from
#'   classify_condition_invariants() -- restricted to rows that clear at
#'   least the lowest threshold (NULL if none do).
condition_invariants_from_pool <- function(pooled, condition_samples,
                                           thresholds = c(0.60, 0.80, 0.90)) {
  sub <- pooled[pooled$sample %in% condition_samples & is.finite(pooled$stands_out), ,
               drop = FALSE]
  if (!nrow(sub)) return(NULL)
  key <- paste(sub$chr, sub$N, sub$k, sep = "|")
  groups <- split(seq_len(nrow(sub)), key)
  rows <- lapply(names(groups), function(g) {
    i <- groups[[g]]
    parts <- strsplit(g, "|", fixed = TRUE)[[1]]
    pl <- phase_locking(sub$phase[i])
    med_p <- stats::median(sub$pnorm[i], na.rm = TRUE)
    prev <- mean(sub$stands_out[i])
    data.frame(chr = parts[1], N = as.integer(parts[2]), k = as.integer(parts[3]),
               period = as.integer(parts[2]) / as.integer(parts[3]),
               n_samples_condition = length(unique(sub$sample[i])),
               prevalence = prev,
               median_power_normalised = med_p,
               plv = unname(pl[["plv"]]),
               plv_rayleigh_p = unname(pl[["rayleigh_p"]]),
               # Same product as consensus_spectrum()'s consensus_score_rank
               # (median power x prevalence x PLV), kept under that name so
               # null_component_pvalues() -- built for that column -- can be
               # reused unchanged to attach a null p-value to these rows too.
               consensus_score_rank = med_p * prev * unname(pl[["plv"]]),
               stringsAsFactors = FALSE)
  })
  cs <- do.call(rbind, rows)
  cs <- classify_condition_invariants(cs, thresholds)
  cs <- cs[cs$condition_invariant_class != "none", , drop = FALSE]
  if (!nrow(cs)) return(NULL)
  cs[order(-cs$prevalence), ]
}

#' Bootstrap confidence intervals for a condition's own invariants, by
#' resampling THAT CONDITION'S OWN samples (with replacement). Descriptive
#' only -- it does not gate condition_invariant_class, which is already
#' decided by the time this runs. Answers a different question from the null
#' p-value: not "is this more than chance", but "if a slightly different set
#' of patients had made up this condition, would prevalence/power/PLV still
#' look like this" -- i.e. how much of the estimate is this specific set of
#' patients versus the condition in general.
#'
#' Only bootstraps the (chr, N, k) rows already in `cond_inv` -- the
#' components that already cleared a prevalence tier -- not every frequency
#' the condition was ever tested at, which is what made the legacy route
#' expensive.
#'
#' @param pooled the table pooled_stands_out() returns
#' @param condition_samples sample ids belonging to ONE condition
#' @param cond_inv the table condition_invariants_from_pool() returned for
#'   this same condition
#' @param n_boot resamples (default 200, the same default consensus_spectrum()
#'   uses)
#' @param seed for reproducibility
#' @return `cond_inv` with `prevalence_ci_lower/upper`,
#'   `median_power_ci_lower/upper`, `plv_ci_lower/upper` and
#'   `consensus_score_ci_lower/upper` added. The last is needed by
#'   null_component_pvalues()'s `beats_global_null` (consensus_spectrum()'s
#'   own signature-selection column), so call this BEFORE
#'   null_component_pvalues(), not after.
condition_invariant_bootstrap <- function(pooled, condition_samples, cond_inv,
                                          n_boot = 200L, seed = 42L) {
  if (is.null(cond_inv) || !nrow(cond_inv)) return(cond_inv)
  sub <- pooled[pooled$sample %in% condition_samples & is.finite(pooled$stands_out), ,
               drop = FALSE]
  key <- paste(sub$chr, sub$N, sub$k, sep = "|")
  want <- paste(cond_inv$chr, cond_inv$N, cond_inv$k, sep = "|")

  set.seed(seed)
  bounds <- lapply(want, function(k0) {
    i <- which(key == k0)
    if (!length(i)) {
      return(c(prev_lo = NA_real_, prev_hi = NA_real_, power_lo = NA_real_,
               power_hi = NA_real_, plv_lo = NA_real_, plv_hi = NA_real_,
               score_lo = NA_real_, score_hi = NA_real_))
    }
    samples_here <- unique(sub$sample[i])
    draws <- vapply(seq_len(n_boot), function(b) {
      pick <- sample(samples_here, replace = TRUE)
      idx <- unlist(lapply(pick, function(s) i[sub$sample[i] == s]), use.names = FALSE)
      if (!length(idx)) return(c(NA_real_, NA_real_, NA_real_, NA_real_))
      prev_b <- mean(sub$stands_out[idx], na.rm = TRUE)
      pow_b <- stats::median(sub$pnorm[idx], na.rm = TRUE)
      plv_b <- Mod(mean(exp(1i * sub$phase[idx])))
      c(prev_b, pow_b, plv_b, prev_b * pow_b * plv_b)
    }, numeric(4))
    qs <- function(v) stats::quantile(v, c(0.025, 0.975), na.rm = TRUE)
    q1 <- qs(draws[1, ]); q2 <- qs(draws[2, ]); q3 <- qs(draws[3, ]); q4 <- qs(draws[4, ])
    c(prev_lo = unname(q1[1]), prev_hi = unname(q1[2]),
      power_lo = unname(q2[1]), power_hi = unname(q2[2]),
      plv_lo = unname(q3[1]), plv_hi = unname(q3[2]),
      score_lo = unname(q4[1]), score_hi = unname(q4[2]))
  })
  bmat <- do.call(rbind, bounds)
  cond_inv$prevalence_ci_lower <- bmat[, "prev_lo"]
  cond_inv$prevalence_ci_upper <- bmat[, "prev_hi"]
  cond_inv$median_power_ci_lower <- bmat[, "power_lo"]
  cond_inv$median_power_ci_upper <- bmat[, "power_hi"]
  cond_inv$plv_ci_lower <- bmat[, "plv_lo"]
  cond_inv$plv_ci_upper <- bmat[, "plv_hi"]
  cond_inv$consensus_score_ci_lower <- bmat[, "score_lo"]
  cond_inv$consensus_score_ci_upper <- bmat[, "score_hi"]
  cond_inv
}
