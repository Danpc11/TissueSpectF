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
prepare_null_consensus_matrices <- function(spectra_all, quantile_cut = 0.95) {
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
  list(pnorm = pn, stands = so, phase = ph, keys = keys, samples = samples)
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
  score[n_valid < 2L | !is.finite(score)] <- NA_real_
  names(score) <- prepared$keys
  score
}

null_consensus_distribution <- function(spectra_all, n_samples, n_null = 50L,
                                        seed = 42L, quantile_cut = 0.95,
                                        q_global = 0.95, blocks = NULL,
                                        n_cores = 1L,
                                        engine = c("matrix", "reference"),
                                        prepared = NULL) {
  engine <- match.arg(engine)
  samples <- unique(spectra_all$sample)
  if (length(samples) <= n_samples || n_null < 10L) return(NULL)
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
    values <- if (engine == "matrix") {
      tryCatch(null_matrix_draw(prepared, pick), error = function(e) NULL)
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
    if (!length(values)) return(NULL)
    list(values = values, best = max(values))
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

  list(per_key = mat, global = unname(stats::quantile(best, q_global)),
       global_max_draws = best, n_null = ncol(mat), n_samples = n_samples)
}

#' Empirical p-value per component against its own null, then BH.
#'
#' The global maximum null controls the error rate over the whole signature and
#' is very strict -- in practice nothing is confirmed under it, which is correct
#' for "is there any component at all" and useless for "which components".
#' Keeping the per-(chr, k) null as well gives a p-value per component; BH over
#' those is the intermediate that lets a signature be localised. Both are
#' reported, and which one a claim rests on is a stated choice.
null_component_pvalues <- function(cs, null_dist, null_q = 0.05) {
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
  if (reachable > plv_q) {
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
  keep <- cs$prevalence >= min_prevalence &
    !is.na(cs$plv_rayleigh_q) & cs$plv_rayleigh_q <= plv_q
  hit <- cs[keep, , drop = FALSE]
  if (!nrow(hit)) return(NULL)
  hit$phase_alignment_testable <- TRUE
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
  hit$signature_class <- ifelse(beats_null & hit$plv_rayleigh_q <= plv_q,
                                "confirmed", "exploratory")
  if (!has_null) {
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
