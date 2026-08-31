# condition_test.R -- significance at the level of a condition.
#
# WHY THIS EXISTS
# ---------------
# The stability criterion (>= stable_frac of samples individually significant)
# is a consistency requirement, not a test with aggregated power. Each sample
# has to reach significance on its own under family-wise control across
# frequencies; evidence is never pooled. A periodicity that is real but moderate
# -- present in every sample, individually short of significance in each -- is
# invisible to it no matter how many samples the condition has. Adding samples
# does not increase power, which is backwards.
#
# Two condition-level tests are computed here. Neither replaces the consistency
# figure: that stays, reported as a descriptor of reproducibility rather than as
# the criterion that decides what exists.
#
#   1. Permutation null on the summary signal. The condition's mean (or median)
#      profile is the signal actually analysed downstream, so it is the thing to
#      test. Values are permuted among the observed grid positions, the spectrum
#      is recomputed, and the maximum power over k gives the null -- the same
#      maxT logic, applied once per condition instead of once per sample. This
#      is the primary test.
#
#   2. Stouffer combination of the per-sample maxT p-values. Aggregates evidence
#      across samples and gains power with n. It is conservative here, because
#      the inputs are already family-wise adjusted within a chromosome, and it
#      assumes independence between samples, which biological replicates only
#      approximately satisfy. Reported as a second opinion, not as the primary.
#
# Cost: one permutation set per (condition, chromosome) instead of one per
# (sample, chromosome). For a 54-sample condition that is ~54x cheaper, which is
# what makes a full run feasible on two cores.

#' Stouffer's combination of independent p-values.
#'
#' Permutation p-values are bounded below by 1/(B+1); the corresponding z is
#' capped accordingly, so a combined value cannot be driven by one sample
#' claiming impossible precision.
stouffer_combine <- function(p, B = NULL) {
  p <- p[is.finite(p)]
  if (!length(p)) return(NA_real_)
  lo <- if (is.null(B)) .Machine$double.eps else 1 / (B + 1)
  p <- pmin(pmax(p, lo), 1 - lo)
  z <- stats::qnorm(1 - p)
  stats::pnorm(sum(z) / sqrt(length(z)), lower.tail = FALSE)
}

#' Permutation test on one condition's summary signal, one chromosome.
condition_permutation_test <- function(signal, terms, B = 1000L, seed = 42L,
                                       block_sizes = c(10L, 20L, 50L),
                                       primary_scheme = "full") {
  permutation_gls_test(signal, terms, B = B, seed = seed,
                       block_sizes = block_sizes,
                       primary_scheme = primary_scheme)
}

#' Condition-level significance for every chromosome of one condition.
#'
#' @param maxt_individual optional per-sample maxT table; when supplied the
#'   Stouffer combination is added alongside the permutation result.
condition_significance <- function(dataset, cond, chrom_idx, maxt_cfg,
                                   branch = "average", maxt_individual = NULL,
                                   n_cores = 1L) {
  sig <- condition_signals(dataset, cond)
  if (is.null(sig)) return(NULL)
  summary_signal <- if (branch == "average") sig$avg_signal else sig$median_signal
  chrom_levels <- names(chrom_idx)
  B <- maxt_cfg$condition_B %||% maxt_cfg$B

  per_chr <- parallel::mclapply(chrom_levels, function(chr_now) {
    ci <- chrom_idx[[chr_now]]
    terms <- gls_prepare(ci$t, ci$N)
    res <- condition_permutation_test(
      summary_signal[ci$rows], terms, B = B,
      seed = maxt_cfg$seed + 100000L * match(chr_now, chrom_levels),
      block_sizes = maxt_cfg$block_sizes,
      primary_scheme = maxt_cfg$primary_scheme %||% "full")
    if (is.null(res) || !nrow(res)) return(NULL)

    out <- data.frame(
      chr = chr_now, N = res$N, k = res$k, freq = res$freq, period = res$period,
      n_observed = res$n_observed, coverage = ci$coverage,
      amplitude = res$amplitude, phase = res$phase, power = res$power,
      power_normalised = res$power_normalised, window_power = res$window_power,
      p_condition = res$p_empirical_maxT,
      p_condition_pointwise = res$p_pointwise %||% NA_real_,
      p_condition_full = res$p_empirical_maxT_full,
      p_condition_all = res$p_empirical_maxT_all,
      stringsAsFactors = FALSE)
    out$window_rank <- rank(-out$window_power)
    out
  }, mc.cores = n_cores, mc.set.seed = FALSE)

  per_chr <- per_chr[!vapply(per_chr, is.null, logical(1))]
  if (!length(per_chr)) return(NULL)
  out <- do.call(rbind, per_chr)

  # CORRECTION FAMILY: chromosomes, not frequencies.
  #
  # maxT already controls the family-wise error rate ACROSS FREQUENCIES within a
  # chromosome -- that is what comparing against the null maximum buys. Applying
  # BH over all ~9,500 frequencies on top of that corrects twice, and fatally:
  # a permutation p-value cannot go below 1/(B+1), so with B = 2000 the smallest
  # attainable q would be 5e-4 * 9500 = 2.4, capped at 1. No peak could ever
  # reach any threshold, whatever the data. The remaining multiplicity is the
  # number of chromosomes tested, so that is what is corrected here.
  #
  # This also sets a floor on B: with n_chr chromosomes, q can only reach 0.05
  # if 1/(B+1) * n_chr <= 0.05, i.e. B >= 20 * n_chr - 1 (about 460 for 23
  # chromosomes). The stage warns when B is too small to be able to conclude.
  n_chr <- length(unique(out$chr))
  out$n_chromosomes_tested <- n_chr

  # TWO MULTIPLICITY REGIMES, BOTH REPORTED, ONE DECLARED AS THE CRITERION.
  #
  # fwer: the maxT p-value already controls the family-wise error across the
  #   frequencies of a chromosome, so only the chromosome multiplicity is left.
  #   It asks whether a frequency beats the strongest frequency of a permuted
  #   spectrum, which on 500-800 frequencies is close to asking whether it is
  #   the dominant component of its chromosome. Almost nothing passes, and what
  #   does is worth believing.
  #
  # fdr: the pointwise p-value asks whether a frequency beats its own null, with
  #   BH across every frequency of every chromosome. It answers the question a
  #   signature is actually about -- which frequencies carry structure -- and
  #   controls the expected proportion of false ones rather than the probability
  #   of any. It selects far more, and a fraction of them are false by design.
  #
  # Neither is more correct in the abstract. Which one a claim rests on is a
  # stated choice, and both columns are always written so the effect of that
  # choice is visible without a rerun.
  out$q_condition <- pmin(1, out$p_condition * n_chr)
  out$q_condition_fdr <- if (all(is.na(out$p_condition_pointwise))) NA_real_ else
    stats::p.adjust(out$p_condition_pointwise, method = "BH")

  n_fwer <- sum(out$q_condition <= 0.05, na.rm = TRUE)
  n_fdr <- sum(out$q_condition_fdr <= 0.05, na.rm = TRUE)
  tsf_log("  ", cond, "/", branch, ": ", n_fwer, " frequencies at FWER q<=0.05, ",
          n_fdr, " at FDR q<=0.05 (of ", nrow(out), ")")
  min_attainable <- n_chr / (B + 1)
  if (min_attainable > 0.05) {
    tsf_warn("B = ", B, " over ", n_chr, " chromosomes: the smallest attainable ",
             "q is ", signif(min_attainable, 2), ", so nothing can reach 0.05. ",
             "Raise maxt$condition_B to at least ", 20 * n_chr, ".")
  }

  if (!is.null(maxt_individual) && nrow(maxt_individual)) {
    key <- paste(maxt_individual$chr, maxt_individual$N, maxt_individual$k, sep = "|")
    agg <- stats::aggregate(maxt_individual$p_empirical_maxT, list(key = key),
                            function(p) stouffer_combine(p, B = maxt_cfg$B))
    out$p_stouffer <- agg$x[match(paste(out$chr, out$N, out$k, sep = "|"), agg$key)]
    # Same family as above: the inputs are already FWER-adjusted within a
    # chromosome, so only the chromosome multiplicity is left to correct.
    out$q_stouffer <- pmin(1, out$p_stouffer * n_chr)
    n_samp <- stats::aggregate(maxt_individual$sample, list(key = key),
                               function(x) length(unique(x)))
    out$n_samples <- n_samp$x[match(paste(out$chr, out$N, out$k, sep = "|"), n_samp$key)]
  }

  out$condition <- cond
  out$branch <- branch
  out[order(out$q_condition, out$p_condition), ]
}

#' Peaks a condition-level test calls significant.
select_condition_peaks <- function(cs, q_threshold = 0.05, use = c("condition", "stouffer")) {
  use <- match.arg(use)
  qcol <- if (use == "condition") "q_condition" else "q_stouffer"
  if (is.null(cs) || !qcol %in% colnames(cs)) return(NULL)
  hit <- cs[!is.na(cs[[qcol]]) & cs[[qcol]] <= q_threshold, , drop = FALSE]
  if (!nrow(hit)) return(NULL)
  hit[order(hit[[qcol]], -hit$power), ]
}
