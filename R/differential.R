# differential.R -- differential expression, on frequencies instead of genes.
#
# WHY THIS IS THE RIGHT QUESTION
# ------------------------------
# Everything else in this pipeline asks DETECTION: is this frequency's power
# larger than a null in which gene positions were shuffled? That null destroys
# the signal, so a frequency has to STAND OUT within a condition to pass.
#
# The question that matters here is COMPARISON, and it has a different null.
# The spectrum is a linear transform of the ordered expression vector, so if two
# conditions differ in expression -- and they demonstrably do, that is what
# differential expression means -- their spectra differ too. The difference may
# be small. It is not zero. A test that can only see components strong enough to
# beat a no-signal null will miss all of it, and report the absence of something
# that is present by construction.
#
# So: for each frequency, compare power ACROSS conditions. Small, consistent
# shifts are exactly what a two-sample test on tens of samples detects well, and
# exactly what a detection test throws away.
#
# `condition_contrast()` was a step in this direction but binarises: it works on
# "stands out", so a frequency whose power rises consistently in F4 without ever
# entering the top of any sample's spectrum contributes nothing. This works on
# the magnitude.
#
# WHICH POWER
# -----------
# `power_normalised`, normalised within (sample, chromosome), so library depth
# and per-chromosome scale are already out. Comparing raw power across cohorts
# would mostly compare sequencing depth.
#
# THREE TESTS, DELIBERATELY
# -------------------------
#   one-vs-rest   per condition: is this frequency's power different HERE than
#                 in every other condition pooled? This is the characteristic
#                 spectrum of a condition, frequency by frequency.
#   omnibus       Kruskal-Wallis across all conditions: does this frequency
#                 differ ANYWHERE? More powerful when the difference is spread
#                 over several stages rather than concentrated in one.
#   trend         Spearman of power against ordinal stage. F0 < F1 < ... < F4 is
#                 ordered, and a frequency that rises monotonically across the
#                 scale is the most interpretable result the design can produce
#                 -- and is invisible to both tests above, which treat the
#                 stages as unordered labels.
#
# Rank-based throughout: spectral power is heavy-tailed and a handful of samples
# per condition is not enough to trust a normal approximation.

#' Differential spectral power between conditions.
#'
#' @param sp per-sample spectra: chr, N, k, sample, power_normalised
#' @param groups named vector, condition per sample id
#' @param stage_order optional ordered condition levels for the trend test
#' @param min_per_condition conditions below this are dropped, not tested
#' @return one row per (frequency, test) with effect size, p and BH q
differential_spectrum <- function(sp, groups, stage_order = NULL,
                                  min_per_condition = 4L,
                                  value_col = "power_normalised") {
  if (!value_col %in% names(sp)) {
    tsf_abort("differential_spectrum: no column '", value_col, "' in the ",
              "per-sample spectra. Re-run the spectra stage.")
  }
  sp$condition <- unname(groups[as.character(sp$sample)])
  sp <- sp[!is.na(sp$condition), , drop = FALSE]

  n_by <- tapply(sp$sample, sp$condition, function(x) length(unique(x)))
  small <- names(n_by)[n_by < min_per_condition]
  if (length(small)) {
    tsf_warn("differential_spectrum: dropping ", paste(small, collapse = ", "),
             " (fewer than ", min_per_condition, " samples). A difference ",
             "estimated from two or three samples is not a difference.")
    sp <- sp[!(sp$condition %in% small), , drop = FALSE]
  }
  conds <- sort(unique(sp$condition))
  if (length(conds) < 2L) {
    tsf_warn("differential_spectrum: fewer than two usable conditions")
    return(NULL)
  }

  sp$key <- paste(sp$chr, sp$N, sp$k, sep = "|")
  by_key <- split(seq_len(nrow(sp)), sp$key)
  tsf_log("  differential spectrum: ", length(by_key), " frequencies over ",
          length(conds), " condition(s)")

  rank_biserial <- function(a, b) {
    # Effect size that matches a rank test: the probability that a random value
    # from `a` exceeds one from `b`, rescaled to [-1, 1]. Reported because a
    # p-value alone cannot say whether a difference is worth anything, and with
    # 50 samples a trivial shift reaches significance.
    if (!length(a) || !length(b)) return(NA_real_)
    r <- rank(c(a, b))
    u <- sum(r[seq_along(a)]) - length(a) * (length(a) + 1) / 2
    2 * u / (length(a) * length(b)) - 1
  }

  rows <- lapply(names(by_key), function(kk) {
    d <- sp[by_key[[kk]], , drop = FALSE]
    v <- suppressWarnings(as.numeric(d[[value_col]]))
    g <- d$condition
    ok <- is.finite(v)
    v <- v[ok]; g <- g[ok]
    if (length(unique(g)) < 2L) return(NULL)

    parts <- strsplit(kk, "|", fixed = TRUE)[[1]]
    base <- data.frame(chr = parts[1], N = as.integer(parts[2]),
                       k = as.integer(parts[3]), key = kk,
                       stringsAsFactors = FALSE)

    out <- list()

    # one-vs-rest, per condition
    for (cnd in intersect(conds, unique(g))) {
      a <- v[g == cnd]; b <- v[g != cnd]
      if (length(a) < 2L || length(b) < 2L) next
      p <- tryCatch(stats::wilcox.test(a, b, exact = FALSE)$p.value,
                    error = function(e) NA_real_)
      out[[length(out) + 1]] <- cbind(base, data.frame(
        test = "one_vs_rest", condition = cnd,
        effect = rank_biserial(a, b),
        median_here = stats::median(a), median_rest = stats::median(b),
        p = p, stringsAsFactors = FALSE))
    }

    # omnibus
    p_om <- tryCatch(stats::kruskal.test(v, factor(g))$p.value,
                     error = function(e) NA_real_)
    out[[length(out) + 1]] <- cbind(base, data.frame(
      test = "omnibus", condition = NA_character_,
      effect = NA_real_, median_here = NA_real_, median_rest = NA_real_,
      p = p_om, stringsAsFactors = FALSE))

    # monotone trend across the ordered scale
    if (!is.null(stage_order)) {
      idx <- match(g, stage_order)
      if (sum(!is.na(idx)) >= 4L && length(unique(idx[!is.na(idx)])) >= 3L) {
        ct <- tryCatch(stats::cor.test(v[!is.na(idx)], idx[!is.na(idx)],
                                       method = "spearman", exact = FALSE),
                       error = function(e) NULL)
        if (!is.null(ct)) {
          out[[length(out) + 1]] <- cbind(base, data.frame(
            test = "trend", condition = NA_character_,
            effect = unname(ct$estimate),
            median_here = NA_real_, median_rest = NA_real_,
            p = ct$p.value, stringsAsFactors = FALSE))
        }
      }
    }
    do.call(rbind, out)
  })

  res <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(res) || !nrow(res)) return(NULL)

  # BH within each test family, and within each condition for one-vs-rest.
  # Pooling them would let a condition with many candidates raise the bar for
  # one with few, and would mix three questions into one family.
  res$q <- NA_real_
  fam <- ifelse(res$test == "one_vs_rest",
                paste(res$test, res$condition), res$test)
  for (f in unique(fam)) {
    i <- fam == f & is.finite(res$p)
    if (any(i)) res$q[i] <- stats::p.adjust(res$p[i], "BH")
  }

  res[order(res$test, res$condition, res$p), , drop = FALSE]
}
