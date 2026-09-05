# contrast.R -- which frequencies are characteristic OF a condition.
#
# THE GAP THIS FILLS
# ------------------
# `consensus_score = power x prevalence x PLV` is computed entirely WITHIN one
# condition. None of the three factors asks whether the component looks any
# different in the other conditions, so the score answers "how strong and how
# common is this frequency here" and never "is this frequency here rather than
# there".
#
# The difference is not academic. On real liver data chr4 k36 had prevalence
#
#   Normal_histology 0.97 | F0 1.00 | F1 1.00 | F2 1.00 | F3 1.00 | F4 1.00
#
# It stands out in essentially every sample of every stage. It is therefore the
# top-scoring component of most conditions and characteristic of none: a
# frequency that does not separate a healthy liver from a cirrhotic one is a
# property of the tissue or of the axis, not a signature.
#
# WHAT IS COMPUTED
# ----------------
# One-vs-rest, because the question is "does this spectrum say F4" rather than
# "does F4 differ from health":
#
#   contrast(f, c) = prevalence(f, c) - max_{c' != c} prevalence(f, c')
#
# The maximum rather than the mean of the others: a frequency equally prevalent
# in F3 cannot identify F4, however rare it is in the rest, and averaging would
# hide exactly that.
#
# Positive and large means characteristic of c. Near zero means invariant --
# still worth reporting, as the constant signature, but never as a marker.
# Negative means characteristic of some other condition.
#
# THE NULL
# --------
# Condition labels are permuted among the samples of the dataset and the
# contrast recomputed. That holds the spectra, the grid, the tissue and the
# per-sample significance calls fixed, and destroys only the association between
# a sample and its stage -- which is precisely the thing being claimed.
#
# It is cheap: prevalence under a permutation is a `rowMeans` over a different
# set of columns of a standing matrix that already exists. No spectrum and no
# maxT is recomputed.

#' Prevalence of every frequency in every condition, from a standing matrix.
#'
#' @param stands frequencies x samples, 1 where the frequency stands out
#' @param groups condition label per column of `stands`
prevalence_by_condition <- function(stands, groups) {
  conds <- sort(unique(groups))
  out <- vapply(conds, function(g) {
    cols <- which(groups == g)
    if (!length(cols)) return(rep(NA_real_, nrow(stands)))
    rowMeans(stands[, cols, drop = FALSE], na.rm = TRUE)
  }, numeric(nrow(stands)))
  if (is.null(dim(out))) out <- matrix(out, nrow = nrow(stands))
  dimnames(out) <- list(rownames(stands), conds)
  out
}

#' One-vs-rest contrast from a prevalence matrix.
contrast_from_prevalence <- function(prev) {
  conds <- colnames(prev)
  out <- prev
  for (j in seq_along(conds)) {
    others <- prev[, -j, drop = FALSE]
    rest <- if (!ncol(others)) rep(0, nrow(prev)) else
      suppressWarnings(apply(others, 1L, max, na.rm = TRUE))
    rest[!is.finite(rest)] <- 0
    out[, j] <- prev[, j] - rest
  }
  out
}

#' Contrast, with a p-value from permuted condition labels.
#'
#' @param stands frequencies x samples standing matrix (maxT-based when maxT
#'   ran, rank-based otherwise -- the same matrix the consensus score used, so
#'   the contrast is tested on the statistic it was computed from)
#' @param groups condition label per column
#' @param B permutations; the floor on any p is 1/(B+1)
#' @return long data frame: key, condition, prevalence, contrast, p, q
condition_contrast <- function(stands, groups, B = 999L, seed = 42L,
                               min_per_condition = 3L) {
  groups <- as.character(groups)
  if (ncol(stands) != length(groups)) {
    tsf_abort("condition_contrast: ", ncol(stands), " sample column(s) but ",
              length(groups), " label(s)")
  }
  n_by <- table(groups)
  too_small <- names(n_by)[n_by < min_per_condition]
  if (length(too_small)) {
    tsf_warn("condition_contrast: dropping condition(s) with fewer than ",
             min_per_condition, " samples: ", paste(too_small, collapse = ", "),
             ". A contrast estimated from one or two samples is not a signature.")
    keep <- !(groups %in% too_small)
    stands <- stands[, keep, drop = FALSE]; groups <- groups[keep]
  }
  if (length(unique(groups)) < 2L) {
    tsf_warn("condition_contrast: fewer than two conditions; nothing to contrast")
    return(NULL)
  }

  obs <- contrast_from_prevalence(prevalence_by_condition(stands, groups))
  prev_obs <- prevalence_by_condition(stands, groups)

  # One-sided: the claim is that the contrast is LARGER than chance, so a
  # frequency prevalent everywhere gets a p near 1 rather than a small one for
  # being unusually flat.
  set.seed(seed)
  ge <- matrix(0L, nrow(obs), ncol(obs), dimnames = dimnames(obs))
  for (b in seq_len(B)) {
    perm <- contrast_from_prevalence(
      prevalence_by_condition(stands, sample(groups)))
    ge <- ge + (perm >= obs)
  }
  p <- (1 + ge) / (B + 1)

  long <- data.frame(
    key        = rep(rownames(obs), times = ncol(obs)),
    condition  = rep(colnames(obs), each = nrow(obs)),
    prevalence = as.numeric(prev_obs),
    contrast   = as.numeric(obs),
    p_contrast = as.numeric(p),
    stringsAsFactors = FALSE)

  # BH within each condition: each is its own family of hypotheses, and pooling
  # them would let a condition with many candidates raise the bar for one with
  # few.
  long$q_contrast <- NA_real_
  for (g in unique(long$condition)) {
    i <- long$condition == g
    long$q_contrast[i] <- stats::p.adjust(long$p_contrast[i], "BH")
  }

  parts <- do.call(rbind, strsplit(long$key, "|", fixed = TRUE))
  long$chr <- parts[, 1]
  long$N <- as.integer(parts[, 2])
  long$k <- as.integer(parts[, 3])

  # Invariant means prevalent in EVERY condition, so it is a property of the
  # minimum across conditions, not of the one in hand. Testing prevalence in
  # this condition alone called a frequency shared by F3 and F4 only -- absent
  # from Control and F0 -- "invariant", because within that pair the contrast
  # is zero and prevalence is high. The minimum catches it.
  min_prev <- apply(prev_obs, 1L, min, na.rm = TRUE)
  long$min_prevalence <- min_prev[match(long$key, rownames(prev_obs))]

  long$class <- ifelse(
    long$q_contrast <= 0.05 & long$contrast > 0, "characteristic",
    ifelse(abs(long$contrast) < 0.05 & long$min_prevalence >= 0.8, "invariant",
           "unremarkable"))

  # Same reachability constraint as every other BH route here, and worth saying
  # out loud: at rank 1 with no ties the smallest q is m/(B+1), so with 200
  # frequencies and 499 draws nothing can clear 0.05 however large the contrast.
  # Reported per condition, since BH is applied per condition.
  m_family <- nrow(obs)
  rank1 <- bh_rank1_diagnostic(m_family, B)
  if (rank1 > 0.05) {
    tsf_log("Contrast: ", m_family, " frequencies at ", B, " draws gives a ",
            "rank-1 BH diagnostic of ", signif(rank1, 3),
            ". Raise --n-contrast to at least ", draws_for_bh(m_family),
            " or read p_contrast with the tie count, not q_contrast alone.")
  }

  long[order(long$condition, -long$contrast), , drop = FALSE]
}
