# clean.R -- greedy spectral deflation (CLEAN / orthogonal matching pursuit)
# with a BIC stopping rule.
#
# WHY
# ---
# On an incomplete grid the sinusoids are not orthogonal over the observed
# positions, so a strong component leaks power into neighbouring frequencies and
# into the aliases the sampling window imposes. Those sidelobes are local maxima
# of the periodogram without being components of the signal. Taking the M
# highest peaks of the raw spectrum therefore returns the real component
# together with its own artefacts, and counts them as independent structure.
#
# CLEAN removes the leakage instead of thresholding around it: take the strongest
# frequency, refit every selected component jointly by least squares, subtract
# the fit, and search the residual. Each component is chosen against what is left
# rather than against what the previous one contaminated.
#
# HOW MANY COMPONENTS
# -------------------
# Not a chosen number. Extraction continues while a component lowers the BIC of
# the joint fit and stops when it does not, so the data decide chromosome by
# chromosome: a chromosome with rich structure yields several components, a flat
# one yields none. Each component costs 2 parameters (a and b, equivalently
# amplitude and phase), plus 1 for the floating mean.
#
# The plain BIC is wrong here, and visibly so: it assumes the model was fixed in
# advance, while each component is chosen as the maximum over ~N/2 candidate
# frequencies. That search is unpaid for, and plain BIC duly selects components
# from pure noise. The extended BIC (Chen & Chen 2008) adds the cost of the
# selection itself:
#
#   EBIC = n*log(RSS/n) + p*log(n) + 2*gamma*log(choose(K, m))
#
#   p = 2*m + 1 parameters (a and b per component, plus the floating mean)
#   K = number of candidate frequencies searched
#   m = number of components selected
#
# gamma = 1 is the strict end of the range and is the default: with a candidate
# set as large as the frequency grid, anything weaker selects noise. BIC rather
# than AIC because the goal is to identify which components are really there,
# not to minimise prediction error. `penalty_factor` scales the whole penalty
# for sensitivity analysis; `ebic_gamma` scales only the selection term.
#
# WHAT THIS IS NOT
# ----------------
# BIC selection is not a significance test and its output carries no error rate.
# A component with a large BIC drop is a component the data prefer to keep; that
# is a different claim from "this periodicity is unlikely under a null". Use
# condition_test.R for the second question. The two disagree by design: CLEAN
# will report components that the permutation test calls non-significant, which
# is exactly what a fingerprint needs and exactly what an inferential claim
# must not use.

#' Design matrix for a set of frequencies at the observed positions.
#'
#' Positions are t0 = grid_index - 1, matching gls_prepare() and
#' reconstruct_signal(); see the phase convention note in grid.R.
clean_design <- function(k_vec, t, N) {
  t0 <- t - 1L
  X <- matrix(1, nrow = length(t), ncol = 1L)
  for (k in k_vec) {
    X <- cbind(X, cos(2 * pi * k * t0 / N), sin(2 * pi * k * t0 / N))
  }
  X
}

#' Joint least-squares fit of the currently selected components.
clean_fit <- function(y, k_vec, t, N) {
  X <- clean_design(k_vec, t, N)
  cf <- tryCatch(qr.solve(X, y), error = function(e) NULL)
  if (is.null(cf)) return(NULL)
  fitted <- as.vector(X %*% cf)
  list(coef = cf, fitted = fitted, residual = y - fitted,
       rss = sum((y - fitted)^2))
}

bic_of <- function(rss, n, n_components, penalty_factor = 1,
                   n_candidates = NULL, ebic_gamma = 1) {
  if (!is.finite(rss) || rss <= 0) return(-Inf)
  p <- 2 * n_components + 1
  selection_cost <- if (is.null(n_candidates) || n_components == 0L) 0 else
    2 * ebic_gamma * lchoose(n_candidates, n_components)
  n * log(rss / n) + penalty_factor * (p * log(n) + selection_cost)
}

#' CLEAN decomposition of one signal on one chromosome.
#'
#' @param y values at the observed positions, in the order of terms$t_index
#' @param terms output of gls_prepare()
#' @param max_components hard ceiling, a guard against pathological input rather
#'   than a modelling choice; BIC normally stops first
#' @return data.frame of components in the order they were extracted, with the
#'   BIC drop each one bought and the variance it explains in the joint fit
clean_decompose <- function(y, terms, max_components = 20L, penalty_factor = 1,
                            min_bic_drop = 0, ebic_gamma = 1) {
  y <- as.numeric(y)
  y[!is.finite(y)] <- 0
  t <- terms$t_index
  N <- terms$N
  n <- length(y)
  if (n < 8L || stats::sd(y) == 0) return(NULL)

  total_ss <- sum((y - mean(y))^2)
  selected <- integer(0)
  residual <- y - mean(y)
  n_candidates <- length(terms$k)
  bic_current <- bic_of(total_ss, n, 0L, penalty_factor, n_candidates, ebic_gamma)
  rows <- list()

  for (step in seq_len(max_components)) {
    sp <- gls_spectrum(residual, terms)
    sp <- sp[!sp$k %in% selected, , drop = FALSE]
    if (!nrow(sp)) break
    k_new <- sp$k[which.max(sp$power)]

    trial <- c(selected, k_new)
    fit <- clean_fit(y, trial, t, N)
    if (is.null(fit)) break
    bic_new <- bic_of(fit$rss, n, length(trial), penalty_factor,
                      n_candidates, ebic_gamma)
    drop <- bic_current - bic_new
    if (!is.finite(drop) || drop <= min_bic_drop) break

    # Amplitude and phase of this component in the JOINT fit, not in the
    # residual spectrum: the joint values are what the model actually says.
    j <- 1L + 2L * length(trial)
    a <- fit$coef[j - 1L]; b <- fit$coef[j]
    rows[[step]] <- data.frame(
      component = step, k = k_new, freq = k_new / N, period = N / k_new,
      amplitude = sqrt(a^2 + b^2), phase = atan2(-b, a),
      power = (a^2 + b^2) / 2,
      bic_drop = drop, bic = bic_new,
      rss = fit$rss,
      var_explained_cumulative = 1 - fit$rss / total_ss,
      window_power = terms$window_power[match(k_new, terms$k)],
      N = N, n_observed = n,
      stringsAsFactors = FALSE)

    selected <- trial
    residual <- fit$residual
    bic_current <- bic_new
  }

  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  out$var_explained <- c(out$var_explained_cumulative[1],
                         diff(out$var_explained_cumulative))
  # Final amplitudes and phases from the fit with every component present.
  final <- clean_fit(y, out$k, t, N)
  if (!is.null(final)) {
    idx <- 1L + 2L * seq_len(nrow(out))
    a <- final$coef[idx - 1L]; b <- final$coef[idx]
    out$amplitude <- sqrt(a^2 + b^2)
    out$phase <- atan2(-b, a)
    out$power <- (a^2 + b^2) / 2
  }
  out$window_rank <- rank(-terms$window_power)[match(out$k, terms$k)]
  out
}

#' CLEAN over every chromosome of one signal set.
clean_condition <- function(dataset, cond, chrom_idx, branch = "average",
                            per_sample = FALSE, clean_cfg = list(),
                            n_cores = 1L) {
  sig <- condition_signals(dataset, cond)
  if (is.null(sig)) return(NULL)
  summary_signal <- if (branch == "average") sig$avg_signal else sig$median_signal
  signals <- if (per_sample) {
    stats::setNames(lapply(sig$samples, function(s) sig$matrix[, s]), sig$samples)
  } else {
    stats::setNames(list(summary_signal), branch)
  }

  per_chr <- parallel::mclapply(names(chrom_idx), function(chr_now) {
    ci <- chrom_idx[[chr_now]]
    shared <- gls_prepare(ci$t, ci$N)
    rows <- list()
    for (nm in names(signals)) {
      # Deflation subtracts a fitted sinusoid from the residual, so a
      # zero-filled hole would be fitted and then removed as though it were
      # signal. Drop the position instead.
      fit <- gls_observed(signals[[nm]][ci$rows], ci$t, ci$N, terms = shared)
      if (is.null(fit)) next
      cp <- clean_decompose(fit$y, fit$terms,
                            max_components = clean_cfg$max_components %||% 20L,
                            penalty_factor = clean_cfg$penalty_factor %||% 1,
                            ebic_gamma = clean_cfg$ebic_gamma %||% 1)
      if (is.null(cp)) next
      cp$chr <- chr_now
      cp$signal <- nm
      cp$coverage <- fit$terms$n / ci$N
      rows[[length(rows) + 1]] <- cp
    }
    if (!length(rows)) NULL else do.call(rbind, rows)
  }, mc.cores = n_cores, mc.set.seed = FALSE)

  per_chr <- per_chr[!vapply(per_chr, is.null, logical(1))]
  if (!length(per_chr)) return(NULL)
  out <- do.call(rbind, per_chr)
  out$condition <- cond
  out$branch <- branch
  out[order(out$chr, out$component), ]
}
