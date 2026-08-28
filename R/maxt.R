# maxt.R -- maxT permutation test for spectral peaks.
#
# Null distribution: for each of B permutations, the maximum power over all k.
# Comparing every observed power against that maximum controls the family-wise
# error rate across frequencies within a chromosome. Two permutation schemes are
# used: full (destroys all positional structure) and block (preserves local
# correlation at the given block sizes). The primary p-value is the full scheme;
# the max across schemes is reported as the conservative alternative.
#
# Direct port of permutation_spectrum_test(); seeds are derived per chromosome
# and per sample exactly as before, so results are reproducible across runs and
# independent of how the work is split across cores.

permutation_spectrum_test <- function(signal, B = 1000L, seed = 42L,
                                      block_sizes = c(10L, 20L, 50L)) {
  signal <- as.numeric(signal)
  signal[!is.finite(signal)] <- 0
  signal <- signal - mean(signal)

  N <- length(signal)
  if (N < 8 || all(signal == 0)) return(NULL)

  k <- seq_len(floor((N - 1) / 2))
  spectral_power <- function(x) {
    Y <- stats::fft(x)
    (Mod(Y[k + 1]) / N)^2 * 2
  }
  power_observed <- spectral_power(signal)

  permute_blocks <- function(x, block_size) {
    block_id <- ceiling(seq_along(x) / block_size)
    blocks <- split(x, block_id)
    unlist(blocks[sample.int(length(blocks))], use.names = FALSE)
  }

  block_sizes <- unique(as.integer(block_sizes[is.finite(block_sizes) & block_sizes >= 2]))
  scheme_names <- c("full", paste0("block", block_sizes))
  null_max <- matrix(NA_real_, nrow = B, ncol = length(scheme_names),
                     dimnames = list(NULL, scheme_names))
  set.seed(seed)
  for (b in seq_len(B)) {
    null_max[b, "full"] <- max(spectral_power(sample(signal, replace = FALSE)))
    for (bs in block_sizes) {
      null_max[b, paste0("block", bs)] <- max(spectral_power(permute_blocks(signal, bs)))
    }
  }

  p_by_scheme <- vapply(seq_len(ncol(null_max)), function(j) {
    vapply(power_observed, function(pk) (1 + sum(null_max[, j] >= pk)) / (B + 1), numeric(1))
  }, numeric(length(power_observed)))
  if (is.null(dim(p_by_scheme))) p_by_scheme <- matrix(p_by_scheme, ncol = 1)
  colnames(p_by_scheme) <- paste0("p_empirical_maxT_", scheme_names)

  p_primary <- p_by_scheme[, "p_empirical_maxT_full"]
  p_all <- apply(p_by_scheme, 1, max, na.rm = TRUE)

  out <- data.frame(
    N = N, k = k, freq = k / N, period = N / k,
    power = power_observed,
    p_empirical_maxT = p_primary,
    significant = p_primary <= 0.05,
    significant_all_schemes = p_all <= 0.05,
    stringsAsFactors = FALSE)
  cbind(out, as.data.frame(p_by_scheme, stringsAsFactors = FALSE))
}

#' Run maxT for every sample of one condition, parallel over chromosomes.
maxt_condition <- function(dataset, cond, chrom_idx, maxt_cfg, n_cores = 1L) {
  sig <- condition_signals(dataset, cond)
  if (is.null(sig)) return(NULL)
  chrom_levels <- names(chrom_idx)

  per_chr <- parallel::mclapply(chrom_levels, function(chr_now) {
    ord <- chrom_idx[[chr_now]]
    chr_seed <- maxt_cfg$seed + 100000L * match(chr_now, chrom_levels)
    rows <- list()
    for (s in sig$samples) {
      res <- permutation_spectrum_test(
        sig$matrix[ord, s],
        B = maxt_cfg$B,
        seed = chr_seed + match(s, sig$samples),
        block_sizes = maxt_cfg$block_sizes)
      if (is.null(res) || !nrow(res)) next
      res$chr <- chr_now
      res$sample <- s
      rows[[length(rows) + 1]] <- res
    }
    if (!length(rows)) NULL else do.call(rbind, rows)
  }, mc.cores = n_cores, mc.set.seed = FALSE)

  per_chr <- per_chr[!vapply(per_chr, is.null, logical(1))]
  if (!length(per_chr)) return(NULL)
  out <- do.call(rbind, per_chr)
  out$condition <- cond
  attr(out, "expected_samples") <- sig$samples
  out
}

#' Cores to use: leave one free, never more than the number of chromosomes.
maxt_cores <- function(chrom_idx) {
  detected <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (is.na(detected) || detected < 1) return(1L)
  max(1L, min(detected - 1L, length(chrom_idx)))
}
