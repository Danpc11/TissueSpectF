# maxt.R -- permutation test for spectral peaks on a gappy reference grid.
#
# Null: the observed values are permuted AMONG THE OBSERVED POSITIONS, with the
# positions themselves held fixed. The sampling pattern is therefore identical
# in the data and in every permutation, so the pattern of missing genes cannot
# by itself produce significance -- which is the guarantee that zero-filling or
# mean-imputation would destroy.
#
# For each permutation the maximum power over all k is recorded; comparing an
# observed power against that maximum controls the family-wise error rate across
# frequencies within a chromosome. This matters more with gaps than without,
# because the frequency bins are no longer orthogonal and a per-bin p-value
# would be badly calibrated.
#
# Two schemes are used: `full` (permute all observed values) and `block`
# (permute contiguous blocks of grid positions, preserving local correlation).
# The primary p-value is `full`; `p_empirical_maxT_all` is the max across
# schemes, reported for the conservative reading.

permutation_gls_test <- function(y, terms, B = 1000L, seed = 42L,
                                 block_sizes = c(10L, 20L, 50L)) {
  y <- as.numeric(y)
  y[!is.finite(y)] <- 0
  n <- terms$n
  if (n < 8L || stats::sd(y) == 0) return(NULL)

  obs <- gls_spectrum(y, terms)
  power_observed <- obs$power

  permute_blocks <- function(v, block_size) {
    block_id <- ceiling(seq_along(v) / block_size)
    blocks <- split(v, block_id)
    unlist(blocks[sample.int(length(blocks))], use.names = FALSE)
  }

  # Every declared scheme keeps a column, even when it cannot run on this
  # chromosome: a block size at or above the number of observed genes would
  # permute a single block, which is not a permutation at all. Those columns
  # stay NA so the table has the same shape for every chromosome -- chromosomes
  # differ in how many genes are observed, and dropping columns per chromosome
  # made the per-chromosome results impossible to rbind.
  block_sizes <- unique(as.integer(block_sizes[is.finite(block_sizes) & block_sizes >= 2]))
  usable <- block_sizes[block_sizes < n / 2]
  skipped <- setdiff(block_sizes, usable)
  scheme_names <- c("full", paste0("block", block_sizes))
  null_max <- matrix(NA_real_, nrow = B, ncol = length(scheme_names),
                     dimnames = list(NULL, scheme_names))
  set.seed(seed)
  for (b in seq_len(B)) {
    null_max[b, "full"] <- max(gls_spectrum(sample(y), terms)$power)
    for (bs in usable) {
      null_max[b, paste0("block", bs)] <-
        max(gls_spectrum(permute_blocks(y, bs), terms)$power)
    }
  }

  p_by_scheme <- vapply(seq_len(ncol(null_max)), function(j) {
    if (all(is.na(null_max[, j]))) return(rep(NA_real_, length(power_observed)))
    vapply(power_observed, function(pk) (1 + sum(null_max[, j] >= pk)) / (B + 1), numeric(1))
  }, numeric(length(power_observed)))
  if (is.null(dim(p_by_scheme))) p_by_scheme <- matrix(p_by_scheme, ncol = 1)
  colnames(p_by_scheme) <- paste0("p_empirical_maxT_", scheme_names)

  p_primary <- p_by_scheme[, "p_empirical_maxT_full"]
  p_all <- apply(p_by_scheme, 1, function(v)
    if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE))

  out <- obs
  out$p_empirical_maxT <- p_primary
  out$p_empirical_maxT_all <- p_all
  out$significant <- p_primary <= 0.05
  out$significant_all_schemes <- !is.na(p_all) & p_all <= 0.05
  out$block_schemes_skipped <- if (length(skipped)) paste(skipped, collapse = ",") else NA_character_
  cbind(out, as.data.frame(p_by_scheme, stringsAsFactors = FALSE))
}

#' maxT for every sample of one condition, parallel over chromosomes.
maxt_condition <- function(dataset, cond, chrom_idx, maxt_cfg, n_cores = 1L) {
  sig <- condition_signals(dataset, cond)
  if (is.null(sig)) return(NULL)
  chrom_levels <- names(chrom_idx)

  per_chr <- parallel::mclapply(chrom_levels, function(chr_now) {
    ci <- chrom_idx[[chr_now]]
    terms <- gls_prepare(ci$t, ci$N)     # window terms reused across samples
    chr_seed <- maxt_cfg$seed + 100000L * match(chr_now, chrom_levels)
    rows <- list()
    for (s in sig$samples) {
      res <- permutation_gls_test(sig$matrix[ci$rows, s], terms,
                                  B = maxt_cfg$B,
                                  seed = chr_seed + match(s, sig$samples),
                                  block_sizes = maxt_cfg$block_sizes)
      if (is.null(res) || !nrow(res)) next
      res$chr <- chr_now
      res$sample <- s
      res$coverage <- ci$coverage
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

maxt_cores <- function(chrom_idx) {
  detected <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (is.na(detected) || detected < 1) return(1L)
  max(1L, min(detected - 1L, length(chrom_idx)))
}
