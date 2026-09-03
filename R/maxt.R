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

MIN_PERMUTATION_BLOCKS <- 10L

permutation_gls_test <- function(y, terms, B = 1000L, seed = 42L,
                                 block_sizes = c(10L, 20L, 50L),
                                 primary_scheme = "full") {
  y <- as.numeric(y)
  y[!is.finite(y)] <- 0
  n <- terms$n
  if (n < 8L || stats::sd(y) == 0) return(NULL)

  obs <- gls_spectrum(y, terms)
  power_observed <- obs$power

  # Blocks are intervals of the REFERENCE GRID, not runs of consecutive entries
  # in the observed list. With gaps the two differ badly: at 20% coverage a run
  # of 10 observed genes spans ~50 grid positions, so grouping by list index
  # would treat genes 50 slots apart as neighbours and the scheme would preserve
  # no local structure at all.
  #
  # Blocks therefore hold a variable number of observed genes. The permutation
  # shuffles the order of the blocks and writes the concatenated values back
  # onto the same observed positions, in order: positions are untouched, values
  # keep their within-block ordering, and long-range structure is destroyed.
  # The compromise is that when block sizes differ, the spacing between a value
  # and its neighbours is not preserved exactly inside a relocated block; there
  # is no permutation that keeps both the position set and every within-block
  # distance, and holding the positions fixed is the property that matters here.
  grid_blocks <- function(block_size) {
    split(seq_along(terms$t_index), ceiling(terms$t_index / block_size))
  }
  permute_blocks <- function(v, blocks) {
    unlist(lapply(blocks[sample.int(length(blocks))], function(i) v[i]),
           use.names = FALSE)
  }

  # Every declared scheme keeps a column, even when it cannot run on this
  # chromosome: a block size at or above the number of observed genes would
  # permute a single block, which is not a permutation at all. Those columns
  # stay NA so the table has the same shape for every chromosome -- chromosomes
  # differ in how many genes are observed, and dropping columns per chromosome
  # made the per-chromosome results impossible to rbind.
  # Pointwise exceedance counts, accumulated alongside the maxima at no extra
  # cost. The maxT p-value asks whether a frequency beats the strongest
  # frequency of an entire permuted spectrum -- family-wise control across the
  # chromosome, and with 500 to 800 frequencies that is a very high bar: it
  # effectively asks whether this is the dominant component of the chromosome.
  # The pointwise p-value asks whether the frequency beats its own null, which
  # is the question a signature is about; the multiplicity is then handled by
  # BH across frequencies rather than by taking a maximum.
  #
  # A per-frequency p-value computed against that frequency's own null has the
  # same floor as any permutation p, 1/(B+1), and BH across ~10,000 genome-wide
  # frequencies then needs B >= n_freq/q -- two hundred thousand draws for
  # q = 0.05. Infeasible, and the fourth time this floor has bitten in this
  # project.
  #
  # The pooled null is the standard way out. Each frequency's null is
  # standardised by its own mean and standard deviation, and the standardised
  # values are pooled across the frequencies of the chromosome, giving an
  # empirical null of size B*K and a floor of 1/(B*K+1) -- 6e-7 with 200 draws
  # over 600 frequencies, small enough for BH to reach anything.
  #
  # The assumption is that the STANDARDISED null is exchangeable across the
  # frequencies of a chromosome. Standardising removes the systematic difference
  # between frequencies; what it does not remove is the effect of the sampling
  # window, which is why window_power travels with every peak and ./tsf window
  # exists. The unpooled family-wise p-value is kept alongside and makes no such
  # assumption.
  null_power_mat <- matrix(NA_real_, nrow = B, ncol = length(power_observed))

  block_sizes <- unique(as.integer(block_sizes[is.finite(block_sizes) & block_sizes >= 2]))
  block_map <- lapply(stats::setNames(block_sizes, block_sizes), grid_blocks)
  # A scheme needs enough blocks for shuffling to have resolution. With only 4
  # blocks there are 24 distinct orderings, so a sizeable fraction of the null
  # draws keep the signal nearly intact and no p-value can get small. Such a
  # scheme would then veto every peak under primary_scheme = "all" -- not
  # because the peak is local structure, but because the test cannot see. 10
  # blocks (3.6 million orderings) is the floor.
  n_blocks <- vapply(block_map, length, integer(1))
  usable <- block_sizes[n_blocks >= MIN_PERMUTATION_BLOCKS]
  skipped <- setdiff(block_sizes, usable)
  scheme_names <- c("full", paste0("block", block_sizes))
  null_max <- matrix(NA_real_, nrow = B, ncol = length(scheme_names),
                     dimnames = list(NULL, scheme_names))
  set.seed(seed)
  for (b in seq_len(B)) {
    null_power <- gls_spectrum(sample(y), terms)$power
    null_max[b, "full"] <- max(null_power)
    null_power_mat[b, ] <- null_power
    for (bs in usable) {
      null_max[b, paste0("block", bs)] <-
        max(gls_spectrum(permute_blocks(y, block_map[[as.character(bs)]]), terms)$power)
    }
  }

  p_by_scheme <- vapply(seq_len(ncol(null_max)), function(j) {
    if (all(is.na(null_max[, j]))) return(rep(NA_real_, length(power_observed)))
    vapply(power_observed, function(pk) (1 + sum(null_max[, j] >= pk)) / (B + 1), numeric(1))
  }, numeric(length(power_observed)))
  if (is.null(dim(p_by_scheme))) p_by_scheme <- matrix(p_by_scheme, ncol = 1)
  colnames(p_by_scheme) <- paste0("p_empirical_maxT_", scheme_names)

  p_full <- p_by_scheme[, "p_empirical_maxT_full"]
  p_all <- apply(p_by_scheme, 1, function(v)
    if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE))
  # `full` destroys every scale of structure, so it is the most permissive null.
  # `all` requires the peak to survive the block schemes too, i.e. to be more
  # than local autocorrelation. Which one is primary is a scientific choice and
  # is declared in config/project.R, never left implicit here.
  p_primary <- switch(primary_scheme,
                      full = p_full,
                      all  = p_all,
                      tsf_abort("maxt$primary_scheme must be 'full' or 'all'"))

  # The pooled pointwise route standardises each frequency by its own null mean
  # and sd, then pools all B x m standardised values into one reference
  # distribution. Its resolution is therefore 1/(B*m + 1), far finer than
  # 1/(B + 1) -- which is legitimate only if the standardised nulls really are
  # exchangeable across frequencies, and only if the sd is estimable at all.
  #
  # It is not estimable from one draw. With B = 1 every sd is NA, the guard
  # below replaced it with 1, the pooled reference collapsed to a vector of
  # zeros, and every frequency above its single null draw received
  # p = 1/(B*m + 1) ~ 1e-4. On a real run that produced 4950 "discoveries" at
  # BH q <= 0.05 out of 10030 -- a spectacular-looking result manufactured
  # entirely by a degenerate standardisation. Substituting 1 for a quantity that
  # could not be computed turns "unknown" into "confident", which is the one
  # substitution a test must never make.
  centre <- colMeans(null_power_mat)
  scale <- apply(null_power_mat, 2, stats::sd)

  MIN_DRAWS_FOR_POOLED <- 10L
  if (B < MIN_DRAWS_FOR_POOLED) {
    tsf_warn("B = ", B, " is too few draws to standardise the pooled pointwise ",
             "null (needs at least ", MIN_DRAWS_FOR_POOLED, "). p_pointwise is ",
             "reported as NA rather than as a number the draws do not support; ",
             "the family-wise maxT route is unaffected.")
    p_pooled <- rep(NA_real_, length(power_observed))
    z_obs <- rep(NA_real_, length(power_observed))
  } else {
    # A frequency whose null has no spread has no scale to divide by. It is
    # dropped from the pointwise route rather than assigned an arbitrary one.
    degenerate <- !is.finite(scale) | scale <= 0
    if (any(degenerate)) {
      tsf_warn(sum(degenerate), " of ", length(scale), " frequencies have a ",
               "null with no spread; p_pointwise is NA for those.")
    }
    scale[degenerate] <- NA_real_
    z_obs <- (power_observed - centre) / scale
    pooled <- as.vector(sweep(sweep(null_power_mat, 2, centre), 2, scale, "/"))
    pooled <- sort(pooled[is.finite(pooled)])
    if (!length(pooled)) {
      p_pooled <- rep(NA_real_, length(z_obs))
    } else {
      # (1 + number of pooled null values at least as extreme) / (1 + pool size)
      p_pooled <- (1 + (length(pooled) -
                          findInterval(z_obs, pooled, left.open = TRUE))) /
        (length(pooled) + 1)
      p_pooled[!is.finite(z_obs)] <- NA_real_
    }
  }

  out <- obs
  out$p_pointwise <- pmin(1, p_pooled)
  out$null_z <- z_obs
  out$p_empirical_maxT <- p_primary
  out$p_empirical_maxT_full <- p_full
  out$p_empirical_maxT_all <- p_all
  out$primary_scheme <- primary_scheme
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
    shared <- gls_prepare(ci$t, ci$N)    # window terms reused across samples
    chr_seed <- maxt_cfg$seed + 100000L * match(chr_now, chrom_levels)
    rows <- list()
    for (s in sig$samples) {
      # The permutation null must permute among the positions this sample was
      # actually measured at. Handing it a zero-filled hole would put a
      # constant in the observed set, and every permutation would carry it.
      fit <- gls_observed(sig$matrix[ci$rows, s], ci$t, ci$N, terms = shared)
      if (is.null(fit)) next
      res <- permutation_gls_test(fit$y, fit$terms,
                                  B = maxt_cfg$B,
                                  seed = chr_seed + match(s, sig$samples),
                                  block_sizes = maxt_cfg$block_sizes,
                                  primary_scheme = maxt_cfg$primary_scheme %||% "full")
      if (is.null(res) || !nrow(res)) next
      res$chr <- chr_now
      res$sample <- s
      res$coverage <- fit$terms$n / ci$N
      res$n_dropped_unmeasured <- fit$n_dropped
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

#' Workers for the chromosome-parallel stages: maxT, condition, CLEAN.
#'
#' This used to read detectCores() directly and ignore both --cores and
#' N_WORKERS, so the documented flag silently did nothing on three stages --
#' including maxT, which the help calls the slow one. It now goes through
#' local_workers() like consensus and compare do.
#'
#' `default_max` is the chromosome count rather than local_workers()' usual 8,
#' which preserves what this function did when no flag is given: one worker per
#' chromosome, minus one core for the parent.
#'
#' The ceiling is structural, not a setting. These stages split on chromosomes,
#' so with 24 chromosomes the 25th core has nothing to do. Going wider means
#' splitting on (chromosome, sample) instead, which is a change to the stage,
#' not a larger number on the command line.
maxt_cores <- function(chrom_idx, opt = NULL) {
  n_chr <- max(1L, length(chrom_idx))
  local_workers(opt %||% list(), n_tasks = n_chr, default_max = n_chr)
}
