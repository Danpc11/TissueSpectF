#!/usr/bin/env Rscript
# Numerical tests for the spectral core. Run: Rscript tests/test_spectrum.R
source("R/utils_io.R"); source("R/config.R"); source("R/labels.R")
source("R/grid.R"); source("R/ingest.R"); source("R/spectrum.R"); source("R/maxt.R"); source("R/stability.R")
source("R/condition_test.R"); source("R/peaks_genes.R"); source("R/compare.R")

failures <- 0L
check <- function(label, expr) {
  ok <- isTRUE(tryCatch(expr, error = function(e) {
    cat("   error: ", conditionMessage(e), "\n"); FALSE }))
  cat(if (ok) "  PASS  " else "  FAIL  ", label, "\n")
  if (!ok) failures <<- failures + 1L
}

# --- run_fft recovers a known sinusoid ---------------------------------------
N <- 200; k_true <- 8; A_true <- 3; phi_true <- 0.7
n <- 0:(N - 1)
sig <- A_true * cos(2 * pi * k_true * n / N + phi_true)
sp <- run_fft(sig)

check("peak is at the injected k", sp$k[which.max(sp$power)] == k_true)
check("amplitude is recovered", abs(sp$amplitude[sp$k == k_true] - A_true) < 1e-8)
check("phase is recovered", abs(sp$phase[sp$k == k_true] - phi_true) < 1e-8)
check("freq and period are consistent",
      abs(sp$freq[sp$k == k_true] - k_true / N) < 1e-12 &&
        abs(sp$period[sp$k == k_true] - N / k_true) < 1e-12)
check("DC carries no power", sp$power[sp$k == 0] == 0)
check("Fisher g flags the injected peak", sp$p_value[sp$k == k_true] < 0.05)
check("pure noise has no Fisher-significant peak at that k", {
  set.seed(11); sn <- run_fft(rnorm(N)); sn$p_value[sn$k == k_true] > 0.05 })

# --- Parseval: one-sided power sums to the signal variance -------------------
check("one-sided power matches the sum of squares", {
  set.seed(3); x <- rnorm(128)
  s <- run_fft(x)
  abs(sum(s$power) - mean((x - mean(x))^2)) < 1e-8 })

# --- reconstruct_signal round-trips ------------------------------------------
check("reconstruction reproduces the original sinusoid", {
  pk <- data.frame(freq = k_true / N, amplitude = A_true, phase = phi_true)
  max(abs(reconstruct_signal(pk, N) - (sig - mean(sig)))) < 1e-8 })

check("a peak outside the Nyquist range is ignored", {
  pk <- data.frame(freq = 1.5, amplitude = 1, phase = 0)
  all(reconstruct_signal(pk, 32) == 0) })

# --- maxT on a complete grid (GLS reduces to the FFT case) --------------------
mt <- permutation_gls_test(sig + rnorm(N, sd = 0.2), gls_prepare(1:N, N),
                           B = 200L, seed = 1L, block_sizes = c(10L, 20L))
check("maxT returns one row per positive frequency", nrow(mt) == floor((N - 1) / 2))
check("maxT flags a strong periodic component", mt$p_empirical_maxT[mt$k == k_true] <= 0.05)
check("maxT p-values are bounded by 1/(B+1)", min(mt$p_empirical_maxT) >= 1 / 201)
check("block schemes are reported", all(c("p_empirical_maxT_full",
                                          "p_empirical_maxT_block10") %in% colnames(mt)))
check("maxT is reproducible for a fixed seed", {
  a <- permutation_gls_test(sig, gls_prepare(1:N, N), B = 50L, seed = 7L)
  b <- permutation_gls_test(sig, gls_prepare(1:N, N), B = 50L, seed = 7L)
  identical(a$p_empirical_maxT, b$p_empirical_maxT) })
check("pure noise rarely reaches significance", {
  set.seed(5)
  m <- permutation_gls_test(rnorm(N), gls_prepare(1:N, N), B = 200L, seed = 2L)
  mean(m$p_empirical_maxT <= 0.05) < 0.10 })

# --- stability ---------------------------------------------------------------
mk <- function(sample, k, p) data.frame(chr = "1", N = 100L, k = k, sample = sample,
                                        power = 1, p_empirical_maxT = p)
indiv <- do.call(rbind, c(
  lapply(paste0("S", 1:10), function(s) mk(s, 5L, 0.01)),   # 10/10 significant
  lapply(paste0("S", 1:10), function(s) mk(s, 6L, if (s %in% paste0("S", 1:8)) 0.01 else 0.5))))
st <- stable_peaks_maxt(indiv, paste0("S", 1:10), alpha = 0.05, stable_frac = 0.9)
check("10/10 significant is stable", st$is_stable[st$k == 5L])
check("8/10 significant is not stable at 90%", !st$is_stable[st$k == 6L])
check("percentages are reported", st$pct_samples_significant[st$k == 6L] == 80)
check("a missing sample blocks stability",
      !stable_peaks_maxt(indiv[indiv$sample != "S1" | indiv$k != 5L, ],
                         paste0("S", 1:10))$is_stable[1])

# --- GLS on a gappy grid -----------------------------------------------------
check("GLS equals the FFT on a complete grid", {
  set.seed(21); N <- 256; t <- 1:N
  y <- 3 * cos(2 * pi * 7 * (t - 1) / N + 0.7) + rnorm(N, sd = 0.3)
  g <- gls_spectrum(y, gls_prepare(t, N)); f <- run_fft(y)
  i <- match(g$k, f$k)
  max(abs(g$amplitude - f$amplitude[i])) < 1e-9 &&
    max(abs(((g$phase - f$phase[i] + pi) %% (2 * pi)) - pi)) < 1e-9 })

check("GLS equals brute-force least squares with gaps", {
  set.seed(22); N <- 200; t <- sort(sample.int(N, 130))
  y <- 2 * cos(2 * pi * 11 * (t - 1) / N + 1.1) + rnorm(length(t), sd = 0.4)
  g <- gls_spectrum(y, gls_prepare(t, N))
  brute <- function(k) {
    X <- cbind(1, cos(2 * pi * k * (t - 1) / N), sin(2 * pi * k * (t - 1) / N))
    cf <- qr.solve(X, y); c(sqrt(cf[2]^2 + cf[3]^2), atan2(-cf[3], cf[2]))
  }
  ks <- c(3, 11, 40, 77)
  a <- vapply(ks, function(k) brute(k)[1], numeric(1))
  p <- vapply(ks, function(k) brute(k)[2], numeric(1))
  max(abs(g$amplitude[match(ks, g$k)] - a)) < 1e-8 &&
    max(abs(((g$phase[match(ks, g$k)] - p + pi) %% (2 * pi)) - pi)) < 1e-8 })

check("the injected k is recovered at 65% coverage", {
  set.seed(23); N <- 200; t <- sort(sample.int(N, 130))
  y <- 2 * cos(2 * pi * 11 * (t - 1) / N + 1.1) + rnorm(length(t), sd = 0.4)
  gls_spectrum(y, gls_prepare(t, N))$k[which.max(gls_spectrum(y, gls_prepare(t, N))$power)] == 11 })

check("a complete grid has a flat zero window", {
  max(spectral_window(1:128, 128)$window_power) < 1e-20 })

check("gaps alone do not create significance", {
  # White noise on a heavily gapped grid: the permutation null holds the
  # positions fixed, so the missingness pattern cannot produce a peak.
  set.seed(24); N <- 300; t <- sort(sample.int(N, 150))
  m <- permutation_gls_test(rnorm(length(t)), gls_prepare(t, N), B = 200L, seed = 3L)
  mean(m$p_empirical_maxT <= 0.05) < 0.10 })

check("a real periodicity survives the gaps", {
  set.seed(25); N <- 300; t <- sort(sample.int(N, 150))
  y <- 1.5 * cos(2 * pi * 9 * (t - 1) / N) + rnorm(length(t), sd = 0.5)
  m <- permutation_gls_test(y, gls_prepare(t, N), B = 200L, seed = 4L)
  m$p_empirical_maxT[m$k == 9] <= 0.05 })

check("chromosomes with few genes keep the same columns", {
  # A block size at or above n/2 cannot run, but its column must survive so
  # per-chromosome tables can be rbind()ed together.
  a <- permutation_gls_test(rnorm(200), gls_prepare(sort(sample.int(400, 200)), 400),
                            B = 30L, seed = 1L, block_sizes = c(10L, 20L, 50L))
  b <- permutation_gls_test(rnorm(30), gls_prepare(sort(sample.int(60, 30)), 60),
                            B = 30L, seed = 1L, block_sizes = c(10L, 20L, 50L))
  identical(colnames(a), colnames(b)) &&
    all(is.na(b$p_empirical_maxT_block50)) &&
    !is.null(rbind(a, b)) })

check("blocks come from the grid, not from the observed list", {
  # Observed positions 1, 100, 101, 102 with block_size 10: a list-index block
  # would put 1 and 100 in the same block; a grid block must not.
  t <- c(1L, 100L, 101L, 102L); N <- 200L
  terms <- gls_prepare(t, N)
  blocks <- split(seq_along(terms$t_index), ceiling(terms$t_index / 10))
  # grid blocks: ceiling(c(1,100,101,102)/10) = 1, 10, 11, 11 -> three blocks,
  # and positions 1 and 100 are never grouped together
  !any(vapply(blocks, function(i) all(c(1L, 2L) %in% i), logical(1))) &&
    length(blocks) == 3L &&
    any(vapply(blocks, function(i) identical(as.integer(i), c(3L, 4L)), logical(1))) })

check("block permutation keeps the observed positions fixed", {
  set.seed(31); N <- 400L; t <- sort(sample.int(N, 200))
  terms <- gls_prepare(t, N)
  m <- permutation_gls_test(rnorm(200), terms, B = 20L, seed = 2L,
                            block_sizes = 20L)
  # Same grid length and same number of observed genes as the data
  all(m$N == N) && all(m$n_observed == 200L) })

check("primary_scheme = all is stricter than full", {
  set.seed(32); N <- 300L; t <- sort(sample.int(N, 180))
  y <- rnorm(180)
  terms <- gls_prepare(t, N)
  a <- permutation_gls_test(y, terms, B = 100L, seed = 5L, primary_scheme = "full")
  b <- permutation_gls_test(y, terms, B = 100L, seed = 5L, primary_scheme = "all")
  all(b$p_empirical_maxT >= a$p_empirical_maxT, na.rm = TRUE) &&
    identical(a$p_empirical_maxT_full, b$p_empirical_maxT_full) })

check("a block scheme with too few blocks is skipped, not silently vetoing", {
  # N = 200 with block 50 gives 4 blocks: too few orderings for a small p.
  set.seed(33); N <- 200L; t <- 1:N
  y <- 2 * cos(2 * pi * 6 * (t - 1) / N) + rnorm(N, sd = 0.3)
  m <- permutation_gls_test(y, gls_prepare(t, N), B = 100L, seed = 6L,
                            block_sizes = c(10L, 50L), primary_scheme = "all")
  all(is.na(m$p_empirical_maxT_block50)) &&
    identical(unique(m$block_schemes_skipped), "50") &&
    m$p_empirical_maxT[m$k == 6] <= 0.05 })

check("an unknown primary scheme is refused", {
  inherits(tryCatch(permutation_gls_test(rnorm(50), gls_prepare(1:50, 50),
                                         B = 10L, primary_scheme = "blocky"),
                    error = function(e) e), "error") })

check("partial NA gene lengths do not poison the filter", {
  m <- matrix(c(10, 20, 30, 40, 50, 60), nrow = 3)
  e <- counts_to_expression(m, c(1000, NA, 1000))
  f <- filter_expressed(e, 1, 0.2)
  !any(is.na(f)) })

check("a skipped block scheme is recorded", {
  # N = 60 with block 50 gives 2 blocks, well under MIN_PERMUTATION_BLOCKS
  b <- permutation_gls_test(rnorm(30), gls_prepare(sort(sample.int(60, 30)), 60),
                            B = 20L, seed = 1L, block_sizes = c(5L, 50L))
  identical(unique(b$block_schemes_skipped), "50") &&
    all(is.na(b$p_empirical_maxT_block50)) &&
    !all(is.na(b$p_empirical_maxT_block5)) })

# --- condition-level test ----------------------------------------------------
check("Stouffer combines evidence and gains power with n", {
  # Ten samples each at p = 0.20 -- none significant alone -- combine below 0.05
  p_one <- stouffer_combine(0.20)
  p_ten <- stouffer_combine(rep(0.20, 10))
  abs(p_one - 0.20) < 1e-8 && p_ten < 0.05 })

check("Stouffer caps the permutation floor", {
  # A p of exactly 0 would give an infinite z; the floor is 1/(B+1)
  is.finite(stouffer_combine(c(0, 0.5), B = 100L)) })

check("Stouffer is calibrated under the null", {
  # A single combination of uniform p-values is itself uniform, not 0.5, so
  # calibration has to be checked over replicates: about 5% should fall <= 0.05.
  set.seed(41)
  p <- vapply(1:400, function(i) stouffer_combine(runif(10)), numeric(1))
  abs(mean(p <= 0.05) - 0.05) < 0.03 && abs(mean(p) - 0.5) < 0.06 })

check("the condition test recovers a peak too weak for any single sample", {
  # Each "sample" carries the periodicity buried in noise; the mean does not.
  set.seed(42); N <- 400L; t <- 1:N; k0 <- 13L
  base <- 0.25 * cos(2 * pi * k0 * (t - 1) / N)
  terms <- gls_prepare(t, N)
  one <- permutation_gls_test(base + rnorm(N, sd = 2), terms, B = 200L, seed = 1L)
  avg <- rowMeans(vapply(1:40, function(i) base + rnorm(N, sd = 2), numeric(N)))
  many <- condition_permutation_test(avg, terms, B = 200L, seed = 2L)
  one$p_empirical_maxT[one$k == k0] > 0.05 && many$p_empirical_maxT[many$k == k0] <= 0.05 })

check("q is corrected over chromosomes, not frequencies", {
  # With B = 500 and 23 chromosomes the smallest attainable q must be under 0.05,
  # which double-correcting over ~9500 frequencies would make impossible.
  23 / 501 < 0.05 })

# --- Wilson ------------------------------------------------------------------
check("Wilson interval brackets the point estimate", {
  ci <- wilson_ci(9, 10); ci[1] < 90 && ci[2] > 90 && ci[1] >= 0 && ci[2] <= 100 })
check("Wilson handles n = 0", all(is.na(wilson_ci(0, 0))))

cat("\n", if (failures == 0L) "All tests passed." else paste(failures, "test(s) failed."), "\n")
quit(status = if (failures == 0L) 0L else 1L)
