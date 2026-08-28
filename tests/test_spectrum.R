#!/usr/bin/env Rscript
# Numerical tests for the spectral core. Run: Rscript tests/test_spectrum.R
source("R/utils_io.R"); source("R/config.R"); source("R/labels.R")
source("R/spectrum.R"); source("R/maxt.R"); source("R/stability.R")
source("R/peaks_genes.R"); source("R/compare.R")

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

# --- maxT --------------------------------------------------------------------
mt <- permutation_spectrum_test(sig + rnorm(N, sd = 0.2), B = 200L, seed = 1L,
                                block_sizes = c(10L, 20L))
check("maxT returns one row per positive frequency", nrow(mt) == floor((N - 1) / 2))
check("maxT flags a strong periodic component", mt$p_empirical_maxT[mt$k == k_true] <= 0.05)
check("maxT p-values are bounded by 1/(B+1)", min(mt$p_empirical_maxT) >= 1 / 201)
check("block schemes are reported", all(c("p_empirical_maxT_full",
                                          "p_empirical_maxT_block10") %in% colnames(mt)))
check("maxT is reproducible for a fixed seed",
      identical(permutation_spectrum_test(sig, B = 50L, seed = 7L)$p_empirical_maxT,
                permutation_spectrum_test(sig, B = 50L, seed = 7L)$p_empirical_maxT))
check("pure noise rarely reaches significance", {
  set.seed(5)
  m <- permutation_spectrum_test(rnorm(N), B = 200L, seed = 2L)
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

# --- Wilson ------------------------------------------------------------------
check("Wilson interval brackets the point estimate", {
  ci <- wilson_ci(9, 10); ci[1] < 90 && ci[2] > 90 && ci[1] >= 0 && ci[2] <= 100 })
check("Wilson handles n = 0", all(is.na(wilson_ci(0, 0))))

cat("\n", if (failures == 0L) "All tests passed." else paste(failures, "test(s) failed."), "\n")
quit(status = if (failures == 0L) 0L else 1L)
