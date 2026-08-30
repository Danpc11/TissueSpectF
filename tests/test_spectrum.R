#!/usr/bin/env Rscript
# Numerical tests for the spectral core. Run: Rscript tests/test_spectrum.R
source("R/utils_io.R"); source("R/config.R"); source("R/labels.R")
source("R/grid.R"); source("R/ingest.R"); source("R/spectrum.R"); source("R/maxt.R"); source("R/stability.R")
source("R/condition_test.R"); source("R/clean.R"); source("R/fingerprint.R"); source("R/reference.R"); source("R/consensus.R"); source("R/peaks_genes.R"); source("R/compare.R")

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

# --- CLEAN + EBIC ------------------------------------------------------------
check("CLEAN recovers two injected components", {
  set.seed(51); N <- 600L; t <- sort(sample.int(N, 400))
  y <- 2 * cos(2 * pi * 7 * (t - 1) / N + 0.3) +
    1.2 * cos(2 * pi * 31 * (t - 1) / N - 1.0) + rnorm(400, sd = 0.4)
  cp <- clean_decompose(y, gls_prepare(t, N))
  all(c(7L, 31L) %in% cp$k) && nrow(cp) <= 4 &&
    abs(cp$amplitude[cp$k == 7] - 2) < 0.2 })

check("CLEAN does not report a sidelobe as a component", {
  # 50% coverage: the raw periodogram's top-5 contains several false maxima
  # alongside the real one. CLEAN must return the real one only.
  set.seed(52); N <- 600L; t <- sort(sample.int(N, 300))
  y <- 3 * cos(2 * pi * 11 * (t - 1) / N + 0.8) + rnorm(300, sd = 0.3)
  terms <- gls_prepare(t, N)
  top5 <- gls_spectrum(y, terms)$k[order(-gls_spectrum(y, terms)$power)][1:5]
  cp <- clean_decompose(y, terms)
  length(setdiff(top5, 11L)) >= 3 && nrow(cp) == 1L && cp$k[1] == 11L })

check("EBIC selects nothing from pure noise", {
  set.seed(53); N <- 600L; t <- sort(sample.int(N, 300))
  is.null(clean_decompose(rnorm(300), gls_prepare(t, N))) })

check("false positive rate under the null stays low", {
  set.seed(54); N <- 600L; terms <- gls_prepare(sort(sample.int(N, 300)), N)
  fp <- vapply(1:30, function(i) {
    c4 <- clean_decompose(rnorm(300), terms); if (is.null(c4)) 0L else nrow(c4)
  }, integer(1))
  mean(fp > 0) <= 0.1 })

check("plain BIC would have selected noise (why EBIC)", {
  # Same null, selection cost switched off: components appear in most draws.
  # Measured over replicates rather than one draw, since the failure is a rate.
  # gamma = 0 selects something in ~65% of null draws; gamma = 1 in none.
  # Kept as a test so nobody 'simplifies' the penalty away.
  set.seed(53); N <- 600L; terms <- gls_prepare(sort(sample.int(N, 300)), N)
  rate <- function(g) mean(vapply(1:20, function(i) {
    c4 <- clean_decompose(rnorm(300), terms, ebic_gamma = g)
    !is.null(c4)
  }, logical(1)))
  rate(0) > 0.4 && rate(1) < 0.1 })

check("variance explained accumulates and BIC drops each step", {
  set.seed(55); N <- 600L; t <- sort(sample.int(N, 400))
  y <- 2 * cos(2 * pi * 7 * (t - 1) / N) + 1.2 * cos(2 * pi * 31 * (t - 1) / N) +
    rnorm(400, sd = 0.4)
  cp <- clean_decompose(y, gls_prepare(t, N))
  all(diff(cp$var_explained_cumulative) > 0) && all(cp$bic_drop > 0) })

# --- query fingerprints: absence is not zero ---------------------------------
make_grid <- function(n_per_chr = 120L, chrs = c("1", "2")) {
  do.call(rbind, lapply(chrs, function(c) data.frame(
    gene_id = paste0("ENSG", c, "_", seq_len(n_per_chr)),
    entrez_id = paste0(c, sprintf("%04d", seq_len(n_per_chr))),
    chr = c, start = seq_len(n_per_chr) * 1000,
    grid_index = seq_len(n_per_chr), grid_N = n_per_chr,
    stringsAsFactors = FALSE)))
}
fake_ref <- function() list(grid = make_grid(),
                            params = list(k_max = 20L, features = "amplitude"))

check("a query keeps only the genes it actually contains", {
  g <- make_grid(); ref <- fake_ref()
  present <- g$gene_id[seq(1, nrow(g), by = 2)]     # half the grid
  qi <- query_grid_index(ref$grid, present)
  abs(qi$coverage - 0.5) < 0.01 &&
    all(vapply(qi$chrom_idx, function(ci) length(ci$t) == ci$N / 2, logical(1))) &&
    all(vapply(qi$chrom_idx, function(ci) ci$N == 120L, logical(1))) })

check("missing genes are absent, not zero-filled", {
  # The same measured values, once with the missing genes simply absent and once
  # zero-filled, must not give the same fingerprint -- if they did, the query
  # path would be reintroducing the very artefact the grid design removes.
  set.seed(61); g <- make_grid()
  ref <- list(grid = g, params = list(k_max = 40L, features = "amplitude"))
  idx <- seq(1, nrow(g), by = 2)
  vals <- round(exp(rnorm(length(idx), 5, 1)))
  fp_absent <- fingerprint_query(vals, g$gene_id[idx], ref)
  vals_full <- rep(0, nrow(g)); vals_full[idx] <- vals
  fp_zeros <- fingerprint_query(vals_full, g$gene_id, ref)
  shared <- intersect(names(fp_absent$vector), names(fp_zeros$vector))
  # Same measurements, different observed sets, therefore different spectra.
  abs(fp_absent$coverage - 0.5) < 0.01 && abs(fp_zeros$coverage - 1) < 0.01 &&
    length(shared) > 10 &&
    max(abs(fp_absent$vector[shared] - fp_zeros$vector[shared])) > 0.05 })

check("query grid N is the reference grid, not the query length", {
  g <- make_grid(); ref <- fake_ref()
  qi <- query_grid_index(ref$grid, g$gene_id[seq(1, nrow(g), by = 3)])
  all(vapply(qi$chrom_idx, function(ci) ci$N == 120L, logical(1))) })

check("entrez identifiers are detected", {
  g <- make_grid(); ref <- fake_ref()
  query_grid_index(ref$grid, g$entrez_id)$id_type == "entrez" })

# --- calibrated rejection ----------------------------------------------------
check("rejection threshold comes from correct held-out matches", {
  pred <- data.frame(
    truth = c(rep("A", 10), rep("B", 10)),
    predicted = c(rep("A", 8), "B", "B", rep("B", 8), "A", "A"),
    similarity = c(runif(8, 0.8, 0.9), 0.3, 0.3, runif(8, 0.8, 0.9), 0.3, 0.3),
    stringsAsFactors = FALSE)
  cal <- calibrate_rejection(pred)
  !is.null(cal) && cal$global_threshold > 0.7 && cal$separability_auc > 0.9 })

check("a low-similarity match is reported UNKNOWN", {
  cal <- list(global_threshold = 0.8, per_class = NULL, quantile = 0.05)
  hi <- apply_rejection(list(best = "A", similarity = 0.9), cal)
  lo <- apply_rejection(list(best = "A", similarity = 0.5), cal)
  identical(hi$decision, "A") && identical(lo$decision, "UNKNOWN") })

check("an uncalibrated reference says so instead of deciding", {
  identical(apply_rejection(list(best = "A", similarity = 0.9), NULL)$decision,
            "UNCALIBRATED") })

check("similarity uses shared features only, never zero padding", {
  set.seed(71)
  feats <- paste0("f", 1:40)
  cent <- rbind(A = rnorm(40), B = rnorm(40)); colnames(cent) <- feats
  model <- list(features = feats, centroids = cent, classes = c("A", "B"))
  q <- stats::setNames(rnorm(40), feats)
  half <- feats[1:20]
  # Scoring a query that only observed half the features must equal scoring the
  # half-length query directly -- i.e. the absent half contributes nothing,
  # rather than contributing a zero.
  a <- score_query(model, q[half], half)
  b <- score_query(model, q, half)
  isTRUE(all.equal(a$similarity, b$similarity)) })

check("padding with zeros would have changed the answer", {
  set.seed(72)
  feats <- paste0("f", 1:40)
  cent <- rbind(A = rnorm(40), B = rnorm(40)); colnames(cent) <- feats
  model <- list(features = feats, centroids = cent, classes = c("A", "B"))
  q <- stats::setNames(rnorm(40), feats); half <- feats[1:20]
  shared_only <- score_query(model, q[half], half)
  padded <- q; padded[feats[21:40]] <- 0
  with_padding <- score_query(model, padded, feats)
  max(abs(shared_only$similarity[match(with_padding$class, shared_only$class)] -
            with_padding$similarity)) > 0.01 })

check("a query with too few shared features is refused", {
  feats <- paste0("f", 1:40)
  cent <- rbind(A = rnorm(40)); colnames(cent) <- feats
  is.null(score_query(list(features = feats, centroids = cent),
                      stats::setNames(rnorm(2), feats[1:2]), feats[1:2])) })

# --- query hygiene -----------------------------------------------------------
check("duplicate identifiers are summed, not silently dropped", {
  r <- collapse_duplicate_ids(c(10, 5, 7), c("g1", "g1", "g2"))
  r$collapsed == 2L && r$values[r$ids == "g1"] == 15 })

check("duplicates cannot be summed for normalised input", {
  inherits(tryCatch(collapse_duplicate_ids(c(1, 2), c("g1", "g1"), unit = "tpm"),
                    error = function(e) e), "error") })

check("negative values are refused", {
  inherits(tryCatch(collapse_duplicate_ids(c(1, -2), c("g1", "g2")),
                    error = function(e) e), "error") })

check("the input unit changes the transform", {
  v <- c(10, 100, 1000)
  !isTRUE(all.equal(query_signal(v, "counts"), query_signal(v, "tpm"))) &&
    identical(query_signal(v, "logged"), v) &&
    inherits(tryCatch(query_signal(v, "rpkm"), error = function(e) e), "error") })

# --- consensus spectrum ------------------------------------------------------
check("PLV is 1 for aligned phases and near 0 for scattered ones", {
  set.seed(81)
  a <- phase_locking(rep(0.7, 20))[["plv"]]
  b <- phase_locking(runif(200, -pi, pi))[["plv"]]
  abs(a - 1) < 1e-12 && b < 0.2 })

check("the Rayleigh p-value is calibrated under uniform phases", {
  set.seed(82)
  p <- vapply(1:300, function(i) phase_locking(runif(12, -pi, pi))[["rayleigh_p"]],
              numeric(1))
  abs(mean(p <= 0.05) - 0.05) < 0.04 })

check("consensus finds what the mean profile cancels", {
  # Every sample carries the same frequency; the phase differs between samples.
  # The mean profile's spectrum loses it to cancellation; the consensus sees a
  # strong, prevalent component with low PLV -- present but not phase-locked.
  set.seed(83); N <- 200L; t <- 1:N; k0 <- 9L
  phases <- runif(12, -pi, pi)
  sig <- vapply(phases, function(ph)
    2 * cos(2 * pi * k0 * (t - 1) / N + ph) + rnorm(N, sd = 0.3), numeric(N))
  mean_spec <- run_fft(rowMeans(sig))
  amp_mean <- mean_spec$amplitude[mean_spec$k == k0]
  per_sample <- do.call(rbind, lapply(seq_along(phases), function(i) {
    r <- run_fft(sig[, i]); r$sample <- paste0("S", i); r$chr <- "1"; r
  }))
  cs <- consensus_spectrum(per_sample, n_boot = 50L)
  row <- cs[cs$chr == "1" & cs$k == k0, ]
  amp_mean < 0.6 && nrow(row) == 1 && row$prevalence == 1 && row$plv < 0.6 &&
    row$median_amplitude > 1.5 })

check("consensus scores an aligned component above a scattered one", {
  set.seed(84); N <- 200L; t <- 1:N
  mk <- function(k, ph) do.call(rbind, lapply(1:12, function(i) {
    r <- run_fft(2 * cos(2 * pi * k * (t - 1) / N + ph[i]) + rnorm(N, sd = 0.3))
    r$sample <- paste0("S", i); r$chr <- "1"; r }))
  aligned <- mk(9L, rep(0.4, 12))
  scattered <- mk(9L, runif(12, -pi, pi))
  ca <- consensus_spectrum(aligned, n_boot = 50L)
  cb <- consensus_spectrum(scattered, n_boot = 50L)
  ca$consensus_score[ca$k == 9] > cb$consensus_score[cb$k == 9] })

check("bootstrap intervals bracket the point estimate", {
  set.seed(85); N <- 200L; t <- 1:N
  per_sample <- do.call(rbind, lapply(1:10, function(i) {
    r <- run_fft(2 * cos(2 * pi * 7 * (t - 1) / N + 0.3) + rnorm(N, sd = 0.5))
    r$sample <- paste0("S", i); r$chr <- "1"; r }))
  cs <- consensus_spectrum(per_sample, n_boot = 100L)
  row <- cs[cs$k == 7, ]
  row$consensus_score >= row$consensus_score_ci_lower &&
    row$consensus_score <= row$consensus_score_ci_upper })

check("prevalence is judged within a chromosome, not across chromosomes", {
  # chr2 has systematically lower normalised power than chr1. Pooling the two
  # would let chr1 set the threshold and no chr2 frequency would ever count as
  # prevalent, even the one that dominates chr2's own spectrum.
  pn <- c(rep(0.020, 39), 0.300,     # chr1: high power throughout
          rep(0.0001, 39), 0.005)     # chr2: much weaker, one clear peak of its own
  chr <- c(rep("1", 40), rep("2", 40))
  smp <- rep("S1", 80)
  per_chr <- prevalence_from_rank(pn, smp, chr, 0.95)
  pooled <- prevalence_from_rank(pn, smp, rep("1", 80), 0.95)
  per_chr[40] && per_chr[80] && !pooled[80] })

# --- coverage-band calibration -----------------------------------------------
check("masking removes genes, not frequencies", {
  # A gene mask leaves the GLS able to estimate nearly every frequency, so
  # feature coverage stays high while gene coverage falls. A feature mask
  # cannot reproduce that, which is why it calibrates the wrong quantity.
  ci <- list("1" = list(rows = 1:200, t = 1:200, N = 200L, coverage = 1))
  masked <- mask_grid_genes(ci, 0.5, "random", seed = 1L)
  gene_cov <- length(masked[["1"]]$t) / length(ci[["1"]]$t)
  set.seed(9); y <- rnorm(200)
  fp_full <- fingerprint_masked(y, ci, 40L, "amplitude")
  fp_part <- fingerprint_masked(y[masked[["1"]]$t], masked, 40L, "amplitude")
  feat_cov <- length(intersect(names(fp_part), names(fp_full))) / length(fp_full)
  abs(gene_cov - 0.5) < 0.02 && feat_cov > 0.9 })

check("a chromosome mask drops whole chromosomes", {
  ci <- stats::setNames(lapply(1:4, function(i)
    list(rows = 1:100, t = 1:100, N = 100L, coverage = 1)), as.character(1:4))
  masked <- mask_grid_genes(ci, 0.5, "chromosome", seed = 3L)
  length(masked) < 4 && all(vapply(masked, function(m) length(m$t) == 100L, logical(1))) })

check("retained_block keeps one contiguous run", {
  ci <- list("1" = list(rows = 1:200, t = 1:200, N = 200L, coverage = 1))
  t <- mask_grid_genes(ci, 0.5, "retained_block", seed = 4L)[["1"]]$t
  all(diff(t) == 1L) })

check("missing_blocks removes several disjoint intervals", {
  ci <- list("1" = list(rows = 1:400, t = 1:400, N = 400L, coverage = 1))
  t <- mask_grid_genes(ci, 0.6, "missing_blocks", seed = 5L, n_blocks = 4L)[["1"]]$t
  gaps <- sum(diff(t) > 1)
  gaps >= 2 && length(t) < 400 })

check("expression_dropout removes the least expressed first", {
  ci <- list("1" = list(rows = 1:200, t = 1:200, N = 200L, coverage = 1))
  v <- seq_len(200)                       # gene i has expression i
  kept <- mask_grid_genes(ci, 0.5, "expression_dropout", seed = 6L,
                          values = v)[["1"]]$t
  mean(v[kept]) > mean(v) })

check("expression_dropout falls back when no values are given", {
  ci <- list("1" = list(rows = 1:200, t = 1:200, N = 200L, coverage = 1))
  length(mask_grid_genes(ci, 0.5, "expression_dropout", seed = 6L)[["1"]]$t) == 100L })

check("a chromosome left too short is dropped", {
  ci <- list("1" = list(rows = 1:20, t = 1:20, N = 20L, coverage = 1))
  length(mask_grid_genes(ci, 0.2, "random", seed = 5L)) == 0L })

check("different mask seeds give different masks, the same seed repeats", {
  ci <- stats::setNames(lapply(1:8, function(i)
    list(rows = 1:100, t = 1:100, N = 100L, coverage = 1)), as.character(1:8))
  masks <- lapply(1:6, function(s) sort(names(mask_grid_genes(ci, 0.5, "chromosome", seed = s))))
  repeated <- sort(names(mask_grid_genes(ci, 0.5, "chromosome", seed = 1L)))
  length(unique(masks)) > 1 && identical(masks[[1]], repeated) })

check("the applied threshold follows the declared policy", {
  bp <- data.frame(band = "50-75%", held_out = "A", mode = "random",
                   mask = rep(1:4, each = 10),
                   truth = "X", predicted = "X",
                   similarity = c(runif(10, .5, .6), runif(10, .7, .8),
                                  runif(10, .8, .9), runif(10, .9, .95)),
                   stringsAsFactors = FALSE)
  pooled <- summarise_bands(bp, policy = "pooled")
  cons <- summarise_bands(bp, policy = "conservative")
  cons$threshold_applied > pooled$threshold_applied &&
    cons$expected_rejection_of_members >= pooled$expected_rejection_of_members })

check("both rejection rates are reported whichever policy is applied", {
  bp <- data.frame(band = "50-75%", held_out = "A", mode = "random",
                   mask = rep(1:4, each = 10), truth = "X", predicted = "X",
                   similarity = runif(40, .5, .95), stringsAsFactors = FALSE)
  b <- summarise_bands(bp, policy = "pooled")
  all(c("unknown_rate_at_threshold", "unknown_rate_conservative",
        "threshold_applied", "policy") %in% colnames(b)) })

check("confirmation requires beating the permuted null, not just zero", {
  cs <- data.frame(chr = "1", N = 200L, k = 6L, freq = 0.03, period = 33,
                   n_samples_valid = 40L, prevalence = 1,
                   plv = 1, plv_rayleigh_p = 1e-12, plv_rayleigh_q = 1e-10,
                   consensus_score = 0.5, consensus_score_ci_lower = 0.3,
                   stringsAsFactors = FALSE)
  beats <- consensus_signature(cs, null_score = 0.1)
  loses <- consensus_signature(cs, null_score = 0.9)
  none  <- consensus_signature(cs, null_score = NA_real_)
  identical(beats$signature_class, "confirmed") &&
    identical(loses$signature_class, "exploratory") &&
    identical(none$signature_class, "exploratory") })

check("signatures are labelled confirmed or exploratory", {
  # Two samples: exp(-2) * n_freq exceeds the q threshold, so phase alignment
  # is not testable and the component can only be exploratory.
  cs <- data.frame(chr = "1", N = 200L, k = 6L, freq = 0.03, period = 33,
                   n_samples_valid = 2L, prevalence = 1,
                   plv = 1, plv_rayleigh_p = 0.0067, plv_rayleigh_q = 1,
                   consensus_score = 0.2, consensus_score_ci_lower = 0.1,
                   stringsAsFactors = FALSE)
  sig <- consensus_signature(cs)
  identical(sig$signature_class, "exploratory") && !sig$phase_alignment_testable })

check("coverage bands are assigned correctly", {
  identical(vapply(c(0.95, 0.8, 0.6, 0.3), coverage_band, character(1)),
            c("90-100%", "75-90%", "50-75%", "<50%")) })

check("a partial query does not reuse the full-coverage threshold", {
  cal <- list(global_threshold = 0.95, per_class = NULL, quantile = 0.05,
              bands = NULL)
  full <- apply_rejection(list(best = "A", similarity = 0.96), cal, coverage = 0.95)
  part <- apply_rejection(list(best = "A", similarity = 0.96), cal, coverage = 0.6)
  identical(full$decision, "A") && identical(part$decision, "UNCALIBRATED_COVERAGE") })

check("a band threshold is used when one exists", {
  cal <- list(global_threshold = 0.95, per_class = NULL, quantile = 0.05,
              bands = data.frame(band = "50-75%", threshold = 0.4,
                                 stringsAsFactors = FALSE))
  r <- apply_rejection(list(best = "A", similarity = 0.5), cal, coverage = 0.6)
  identical(r$decision, "A") && startsWith(r$threshold_source, "band") })

check("below 50% coverage nothing is classified", {
  cal <- list(global_threshold = 0.2, per_class = NULL, quantile = 0.05,
              bands = data.frame(band = "<50%", threshold = 0.1,
                                 stringsAsFactors = FALSE))
  identical(apply_rejection(list(best = "A", similarity = 0.9), cal,
                            coverage = 0.3)$decision, "LOW_COVERAGE") })

# --- duplicates that are all NA, and units -----------------------------------
check("duplicate rows that are all NA stay unmeasured, not zero", {
  r <- collapse_duplicate_ids(c(NA, NA, 4), c("g1", "g1", "g2"))
  is.na(r$values[r$ids == "g1"]) && r$values[r$ids == "g2"] == 4 && r$all_na == 1L })

check("a partly-NA duplicate group sums the finite values", {
  r <- collapse_duplicate_ids(c(NA, 3, 4), c("g1", "g1", "g2"))
  r$values[r$ids == "g1"] == 3 })

check("a TPM query is refused against a CPM reference", {
  ref <- list(params = list(expression_unit = "asinh(CPM)"))
  inherits(tryCatch(assert_unit_compatible(ref, "tpm"), error = function(e) e),
           "error") &&
    !inherits(tryCatch(assert_unit_compatible(ref, "cpm"), error = function(e) e),
              "error") })

check("counts are accepted against a CPM reference", {
  ref <- list(params = list(expression_unit = "asinh(CPM)"))
  !inherits(tryCatch(assert_unit_compatible(ref, "counts"), error = function(e) e),
            "error") })

check("a TPM reference accepts a TPM query only", {
  ref <- list(params = list(expression_unit = "asinh(TPM)"))
  !inherits(tryCatch(assert_unit_compatible(ref, "tpm"), error = function(e) e), "error") &&
    inherits(tryCatch(assert_unit_compatible(ref, "cpm"), error = function(e) e), "error") })

# --- Wilson ------------------------------------------------------------------
check("Wilson interval brackets the point estimate", {
  ci <- wilson_ci(9, 10); ci[1] < 90 && ci[2] > 90 && ci[1] >= 0 && ci[2] <= 100 })
check("Wilson handles n = 0", all(is.na(wilson_ci(0, 0))))

cat("\n", if (failures == 0L) "All tests passed." else paste(failures, "test(s) failed."), "\n")
quit(status = if (failures == 0L) 0L else 1L)
