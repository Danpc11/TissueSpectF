#!/usr/bin/env Rscript
# Numerical tests for the spectral core. Run: Rscript tests/test_spectrum.R
source("R/utils_io.R"); source("R/config.R"); source("R/labels.R")
source("R/grid.R"); source("R/period_floor.R"); source("R/ingest.R"); source("R/spectrum.R"); source("R/maxt.R"); source("R/stability.R")
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

check("the pooled null reaches a floor BH can use", {
  # A per-frequency p against its own null floors at 1/(B+1); BH over ~10,000
  # genome-wide frequencies would then need B >= n_freq/q, two hundred thousand
  # draws. Pooling the standardised nulls across frequencies floors at 1/(B*K+1).
  set.seed(41); N <- 400L; tt <- 1:N
  m <- permutation_gls_test(rnorm(N), gls_prepare(tt, N), B = 100L, seed = 1L)
  # The p-values go below the unpooled floor, so BH is no longer capped by the
  # number of draws. On this input BH still finds nothing, which is correct --
  # it is pure noise; what is being checked is the floor, not a detection.
  min(m$p_pointwise) < 1 / (100 + 1) && min(p.adjust(m$p_pointwise, "BH")) < 1 })

check("FDR finds components that do not dominate their chromosome", {
  # Two real components, neither the strongest thing in a noisy spectrum. The
  # family-wise rule asks whether a frequency beats the maximum of a permuted
  # spectrum, which is close to asking whether it dominates the chromosome; the
  # pointwise rule asks whether it beats its own null.
  set.seed(42); N <- 400L; tt <- 1:N
  y <- 0.6 * cos(2 * pi * 7 * (tt - 1) / N) +
       0.5 * cos(2 * pi * 23 * (tt - 1) / N + 1) + rnorm(N, sd = 1)
  m <- permutation_gls_test(y, gls_prepare(tt, N), B = 200L, seed = 2L)
  q <- p.adjust(m$p_pointwise, "BH")
  all(c(7L, 23L) %in% m$k[q <= 0.05]) && sum(q <= 0.05) <= 6 })

check("the pooled null is calibrated on noise", {
  set.seed(43); N <- 400L; tt <- 1:N
  m <- permutation_gls_test(rnorm(N), gls_prepare(tt, N), B = 200L, seed = 3L)
  mean(m$p_pointwise <= 0.05) < 0.12 })

check("the no-maxT path honours the declared criterion", {
  # This path used to select on q_condition whatever the caller asked for, so
  # --criterion condition_fdr silently returned the family-wise result.
  d <- file.path(tempdir(), "crit", "condition")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  write.table(data.frame(chr = "1", N = 200L, k = 1:3, freq = (1:3) / 200,
                         period = 200 / (1:3), power = 1,
                         window_power = 0, window_rank = 1,
                         p_condition = c(0.001, 0.20, 0.30),
                         q_condition = c(0.02, 0.90, 0.99),
                         q_condition_fdr = c(0.001, 0.01, 0.60)),
              file.path(d, "condition_significance_average_F0.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  p <- list(base = dirname(d))
  fwer <- stability_from_condition(p, "F0", "average", "condition")
  fdr <- stability_from_condition(p, "F0", "average", "condition_fdr")
  sum(fwer$is_stable) == 1 && sum(fdr$is_stable) == 2 })

check("asking for FDR against a file that predates it aborts", {
  d <- file.path(tempdir(), "old", "condition")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  write.table(data.frame(chr = "1", N = 200L, k = 1L, freq = 0.005, period = 200,
                         power = 1, window_power = 0, window_rank = 1,
                         p_condition = 0.001, q_condition = 0.02),
              file.path(d, "condition_significance_average_F0.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  inherits(tryCatch(stability_from_condition(list(base = dirname(d)), "F0",
                                             "average", "condition_fdr"),
                    error = function(e) e), "error") })

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

check("missing_blocks removes several intervals and hits its target", {
  ci <- list("1" = list(rows = 1:400, t = 1:400, N = 400L, coverage = 1))
  kept <- vapply(1:20, function(s)
    length(mask_grid_genes(ci, 0.6, "missing_blocks", seed = s, n_blocks = 4L)[["1"]]$t),
    integer(1))
  t <- mask_grid_genes(ci, 0.6, "missing_blocks", seed = 5L, n_blocks = 4L)[["1"]]$t
  # Overlapping draws used to remove fewer genes than asked; the retained count
  # must now match the target rather than drifting above it.
  sum(diff(t) > 1) >= 2 && all(abs(kept - 240) <= 2) })

check("the two consensus scores are both reported and differ in meaning", {
  set.seed(91); N <- 200L; tt <- 1:N
  per_sample <- do.call(rbind, lapply(1:8, function(i) {
    r <- run_fft(2 * cos(2 * pi * 7 * (tt - 1) / N + 0.3) + rnorm(N, sd = 0.5))
    r$sample <- paste0("S", i); r$chr <- "1"; r }))
  no_maxt <- consensus_spectrum(per_sample, n_boot = 20L)
  # Without maxT the maxT-based score is NA and the rank score is populated;
  # the permuted comparison must use the rank one on both sides.
  all(is.na(no_maxt$consensus_score_maxt)) &&
    all(is.finite(no_maxt$consensus_score_rank)) &&
    isTRUE(all.equal(no_maxt$consensus_score, no_maxt$consensus_score_rank)) })

check("confirmation uses a family-wise p with a reachable floor", {
  # 50 draws over 297 frequencies: BH on the pointwise p cannot go below 5.8,
  # so nothing could ever be confirmed that way. The family-wise p against the
  # distribution of the null's maximum has a floor of 1/51 and can.
  cs <- data.frame(chr = "1", N = 200L, k = seq_len(297),
                   consensus_score_rank = c(0.95, runif(296, 0, 0.3)),
                   consensus_score_ci_lower = 0.9, stringsAsFactors = FALSE)
  nd <- list(per_key = matrix(runif(297 * 50, 0, 0.3), nrow = 297,
                              dimnames = list(paste("1", 200L, seq_len(297), sep = "|"), NULL)),
             global = 0.4, global_max_draws = runif(50, 0.2, 0.4),
             n_null = 50L, n_samples = 10L)
  out <- null_component_pvalues(cs, nd)
  min(out$q_null, na.rm = TRUE) > 0.05 &&      # the pointwise route cannot confirm
    out$p_null_fwer[1] <= 0.05 &&              # the family-wise route can
    out$p_null_fwer[1] >= 1 / 51 })

check("a blocked null draws whole blocks", {
  set.seed(95); N <- 120L; tt <- 1:N
  per_sample <- do.call(rbind, lapply(1:12, function(i) {
    r <- run_fft(cos(2 * pi * 5 * (tt - 1) / N) + rnorm(N, sd = 0.5))
    r$sample <- paste0("S", i); r$chr <- "1"; r }))
  blocks <- stats::setNames(rep(paste0("subj", 1:6), each = 2), paste0("S", 1:12))
  nd <- null_consensus_distribution(per_sample, n_samples = 4L, n_null = 12L,
                                    blocks = blocks)
  # Blocks of two, four samples per draw: a draw can never split a subject.
  !is.null(nd) && nd$n_null > 0 })

check("the per-component null gives a p-value and a BH q-value", {
  cs <- data.frame(chr = "1", N = 200L, k = c(6L, 7L), freq = c(.03, .035),
                   period = c(33, 28), n_samples_valid = 10L, prevalence = 1,
                   plv = 1, plv_rayleigh_p = 1e-6, plv_rayleigh_q = 1e-5,
                   consensus_score_rank = c(0.9, 0.1),
                   consensus_score_ci_lower = c(0.8, 0.05),
                   stringsAsFactors = FALSE)
  nd <- list(per_key = matrix(runif(2 * 50, 0, 0.3), nrow = 2,
                              dimnames = list(c("1|200|6", "1|200|7"), NULL)),
             global = 0.95, n_null = 50L, n_samples = 10L)
  out <- null_component_pvalues(cs, nd)
  out$p_null[1] < 0.05 && out$p_null[2] > 0.05 &&
    all(out$q_null >= out$p_null) && !out$beats_global_null[1] })

check("the policy reaches the band+class level too", {
  bp <- data.frame(band = "50-75%", held_out = "A", mode = "random",
                   mask = rep(1:4, each = 10), truth = "X", predicted = "X",
                   similarity = c(runif(10, .5, .6), runif(10, .7, .8),
                                  runif(10, .8, .9), runif(10, .9, .95)),
                   stringsAsFactors = FALSE)
  po <- summarise_bands_by_class(bp, policy = "pooled", min_correct = 8L)
  co <- summarise_bands_by_class(bp, policy = "conservative", min_correct = 8L)
  co$threshold_applied > po$threshold_applied &&
    identical(po$threshold_applied, po$threshold_pooled) })

check("band and class thresholds beat the band-only threshold", {
  cal <- list(global_threshold = 0.9, per_class = NULL, quantile = 0.05,
              bands = data.frame(band = "50-75%", threshold = 0.4,
                                 threshold_applied = 0.4, policy = "pooled",
                                 stringsAsFactors = FALSE),
              bands_by_class = data.frame(band = "50-75%", class = "A",
                                          threshold = 0.7, threshold_applied = 0.7,
                                          policy = "pooled", stringsAsFactors = FALSE))
  a <- apply_rejection(list(best = "A", similarity = 0.5), cal, coverage = 0.6)
  b <- apply_rejection(list(best = "B", similarity = 0.5), cal, coverage = 0.6)
  identical(a$decision, "UNKNOWN") && startsWith(a$threshold_source, "band+class") &&
    identical(b$decision, "B") && startsWith(b$threshold_source, "band:") })

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
  base <- data.frame(chr = "1", N = 200L, k = 6L, freq = 0.03, period = 33,
                     n_samples_valid = 40L, prevalence = 1,
                     plv = 1, plv_rayleigh_p = 1e-12, plv_rayleigh_q = 1e-10,
                     consensus_score = 0.5, consensus_score_rank = 0.5,
                     consensus_score_ci_lower = 0.3,
                     stringsAsFactors = FALSE)
  beats <- consensus_signature(within(base, { p_null_fwer <- 0.02 }))
  loses <- consensus_signature(within(base, { p_null_fwer <- 0.4 }))
  none  <- consensus_signature(within(base, { p_null_fwer <- NA_real_ }))
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

# --- unmeasured positions are dropped, never zero-filled ---------------------
#
# The rule the whole grid design rests on: an unmeasured gene keeps its slot on
# the axis and leaves the observed set. It never becomes a zero at a fixed
# position. A zero-fill here is not a rounding difference -- it injects a
# deterministic pattern whose own spectrum can manufacture a peak -- so the
# invariant is asserted directly rather than trusted to code comments.

check("gls_observed drops an unmeasured position and rebuilds the window", {
  N <- 128L; t <- as.integer(seq(1, N, by = 2))
  y <- cos(2 * pi * 7 * (t - 1) / N)
  y[5] <- NA_real_
  fit <- gls_observed(y, t, N)
  fit$n_dropped == 1L &&
    length(fit$y) == length(t) - 1L &&
    fit$terms$n == length(t) - 1L &&
    all(is.finite(fit$y)) })

check("a complete signal reuses the shared window terms untouched", {
  N <- 128L; t <- as.integer(seq(1, N, by = 2))
  shared <- gls_prepare(t, N)
  fit <- gls_observed(rnorm(length(t)), t, N, terms = shared)
  fit$n_dropped == 0L && identical(fit$terms, shared) })

check("gls_observed refuses to fit when too little survives", {
  N <- 128L; t <- as.integer(seq(1, N, by = 2))
  y <- rep(NA_real_, length(t)); y[1:4] <- 1:4
  is.null(gls_observed(y, t, N, min_observed = 8L)) })

check("gls_spectrum refuses a non-finite value at an observed position", {
  N <- 64L; t <- seq_len(N)
  y <- rnorm(N); y[3] <- NA_real_
  inherits(tryCatch(gls_spectrum(y, gls_prepare(t, N)), error = function(e) e),
           "error") })

check("dropping a hole is not the same as zero-filling it", {
  # A hole filled with zero pulls the fit towards a spike at a fixed position;
  # dropping it leaves the underlying sinusoid's amplitude essentially intact.
  # If these two ever agree, the distinction this pipeline is built on has
  # stopped being enforced somewhere.
  set.seed(4)
  N <- 256L; t <- seq_len(N); k <- 9L
  y <- 3 * cos(2 * pi * k * (t - 1) / N)
  holes <- seq(10, 120, by = 10)

  y_na <- y; y_na[holes] <- NA_real_
  dropped <- gls_observed(y_na, t, N)
  amp_dropped <- gls_spectrum(dropped$y, dropped$terms)
  amp_dropped <- amp_dropped$amplitude[amp_dropped$k == k]

  y_zero <- y; y_zero[holes] <- 0
  amp_zero <- gls_spectrum(y_zero, gls_prepare(t, N))
  amp_zero <- amp_zero$amplitude[amp_zero$k == k]

  abs(amp_dropped - 3) < 0.05 && abs(amp_zero - 3) > 0.1 })

check("condition_signals keeps a hole as NA instead of turning it into zero", {
  expr <- matrix(rnorm(30), nrow = 10,
                 dimnames = list(paste0("g", 1:10), paste0("S", 1:3)))
  expr[4, 2] <- NA_real_
  ds <- list(samples = data.frame(sample_id = paste0("S", 1:3),
                                  condition = rep("F1", 3),
                                  stringsAsFactors = FALSE),
             expression = expr)
  sig <- suppressWarnings(condition_signals(ds, "F1"))
  # The hole survives in the matrix, and the summary signal averages the two
  # samples that do carry a value rather than averaging in a zero.
  is.na(sig$matrix[4, 2]) &&
    isTRUE(all.equal(sig$avg_signal[[4]], mean(expr[4, c(1, 3)]))) })

check("a grid position measured in no sample stays NA in the summary signal", {
  expr <- matrix(rnorm(30), nrow = 10,
                 dimnames = list(paste0("g", 1:10), paste0("S", 1:3)))
  expr[7, ] <- NA_real_
  ds <- list(samples = data.frame(sample_id = paste0("S", 1:3),
                                  condition = rep("F1", 3),
                                  stringsAsFactors = FALSE),
             expression = expr)
  sig <- suppressWarnings(condition_signals(ds, "F1"))
  is.na(sig$avg_signal[[7]]) && is.na(sig$median_signal[[7]]) })

# --- period floor, applied before the null ---------------------------------

check("the technical floor is the Nyquist of the sampling that happened", {
  # coverage 0.5 means observed genes sit ~2 apart, so 2/0.5 = 4 is the
  # shortest resolvable period; +2 margin makes 6.
  long <- data.frame(period = c(3, 5, 7, 40), coverage = 0.5)
  kept <- apply_period_floor(long, min_period = "auto", period_margin = 2,
                             quiet = TRUE)
  identical(kept$period, c(7, 40)) })

check("the biological floor applies even with the technical floor off", {
  long <- data.frame(period = c(2, 4, 12, 90), coverage = 1)
  kept <- apply_period_floor(long, min_period = "off",
                             min_period_biological = 10, quiet = TRUE)
  identical(kept$period, c(12, 90)) })

check("the applied floor is the larger of the two", {
  # technical = 2/0.1 + 2 = 22, biological = 10 -> 22 wins.
  long <- data.frame(period = c(15, 25), coverage = 0.1)
  kept <- apply_period_floor(long, min_period = "auto", period_margin = 2,
                             min_period_biological = 10, quiet = TRUE)
  identical(kept$period, 25) })

check("a reachable BH q is never above 1", {
  # m/(B+1) is the uncapped adjusted p at rank 1; a q-value is capped at 1, and
  # reporting 10 reads like a threshold narrowly missed rather than one that
  # cannot be reached at all.
  bh_rank1_diagnostic(10030, 999) == 1 &&
    abs(bh_rank1_diagnostic(500, 999) - 0.5) < 1e-9 })

check("draws_for_bh inverts the rank-1 figure exactly", {
  # B is the count of draws, so the family of B+1 values includes the observed
  # one: m/(B+1) <= target at B, and not at B-1. Off by one here would either
  # promise a threshold that is not reached or demand a draw that is not needed.
  m <- 1000L
  B <- draws_for_bh(m, 0.05)
  B == 19999L &&
    bh_rank1_diagnostic(m, B) <= 0.05 &&
    bh_rank1_diagnostic(m, B - 1L) > 0.05 })

check("shrinking the family is what makes the pointwise route reachable", {
  # The whole reason the floor goes before the null: at 999 draws a family of
  # 10030 cannot reach 0.05, and a family of 500 can.
  bh_rank1_diagnostic(10030, 999) > 0.05 &&
    bh_rank1_diagnostic(500, 999) <= 0.5 })

check("the rank-1 BH figure is a conservative diagnostic, not a bound", {
  # This was documented and reported as "the smallest reachable q", which is
  # wrong: BH takes min over j>=i of m*p_j/j, so ties at the permutation floor
  # lower the achievable value by up to a factor of m. Reporting the rank-1
  # value as a bound concluded "nothing can pass" when 0.001 was reachable.
  m <- 10030L; B <- 999L
  rank1 <- bh_rank1_diagnostic(m, B)
  all_tied <- bh_achievable_q(m, B, m)
  some_tied <- bh_achievable_q(m, B, 100L)
  rank1 == 1 &&
    abs(all_tied - 1 / (B + 1)) < 1e-12 &&
    some_tied < rank1 && all_tied < some_tied })

check("the achievable q equals the rank-1 figure when nothing ties", {
  identical(bh_achievable_q(500L, 999L, 1L), bh_rank1_diagnostic(500L, 999L)) })

check("a q-value is never reported above 1", {
  all(vapply(list(c(10030, 999), c(1e6, 10), c(5, 1)),
             function(x) bh_rank1_diagnostic(x[1], x[2]) <= 1, logical(1))) })

check("the frequency-family digest is order-independent", {
  identical(digest_keys(c("1|100|3", "2|50|7")),
            digest_keys(c("2|50|7", "1|100|3"))) })

check("the digest changes when the retained family changes", {
  # The null cache was keyed on sample count alone, so a run with a period
  # floor reused the null computed without one -- a different family, a
  # different global maximum, and no error to show for it.
  base <- digest_keys(c("21|500|6", "21|500|7", "1|900|3"))
  base != digest_keys(c("21|500|6", "1|900|3")) &&
    base != digest_keys(c("21|500|6", "21|500|8", "1|900|3")) })

check("the null cache key records the period configuration", {
  src <- paste(readLines("R/stages.R", warn = FALSE), collapse = " ")
  all(vapply(c("min_period", "period_margin", "margin_mode",
               "min_period_biological", "retained_digest"),
             function(k) grepl(k, src, fixed = TRUE), logical(1))) })

# --- peak reconstruction needs one shared axis -------------------------------
#
# A frequency is identified by (chr, N, k). Reconstruction lays that sinusoid on
# a dataset's axis of length ci$N and reads off which genes sit at the crests.
# If the peak's N and the axis N disagree the sinusoid is placed at the wrong
# frequency, every crest lands on the wrong gene, and the table is still written
# -- carrying `N` and `grid_N` as two different numbers in adjacent columns.

check("reconstruction refuses a peak measured on a different grid", {
  genes <- data.frame(gene_id = paste0("g", 1:20),
                      gene_name = paste0("G", 1:20),
                      stringsAsFactors = FALSE)
  chrom_idx <- list("21" = list(rows = 1:20, t = 1:20, N = 100L,
                                coverage = 0.2))
  # N = 120 on the peak against N = 100 on the axis.
  peaks <- data.frame(chr = "21", N = 120L, k = 6L, freq = 6 / 120,
                      period = 20, phase = 0.5, amplitude_mean = 1,
                      stringsAsFactors = FALSE)
  out <- tempfile(); dir.create(out)
  err <- tryCatch(
    write_peak_gene_tables(peaks, genes, chrom_idx, out, "average", "F2"),
    error = function(e) e)
  unlink(out, recursive = TRUE)
  inherits(err, "error") &&
    grepl("wrong frequency", conditionMessage(err), fixed = TRUE) })

check("reconstruction proceeds when the axes agree", {
  genes <- data.frame(gene_id = paste0("g", 1:20),
                      gene_name = paste0("G", 1:20),
                      stringsAsFactors = FALSE)
  chrom_idx <- list("21" = list(rows = 1:20, t = 1:20, N = 100L,
                                coverage = 0.2))
  peaks <- data.frame(chr = "21", N = 100L, k = 6L, freq = 6 / 100,
                      period = 100 / 6, phase = 0.5, amplitude_mean = 1,
                      stringsAsFactors = FALSE)
  out <- tempfile(); dir.create(out)
  n <- write_peak_gene_tables(peaks, genes, chrom_idx, out, "average", "F2")
  files <- list.files(out, recursive = TRUE)
  unlink(out, recursive = TRUE)
  n == 1L && length(files) == 1L })

check("the reconstructed sinusoid uses the reference grid length, not the observed count", {
  # 20 observed genes on a grid of 100. The sinusoid must complete k cycles over
  # 100 positions and be sampled at the observed ones -- not complete k cycles
  # over the 20 genes that happen to be present.
  N <- 100L; k <- 4L
  t_obs <- as.integer(seq(1, 100, by = 5))
  full <- reconstruct_signal(
    data.frame(freq = k / N, amplitude = 1, phase = 0), N)
  expected <- cos(2 * pi * k * (t_obs - 1) / N)
  max(abs(full[t_obs] - expected)) < 1e-9 })

check("peaks and compare check grid compatibility", {
  # assert_compatible_grids() was called only by stage_reference, so the two
  # stages whose purpose is lining cohorts up against each other never
  # confirmed they shared an axis.
  src <- paste(readLines("R/stages.R", warn = FALSE), collapse = "\n")
  peaks_body <- sub(".*stage_peaks <- function.*?\\{", "", src)
  compare_body <- sub(".*stage_compare <- function.*?\\{", "", src)
  grepl("assert_stage_grids", substr(peaks_body, 1, 400), fixed = TRUE) &&
    grepl("assert_stage_grids", substr(compare_body, 1, 600), fixed = TRUE) })

check("the consensus table gets coverage joined from the per-sample spectra", {
  # consensus_spectrum() builds its table from the (chr, N, k) key alone and has
  # no coverage column, so `--min-period auto` -- whose technical floor is
  # 2/coverage -- aborted the whole stage. The join happens in stage_consensus
  # because the coverage that matters is the one the consensus was built from:
  # gls_observed() drops unmeasured positions per signal, so per-signal and
  # nominal chromosome coverage genuinely differ now.
  sp <- data.frame(chr = c("21", "21", "21", "21"),
                   N = c(100L, 100L, 100L, 100L),
                   k = c(6L, 6L, 7L, 7L),
                   coverage = c(0.40, 0.60, 0.50, 0.50),
                   stringsAsFactors = FALSE)
  cs <- data.frame(chr = "21", N = c(100L, 100L), k = c(6L, 7L),
                   stringsAsFactors = FALSE)
  sp_key <- paste(sp$chr, sp$N, sp$k, sep = "|")
  cov_by_key <- tapply(sp$coverage, sp_key, stats::median, na.rm = TRUE)
  # as.numeric(), not unname(): tapply returns a 1-d array and unname() leaves
  # the dim attribute on, which would make this an array column.
  cs$coverage <- as.numeric(cov_by_key[paste(cs$chr, cs$N, cs$k, sep = "|")])
  # Median across contributing samples, not the worst or the best.
  identical(cs$coverage, c(0.5, 0.5)) && is.null(dim(cs$coverage)) })

check("auto floor fails with an actionable message, not a bare abort", {
  bad <- data.frame(period = c(3, 40))
  msg <- tryCatch(apply_period_floor(bad, label = "F2", min_period = "auto"),
                  error = function(e) conditionMessage(e))
  is.character(msg) &&
    grepl("--min-period <genes>", msg, fixed = TRUE) &&
    grepl("spectra", msg, fixed = TRUE) })

# --- the pooled pointwise null cannot be conjured from too few draws ---------

check("too few draws yields NA, not a confident p-value", {
  # With B = 1 the per-frequency sd is NA, the old guard replaced it with 1, the
  # pooled reference collapsed to zeros, and every frequency above its single
  # draw got p ~ 1/(B*m+1). On a real run that was 4950 "discoveries" at BH
  # q <= 0.05 out of 10030, manufactured entirely by a degenerate scale.
  set.seed(7)
  N <- 128L; t <- seq_len(N)
  y <- rnorm(N)
  res <- permutation_gls_test(y, gls_prepare(t, N), B = 1L, seed = 1L)
  all(is.na(res$p_pointwise)) })

check("a substantive B still produces pointwise p-values", {
  set.seed(7)
  N <- 128L; t <- seq_len(N)
  res <- permutation_gls_test(rnorm(N), gls_prepare(t, N), B = 50L, seed = 1L)
  any(is.finite(res$p_pointwise)) &&
    all(res$p_pointwise >= 0 & res$p_pointwise <= 1, na.rm = TRUE) })

check("the family-wise route survives a B too small for the pointwise one", {
  # The two routes fail independently: maxT needs only the per-draw maximum,
  # which is defined however few draws there are.
  set.seed(7)
  N <- 128L; t <- seq_len(N)
  res <- permutation_gls_test(rnorm(N), gls_prepare(t, N), B = 1L, seed = 1L)
  all(is.finite(res$p_empirical_maxT)) })

check("an unknown dataset name is refused", {
  out <- suppressWarnings(system2("./tsf",
    c("condition", "GSE130970", "499", "--config", "config/project.R",
      "--geo-dir", "d", "--interim-dir", "i", "--results-dir", "r",
      "--dry-run"), stdout = TRUE, stderr = TRUE))
  txt <- paste(out, collapse = " ")
  !is.null(attr(out, "status")) &&
    grepl("No config for dataset", txt, fixed = TRUE) &&
    grepl("499", txt, fixed = TRUE) })

# --- PLV must be calibrated, not Rayleigh-tested ------------------------------
#
# plv_rayleigh_p assumes phases iid uniform across samples. Samples here share
# one grid and one tissue, so phases agree structurally: on a real run the PLV
# distribution had median 0.971, and 98.9% of 1995 frequencies cleared
# plv_rayleigh_q <= 0.05. The permutation null carries the same shared structure,
# so it is the right reference.

check("the null draw returns its PLV instead of discarding it", {
  set.seed(3)
  m <- 40L; ns <- 8L
  keys <- paste0("1|100|", seq_len(m))
  prepared <- list(
    keys = keys, samples = paste0("S", seq_len(ns)),
    pnorm = matrix(runif(m * ns), m, dimnames = list(keys, paste0("S", 1:ns))),
    stands = matrix(rbinom(m * ns, 1, 0.5), m,
                    dimnames = list(keys, paste0("S", 1:ns))),
    phase = matrix(exp(1i * runif(m * ns, 0, 2 * pi)), m,
                   dimnames = list(keys, paste0("S", 1:ns))))
  d <- null_matrix_draw(prepared, paste0("S", 1:ns))
  is.list(d) && all(c("score", "plv") %in% names(d)) &&
    length(d$plv) == m && all(d$plv >= 0 & d$plv <= 1, na.rm = TRUE) })

check("a PLV typical of the null gets a large p, an extreme one a small p", {
  # The point of the calibration: being at 0.97 means nothing when the null
  # sits at 0.97 too.
  cs <- data.frame(chr = "1", N = 100L, k = c(1L, 2L),
                   plv = c(0.97, 0.9999), stringsAsFactors = FALSE)
  draws <- matrix(rep(rnorm(200, 0.97, 0.01), 2), nrow = 2, byrow = TRUE)
  rownames(draws) <- c("1|100|1", "1|100|2")
  out <- plv_null_pvalues(cs, list(per_key_plv = draws))
  out$p_plv_null[1] > 0.2 && out$p_plv_null[2] < 0.05 })

check("plv_rayleigh_p is left in place, not removed", {
  # Deleting a column silently changes what an existing results tree means.
  src <- paste(readLines("R/consensus.R", warn = FALSE), collapse = " ")
  grepl("rayleigh_p", src, fixed = TRUE) })

check("too few null draws leaves the calibrated p as NA", {
  cs <- data.frame(chr = "1", N = 100L, k = 1L, plv = 0.99,
                   stringsAsFactors = FALSE)
  draws <- matrix(runif(5), nrow = 1, dimnames = list("1|100|1", NULL))
  is.na(plv_null_pvalues(cs, list(per_key_plv = draws))$p_plv_null[1]) })

# --- period-indexed features --------------------------------------------------
#
# "amplitude" indexes by (chromosome, k). k is cycles per chromosome, so the
# same k is a different physical scale on each: chr1_k64 is a period of 32 genes
# and chr21_k64 is 12, and the model treats them as parallel columns.

mk_ci <- function(N) { t <- as.integer(seq(1, N, by = 2))
  list(rows = seq_along(t), t = t, N = as.integer(N), coverage = length(t) / N) }

check("a period bin means the same thing on every chromosome", {
  ci <- list("1" = mk_ci(2066), "21" = mk_ci(759))
  v <- fingerprint_vector(rnorm(2066), ci, fingerprint_terms(ci),
                          features = "period_bins")
  suffix <- sub("^.*chr[0-9XY]+_", "", names(v))
  # Both chromosomes contribute the same bin labels, so a column pair
  # chr1_p021 / chr21_p021 describes one period, not two.
  length(unique(suffix)) == length(FINGERPRINT_PERIOD_BREAKS) - 1L &&
    sum(grepl("^1[.]", names(v))) == sum(grepl("^21[.]", names(v))) })

check("the genomic spectrum is one curve over period", {
  ci <- list("1" = mk_ci(2066), "21" = mk_ci(759))
  v <- fingerprint_vector(rnorm(2066), ci, fingerprint_terms(ci),
                          features = "period_bins_genomic")
  length(v) == length(FINGERPRINT_PERIOD_BREAKS) - 1L &&
    !any(grepl("chr", names(v))) })

check("k_max does not truncate the period representations", {
  # k_max caps cycles per chromosome. On chr1 with N = 2066 that puts k <= 64
  # at period >= 32, so every bin below 32 would be empty and the short-period
  # half of the grid would exist only on the short chromosomes.
  ci <- list("1" = mk_ci(2066))
  tc <- fingerprint_terms(ci); y <- rnorm(2066)
  lo <- fingerprint_vector(y, ci, tc, k_max = 8L, features = "period_bins")
  hi <- fingerprint_vector(y, ci, tc, k_max = 512L, features = "period_bins")
  identical(lo, hi) && sum(is.finite(lo)) > 20L })

check("an empty period bin is NA, never zero", {
  # chr21 cannot hold a 500-gene period many times. Zero would claim the
  # spectrum has no power there, which is a measurement the data did not make.
  ci <- list("21" = mk_ci(759))
  v <- fingerprint_vector(rnorm(759), ci, fingerprint_terms(ci),
                          features = "period_bins")
  any(is.na(v)) && !any(v[is.finite(v)] == 0 & is.na(v)) })

check("band ratios are invariant to a global scale factor", {
  # The premise of the representation: what identifies a spectrum is the
  # relation between its landmarks, not their magnitudes, so library size and
  # sequencing depth cancel. This failed on the first attempt because log1p was
  # used -- log1p(5a) - log1p(5b) is not log1p(a) - log1p(b) -- and a feature
  # that does not have the property it exists for is worse than no feature.
  ci <- list("1" = mk_ci(2066), "21" = mk_ci(759))
  tc <- fingerprint_terms(ci); y <- abs(rnorm(2066)) + 1
  base <- fingerprint_vector(y, ci, tc, features = "band_ratios")
  all(vapply(c(5, 100, 0.01), function(k)
    isTRUE(all.equal(base, fingerprint_vector(y * k, ci, tc,
                                              features = "band_ratios"))),
    logical(1))) })

check("the amplitude representations are NOT scale invariant", {
  # Stated as a test so the distinction is deliberate: period_bins keeps log1p
  # and therefore keeps absolute scale, which is what makes band_ratios a
  # different feature rather than a renaming.
  ci <- list("1" = mk_ci(2066))
  tc <- fingerprint_terms(ci); y <- abs(rnorm(2066)) + 1
  !isTRUE(all.equal(fingerprint_vector(y, ci, tc, features = "period_bins"),
                    fingerprint_vector(y * 5, ci, tc, features = "period_bins"))) })

check("the expression baseline returns the expression, in grid order", {
  # The control the pipeline lacked: without it the accuracy of a spectral
  # fingerprint cannot be read as good or bad.
  ci <- list("1" = mk_ci(20), "21" = mk_ci(10))
  y <- seq_len(20) * 1.0
  v <- fingerprint_vector(y, ci, fingerprint_terms(ci),
                          features = "expression_baseline")
  n1 <- length(ci[["1"]]$t); n2 <- length(ci[["21"]]$t)
  length(v) == n1 + n2 &&
    isTRUE(all.equal(unname(v[seq_len(n1)]), y[ci[["1"]]$rows])) })

check("every declared feature representation builds a vector", {
  ci <- list("1" = mk_ci(2066), "21" = mk_ci(759))
  tc <- fingerprint_terms(ci); y <- rnorm(2066)
  all(vapply(FINGERPRINT_FEATURES, function(f) {
    v <- fingerprint_vector(y, ci, tc, features = f)
    !is.null(v) && length(v) > 0L && !is.null(names(v))
  }, logical(1))) })

# --- Wilson ------------------------------------------------------------------
check("Wilson interval brackets the point estimate", {
  ci <- wilson_ci(9, 10); ci[1] < 90 && ci[2] > 90 && ci[1] >= 0 && ci[2] <= 100 })
check("Wilson handles n = 0", all(is.na(wilson_ci(0, 0))))

cat("\n", if (failures == 0L) "All tests passed." else paste(failures, "test(s) failed."), "\n")
quit(status = if (failures == 0L) 0L else 1L)
