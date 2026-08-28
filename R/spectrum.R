# spectrum.R -- FFT of chromosome-ordered expression.
#
# The numerical core is a direct port of run_fft() / poly_equation() from the
# original scripts: same one-sided scaling, same sequential Fisher g-test, same
# reported columns. Verified against the legacy output by
# scripts/99_validate_against_legacy.R -- do not "clean up" the formulas without
# re-running that comparison.

#' Fisher's g-test tail probability (Siegel's sequential form).
poly_equation <- function(x, alpha) {
  if (!is.finite(x) || x <= 0) return(1)
  if (x >= 1) return(0)
  r_max <- min(floor(1 / x), alpha)
  if (!is.finite(r_max) || r_max < 1) return(1)
  log_terms <- rep(NA_real_, r_max + 1)
  signs     <- rep(NA_real_, r_max + 1)
  for (r in 1:r_max) {
    val <- 1 - r * x
    if (val <= 0 || !is.finite(val)) break
    log_terms[r + 1] <- lchoose(alpha, r) + (alpha - 1) * log(val)
    signs[r + 1] <- (-1)^(r - 1)
  }
  valid <- is.finite(log_terms) & is.finite(signs)
  if (!any(valid)) return(1)
  max_log <- max(log_terms[valid])
  s <- sum(signs[valid] * exp(log_terms[valid] - max_log)) * exp(max_log)
  max(min(s, 1), 0)
}

#' One-sided spectrum of a mean-centred signal, with sequential Fisher g p-values.
#'
#' Returns k = 0..N/2 with freq = k/N (cycles per gene) and period = N/k
#' (genes per cycle). k = 0 carries no power by construction.
run_fft <- function(signal) {
  signal <- suppressWarnings(as.numeric(signal))
  signal[!is.finite(signal)] <- 0
  if (length(signal) < 4) return(NULL)
  N <- length(signal)
  signal <- signal - mean(signal)
  if (all(signal == 0)) return(NULL)
  Y <- tryCatch(stats::fft(signal), error = function(e) NULL)
  if (is.null(Y)) return(NULL)

  k  <- 0:(N %/% 2)
  Yh <- Y[k + 1]
  amp <- Mod(Yh) / N
  power_raw <- (Mod(Yh) / N)^2
  # Fold the negative frequencies onto the positive ones, leaving DC (and the
  # Nyquist bin for even N) unscaled.
  if (N %% 2 == 0) {
    if (length(amp) > 2) {
      amp[2:(length(amp) - 1)]             <- 2 * amp[2:(length(amp) - 1)]
      power_raw[2:(length(power_raw) - 1)] <- 2 * power_raw[2:(length(power_raw) - 1)]
    }
  } else {
    if (length(amp) > 1) {
      amp[2:length(amp)]             <- 2 * amp[2:length(amp)]
      power_raw[2:length(power_raw)] <- 2 * power_raw[2:length(power_raw)]
    }
  }

  phase <- Arg(Yh)
  total_power_original <- sum(power_raw[-1])
  power_norm <- rep(NA_real_, length(power_raw))
  if (is.finite(total_power_original) && total_power_original > 0) {
    power_norm[-1] <- power_raw[-1] / total_power_original
  }

  power_test <- power_raw
  power_test[1] <- 0
  if (N %% 2 == 0) power_test[length(power_test)] <- 0
  if (sum(power_test > 0) == 0) return(NULL)

  pvals <- rep(1, length(power_test))
  remaining_idx <- which(power_test > 0)
  repeat {
    if (length(remaining_idx) < 2) break
    total_remaining <- sum(power_test[remaining_idx])
    if (!is.finite(total_remaining) || total_remaining <= 0) break
    x_norm <- power_test[remaining_idx] / total_remaining
    max_pos <- which.max(x_norm)
    max_idx <- remaining_idx[max_pos]
    p_max <- poly_equation(x_norm[max_pos], length(remaining_idx))
    pvals[max_idx] <- p_max
    if (!is.finite(p_max) || p_max > 0.05) break
    remaining_idx <- remaining_idx[-max_pos]
  }

  power_reported <- power_raw
  power_reported[1] <- 0

  data.frame(
    k = k, freq = k / N, amplitude = amp,
    power = power_reported, power_norm = power_norm, phase = phase,
    real = Re(Yh), imag = Im(Yh), p_value = pvals,
    period = ifelse(k == 0, NA_real_, N / k), N = N,
    stringsAsFactors = FALSE)
}

#' Gene indices per chromosome, in gene_order, for chromosomes long enough.
chromosome_index <- function(genes, chrom_levels, min_genes_per_chr = 8L) {
  out <- list()
  for (chr_now in chrom_levels) {
    idx <- which(genes$chr == chr_now)
    if (length(idx) < min_genes_per_chr) next
    out[[chr_now]] <- idx[order(genes$gene_order[idx])]
  }
  out
}

#' Per-sample and per-condition signals for one condition.
#'
#' avg_signal / median_signal are the summary signals the condition-level
#' spectra are computed from; they are NOT averages of per-sample spectra.
condition_signals <- function(dataset, cond) {
  samps <- dataset$samples$sample_id[!is.na(dataset$samples$condition) &
                                       dataset$samples$condition == cond]
  samps <- intersect(samps, colnames(dataset$expression))
  if (length(samps) < 2) return(NULL)
  mat <- dataset$expression[, samps, drop = FALSE]
  mat[!is.finite(mat)] <- 0
  list(samples = samps,
       matrix = mat,
       avg_signal = rowMeans(mat),
       median_signal = apply(mat, 1, stats::median))
}

#' Spectra of every signal of one condition, over every usable chromosome.
compute_condition_spectra <- function(dataset, cond, chrom_idx) {
  sig <- condition_signals(dataset, cond)
  if (is.null(sig)) {
    tsf_warn("Condition ", cond, ": fewer than 2 samples, skipped")
    return(NULL)
  }
  tsf_log("Condition ", cond, ": ", length(sig$samples), " samples, ",
          length(chrom_idx), " chromosomes")

  all_signals <- c(list(avg_signal = sig$avg_signal, median_signal = sig$median_signal),
                   stats::setNames(lapply(sig$samples, function(s) sig$matrix[, s]), sig$samples))

  rows <- list()
  for (chr_now in names(chrom_idx)) {
    ord <- chrom_idx[[chr_now]]
    for (nm in names(all_signals)) {
      res <- run_fft(all_signals[[nm]][ord])
      if (is.null(res) || !nrow(res)) next
      res$chr <- chr_now
      res$sample <- nm
      rows[[length(rows) + 1]] <- res
    }
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  out$condition <- cond
  attr(out, "samples") <- sig$samples
  out
}

#' Split a spectrum table into the two condition-level branches.
branch_spectrum <- function(spectra, branch = c("average", "median")) {
  branch <- match.arg(branch)
  nm <- if (branch == "average") "avg_signal" else "median_signal"
  out <- spectra[spectra$sample == nm, , drop = FALSE]
  colnames(out)[colnames(out) == "amplitude"] <-
    if (branch == "average") "amplitude_mean" else "amplitude_median"
  out$sample <- NULL
  out
}
