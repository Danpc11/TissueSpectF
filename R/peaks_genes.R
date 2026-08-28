# peaks_genes.R -- map a spectral peak back onto the genes it runs over.
#
# For a peak (chr, N, k) the reconstructed signal is a single sinusoid evaluated
# at every gene position along the chromosome. Genes where it is positive sit on
# a crest of that periodic component, genes where it is negative sit in a
# trough. This is a projection, not a statement about the individual gene.
#
# This stage reads everything it needs from disk (stability tables, condition
# spectra, expression, genes) and holds no session state, which is what the old
# regenerar_picos_genes_condicion.R had to work around: re-running it after a
# threshold change is now just `Rscript scripts/05_peaks_genes.R`, with no need
# to keep the previous session alive or to recompute FFT and maxT.

#' Sum of sinusoids for a set of peaks over N positions.
reconstruct_signal <- function(peaks_chr, N) {
  amp_col <- if ("amplitude_mean" %in% colnames(peaks_chr)) "amplitude_mean"
    else if ("amplitude_median" %in% colnames(peaks_chr)) "amplitude_median"
    else if ("amplitude" %in% colnames(peaks_chr)) "amplitude"
    else NULL
  if (is.null(amp_col)) return(rep(0, N))

  n <- 0:(N - 1)
  sig <- rep(0, N)
  for (i in seq_len(nrow(peaks_chr))) {
    A   <- peaks_chr[[amp_col]][i]
    phi <- peaks_chr$phase[i]
    k   <- as.integer(round(peaks_chr$freq[i] * N))
    if (!is.finite(k) || k <= 0 || k >= N) next
    if (!is.finite(A) || !is.finite(phi)) next
    sig <- sig + A * cos(2 * pi * k * n / N + phi)
  }
  sig
}

peak_file_name <- function(chr, N, k) sprintf("pico_chr%s_N%s_k%s.tsv", chr, N, k)

#' One TSV per stable peak, listing every gene of that chromosome in order.
write_peak_gene_tables <- function(peaks, genes, chrom_idx, out_dir, branch, cond) {
  if (is.null(peaks) || !nrow(peaks)) return(invisible(0L))
  ensure_dir(out_dir)
  amp_col <- if (branch == "average") "amplitude_mean" else "amplitude_median"
  written <- 0L

  for (i in seq_len(nrow(peaks))) {
    chr_now <- as.character(peaks$chr[i])
    ord <- chrom_idx[[chr_now]]
    if (is.null(ord)) next
    N_genes <- length(ord)

    amp <- peaks[[amp_col]][i]
    phi <- peaks$phase[i]
    if (!is.finite(amp) || !is.finite(phi)) next

    sig <- reconstruct_signal(
      data.frame(freq = peaks$freq[i], amplitude = amp, phase = phi), N_genes)

    df <- data.frame(
      condition = cond, branch = branch,
      chr = chr_now, N = peaks$N[i], k = peaks$k[i],
      freq = peaks$freq[i], period = peaks$period[i],
      phase = phi, amplitude = amp,
      power = peaks$power[i],
      power_norm = peaks$power_norm[i],
      p_value_fisher = peaks$p_value_fisher[i],
      p_fdr_fisher = peaks$p_fdr_fisher[i],
      n_samples_expected = peaks$n_samples_expected[i],
      n_samples_significant = peaks$n_samples_significant[i],
      pct_samples_significant = peaks$pct_samples_significant[i],
      gene_id = genes$gene_id[ord],
      gene_name = genes$gene_name[ord],
      gene_position = seq_len(N_genes),
      reconstructed_signal = sig,
      sign = ifelse(sig > 0, "positive", ifelse(sig < 0, "negative", "zero")),
      magnitude = abs(sig),
      stringsAsFactors = FALSE)

    write_tsv_tsf(df, file.path(out_dir, peak_file_name(chr_now, peaks$N[i], peaks$k[i])))
    written <- written + 1L
  }
  invisible(written)
}
