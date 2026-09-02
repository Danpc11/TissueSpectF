# peaks_genes.R -- map a spectral peak back onto the genes it runs over.

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
    A <- peaks_chr[[amp_col]][i]
    phi <- peaks_chr$phase[i]
    k <- as.integer(round(peaks_chr$freq[i] * N))
    if (!is.finite(k) || k <= 0 || k >= N) next
    if (!is.finite(A) || !is.finite(phi)) next
    sig <- sig + A * cos(2 * pi * k * n / N + phi)
  }
  sig
}

peak_file_name <- function(chr, N, k) sprintf("pico_chr%s_N%s_k%s.tsv", chr, N, k)

#' One TSV per stable peak, listing every observed gene in grid order.
write_peak_gene_tables <- function(peaks, genes, chrom_idx, out_dir, branch, cond) {
  if (is.null(peaks) || !nrow(peaks)) return(invisible(0L))
  ensure_dir(out_dir)
  amp_col <- if (branch == "average") "amplitude_mean" else "amplitude_median"
  written <- 0L

  for (i in seq_len(nrow(peaks))) {
    chr_now <- as.character(peaks$chr[i])
    ci <- chrom_idx[[chr_now]]
    if (is.null(ci)) next
    ord <- ci$rows
    N_genes <- ci$N

    # The peak's frequency is k/N measured on the grid the spectrum was
    # computed on. Reconstruction places that sinusoid on THIS dataset's axis,
    # of length ci$N. If the two N disagree, the sinusoid is laid down at the
    # wrong frequency and every crest lands on the wrong gene -- and the table
    # would still be written, carrying `N` and `grid_N` as two different
    # numbers in adjacent columns with nothing to say they should have matched.
    #
    # They can only disagree if the grids differ, which means the datasets were
    # ingested against different annotations (project$annotation_format and
    # gene_universe both change the grid). That is a configuration error, not a
    # data property, so it aborts rather than skipping: a partial peak-gene
    # library that looks complete is worse than no library.
    if (!identical(as.integer(peaks$N[i]), as.integer(ci$N))) {
      tsf_abort("chr", chr_now, ": the peak was measured on a grid of N = ",
                peaks$N[i], " but this dataset's grid has N = ", ci$N,
                ". Reconstruction would place the sinusoid at the wrong ",
                "frequency. The datasets were ingested against different ",
                "annotations or gene universes; re-run ingest for all of them ",
                "with the same project$annotation_format and gene_universe.")
    }

    amp <- peaks[[amp_col]][i]
    phi <- peaks$phase[i]
    if (!is.finite(amp) || !is.finite(phi)) next

    full <- reconstruct_signal(
      data.frame(freq = peaks$freq[i], amplitude = amp, phase = phi), N_genes)
    sig <- full[ci$t]

    col <- function(nm) if (nm %in% colnames(peaks)) peaks[[nm]][i] else NA
    # Some annotation sources provide no gene_name column.  A zero-length
    # vector here used to make data.frame() fail for otherwise valid peaks.
    gene_id <- as.character(genes$gene_id[ord])
    gene_name <- if ("gene_name" %in% colnames(genes))
      as.character(genes$gene_name[ord]) else rep(NA_character_, length(ord))
    if (length(gene_name) != length(ord)) {
      gene_name <- rep(NA_character_, length(ord))
    }

    df <- data.frame(
      condition = cond, branch = branch,
      chr = chr_now, N = peaks$N[i], k = peaks$k[i],
      freq = peaks$freq[i], period = peaks$period[i],
      phase = phi, amplitude = amp,
      power = col("power"),
      power_normalised = col("power_normalised"),
      window_power = col("window_power"),
      window_rank = col("window_rank"),
      n_samples_expected = col("n_samples_expected"),
      n_samples_significant = col("n_samples_significant"),
      pct_samples_significant = col("pct_samples_significant"),
      gene_id = gene_id,
      gene_name = gene_name,
      grid_position = ci$t,
      grid_N = ci$N,
      coverage = ci$coverage,
      reconstructed_signal = sig,
      sign = ifelse(sig > 0, "positive", ifelse(sig < 0, "negative", "zero")),
      magnitude = abs(sig),
      stringsAsFactors = FALSE)

    write_tsv_tsf(df, file.path(out_dir,
      peak_file_name(chr_now, peaks$N[i], peaks$k[i])))
    written <- written + 1L
  }
  invisible(written)
}
