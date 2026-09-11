#!/usr/bin/env Rscript
# crest_genes.R -- which genes carry a characteristic peak.
#
# Usage:
#   Rscript scripts/crest_genes.R --signature <condition_signature_F3.tsv> \
#       --datasets GSE135251,GSE130970 --condition F3 [--top 30] [--out <dir>]
#
# A peak at frequency k of chromosome c is a property of the WHOLE profile, not
# of a region; listing the genes "under the peak" answers a different question.
# Two quantities per gene, computed on the condition's mean profile x (the
# signal the condition test used; with a reference profile applied this is the
# mean deviation from healthy tissue):
#
#   contribution  c_j = x_j * u_j     where u is the unit-norm fitted component
#                 a cos(wt_j) + b sin(wt_j) / sqrt(a^2+b^2), fitted by least
#                 squares on the observed positions (the same model GLS fits).
#                 c_j > 0: the gene sits where the component is high and is
#                 above the mean (a crest), or sits low and is below (a trough
#                 gene) -- either way it SUPPORTS the periodicity. Sum_j c_j is
#                 the component's amplitude times n/2. `share` is c_j over the
#                 sum of positive contributions.
#
#   delta_power   fraction of the periodogram power at k lost when gene j is
#                 removed: 1 - |W - x_j e^{-iwt_j}|^2 / |W|^2. Exact for the
#                 periodogram, a close approximation for the GLS power the
#                 pipeline reports; it is an ablation, so it is the number to
#                 quote when saying "this gene carries the peak".
#
# Significance: x is permuted across the observed positions of the chromosome
# B times and the distribution of the LARGEST |c_j| is recorded; each gene gets
# p_maxperm = fraction of permutations whose largest contribution exceeds its
# own. That is a family-wise p over the chromosome's genes, like maxT.
#
# Output: <out>/crest_genes_<condition>.tsv (all genes of all peaks), and
# <out>/crest_genes_<condition>_top.tsv (top N per peak), aggregated over the
# datasets given (mean contribution, and in how many datasets the gene is in
# the peak's top decile).
if (any(commandArgs(TRUE) %in% c("-h", "--help"))) {
  cat(paste(sub("^# ?", "", grep("^#", readLines(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))[-1], value = TRUE)), collapse = "\n"), "\n")
  quit(save = "no")
}
suppressWarnings({ source("R/utils_io.R"); tsf_load_all("R") })
args <- commandArgs(trailingOnly = TRUE)
flag <- function(n, d = NULL) {
  h <- grep(paste0("^", n, "="), args, value = TRUE)
  if (length(h)) return(sub(paste0("^", n, "="), "", h[1]))
  i <- match(n, args); if (!is.na(i) && length(args) > i) args[i + 1] else d
}
# --geo-dir / --interim-dir / --results-dir override the environment, same as ./tsf
for (pair in list(c("--geo-dir", "TSF_GEO_DIR"), c("--interim-dir", "TSF_INTERIM_DIR"),
                  c("--results-dir", "TSF_RESULTS_DIR"), c("--config", "TSF_CONFIG"))) {
  v <- flag(pair[1], NULL); if (!is.null(v)) do.call(Sys.setenv, stats::setNames(list(v), pair[2]))
}
project <- load_project_config(Sys.getenv("TSF_CONFIG", "config/project.R"))
sig_path <- flag("--signature", ""); if (!nzchar(sig_path)) tsf_abort("Pasa --signature <tsv>")
datasets <- trimws(strsplit(flag("--datasets", ""), ",")[[1]])
if (!length(datasets)) tsf_abort("Pasa --datasets A,B")
condition <- flag("--condition", ""); if (!nzchar(condition)) tsf_abort("Pasa --condition F3")
top_n <- as.integer(flag("--top", "30"))
B <- as.integer(flag("--perm", "500"))
out_dir <- flag("--out", file.path(project$results_dir, "crest_genes"))
ensure_dir(out_dir)

sig <- read_tsv_tsf(sig_path)
need <- c("chr", "N", "k")
if (!all(need %in% colnames(sig))) tsf_abort("signature lacks ", paste(setdiff(need, colnames(sig)), collapse = ","))
sig <- sig[!duplicated(paste(sig$chr, sig$N, sig$k)), , drop = FALSE]
tsf_log(nrow(sig), " peak(s) in ", basename(sig_path))

crest_one <- function(x, t0, N, k, B) {
  # x: mean profile on observed positions (centred here), t0: 0-based positions
  x <- x - mean(x)
  w <- 2 * pi * k / N
  C <- cos(w * t0); S <- sin(w * t0)
  X <- cbind(C, S)
  ab <- tryCatch(stats::lm.fit(X, x)$coefficients, error = function(e) c(NA, NA))
  if (any(!is.finite(ab))) return(NULL)
  u <- (ab[1] * C + ab[2] * S); u <- u / sqrt(sum(u^2))
  contrib <- x * u
  W <- sum(x * exp(-1i * w * t0))
  P <- Mod(W)^2
  W_minus <- W - x * exp(-1i * w * t0)
  delta <- 1 - Mod(W_minus)^2 / P
  # permutation of x across positions: distribution of the largest |c_j|
  set.seed(k * 7919L + N)
  null_max <- vapply(seq_len(B), function(b) {
    xp <- sample(x)
    abp <- stats::lm.fit(X, xp)$coefficients
    up <- abp[1] * C + abp[2] * S; up <- up / sqrt(sum(up^2))
    max(abs(xp * up))
  }, numeric(1))
  p_maxperm <- vapply(abs(contrib), function(c) (sum(null_max >= c) + 1) / (B + 1), numeric(1))
  list(contrib = contrib, share = pmax(contrib, 0) / sum(pmax(contrib, 0)),
       delta_power = delta, p_maxperm = p_maxperm, amplitude = 2 * sqrt(sum(ab^2)) / 2,
       phase = atan2(-ab[2], ab[1]), power = P)
}

rows <- list()
for (ds in datasets) {
  d <- tryCatch(load_dataset(ds, project), error = function(e) NULL)
  if (is.null(d)) { tsf_warn(ds, ": not ingested, skipped"); next }
  keep <- d$samples$sample_id[as.character(d$samples$condition) == condition]
  keep <- intersect(keep, colnames(d$expression))
  if (length(keep) < 3L) { tsf_log(ds, ": ", length(keep), " sample(s) of ", condition, ", skipped"); next }
  xbar <- rowMeans(d$expression[, keep, drop = FALSE], na.rm = TRUE)
  g <- d$genes
  g$gene_id <- sub("\\..*$", "", g$gene_id)
  for (i in seq_len(nrow(sig))) {
    chr <- as.character(sig$chr[i]); N <- as.integer(sig$N[i]); k <- as.integer(sig$k[i])
    gi <- g[as.character(g$chr) == chr, , drop = FALSE]
    if (!nrow(gi) || !all(gi$grid_N == N)) next   # a different grid: not comparable
    x <- xbar[match(gi$gene_id, sub("\\..*$", "", names(xbar)))]
    ok <- is.finite(x)
    gi <- gi[ok, , drop = FALSE]; x <- x[ok]
    if (length(x) < 16L) next
    r <- crest_one(x, as.numeric(gi$grid_index) - 1, N, k, B)
    if (is.null(r)) next
    rows[[length(rows) + 1L]] <- data.frame(
      dataset = ds, condition = condition, chr = chr, N = N, k = k, period = N / k,
      gene_id = gi$gene_id,
      gene_name = if ("gene_name" %in% colnames(gi)) gi$gene_name else NA_character_,
      grid_index = gi$grid_index,
      mean_signal = x, contribution = r$contrib, share = r$share,
      delta_power = r$delta_power, p_maxperm = r$p_maxperm,
      component_amplitude = r$amplitude, component_phase = r$phase,
      stringsAsFactors = FALSE)
  }
}
if (!length(rows)) tsf_abort("nothing computed: check --datasets/--condition and that the grid N matches the signature")
all <- do.call(rbind, rows)

# aggregate across datasets
key <- paste(all$chr, all$N, all$k, all$gene_id)
all$top_decile <- ave(all$contribution, paste(all$dataset, all$chr, all$k),
                      FUN = function(v) v >= stats::quantile(v, 0.9, na.rm = TRUE))
agg <- stats::aggregate(cbind(contribution, share, delta_power, top_decile) ~
                          chr + N + k + period + gene_id + gene_name + grid_index,
                        data = all, FUN = mean, na.action = na.pass)
n_ds <- stats::aggregate(dataset ~ chr + N + k + gene_id, data = all,
                         FUN = function(v) length(unique(v)))
agg <- merge(agg, n_ds, by = c("chr", "N", "k", "gene_id"))
names(agg)[names(agg) == "dataset"] <- "n_datasets"
names(agg)[names(agg) == "top_decile"] <- "frac_datasets_top_decile"
pmin_agg <- stats::aggregate(p_maxperm ~ chr + N + k + gene_id, data = all, FUN = min)
names(pmin_agg)[5] <- "p_maxperm_min"
agg <- merge(agg, pmin_agg, by = c("chr", "N", "k", "gene_id"))
agg <- agg[order(agg$chr, agg$k, -agg$contribution), , drop = FALSE]

top <- do.call(rbind, lapply(split(agg, paste(agg$chr, agg$k)), function(t) {
  t <- t[order(-t$contribution), , drop = FALSE]
  utils::head(t, top_n)
}))
rownames(top) <- NULL

write_tsv_tsf(all, file.path(out_dir, sprintf("crest_genes_%s_by_dataset.tsv", condition)))
write_tsv_tsf(agg, file.path(out_dir, sprintf("crest_genes_%s.tsv", condition)))
write_tsv_tsf(top, file.path(out_dir, sprintf("crest_genes_%s_top.tsv", condition)))
tsf_log(nrow(agg), " gene x peak rows, ", nrow(top), " in the top table -> ", out_dir)
