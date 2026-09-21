#!/usr/bin/env Rscript
# sensitivity_quantile_cut.R -- how many candidate/robust condition components
# survive the meta-analysis (R/stages.R, condition_library/) as the "stands
# out" quantile (TSF_QUANTILE_CUT, R/consensus.R::prevalence_from_rank) is
# relaxed from the default top 5% (quantile_cut = 0.95) to top 10% and top 20%.
#
# Reads condition_library/condition_signature_*.tsv from each results tree
# passed on the command line -- these are the files run_differential.sh's own
# crest_genes step already reads to decide whether a condition has anything to
# extract genes from -- and counts, per condition and in total, how many
# components made it in (n_candidates, every row of the file) and how many
# carry signature_class == "robust" (n_robust). Base R only, no packages.
#
# Usage:
#   Rscript sensitivity_quantile_cut.R \
#     0.95=/path/to/results_diff_q5/combined \
#     0.90=/path/to/results_diff_q10/combined \
#     0.80=/path/to/results_diff_q20/combined \
#     --out /path/to/output_dir
#
# Each "QUANTILE=PATH" argument is one results tree's *library* dir (the one
# that directly contains condition_library/), i.e. the same path you passed
# to run_differential.sh's --results-dir plus the library name (here: only
# "combined", since that is what --only combined built for the three runs).

# --help prints this file's own header comment and exits with 0. Before this
# guard, --help was taken as the VALUE of the previous flag, or the script ran
# with an empty configuration. A script that cannot explain itself is a script
# nobody uses correctly (tests/test_labels.R: "todo script ejecutable responde
# a --help y sale con 0").
if (any(commandArgs(TRUE) %in% c("-h", "--help"))) {
  self <- sub("^--file=", "",
              grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  if (!is.na(self) && file.exists(self)) {
    hdr <- readLines(self, warn = FALSE)
    hdr <- hdr[!grepl("^#!", hdr)]                    # drop the shebang
    stop_at <- which(!grepl("^#", hdr) & nzchar(hdr))[1]
    if (is.na(stop_at)) stop_at <- length(hdr) + 1L
    hdr <- hdr[seq_len(stop_at - 1L)]
    cat(paste(sub("^#[ ]?", "", hdr[grepl("^#", hdr)]), collapse = "\n"), "\n")
  }
  quit(save = "no", status = 0)
}

args <- commandArgs(trailingOnly = TRUE)
out_dir <- "."
pairs <- list()
i <- 1
while (i <= length(args)) {
  a <- args[i]
  if (a == "--out") { out_dir <- args[i + 1]; i <- i + 2; next }
  kv <- strsplit(a, "=", fixed = TRUE)[[1]]
  if (length(kv) != 2) stop("bad argument, expected QUANTILE=PATH: ", a)
  pairs[[kv[1]]] <- kv[2]
  i <- i + 1
}
if (!length(pairs)) {
  stop("usage: Rscript sensitivity_quantile_cut.R Q1=PATH1 Q2=PATH2 ... [--out DIR]")
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

count_one <- function(quantile_cut, lib_dir) {
  sig_dir <- file.path(lib_dir, "condition_library")
  files <- list.files(sig_dir, pattern = "^condition_signature_.*\\.tsv$",
                      full.names = TRUE)
  if (!length(files)) {
    warning("no condition_signature_*.tsv under ", sig_dir)
    return(data.frame(quantile_cut = numeric(0), condition = character(0),
                       n_candidates = integer(0), n_robust = integer(0)))
  }
  rows <- lapply(files, function(f) {
    cond <- sub("^condition_signature_", "", sub("\\.tsv$", "", basename(f)))
    # a header-only file (no components at all) is a valid, common outcome --
    # not an error -- so this must not abort the whole comparison
    d <- tryCatch(utils::read.delim(f, sep = "\t", stringsAsFactors = FALSE),
                  error = function(e) NULL)
    n_cand <- if (is.null(d)) 0L else nrow(d)
    n_rob <- if (is.null(d) || !("signature_class" %in% names(d))) 0L
             else sum(d$signature_class == "robust", na.rm = TRUE)
    data.frame(quantile_cut = as.numeric(quantile_cut), condition = cond,
               n_candidates = n_cand, n_robust = n_rob)
  })
  do.call(rbind, rows)
}

summary_df <- do.call(rbind, Map(count_one, names(pairs), pairs))
summary_df <- summary_df[order(-summary_df$quantile_cut, summary_df$condition), ]
rownames(summary_df) <- NULL

write.table(summary_df, file.path(out_dir, "sensitivity_quantile_cut.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

totals <- stats::aggregate(cbind(n_candidates, n_robust) ~ quantile_cut,
                           summary_df, sum)
totals <- totals[order(-totals$quantile_cut), ]
write.table(totals, file.path(out_dir, "sensitivity_quantile_cut_totals.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Por condicion y cuantil:\n"); print(summary_df)
cat("\nTotales (suma de las ocho condiciones), por cuantil:\n"); print(totals)

# The x axis is "top X%" = (1 - quantile_cut) * 100, the quantity that was
# actually declared per run (5, 10, 20) -- reads left-to-right as "how
# permissive the cut is", the opposite order of quantile_cut itself.
totals$top_pct <- round((1 - totals$quantile_cut) * 100, 1)
totals <- totals[order(totals$top_pct), ]

png(file.path(out_dir, "sensitivity_quantile_cut.png"),
    width = 900, height = 650, res = 120)
ylim <- c(0, max(1, totals$n_candidates, totals$n_robust) * 1.15)
plot(totals$top_pct, totals$n_candidates, type = "b", pch = 16, col = "#1f6feb",
     lwd = 2, cex = 1.3,
     xlab = "Cuantil de prevalencia (\"top X%\" de potencia normalizada)",
     ylab = "Componentes (suma de las ocho condiciones)",
     ylim = ylim, xaxt = "n",
     main = "Sensibilidad de candidatos/robustos a quantile_cut (libreria combined)")
axis(1, at = totals$top_pct, labels = paste0(totals$top_pct, "%"))
lines(totals$top_pct, totals$n_robust, type = "b", pch = 17, col = "#d1242f",
      lwd = 2, cex = 1.3)
legend("topleft", legend = c("candidatos", "robustos"),
       col = c("#1f6feb", "#d1242f"), pch = c(16, 17), lty = 1, bty = "n")
grid(nx = NA, ny = NULL, lty = "dotted")
invisible(dev.off())

cat("\nEscritos:\n  ", file.path(out_dir, "sensitivity_quantile_cut.tsv"),
    "\n  ", file.path(out_dir, "sensitivity_quantile_cut_totals.tsv"),
    "\n  ", file.path(out_dir, "sensitivity_quantile_cut.png"), "\n", sep = "")
