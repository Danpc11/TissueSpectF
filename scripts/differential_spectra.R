#!/usr/bin/env Rscript
# differential_spectra.R -- the characteristic spectrum of each condition.
#
#   Rscript scripts/differential_spectra.R \
#     --results-dir results_gencode_2026_09_02 \
#     --datasets GSE130970,GSE135251,GSE162694,GSE276114 \
#     [--period-bins] [--min-period 10] [--out-dir <dir>]
#
# WHY THIS IS SEPARATE FROM THE CONSENSUS STAGE
# ---------------------------------------------
# `consensus` asks DETECTION: is this frequency stronger than a null in which
# gene positions were shuffled? That requires a component to stand out within a
# condition, and a difference between conditions does not have to stand out
# anywhere to exist.
#
# This asks COMPARISON. The spectrum is a linear transform of the ordered
# expression vector, so if two conditions differ in expression -- which is what
# differential expression means -- their spectra differ too. The difference may
# be small; it is not zero. This is that test, one frequency at a time.
#
# WHY --period-bins IS THE DEFAULT RECOMMENDATION
# -----------------------------------------------
# Power, not elegance. The minimum detectable effect at 80% power, one condition
# against the rest, was measured on this design:
#
#   n = 54, m = 1995 frequencies   ->  0.85 sd
#   n = 54, m =  200               ->  0.75 sd
#
# Multiplicity is the binding constraint, and the only lever on it is m.
# Collapsing to ~40 period bands cuts m by fifty and averages across
# chromosomes, which lowers the variance of the estimate as well. A difference
# that is real but modest is detectable there and buried at m = 1995.
#
# The bands are the common log-period grid from R/fingerprint.R, so a band means
# the same thing on every chromosome -- which (chromosome, k) does not, since k
# is cycles per chromosome.

suppressWarnings({
  source("R/utils_io.R"); source("R/config.R")
  source("R/grid.R"); source("R/period_floor.R")
  source("R/fingerprint.R"); source("R/differential.R")
  source("R/paths.R")
})

args <- commandArgs(trailingOnly = TRUE)
flag <- function(n, d = NULL) {
  h <- grep(paste0("^", n, "="), args, value = TRUE)
  if (length(h)) return(sub(paste0("^", n, "="), "", h[1]))
  i <- match(n, args); if (!is.na(i) && length(args) > i) args[i + 1] else d
}
has <- function(n) n %in% args

results_dir <- flag("--results-dir", Sys.getenv("TSF_RESULTS_DIR", ""))
datasets    <- strsplit(flag("--datasets", ""), ",")[[1]]
out_dir     <- flag("--out-dir", file.path(results_dir, "differential"))
period_bins <- has("--period-bins")
min_period  <- as.numeric(flag("--min-period", "10"))
stage_order <- strsplit(flag("--stage-order",
                             "Controles,F0,F1,F2,F3,F4"), ",")[[1]]

if (!nzchar(results_dir)) tsf_abort("Pass --results-dir <dir>.")
if (!length(datasets)) tsf_abort("Pass --datasets A,B,C.")

ensure_dir(out_dir)

for (id in datasets) {
  base <- file.path(results_dir, id, "spectra")
  files <- list.files(base, pattern = "^spectra_samples_.*\\.tsv$",
                      full.names = TRUE)
  if (!length(files)) {
    tsf_warn(id, ": no per-sample spectra in ", base, "; skipped")
    next
  }

  parts <- lapply(files, function(f) {
    d <- read_tsv_tsf(f)
    d$condition <- sub("^spectra_samples_(.*)\\.tsv$", "\\1", basename(f))
    d
  })
  sp <- do.call(rbind, parts)

  need <- c("chr", "N", "k", "sample", "power_normalised", "period")
  miss <- setdiff(need, names(sp))
  if (length(miss)) {
    tsf_warn(id, ": per-sample spectra lack ", paste(miss, collapse = ", "),
             "; re-run the spectra stage. Skipped.")
    next
  }

  groups <- stats::setNames(sp$condition[!duplicated(sp$sample)],
                            sp$sample[!duplicated(sp$sample)])
  tsf_log("=== ", id, ": ", length(groups), " sample(s), ",
          length(unique(groups)), " condition(s) ===")

  # A period floor before anything, so a frequency the sampling cannot resolve
  # never enters the family and never costs multiplicity.
  sp <- sp[is.finite(sp$period) & sp$period >= min_period, , drop = FALSE]

  if (period_bins) {
    br <- FINGERPRINT_PERIOD_BREAKS
    sp <- sp[sp$period >= min(br) & sp$period <= max(br), , drop = FALSE]
    idx <- cut(sp$period, breaks = br, include.lowest = TRUE, labels = FALSE)
    lab <- sprintf("p%05.1f", sqrt(br[-1] * br[-length(br)]))
    # Average power within (sample, band) ACROSS chromosomes: one genome-wide
    # curve over period per sample. chr and N are set to constants so the
    # downstream key stays well formed, and k becomes the band index.
    agg <- stats::aggregate(
      sp$power_normalised,
      by = list(sample = sp$sample, band = idx),
      FUN = mean, na.rm = TRUE)
    sp <- data.frame(chr = "genome", N = length(lab),
                     k = agg$band, sample = agg$sample,
                     power_normalised = agg$x,
                     band = lab[agg$band], stringsAsFactors = FALSE)
    tsf_log("  collapsed to ", length(unique(sp$k)), " period band(s): ",
            "multiplicity is the binding constraint and m is the only lever")
  }

  res <- differential_spectrum(sp, groups, stage_order = stage_order)
  if (is.null(res)) { tsf_warn(id, ": nothing to test"); next }

  if (period_bins) {
    lab <- sprintf("p%05.1f", sqrt(FINGERPRINT_PERIOD_BREAKS[-1] *
      FINGERPRINT_PERIOD_BREAKS[-length(FINGERPRINT_PERIOD_BREAKS)]))
    res$band <- lab[res$k]
  }

  path <- file.path(out_dir, paste0(id, "_differential.tsv"))
  write_tsv_tsf(res, path)

  sig <- res[is.finite(res$q) & res$q <= 0.05, , drop = FALSE]
  tsf_log("  ", nrow(sig), " result(s) at q <= 0.05 of ", nrow(res), " test(s)")
  for (t in unique(res$test)) {
    n <- sum(sig$test == t)
    tsf_log("    ", t, ": ", n)
  }
  tsf_log("  wrote ", path)
}

tsf_log("Read `effect` beside `q`: with tens of samples a trivial shift ",
        "reaches significance, and the effect size is what says whether it ",
        "is worth anything.")
