#!/usr/bin/env Rscript

# Build one cross-cohort power spectrum per biological condition.
#
# Inputs are the cohort-level consensus_spectrum_<condition>.tsv files.  Power
# is combined with a median (robust to one cohort/platform), never by pooling
# samples.  A within-cohort log2 contrast against that cohort's other conditions
# is carried alongside the absolute power so tissue-wide structure is not
# mistaken for condition specificity.

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

abort <- function(...) stop(paste0(...), call. = FALSE)
messagef <- function(...) message(sprintf(...))

parse_args <- function(x) {
  out <- list(results_dir = "results_gencode_v2", out_dir = NULL,
              datasets = NULL, conditions = NULL, min_cohorts = 2L,
              top = NA_integer_, signature_q = 0.05,
              window_cut = 1, min_prevalence = 0.5,
              exclude_chromosomes = "Y,MT", healthy_superclass = "no",
              min_period = "off", period_margin = 2, margin_mode = "add",
              min_period_biological = 0)
  i <- 1L
  while (i <= length(x)) {
    a <- x[[i]]
    if (!startsWith(a, "--")) abort("Unexpected argument: ", a)
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    key <- gsub("-", "_", kv[[1]])
    if (length(kv) > 1L) value <- paste(kv[-1], collapse = "=") else {
      i <- i + 1L
      if (i > length(x)) abort("Missing value for --", key)
      value <- x[[i]]
    }
    if (!key %in% names(out)) abort("Unknown option --", key)
    out[[key]] <- value
    i <- i + 1L
  }
  out$out_dir <- out$out_dir %||% file.path(out$results_dir, "condition_library")
  if (!is.null(out$datasets)) out$datasets <- strsplit(out$datasets, ",", fixed = TRUE)[[1]]
  if (!is.null(out$conditions)) out$conditions <- strsplit(out$conditions, ",", fixed = TRUE)[[1]]
  out$min_cohorts <- as.integer(out$min_cohorts)
  out$min_prevalence <- as.numeric(out$min_prevalence)
  out$exclude_chromosomes <- if (nzchar(out$exclude_chromosomes))
    strsplit(out$exclude_chromosomes, ",", fixed = TRUE)[[1]] else character(0)
  out$healthy_superclass <- tolower(out$healthy_superclass) %in% c("yes", "true", "1")
  out$top <- suppressWarnings(as.integer(out$top))
  out$signature_q <- as.numeric(out$signature_q)
  out$period_margin <- as.numeric(out$period_margin)
  out$min_period_biological <- as.numeric(out$min_period_biological)
  out$margin_mode <- tolower(out$margin_mode)
  if (!out$margin_mode %in% c("add", "mult")) {
    abort("--margin-mode must be add or mult")
  }
  if (!identical(tolower(out$min_period), "off") &&
      !identical(tolower(out$min_period), "auto")) {
    v <- suppressWarnings(as.numeric(out$min_period))
    if (!is.finite(v) || v <= 0) abort("--min-period must be off, auto or a number")
  }
  out$window_cut <- as.numeric(out$window_cut)
  if (is.na(out$min_cohorts) || out$min_cohorts < 2L) abort("--min-cohorts must be >= 2")
  if (!is.na(out$top) && out$top < 1L) abort("--top must be >= 1 when given")
  if (!is.finite(out$signature_q) || out$signature_q <= 0 || out$signature_q > 1) {
    abort("--signature-q must be in (0, 1]")
  }
  out
}

read_tsv <- function(path) {
  utils::read.delim(path, sep = "\t", check.names = FALSE,
                    stringsAsFactors = FALSE, quote = "", comment.char = "")
}

write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", row.names = FALSE,
                     col.names = TRUE, quote = FALSE, na = "NA")
}

# Base rbind refuses data frames whose optional inferential columns differ.
# Align by column name and represent unavailable quantities as NA; this is
# essential for healthy controls such as GSE142530, whose within-cohort null is
# not estimable and therefore may not have every null-derived column.
bind_rows_fill <- function(xs) {
  xs <- xs[!vapply(xs, is.null, logical(1))]
  if (!length(xs)) return(NULL)
  cols <- unique(unlist(lapply(xs, names), use.names = FALSE))
  aligned <- lapply(xs, function(x) {
    missing <- setdiff(cols, names(x))
    for (nm in missing) x[[nm]] <- NA
    x[, cols, drop = FALSE]
  })
  out <- do.call(rbind, aligned)
  rownames(out) <- NULL
  out
}

feature_key <- function(x) paste(x$chr, x$N, x$k, sep = "|")

median_finite <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) stats::median(x) else NA_real_
}

min_finite <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) min(x) else NA_real_
}

max_finite <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) max(x) else NA_real_
}

# Read window metrics independently of the stability-selected peak tables.
read_window_map <- function(results_dir, dataset) {
  wdir <- file.path(results_dir, dataset, "window")
  files <- list.files(wdir, pattern = "^window_chr.*\\.tsv$", full.names = TRUE)
  if (!length(files)) return(NULL)
  rows <- lapply(files, function(path) {
    w <- read_tsv(path)
    if (!all(c("k", "window_power") %in% names(w))) return(NULL)
    if (!"chr" %in% names(w)) {
      w$chr <- sub("^window_chr(.*)\\.tsv$", "\\1", basename(path))
    }
    w$window_rank <- rank(-w$window_power, ties.method = "min", na.last = "keep")
    w$window_pct <- 100 * w$window_rank / sum(is.finite(w$window_power))
    # Coverage travels with the window map: it is the fraction of this
    # chromosome's grid the cohort observed, and the technical period floor is
    # derived from it. The consensus spectrum does not carry it.
    if (!"coverage" %in% names(w)) w$coverage <- NA_real_
    w[, c("chr", "k", "window_power", "window_rank", "window_pct", "coverage")]
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) NULL else do.call(rbind, rows)
}

discover_inputs <- function(opt) {
  ds <- opt$datasets
  if (is.null(ds)) {
    ds <- list.dirs(opt$results_dir, recursive = FALSE, full.names = FALSE)
    ds <- ds[dir.exists(file.path(opt$results_dir, ds, "consensus"))]
  }
  if (!length(ds)) abort("No dataset consensus directories under ", opt$results_dir)

  out <- list()
  for (id in ds) {
    cdir <- file.path(opt$results_dir, id, "consensus")
    files <- list.files(cdir, pattern = "^consensus_spectrum_.*\\.tsv$", full.names = TRUE)
    for (path in files) {
      cond <- sub("^consensus_spectrum_(.*)\\.tsv$", "\\1", basename(path))
      # Load every available condition even when only one output is requested.
      # Non-target conditions provide the within-cohort enrichment background.
      d <- read_tsv(path)
      needed <- c("chr", "N", "k", "period", "median_power_normalised",
                  "prevalence_rank", "plv", "mean_phase")
      miss <- setdiff(needed, names(d))
      if (length(miss)) abort(path, " lacks: ", paste(miss, collapse = ", "))
      d$dataset <- id
      d$condition <- cond
      d$source_condition <- cond

      # Current consensus files may not yet carry window annotations. Join the
      # raw chromosome window files when needed.
      if (!all(c("window_power", "window_rank", "window_pct") %in% names(d))) {
        wm <- read_window_map(opt$results_dir, id)
        if (!is.null(wm)) {
          m <- match(paste(d$chr, d$k), paste(wm$chr, wm$k))
          d$window_power <- wm$window_power[m]
          d$window_rank <- wm$window_rank[m]
          d$window_pct <- wm$window_pct[m]
          d$coverage <- wm$coverage[m]
        } else {
          d$window_power <- d$window_rank <- d$window_pct <- NA_real_
          d$coverage <- NA_real_
        }
      }

      sig_path <- file.path(cdir, paste0("signature_", cond, ".tsv"))
      d$selected_in_signature <- FALSE
      d$confirmed_in_cohort <- FALSE
      if (file.exists(sig_path)) {
        s <- read_tsv(sig_path)
        hit <- match(feature_key(d), feature_key(s))
        d$selected_in_signature <- !is.na(hit)
        if ("signature_class" %in% names(s)) {
          d$confirmed_in_cohort <- !is.na(hit) &
            s$signature_class[hit] %in% "confirmed"
        }
      }
      out[[paste(id, cond, sep = "\r")]] <- d
    }
  }
  if (!length(out)) abort("No consensus spectra matched the requested inputs")
  out
}

#' Add a broad healthy-liver view without erasing control provenance.
#'
#' OFF BY DEFAULT (--healthy-superclass yes to enable), because it merges three
#' classes this project deliberately keeps apart:
#'
#'   Control_disease_cohort   subjects outside a NAFLD cohort -- two of the ten
#'                            carry incidental fibrosis of stage 1 and 2, so it
#'                            is not a strictly healthy extreme
#'   Normal_histology         normal-looking biopsies within a disease series
#'   Control_external_study   controls of a different disease study, recruited,
#'                            sampled and processed under another protocol
#'
#' Taking a median across them treats them as cohorts of one class, which is
#' precisely the equivalence worth testing rather than assuming. Compare the
#' three spectra first; if they agree, the merge is a result and this flag
#' records that it was a decision.
add_healthy_superclass <- function(inputs) {
  healthy_labels <- c("Normal_histology", "Control_disease_cohort",
                      "Control_external_study")
  source_names <- names(inputs)[vapply(inputs, function(x)
    x$condition[1] %in% healthy_labels, logical(1))]
  for (nm in source_names) {
    d <- inputs[[nm]]
    d$source_condition <- d$condition
    d$condition <- "Healthy"
    inputs[[paste0(nm, "\rHealthy")]] <- d
  }
  inputs
}

add_within_cohort_effect <- function(inputs) {
  by_dataset <- split(names(inputs), vapply(inputs, function(x) x$dataset[1], character(1)))
  for (id in names(by_dataset)) {
    nms <- by_dataset[[id]]
    all_power <- lapply(nms, function(nm) {
      d <- inputs[[nm]]
      stats::setNames(d$median_power_normalised, feature_key(d))
    })
    names(all_power) <- vapply(nms, function(nm) inputs[[nm]]$condition[1], character(1))
    for (nm in nms) {
      d <- inputs[[nm]]
      cond <- d$condition[1]
      others <- all_power[names(all_power) != cond]
      key <- feature_key(d)
      if (!length(others)) {
        d$cohort_log2_enrichment <- NA_real_
      } else {
        bg <- vapply(key, function(k) median_finite(vapply(others, function(v)
          unname(v[k]), numeric(1))), numeric(1))
        eps <- 1e-12
        d$cohort_log2_enrichment <- log2((d$median_power_normalised + eps) / (bg + eps))
      }
      inputs[[nm]] <- d
    }
  }
  inputs
}

aggregate_condition <- function(inputs, condition, min_cohorts, window_cut,
                                min_prevalence = 0.5,
                                exclude_chromosomes = character(0),
                                signature_q = 0.05, min_period = "off",
                                period_margin = 2, margin_mode = "add",
                                min_period_biological = 0) {
  pieces <- inputs[vapply(inputs, function(x) identical(x$condition[1], condition), logical(1))]
  if (!length(pieces)) return(NULL)
  long <- bind_rows_fill(pieces)
  # Sex chromosomes and the mitochondrion are excluded by default. chrY carries
  # a handful of genes, so its grid is short and its spectrum unstable, and its
  # expression tracks the sex composition of a group -- which differs between
  # conditions and between cohorts. A component there reports who was recruited,
  # not what the disease does. (In this dataset the only "robust" component of
  # F0 was chrY with a period of 3.6 genes.)
  if (length(exclude_chromosomes)) {
    drop <- as.character(long$chr) %in% exclude_chromosomes
    if (any(drop)) long <- long[!drop, , drop = FALSE]
  }
  if (!nrow(long)) return(NULL)

  # --- period floor ------------------------------------------------------------
  #
  # Frequencies below a stated period are removed from the analysis AND from the
  # multiple-testing family. That is legitimate here for one specific reason:
  # the period of a frequency is N/k, a property of the annotation grid alone.
  # It does not depend on expression, on condition labels or on the null, so the
  # filter can be decided before a single spectrum is seen. Filtering by
  # enrichment or prevalence instead would invalidate the error control, since
  # those feed the same score the p-value comes from.
  #
  # Two floors, and the effective one is whichever is higher:
  #
  #   technical    2 / coverage, the Nyquist limit corrected for the gaps this
  #                chromosome actually has. Below it a "periodicity" is the
  #                sampling pattern, not the expression. A margin is added
  #                because the estimator is noisiest and window effects
  #                concentrate right at that boundary, so the useful floor sits
  #                above it rather than on it. Additive by default; multiplicative
  #                scales with each chromosome's own resolution, which matters
  #                here because coverage runs from 16% to 78%.
  #
  #   biological   a scale the mechanisms of interest could plausibly produce.
  #                Co-regulation domains, chromatin loops and replication timing
  #                act over tens to hundreds of genes; nothing known produces
  #                alternation gene by gene across a genome. This one is a scope
  #                decision and must be declared before looking at results, not
  #                chosen once the answer is visible.
  if (!identical(tolower(min_period), "off")) {
    cov <- if ("coverage" %in% names(long))
      suppressWarnings(as.numeric(long$coverage)) else rep(NA_real_, nrow(long))
    if (identical(tolower(min_period), "auto") && all(is.na(cov))) {
      abort("--min-period auto needs per-chromosome coverage, which comes from ",
            "the window maps. Run ./tsf window first, or give --min-period a ",
            "number.")
    }
    cov[!is.finite(cov) | cov <= 0] <- NA_real_
    technical <- if (identical(tolower(min_period), "auto")) {
      base <- 2 / cov
      if (identical(margin_mode, "mult")) base * (1 + period_margin)
      else base + period_margin
    } else {
      rep(as.numeric(min_period), nrow(long))
    }
    floor_period <- pmax(technical, min_period_biological, na.rm = TRUE)
    keep <- !is.finite(long$period) | long$period >= floor_period
    n_drop <- sum(!keep)
    if (n_drop) {
      message(sprintf("  period floor (%s%s, biological >= %g): removed %d of %d frequencies (%.0f%%)",
                      min_period,
                      if (identical(tolower(min_period), "auto"))
                        paste0(" ", margin_mode, " ", period_margin) else "",
                      min_period_biological, n_drop, nrow(long),
                      100 * n_drop / nrow(long)))
      long <- long[keep, , drop = FALSE]
    }
    if (!nrow(long)) {
      message("  every frequency fell below the period floor for ", condition)
      return(NULL)
    }
  }
  groups <- split(seq_len(nrow(long)), feature_key(long))

  rows <- lapply(groups, function(i) {
    x <- long[i, , drop = FALSE]
    phases <- x$mean_phase[is.finite(x$mean_phase)]
    phase_vec <- if (length(phases)) mean(exp(1i * phases)) else NA_complex_
    effect <- median_finite(x$cohort_log2_enrichment)
    data.frame(
      condition = condition,
      chr = as.character(x$chr[1]), N = as.integer(x$N[1]),
      k = as.integer(x$k[1]), freq = x$k[1] / x$N[1],
      period = x$N[1] / x$k[1],
      n_cohorts = length(unique(x$dataset)),
      cohorts = paste(sort(unique(x$dataset)), collapse = ","),
      source_conditions = paste(sort(unique(if ("source_condition" %in% names(x))
        as.character(x$source_condition) else as.character(x$condition))),
        collapse = ","),
      final_power = median_finite(x$median_power_normalised),
      final_power_min = min_finite(x$median_power_normalised),
      final_power_max = max_finite(x$median_power_normalised),
      condition_log2_enrichment = effect,
      condition_specific_power = if (is.finite(effect)) max(2^effect - 1, 0) else NA_real_,
      median_prevalence = median_finite(x$prevalence_rank),
      median_plv_within_cohort = median_finite(x$plv),
      phase_coherence_between_cohorts = if (is.finite(Re(phase_vec))) Mod(phase_vec) else NA_real_,
      mean_phase_between_cohorts = if (is.finite(Re(phase_vec))) Arg(phase_vec) else NA_real_,
      n_signature_cohorts = sum(x$selected_in_signature, na.rm = TRUE),
      n_confirmed_cohorts = sum(x$confirmed_in_cohort, na.rm = TRUE),
      signature_replication_fraction = mean(x$selected_in_signature, na.rm = TRUE),
      worst_window_pct = min_finite(x$window_pct),
      median_window_pct = median_finite(x$window_pct),
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$window_suspect <- is.finite(out$worst_window_pct) & out$worst_window_pct <= window_cut
  # --- a cut chosen by the data, not by a rank ---------------------------------
  #
  # `--top N` selects the N highest-scoring components, which is a display
  # convention pretending to be a criterion: it returns N whether the evidence
  # supports N, none, or a thousand.
  #
  # Combine the *pointwise* permutation p-values across independent cohorts and
  # apply BH once across frequencies.  p_null_fwer is already family-wise
  # adjusted inside each cohort; combining it and then applying BH again would
  # correct twice for the same frequency search and drive q_meta_null toward 1.
  # Keep p_null_fwer only as the separate single-cohort confirmation field.
  stouffer <- function(p, floor_p) {
    p <- p[is.finite(p)]
    if (!length(p)) return(NA_real_)
    p <- pmin(pmax(p, floor_p), 1 - floor_p)
    stats::pnorm(sum(stats::qnorm(1 - p)) / sqrt(length(p)), lower.tail = FALSE)
  }
  if ("p_null" %in% names(long)) {
    # The permutation floor is 1/(B+1); capping z at it stops one cohort's
    # impossible precision from carrying the combination.
    n_draws <- suppressWarnings(max(long$n_null, na.rm = TRUE))
    floor_p <- if (is.finite(n_draws) && n_draws > 0) 1 / (n_draws + 1) else 1e-4
    meta_p <- vapply(split(long$p_null, feature_key(long)),
                     stouffer, numeric(1), floor_p = floor_p)
    out$p_meta_null <- unname(meta_p[feature_key(out)])
    out$q_meta_null <- stats::p.adjust(out$p_meta_null, method = "BH")
    out$n_null_draws <- n_draws
    # Reachability, computed on the COMBINED p rather than on one cohort's.
    # Stouffer over k cohorts each at the permutation floor gives a far smaller
    # value than the floor itself -- with four cohorts at 1/100 the combination
    # reaches 1.6e-6, which BH over ten thousand frequencies can still use. The
    # obvious calculation (floor x n_frequencies) is the single-cohort one and
    # would have declared the cut unusable when it is not.
    k_min <- suppressWarnings(min(out$n_cohorts[out$n_cohorts >= min_cohorts],
                                  na.rm = TRUE))
    if (is.finite(k_min) && k_min >= 1) {
      best_meta <- stouffer(rep(floor_p, k_min), floor_p)
      reachable <- best_meta * nrow(out)
      msg <- paste0("  calibrated cut: ", n_draws, " draws, ", k_min,
                    " cohort(s), ", nrow(out), " frequencies -> smallest ",
                    "reachable BH q = ", signif(reachable, 2))
      if (reachable > signature_q) {
        message(msg, ". Nothing can pass at q <= ", signature_q,
                "; raise --n-null in the consensus stage.")
      } else {
        message(msg, " (usable).")
      }
    }
  } else {
    out$p_meta_null <- NA_real_
    out$q_meta_null <- NA_real_
    message("  NOTE: the consensus spectra carry no pointwise p_null, so the ",
            "signature falls back to the rank cut. Re-run ./tsf consensus ",
            "with a permutation null to get a calibrated one.")
  }

  out$eligible <- out$n_cohorts >= min_cohorts

  # With k cohorts, the mean of k unit phase vectors has expected modulus about
  # sqrt(pi)/(2*sqrt(k)) under random phases: 0.63 for two cohorts, 0.44 for
  # four. A fixed 0.8 therefore filters almost nothing when k is small -- here
  # 84% of frequencies clear it. The bar scales with k instead.
  k_cohorts <- max(out$n_cohorts, na.rm = TRUE)
  coherence_floor <- max(0.8, min(0.95, 2 * sqrt(pi) / (2 * sqrt(k_cohorts))))
  out$coherence_floor <- coherence_floor
  out$meta_score <- with(out,
    condition_specific_power * median_prevalence * median_plv_within_cohort *
      phase_coherence_between_cohorts)
  out$meta_score[!is.finite(out$meta_score)] <- 0
  out$signature_class <- "background"
  # Prevalence enters the rule, not only the score. Its median across
  # frequencies is 0 here, so without a floor a component can be called a
  # candidate while standing out in no sample at all.
  candidate <- out$eligible & out$condition_log2_enrichment > 0 &
    out$signature_replication_fraction >= 0.5 & !out$window_suspect &
    is.finite(out$median_prevalence) & out$median_prevalence >= min_prevalence
  out$signature_class[candidate] <- "candidate"

  # Calibrated membership of the final cross-cohort signature. A meta-analysis
  # exists precisely to recover effects that are individually underpowered but
  # reproducible. Requiring one cohort to be confirmed here would discard that
  # gain and make robust empty whenever every single study falls just short.
  out$signature_selected <- if (all(is.na(out$q_meta_null))) candidate else
    candidate & !is.na(out$q_meta_null) & out$q_meta_null <= signature_q
  robust <- candidate & out$signature_selected &
    out$n_signature_cohorts >= 2L &
    out$phase_coherence_between_cohorts >= coherence_floor
  out$signature_class[robust] <- "robust"
  out$signature_class[out$window_suspect & out$eligible] <- "window_suspect"
  out[order(match(out$chr, c(as.character(1:22), "X", "Y")), out$k), ]
}

open_plot_device <- function(path, width = 2400, height = 1800, res = 180) {
  # Compute nodes commonly have no DISPLAY. The default png(type = "Xlib")
  # therefore fails even though no interactive graphics are needed.
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(path, width = width, height = height, units = "px", res = res)
    return(path)
  }
  if (isTRUE(capabilities("cairo"))) {
    grDevices::png(path, width = width, height = height, res = res,
                   type = "cairo")
    return(path)
  }

  pdf_path <- sub("\\.png$", ".pdf", path, ignore.case = TRUE)
  warning("Neither ragg nor Cairo is available; writing PDF instead: ", pdf_path,
          call. = FALSE)
  grDevices::pdf(pdf_path, width = width / res, height = height / res,
                 onefile = TRUE)
  pdf_path
}

#' Signature-only view, one panel per chromosome.
#'
#' The full spectrum panel shows every frequency with the selected ones marked,
#' which is right for judging whether a component stands out from its
#' neighbourhood. This one drops the background entirely and plots only what was
#' selected, against period rather than k: period is comparable across
#' chromosomes, k is not, since a given k means a different scale on a
#' chromosome with a different N.
#'
#' A chromosome with nothing selected is drawn empty rather than skipped. That
#' absence is a result, and a grid with holes in it says so more plainly than a
#' figure that quietly renumbers the panels it has.
plot_signature_by_chromosome <- function(signature, path, condition,
                                         chr_order = c(as.character(1:22), "X", "Y")) {
  actual_path <- open_plot_device(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old), add = TRUE)
  graphics::par(mfrow = c(6, 4), mar = c(2.6, 3.2, 2.1, 0.8), oma = c(4, 4, 4, 1))

  amp <- if (nrow(signature)) pmax(signature$final_power, 1e-15) else numeric(0)
  ylim <- if (length(amp)) range(amp) * c(0.5, 2) else c(1e-6, 1)
  xlim <- if (nrow(signature)) range(signature$period) * c(0.8, 1.2) else c(1, 100)

  for (chr in chr_order) {
    s <- signature[as.character(signature$chr) == chr, , drop = FALSE]
    graphics::plot(NA, xlim = xlim, ylim = ylim, log = "xy",
                   xlab = "period (genes)", ylab = "power",
                   main = paste0("chr", chr, if (nrow(s)) paste0(" (", nrow(s), ")") else ""),
                   col.main = if (nrow(s)) "black" else "grey60")
    if (!nrow(s)) next
    col <- ifelse(s$signature_class == "robust", "#B2182B", "#2166AC")
    graphics::segments(s$period, ylim[1], s$period, pmax(s$final_power, 1e-15),
                       col = col, lwd = 1.2)
    graphics::points(s$period, pmax(s$final_power, 1e-15), pch = 19, cex = 0.8, col = col)
  }
  graphics::mtext(paste("Selected components:", condition), outer = TRUE,
                  side = 3, line = 1.5, cex = 1.35, font = 2)
  graphics::mtext("period (genes per cycle, log scale)", outer = TRUE,
                  side = 1, line = 1.4)
  graphics::mtext("median normalised power across cohorts", outer = TRUE,
                  side = 2, line = 2)
  graphics::mtext("red: robust    blue: candidate    grey panel: no components selected",
                  outer = TRUE, side = 1, line = 2.6, cex = 0.85)
  invisible(actual_path)
}

#' The whole genome on one axis, chromosomes concatenated in order.
#'
#' Useful for seeing where the selected components sit relative to everything
#' else at a glance. The x axis is a running index, NOT a genomic coordinate and
#' not comparable between chromosomes: each contributes as many positions as it
#' has frequencies, so a long chromosome occupies more of the axis than a short
#' one regardless of its physical size.
plot_genome_wide <- function(tab, signature, path, condition,
                             chr_order = c(as.character(1:22), "X", "Y")) {
  actual_path <- open_plot_device(path, width = 3600, height = 1200)
  on.exit(grDevices::dev.off(), add = TRUE)
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old), add = TRUE)
  graphics::par(mar = c(4.2, 4.5, 3.2, 1))

  present <- chr_order[chr_order %in% as.character(tab$chr)]
  tab <- tab[order(match(as.character(tab$chr), present), tab$k), , drop = FALSE]
  tab$pos <- seq_len(nrow(tab))
  bounds <- cumsum(as.integer(table(factor(as.character(tab$chr), levels = present))))
  mids <- c(0, utils::head(bounds, -1)) + diff(c(0, bounds)) / 2

  y <- pmax(tab$final_power, 1e-15)
  graphics::plot(tab$pos, y, type = "h", log = "y", lwd = 0.4, col = "grey78",
                 xlab = "", ylab = "power (median across cohorts)",
                 main = paste("Full spectrum:", condition), xaxt = "n")
  graphics::abline(v = utils::head(bounds, -1), col = "grey90", lwd = 0.6)
  graphics::axis(1, at = mids, labels = present, tick = FALSE, las = 1,
                 cex.axis = 0.8, line = -0.5)

  if (nrow(signature)) {
    m <- match(feature_key(signature), feature_key(tab))
    ok <- !is.na(m)
    col <- ifelse(signature$signature_class[ok] == "robust", "#B2182B", "#2166AC")
    graphics::points(tab$pos[m[ok]], pmax(signature$final_power[ok], 1e-15),
                     pch = 19, cex = 0.7, col = col)
  }
  graphics::mtext("chromosome (axis: running frequency index, not a genomic coordinate)",
                  side = 1, line = 3, cex = 0.9)
  invisible(actual_path)
}

plot_condition <- function(tab, signature, path, condition) {
  actual_path <- open_plot_device(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old), add = TRUE)
  graphics::par(mfrow = c(6, 4), mar = c(2.5, 3.2, 2.1, 0.8),
                oma = c(3, 4, 4, 1))
  chr_order <- c(as.character(1:22), "X", "Y")
  for (chr in chr_order) {
    x <- tab[as.character(tab$chr) == chr, , drop = FALSE]
    if (!nrow(x)) { graphics::plot.new(); next }
    y <- pmax(x$final_power, 1e-15)
    graphics::plot(x$k, y, type = "l", log = "y", lwd = 0.7,
                   xlab = "k", ylab = "power", main = paste0("chr", chr),
                   col = "grey25")
    s <- signature[as.character(signature$chr) == chr, , drop = FALSE]
    if (nrow(s)) {
      graphics::points(s$k, pmax(s$final_power, 1e-15), pch = 19, cex = 0.65,
                       col = ifelse(s$signature_class == "robust", "#B2182B", "#2166AC"))
    }
  }
  graphics::mtext(paste("Espectro final de potencia —", condition), outer = TRUE,
                  side = 3, line = 1.5, cex = 1.35, font = 2)
  graphics::mtext("Frecuencia discreta dentro de cada cromosoma", outer = TRUE,
                  side = 1, line = 1.2)
  graphics::mtext("Mediana intercohorte de potencia normalizada (escala log)",
                  outer = TRUE, side = 2, line = 2)
  invisible(actual_path)
}

main <- function() {
  opt <- parse_args(commandArgs(trailingOnly = TRUE))
  dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
  inputs <- add_within_cohort_effect(discover_inputs(opt))
  if (opt$healthy_superclass) {
    message("NOTE: --healthy-superclass merges three healthy classes whose ",
            "equivalence has not been demonstrated. source_conditions records ",
            "what went in.")
    inputs <- add_healthy_superclass(inputs)
  }
  conditions <- unique(vapply(inputs, function(x) x$condition[1], character(1)))
  preferred <- c("Healthy", "Normal_histology", "Control_disease_cohort",
                 "Control_external_study", "F0", "F1", "F2", "F3", "F4")
  conditions <- c(intersect(preferred, conditions), setdiff(sort(conditions), preferred))
  if (!is.null(opt$conditions)) {
    conditions <- opt$conditions[opt$conditions %in% conditions]
  }

  manifest <- list()
  for (cond in conditions) {
    tab <- aggregate_condition(inputs, cond, opt$min_cohorts, opt$window_cut,
                               opt$min_prevalence, opt$exclude_chromosomes,
                               opt$signature_q, opt$min_period,
                               opt$period_margin, opt$margin_mode,
                               opt$min_period_biological)
    if (is.null(tab)) next
    full_path <- file.path(opt$out_dir, paste0("condition_spectrum_", cond, ".tsv"))
    write_tsv(tab, full_path)

    # Avoid NA indices creating all-NA rows and invalid plotting ranges.
    sig <- tab[tab$eligible %in% TRUE &
                 tab$signature_selected %in% TRUE, , drop = FALSE]
    sig <- sig[order(-sig$meta_score, -sig$final_power), , drop = FALSE]
    # --top is a safety cap on file size, not the criterion. When it bites, say
    # so, because a truncated signature and a signature are different objects.
    if (!is.na(opt$top) && nrow(sig) > opt$top) {
      message("  NOTE: ", nrow(sig), " components passed the calibrated cut for ",
              cond, "; writing the ", opt$top, " highest-scoring. Raise --top ",
              "to keep them all.")
      sig <- utils::head(sig, opt$top)
    }
    sig_path <- file.path(opt$out_dir, paste0("condition_signature_", cond, ".tsv"))
    write_tsv(sig, sig_path)
    png_path <- file.path(opt$out_dir, paste0("power_spectrum_", cond, ".png"))
    plot_path <- plot_condition(tab, sig, png_path, cond)
    sig_png <- plot_signature_by_chromosome(
      sig, file.path(opt$out_dir, paste0("signature_peaks_", cond, ".png")), cond)
    genome_png <- plot_genome_wide(
      tab, sig, file.path(opt$out_dir, paste0("genome_spectrum_", cond, ".png")), cond)

    if (sum(tab$n_confirmed_cohorts, na.rm = TRUE) == 0L) {
      has_meta_null <- any(is.finite(tab$q_meta_null))
      n_draws_here <- suppressWarnings(max(tab$n_null_draws, na.rm = TRUE))
      if (has_meta_null && is.finite(n_draws_here) && n_draws_here > 0) {
        message("  NOTE: no component retained after chromosome/period filters ",
                "was confirmed by one cohort alone for ", cond, ". The final ",
                "robust call uses the calibrated cross-cohort q_meta_null; ",
                n_draws_here, " null draws are already present.")
      } else {
        message("  NOTE: no usable permutation null is present for ", cond,
                ". Re-run ./tsf consensus with --n-null 199 or more.")
      }
    }
    manifest[[cond]] <- data.frame(
      condition = cond, n_cohorts_max = max(tab$n_cohorts),
      coherence_floor = tab$coherence_floor[1],
      min_prevalence = opt$min_prevalence,
      window_cut_pct = opt$window_cut,
      excluded_chromosomes = paste(opt$exclude_chromosomes, collapse = ","),
      healthy_superclass = opt$healthy_superclass,
      n_frequencies = nrow(tab), n_robust = sum(tab$signature_class == "robust"),
      n_candidates = sum(tab$signature_class == "candidate"),
      n_signature_calibrated = sum(tab$signature_selected, na.rm = TRUE),
      plot_signature = basename(sig_png),
      plot_genome = basename(genome_png),
      signature_q = opt$signature_q,
      min_period = opt$min_period,
      period_margin = opt$period_margin,
      margin_mode = opt$margin_mode,
      min_period_biological = opt$min_period_biological,
      n_window_suspect = sum(tab$signature_class == "window_suspect"),
      plot_file = basename(plot_path),
      stringsAsFactors = FALSE)
    messagef("%-24s cohorts=%d frequencies=%d robust=%d candidates=%d",
             cond, max(tab$n_cohorts), nrow(tab),
             sum(tab$signature_class == "robust"),
             sum(tab$signature_class == "candidate"))
  }
  write_tsv(do.call(rbind, manifest), file.path(opt$out_dir, "condition_library_manifest.tsv"))
  message("Wrote final condition spectra to: ", normalizePath(opt$out_dir))
}

main()
