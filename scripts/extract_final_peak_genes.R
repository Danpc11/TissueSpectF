#!/usr/bin/env Rscript
# extract_final_peak_genes.R
#
# Map final TissueSpectF condition peaks and invariant peaks back to genes on
# the canonical chromosome gene grid, then build chromosome-normalised gene
# participation summaries.
#
# A Fourier peak is a sinusoidal mode spanning the chromosome, not a discrete
# genomic interval. For each gene j:
#
#   loading_j = cos(2*pi*k*(grid_index_j - 1)/N + phase)
#
# The sign indicates phase direction; magnitude=abs(loading) indicates how close
# the gene lies to a crest/trough of that mode.
#
# New normalised gene metrics:
#   chromosome_peak_count
#   peak_participation_fraction
#   mean_abs_loading
#   phase_consistency
#   normalized_weighted_score
#
# Outputs:
#   peak_genes/
#     manifest.tsv
#     all_peak_genes_long.tsv
#     conditions/
#       condition_peak_genes_long.tsv
#       condition_gene_summary_all.tsv
#       condition_gene_summary_robust.tsv
#       condition_gene_summary_candidate.tsv
#       by_condition/
#         <condition>_all_genes.tsv
#         <condition>_robust_genes.tsv
#         <condition>_candidate_genes.tsv
#       <condition>/
#         peak_chr..._N..._k....tsv
#         peak_chr..._N..._k..._top.tsv
#     invariants/
#       invariant_peak_genes_long.tsv
#       invariant_gene_summary.tsv
#       invariant_core_gene_summary.tsv
#       invariant_core_top_genes.tsv
#       <class>/
#         peak_chr..._N..._k....tsv
#         peak_chr..._N..._k..._top.tsv

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
abort <- function(...) stop(paste0(...), call. = FALSE)
messagef <- function(...) message(sprintf(...))

read_tsv <- function(path) {
  utils::read.delim(
    path,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    x,
    path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    na = "NA"
  )
}

bind_rows_fill <- function(xs) {
  xs <- xs[!vapply(xs, is.null, logical(1))]
  if (!length(xs)) return(NULL)

  cols <- unique(unlist(lapply(xs, names), use.names = FALSE))
  aligned <- lapply(xs, function(x) {
    miss <- setdiff(cols, names(x))
    for (nm in miss) x[[nm]] <- NA
    x[, cols, drop = FALSE]
  })

  out <- do.call(rbind, aligned)
  rownames(out) <- NULL
  out
}

parse_args <- function(x) {
  out <- list(
    library_dir = "results_gencode_v3/library_domains",
    out_dir = NULL,
    genes_file = NULL,
    interim_dir = NULL,
    datasets = NULL,
    conditions = NULL,
    invariant_classes = "core_invariant",
    top_genes_per_peak = 50L,
    top_genes_summary = 200L,
    min_abs_loading = 0,
    include_candidates = "yes"
  )

  i <- 1L
  while (i <= length(x)) {
    a <- x[[i]]
    if (!startsWith(a, "--")) abort("Unexpected argument: ", a)

    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    key <- gsub("-", "_", kv[[1]])

    if (length(kv) > 1L) {
      value <- paste(kv[-1], collapse = "=")
    } else {
      i <- i + 1L
      if (i > length(x)) abort("Missing value for --", key)
      value <- x[[i]]
    }

    if (!key %in% names(out)) abort("Unknown option --", key)
    out[[key]] <- value
    i <- i + 1L
  }

  split_csv <- function(z) {
    if (is.null(z) || !length(z) || !nzchar(z)) return(NULL)
    trimws(strsplit(z, ",", fixed = TRUE)[[1]])
  }

  out$out_dir <- out$out_dir %||%
    file.path(out$library_dir, "peak_genes")

  out$datasets <- split_csv(out$datasets)
  out$conditions <- split_csv(out$conditions)
  out$invariant_classes <- split_csv(out$invariant_classes) %||%
    "core_invariant"

  out$top_genes_per_peak <- as.integer(out$top_genes_per_peak)
  out$top_genes_summary <- as.integer(out$top_genes_summary)
  out$min_abs_loading <- as.numeric(out$min_abs_loading)
  out$include_candidates <- tolower(out$include_candidates) %in%
    c("yes", "true", "1")

  if (!is.finite(out$top_genes_per_peak) || out$top_genes_per_peak < 1L) {
    abort("--top-genes-per-peak must be >= 1")
  }
  if (!is.finite(out$top_genes_summary) || out$top_genes_summary < 1L) {
    abort("--top-genes-summary must be >= 1")
  }
  if (!is.finite(out$min_abs_loading) ||
      out$min_abs_loading < 0 ||
      out$min_abs_loading > 1) {
    abort("--min-abs-loading must be in [0,1]")
  }

  out
}

normalise_chr <- function(x) {
  x <- as.character(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x[x %in% c("M", "Mt", "mt")] <- "MT"
  x
}

discover_genes_file <- function(opt) {
  if (!is.null(opt$genes_file)) {
    if (!file.exists(opt$genes_file)) {
      abort("genes file does not exist: ", opt$genes_file)
    }
    return(opt$genes_file)
  }

  candidates <- character(0)

  if (!is.null(opt$interim_dir) && dir.exists(opt$interim_dir)) {
    if (!is.null(opt$datasets)) {
      candidates <- file.path(opt$interim_dir, opt$datasets, "genes.tsv")
      candidates <- candidates[file.exists(candidates)]
    }

    if (!length(candidates)) {
      candidates <- list.files(
        opt$interim_dir,
        pattern = "^genes\\.tsv$",
        recursive = TRUE,
        full.names = TRUE
      )
    }
  }

  if (!length(candidates)) {
    roots <- c("interim", "interim_gencode_v2", "interim_gencode_v3")
    for (root in roots) {
      if (!dir.exists(root)) next
      hit <- list.files(
        root,
        pattern = "^genes\\.tsv$",
        recursive = TRUE,
        full.names = TRUE
      )
      candidates <- c(candidates, hit)
    }
  }

  candidates <- unique(candidates[file.exists(candidates)])

  if (!length(candidates)) {
    abort(
      "Could not find genes.tsv. Supply --genes-file PATH or ",
      "--interim-dir PATH."
    )
  }

  message("Using canonical gene grid: ", candidates[[1]])
  candidates[[1]]
}

prepare_genes <- function(path) {
  g <- read_tsv(path)

  required <- c("gene_id", "chr", "grid_index", "grid_N")
  miss <- setdiff(required, names(g))
  if (length(miss)) {
    abort(path, " lacks required columns: ", paste(miss, collapse = ", "))
  }

  if (!"gene_name" %in% names(g)) g$gene_name <- NA_character_
  if (!"start" %in% names(g)) g$start <- NA_real_

  g$chr <- normalise_chr(g$chr)
  g$grid_index <- suppressWarnings(as.integer(g$grid_index))
  g$grid_N <- suppressWarnings(as.integer(g$grid_N))

  bad <- !is.finite(g$grid_index) | !is.finite(g$grid_N) |
    g$grid_index < 1 | g$grid_N < 1
  if (any(bad)) abort("Invalid grid_index/grid_N values in ", path)

  n_by_chr <- split(g$grid_N, g$chr)
  bad_chr <- names(n_by_chr)[vapply(
    n_by_chr,
    function(z) length(unique(z)) != 1L,
    logical(1)
  )]
  if (length(bad_chr)) {
    abort(
      "Multiple grid_N values found within chromosome(s): ",
      paste(bad_chr, collapse = ", ")
    )
  }

  g <- g[
    order(match(g$chr, c(as.character(1:22), "X", "Y", "MT")),
          g$grid_index),
    ,
    drop = FALSE
  ]
  rownames(g) <- NULL
  g
}

feature_key <- function(chr, N, k) {
  paste(normalise_chr(chr), as.integer(N), as.integer(k), sep = "|")
}

safe_file_label <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

choose_phase_column <- function(peaks, peak_type) {
  preferred <- if (peak_type == "condition") {
    c("mean_phase_between_cohorts", "mean_phase", "phase")
  } else {
    c("mean_phase", "mean_phase_between_cohorts", "phase")
  }

  hit <- preferred[preferred %in% names(peaks)]
  if (!length(hit)) {
    abort(
      "No usable phase column found for ", peak_type,
      " peaks. Expected one of: ",
      paste(preferred, collapse = ", ")
    )
  }
  hit[[1]]
}

choose_peak_score_column <- function(peaks, peak_type) {
  preferred <- if (peak_type == "condition") {
    c("meta_score", "condition_specific_power", "final_power")
  } else {
    c("invariant_score", "median_power")
  }

  hit <- preferred[preferred %in% names(peaks)]
  if (!length(hit)) return(NULL)
  hit[[1]]
}

choose_power_column <- function(peaks, peak_type) {
  preferred <- if (peak_type == "condition") {
    c("final_power", "power", "median_power")
  } else {
    c("median_power", "final_power", "power")
  }

  hit <- preferred[preferred %in% names(peaks)]
  if (!length(hit)) return(NULL)
  hit[[1]]
}

project_peak_to_genes <- function(
  peak_row,
  genes,
  peak_type,
  source_label,
  class_label,
  min_abs_loading = 0
) {
  chr_now <- normalise_chr(peak_row$chr[[1]])
  N <- as.integer(peak_row$N[[1]])
  k <- as.integer(peak_row$k[[1]])

  phase_col <- choose_phase_column(peak_row, peak_type)
  phase <- suppressWarnings(as.numeric(peak_row[[phase_col]][[1]]))

  if (!is.finite(phase)) {
    warning(
      "Skipping peak ", chr_now, " N=", N, " k=", k,
      " because phase is not finite.",
      call. = FALSE
    )
    return(NULL)
  }

  g <- genes[genes$chr == chr_now, , drop = FALSE]

  if (!nrow(g)) {
    warning("No genes on chromosome ", chr_now, call. = FALSE)
    return(NULL)
  }

  grid_N <- unique(g$grid_N)
  if (length(grid_N) != 1L || grid_N[[1]] != N) {
    abort(
      "Canonical grid mismatch for chr", chr_now,
      ": peak N=", N,
      ", genes.tsv grid_N=", paste(grid_N, collapse = ",")
    )
  }

  n <- g$grid_index - 1L
  loading <- cos(2 * pi * k * n / N + phase)
  magnitude <- abs(loading)

  score_col <- choose_peak_score_column(peak_row, peak_type)
  peak_score <- if (!is.null(score_col)) {
    suppressWarnings(as.numeric(peak_row[[score_col]][[1]]))
  } else {
    1
  }
  if (!is.finite(peak_score)) peak_score <- 0

  power_col <- choose_power_column(peak_row, peak_type)
  peak_power <- if (!is.null(power_col)) {
    suppressWarnings(as.numeric(peak_row[[power_col]][[1]]))
  } else {
    NA_real_
  }

  weighted_score <- magnitude * max(peak_score, 0)

  out <- data.frame(
    peak_type = peak_type,
    source = source_label,
    peak_class = class_label,
    chr = chr_now,
    N = N,
    k = k,
    freq = if ("freq" %in% names(peak_row)) peak_row$freq[[1]] else k / N,
    period = if ("period" %in% names(peak_row)) peak_row$period[[1]] else N / k,
    phase = phase,
    phase_column = phase_col,
    peak_power = peak_power,
    peak_score = peak_score,
    gene_id = as.character(g$gene_id),
    gene_name = as.character(g$gene_name),
    gene_start = suppressWarnings(as.numeric(g$start)),
    grid_index = g$grid_index,
    grid_N = g$grid_N,
    loading = loading,
    sign = ifelse(
      loading > 0,
      "positive",
      ifelse(loading < 0, "negative", "zero")
    ),
    magnitude = magnitude,
    weighted_score = weighted_score,
    stringsAsFactors = FALSE
  )

  out$rank_magnitude <- rank(
    -out$magnitude,
    ties.method = "min",
    na.last = "keep"
  )

  out$rank_weighted <- rank(
    -out$weighted_score,
    ties.method = "min",
    na.last = "keep"
  )

  if (min_abs_loading > 0) {
    out <- out[
      is.finite(out$magnitude) &
        out$magnitude >= min_abs_loading,
      ,
      drop = FALSE
    ]
  }

  out
}

condition_files <- function(library_dir, requested = NULL) {
  files <- list.files(
    library_dir,
    pattern = "^condition_signature_.*\\.tsv$",
    full.names = TRUE
  )

  if (!length(files)) return(character(0))

  labels <- sub(
    "^condition_signature_(.*)\\.tsv$",
    "\\1",
    basename(files)
  )
  names(files) <- labels

  if (!is.null(requested)) {
    files <- files[names(files) %in% requested]
  }

  files
}

read_invariant_peaks <- function(library_dir, invariant_classes) {
  inv_dir <- file.path(library_dir, "invariants")

  core_path <- file.path(inv_dir, "invariant_core.tsv")
  components_path <- file.path(inv_dir, "invariant_components.tsv")
  spectrum_path <- file.path(inv_dir, "invariant_spectrum.tsv")

  if (file.exists(core_path) &&
      length(invariant_classes) == 1L &&
      identical(invariant_classes, "core_invariant")) {
    x <- read_tsv(core_path)
    if (!"invariant_class" %in% names(x)) {
      x$invariant_class <- "core_invariant"
    }
    return(x)
  }

  path <- if (file.exists(components_path)) {
    components_path
  } else if (file.exists(spectrum_path)) {
    spectrum_path
  } else {
    return(NULL)
  }

  x <- read_tsv(path)

  if ("invariant_class" %in% names(x)) {
    x <- x[
      x$invariant_class %in% invariant_classes,
      ,
      drop = FALSE
    ]
  }

  x
}

write_peak_outputs <- function(
  peaks,
  genes,
  peak_type,
  source_label,
  out_dir,
  top_n,
  min_abs_loading
) {
  if (is.null(peaks) || !nrow(peaks)) return(NULL)

  needed <- c("chr", "N", "k")
  miss <- setdiff(needed, names(peaks))
  if (length(miss)) {
    abort(
      source_label, " peak table lacks: ",
      paste(miss, collapse = ", ")
    )
  }

  rows <- vector("list", nrow(peaks))

  for (i in seq_len(nrow(peaks))) {
    p <- peaks[i, , drop = FALSE]

    class_label <- if (peak_type == "condition") {
      if ("signature_class" %in% names(p)) {
        as.character(p$signature_class[[1]])
      } else {
        "selected"
      }
    } else {
      if ("invariant_class" %in% names(p)) {
        as.character(p$invariant_class[[1]])
      } else {
        "core_invariant"
      }
    }

    mapped <- project_peak_to_genes(
      peak_row = p,
      genes = genes,
      peak_type = peak_type,
      source_label = source_label,
      class_label = class_label,
      min_abs_loading = min_abs_loading
    )

    if (is.null(mapped) || !nrow(mapped)) next

    mapped <- mapped[
      order(-mapped$magnitude, mapped$grid_index),
      ,
      drop = FALSE
    ]

    base <- sprintf(
      "peak_chr%s_N%s_k%s",
      safe_file_label(mapped$chr[[1]]),
      mapped$N[[1]],
      mapped$k[[1]]
    )

    write_tsv(
      mapped,
      file.path(out_dir, paste0(base, ".tsv"))
    )

    write_tsv(
      utils::head(mapped, top_n),
      file.path(out_dir, paste0(base, "_top.tsv"))
    )

    rows[[i]] <- mapped
  }

  bind_rows_fill(rows)
}

build_peak_inventory <- function(long) {
  if (is.null(long) || !nrow(long)) return(NULL)

  unique(
    long[
      ,
      c(
        "peak_type", "source", "peak_class",
        "chr", "N", "k", "period", "peak_score"
      ),
      drop = FALSE
    ]
  )
}

summarise_genes <- function(long) {
  if (is.null(long) || !nrow(long)) return(NULL)

  peak_inventory <- build_peak_inventory(long)

  chr_peak_counts <- aggregate(
    k ~ source + peak_class + chr,
    data = peak_inventory,
    FUN = length
  )
  names(chr_peak_counts)[names(chr_peak_counts) == "k"] <-
    "chromosome_peak_count"

  key <- paste(long$source, long$peak_class, long$gene_id, sep = "\r")
  groups <- split(seq_len(nrow(long)), key)

  rows <- lapply(groups, function(i) {
    x <- long[i, , drop = FALSE]

    pos <- sum(x$sign == "positive", na.rm = TRUE)
    neg <- sum(x$sign == "negative", na.rm = TRUE)
    signed_n <- pos + neg

    phase_consistency <- if (signed_n > 0) {
      max(pos, neg) / signed_n
    } else {
      NA_real_
    }

    best_i <- which.max(x$weighted_score)
    if (!length(best_i) || !is.finite(x$weighted_score[best_i])) {
      best_i <- which.max(x$magnitude)
    }

    data.frame(
      peak_type = x$peak_type[[1]],
      source = x$source[[1]],
      peak_class = x$peak_class[[1]],
      gene_id = x$gene_id[[1]],
      gene_name = x$gene_name[[1]],
      chr = x$chr[[1]],
      n_peaks = nrow(x),
      n_positive = pos,
      n_negative = neg,
      sign_balance = if (signed_n > 0) (pos - neg) / signed_n else NA_real_,
      phase_consistency = phase_consistency,
      max_magnitude = max(x$magnitude, na.rm = TRUE),
      mean_abs_loading = mean(x$magnitude, na.rm = TRUE),
      sum_weighted_score = sum(x$weighted_score, na.rm = TRUE),
      mean_weighted_score = mean(x$weighted_score, na.rm = TRUE),
      max_weighted_score = max(x$weighted_score, na.rm = TRUE),
      best_peak_k = x$k[best_i],
      best_peak_period = x$period[best_i],
      best_peak_class = x$peak_class[best_i],
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)

  m <- match(
    paste(out$source, out$peak_class, out$chr, sep = "\r"),
    paste(
      chr_peak_counts$source,
      chr_peak_counts$peak_class,
      chr_peak_counts$chr,
      sep = "\r"
    )
  )

  out$chromosome_peak_count <- chr_peak_counts$chromosome_peak_count[m]

  out$peak_participation_fraction <- ifelse(
    is.finite(out$chromosome_peak_count) &
      out$chromosome_peak_count > 0,
    out$n_peaks / out$chromosome_peak_count,
    NA_real_
  )

  out$normalized_weighted_score <- ifelse(
    is.finite(out$chromosome_peak_count) &
      out$chromosome_peak_count > 0,
    out$sum_weighted_score / out$chromosome_peak_count,
    NA_real_
  )

  out <- out[
    order(
      out$source,
      out$peak_class,
      -out$normalized_weighted_score,
      -out$phase_consistency,
      -out$mean_abs_loading
    ),
    ,
    drop = FALSE
  ]

  rownames(out) <- NULL
  out
}

summarise_across_classes <- function(summary_table) {
  if (is.null(summary_table) || !nrow(summary_table)) return(NULL)

  key <- paste(summary_table$source, summary_table$gene_id, sep = "\r")
  groups <- split(seq_len(nrow(summary_table)), key)

  rows <- lapply(groups, function(i) {
    x <- summary_table[i, , drop = FALSE]

    data.frame(
      peak_type = x$peak_type[[1]],
      source = x$source[[1]],
      gene_id = x$gene_id[[1]],
      gene_name = x$gene_name[[1]],
      chr = x$chr[[1]],
      classes = paste(sort(unique(x$peak_class)), collapse = ","),
      n_classes = length(unique(x$peak_class)),
      n_peaks = sum(x$n_peaks, na.rm = TRUE),
      chromosome_peak_count = sum(x$chromosome_peak_count, na.rm = TRUE),
      peak_participation_fraction = if (
        sum(x$chromosome_peak_count, na.rm = TRUE) > 0
      ) {
        sum(x$n_peaks, na.rm = TRUE) /
          sum(x$chromosome_peak_count, na.rm = TRUE)
      } else {
        NA_real_
      },
      mean_abs_loading = weighted.mean(
        x$mean_abs_loading,
        w = pmax(x$n_peaks, 1),
        na.rm = TRUE
      ),
      phase_consistency = weighted.mean(
        x$phase_consistency,
        w = pmax(x$n_peaks, 1),
        na.rm = TRUE
      ),
      normalized_weighted_score = sum(
        x$normalized_weighted_score,
        na.rm = TRUE
      ),
      max_weighted_score = max(
        x$max_weighted_score,
        na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  out <- out[
    order(
      out$source,
      -out$normalized_weighted_score,
      -out$phase_consistency,
      -out$mean_abs_loading
    ),
    ,
    drop = FALSE
  ]
  rownames(out) <- NULL
  out
}

write_condition_summaries <- function(condition_long, out_dir) {
  if (is.null(condition_long) || !nrow(condition_long)) return(NULL)

  summary_by_class <- summarise_genes(condition_long)
  summary_all <- summarise_across_classes(summary_by_class)

  write_tsv(
    condition_long,
    file.path(out_dir, "condition_peak_genes_long.tsv")
  )

  write_tsv(
    summary_all,
    file.path(out_dir, "condition_gene_summary_all.tsv")
  )

  for (cl in c("robust", "candidate")) {
    sub <- summary_by_class[
      summary_by_class$peak_class == cl,
      ,
      drop = FALSE
    ]
    write_tsv(
      sub,
      file.path(out_dir, paste0("condition_gene_summary_", cl, ".tsv"))
    )
  }

  by_condition_dir <- file.path(out_dir, "by_condition")
  dir.create(by_condition_dir, recursive = TRUE, showWarnings = FALSE)

  for (cond in unique(summary_all$source)) {
    zall <- summary_all[
      summary_all$source == cond,
      ,
      drop = FALSE
    ]
    write_tsv(
      zall,
      file.path(
        by_condition_dir,
        paste0(safe_file_label(cond), "_all_genes.tsv")
      )
    )

    for (cl in c("robust", "candidate")) {
      z <- summary_by_class[
        summary_by_class$source == cond &
          summary_by_class$peak_class == cl,
        ,
        drop = FALSE
      ]
      write_tsv(
        z,
        file.path(
          by_condition_dir,
          paste0(
            safe_file_label(cond),
            "_",
            cl,
            "_genes.tsv"
          )
        )
      )
    }
  }

  invisible(
    list(
      by_class = summary_by_class,
      all = summary_all
    )
  )
}

write_invariant_summaries <- function(
  invariant_long,
  out_dir,
  top_n
) {
  if (is.null(invariant_long) || !nrow(invariant_long)) return(NULL)

  summary_by_class <- summarise_genes(invariant_long)
  summary_all <- summarise_across_classes(summary_by_class)

  write_tsv(
    invariant_long,
    file.path(out_dir, "invariant_peak_genes_long.tsv")
  )

  write_tsv(
    summary_all,
    file.path(out_dir, "invariant_gene_summary.tsv")
  )

  core <- summary_by_class[
    summary_by_class$peak_class == "core_invariant",
    ,
    drop = FALSE
  ]

  write_tsv(
    core,
    file.path(out_dir, "invariant_core_gene_summary.tsv")
  )

  if (nrow(core)) {
    core_top <- core[
      order(
        -core$normalized_weighted_score,
        -core$phase_consistency,
        -core$mean_abs_loading
      ),
      ,
      drop = FALSE
    ]

    write_tsv(
      utils::head(core_top, top_n),
      file.path(out_dir, "invariant_core_top_genes.tsv")
    )
  } else {
    write_tsv(
      core,
      file.path(out_dir, "invariant_core_top_genes.tsv")
    )
  }

  invisible(
    list(
      by_class = summary_by_class,
      all = summary_all,
      core = core
    )
  )
}

main <- function() {
  opt <- parse_args(commandArgs(trailingOnly = TRUE))

  if (!dir.exists(opt$library_dir)) {
    abort("Library directory does not exist: ", opt$library_dir)
  }

  dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

  genes_path <- discover_genes_file(opt)
  genes <- prepare_genes(genes_path)

  condition_long <- list()

  cfiles <- condition_files(
    opt$library_dir,
    opt$conditions
  )

  if (length(cfiles)) {
    for (cond in names(cfiles)) {
      peaks <- read_tsv(cfiles[[cond]])

      if (!opt$include_candidates &&
          "signature_class" %in% names(peaks)) {
        peaks <- peaks[
          peaks$signature_class == "robust",
          ,
          drop = FALSE
        ]
      }

      messagef("CONDITION %-20s peaks=%d", cond, nrow(peaks))

      out_dir <- file.path(
        opt$out_dir,
        "conditions",
        safe_file_label(cond)
      )

      condition_long[[cond]] <- write_peak_outputs(
        peaks = peaks,
        genes = genes,
        peak_type = "condition",
        source_label = cond,
        out_dir = out_dir,
        top_n = opt$top_genes_per_peak,
        min_abs_loading = opt$min_abs_loading
      )
    }
  } else {
    message("No condition_signature_*.tsv files found.")
  }

  condition_long <- bind_rows_fill(condition_long)

  condition_summary <- NULL
  if (!is.null(condition_long) && nrow(condition_long)) {
    condition_summary <- write_condition_summaries(
      condition_long,
      file.path(opt$out_dir, "conditions")
    )
  }

  invariant_peaks <- read_invariant_peaks(
    opt$library_dir,
    opt$invariant_classes
  )

  invariant_long <- NULL
  invariant_summary <- NULL

  if (!is.null(invariant_peaks) && nrow(invariant_peaks)) {
    messagef(
      "INVARIANTS peaks=%d classes=%s",
      nrow(invariant_peaks),
      paste(opt$invariant_classes, collapse = ",")
    )

    inv_out <- file.path(
      opt$out_dir,
      "invariants",
      paste(
        safe_file_label(opt$invariant_classes),
        collapse = "_"
      )
    )

    invariant_long <- write_peak_outputs(
      peaks = invariant_peaks,
      genes = genes,
      peak_type = "invariant",
      source_label = paste(opt$invariant_classes, collapse = ","),
      out_dir = inv_out,
      top_n = opt$top_genes_per_peak,
      min_abs_loading = opt$min_abs_loading
    )

    if (!is.null(invariant_long) && nrow(invariant_long)) {
      invariant_summary <- write_invariant_summaries(
        invariant_long,
        file.path(opt$out_dir, "invariants"),
        top_n = opt$top_genes_summary
      )
    }
  } else {
    message(
      "No invariant peaks found for class(es): ",
      paste(opt$invariant_classes, collapse = ", ")
    )
  }

  all_long <- bind_rows_fill(
    list(condition_long, invariant_long)
  )

  if (!is.null(all_long) && nrow(all_long)) {
    write_tsv(
      all_long,
      file.path(opt$out_dir, "all_peak_genes_long.tsv")
    )
  }

  manifest <- data.frame(
    library_dir = normalizePath(opt$library_dir),
    genes_file = normalizePath(genes_path),
    n_grid_genes = nrow(genes),
    conditions = if (length(cfiles)) paste(names(cfiles), collapse = ",") else "",
    include_candidates = opt$include_candidates,
    invariant_classes = paste(opt$invariant_classes, collapse = ","),
    top_genes_per_peak = opt$top_genes_per_peak,
    top_genes_summary = opt$top_genes_summary,
    min_abs_loading = opt$min_abs_loading,
    n_condition_gene_rows = if (!is.null(condition_long)) nrow(condition_long) else 0L,
    n_invariant_gene_rows = if (!is.null(invariant_long)) nrow(invariant_long) else 0L,
    stringsAsFactors = FALSE
  )

  write_tsv(
    manifest,
    file.path(opt$out_dir, "manifest.tsv")
  )

  message("Wrote peak-gene maps to: ", normalizePath(opt$out_dir))
}

if (sys.nframe() == 0L) {
  main()
}
