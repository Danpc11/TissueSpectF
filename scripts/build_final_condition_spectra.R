#!/usr/bin/env Rscript
# Build one cross-cohort power spectrum per biological condition and a
# cross-condition invariant spectral architecture.
#
# Condition axis:
#   condition_spectrum_<condition>.tsv  -> complete continuous spectrum
#   condition_signature_<condition>.tsv -> selected candidate/robust components
#
# Invariant axis:
#   invariants/invariant_spectrum.tsv   -> complete invariant spectrum
#   invariants/invariant_components.tsv -> shared/core invariant components
#   invariants/invariant_manifest.tsv   -> thresholds and summary
#
# The invariant layer is intentionally separate from condition signatures:
#
#   S_c(f) = I_tissue(f) + D_c(f)
#
# where I_tissue is shared architecture and D_c is condition-specific deviation.
# Existing per-cohort constant_signature() outputs are NOT replaced.

TSF_BUILD_VERSION <- "2026-09-01-core-invariant-v4"

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
abort <- function(...) stop(paste0(...), call. = FALSE)
messagef <- function(...) message(sprintf(...))

# The results tree defaults the same way the rest of the pipeline does: the
# environment first, then <repo>/results. It used to default to one particular
# run directory ("results_gencode_v2"), which meant the script silently read
# somebody else's tree on any other machine, or aborted on a path the user had
# never created.
default_results_dir <- function() {
  v <- Sys.getenv("TSF_RESULTS_DIR", unset = "")
  if (nzchar(v)) return(v)
  root <- Sys.getenv("TSF_ROOT", unset = getwd())
  file.path(root, "results")
}

parse_args <- function(x) {
  out <- list(
    results_dir = default_results_dir(),
    # Workers for the per-condition aggregation. Defaults to N_WORKERS if set,
    # the same variable the R pipeline honours, then to serial. Capped at the
    # number of conditions, because that is the only axis being split.
    cores = as.integer(Sys.getenv("N_WORKERS", unset = "1")),
    out_dir = NULL,
    datasets = NULL,
    conditions = NULL,
    min_cohorts = 2L,
    top = NA_integer_,
    signature_q = 0.05,
    window_cut = 1,
    min_prevalence = 0.5,
    exclude_chromosomes = "Y,MT",
    healthy_superclass = "no",
    min_period = "off",
    period_margin = 2,
    margin_mode = "add",
    min_period_biological = 0,

    # Invariant-layer parameters
    build_invariants = "yes",
    invariant_conditions = NULL,
    invariant_min_condition_fraction = 0.70,
    invariant_min_prevalence = 0.50,
    invariant_min_phase_coherence = 0.80,
    invariant_max_power_cv = 1.00,
    invariant_max_abs_log2_enrichment = 0.50,
    invariant_core_phase_coherence = 0.90,
    invariant_core_max_power_cv = 0.30,
    invariant_core_max_abs_log2_enrichment = 0.50
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

  out$out_dir <- out$out_dir %||%
    file.path(out$results_dir, "condition_library")

  split_csv <- function(z) {
    if (is.null(z) || !length(z) || !nzchar(z)) return(NULL)
    trimws(strsplit(z, ",", fixed = TRUE)[[1]])
  }

  out$datasets <- split_csv(out$datasets)
  out$conditions <- split_csv(out$conditions)
  out$invariant_conditions <- split_csv(out$invariant_conditions)

  out$min_cohorts <- as.integer(out$min_cohorts)
  out$min_prevalence <- as.numeric(out$min_prevalence)
  out$window_cut <- as.numeric(out$window_cut)
  out$top <- suppressWarnings(as.integer(out$top))
  out$signature_q <- as.numeric(out$signature_q)
  out$period_margin <- as.numeric(out$period_margin)
  out$min_period_biological <- as.numeric(out$min_period_biological)

  out$invariant_min_condition_fraction <-
    as.numeric(out$invariant_min_condition_fraction)
  out$invariant_min_prevalence <-
    as.numeric(out$invariant_min_prevalence)
  out$invariant_min_phase_coherence <-
    as.numeric(out$invariant_min_phase_coherence)
  out$invariant_max_power_cv <-
    as.numeric(out$invariant_max_power_cv)
  out$invariant_max_abs_log2_enrichment <-
    as.numeric(out$invariant_max_abs_log2_enrichment)
  out$invariant_core_phase_coherence <-
    as.numeric(out$invariant_core_phase_coherence)
  out$invariant_core_max_power_cv <-
    as.numeric(out$invariant_core_max_power_cv)
  out$invariant_core_max_abs_log2_enrichment <-
    as.numeric(out$invariant_core_max_abs_log2_enrichment)

  out$exclude_chromosomes <- split_csv(out$exclude_chromosomes) %||%
    character(0)

  out$healthy_superclass <-
    tolower(out$healthy_superclass) %in% c("yes", "true", "1")

  out$build_invariants <-
    tolower(out$build_invariants) %in% c("yes", "true", "1")

  out$margin_mode <- tolower(out$margin_mode)
  if (!out$margin_mode %in% c("add", "mult")) {
    abort("--margin-mode must be add or mult")
  }

  if (!identical(tolower(out$min_period), "off") &&
      !identical(tolower(out$min_period), "auto")) {
    v <- suppressWarnings(as.numeric(out$min_period))
    if (!is.finite(v) || v <= 0) {
      abort("--min-period must be off, auto or a positive number")
    }
  }

  if (is.na(out$min_cohorts) || out$min_cohorts < 2L) {
    abort("--min-cohorts must be >= 2")
  }

  if (!is.na(out$top) && out$top < 1L) {
    abort("--top must be >= 1 when given")
  }

  if (!is.finite(out$signature_q) ||
      out$signature_q <= 0 ||
      out$signature_q > 1) {
    abort("--signature-q must be in (0, 1]")
  }

  if (!is.finite(out$invariant_min_condition_fraction) ||
      out$invariant_min_condition_fraction <= 0 ||
      out$invariant_min_condition_fraction > 1) {
    abort("--invariant-min-condition-fraction must be in (0, 1]")
  }

  if (!is.finite(out$invariant_min_prevalence) ||
      out$invariant_min_prevalence < 0 ||
      out$invariant_min_prevalence > 1) {
    abort("--invariant-min-prevalence must be in [0, 1]")
  }

  out$cores <- suppressWarnings(as.integer(out$cores))
  if (!is.finite(out$cores) || out$cores < 1L) {
    abort("--cores must be an integer >= 1")
  }

  if (!is.finite(out$invariant_min_phase_coherence) ||
      out$invariant_min_phase_coherence < 0 ||
      out$invariant_min_phase_coherence > 1) {
    abort("--invariant-min-phase-coherence must be in [0, 1]")
  }

  if (!is.finite(out$invariant_max_power_cv) ||
      out$invariant_max_power_cv < 0) {
    abort("--invariant-max-power-cv must be >= 0")
  }

  if (!is.finite(out$invariant_max_abs_log2_enrichment) ||
      out$invariant_max_abs_log2_enrichment < 0) {
    abort("--invariant-max-abs-log2-enrichment must be >= 0")
  }

  if (!is.finite(out$invariant_core_phase_coherence) ||
      out$invariant_core_phase_coherence < 0 ||
      out$invariant_core_phase_coherence > 1) {
    abort("--invariant-core-phase-coherence must be in [0, 1]")
  }

  if (!is.finite(out$invariant_core_max_power_cv) ||
      out$invariant_core_max_power_cv < 0) {
    abort("--invariant-core-max-power-cv must be >= 0")
  }

  if (!is.finite(out$invariant_core_max_abs_log2_enrichment) ||
      out$invariant_core_max_abs_log2_enrichment < 0) {
    abort("--invariant-core-max-abs-log2-enrichment must be >= 0")
  }

  out
}

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
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x)) stats::median(x) else NA_real_
}

mean_finite <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

sd_finite <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) > 1L) stats::sd(x) else NA_real_
}

min_finite <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x)) min(x) else NA_real_
}

max_finite <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x)) max(x) else NA_real_
}

chr_order <- c(as.character(1:22), "X", "Y")

order_spectrum <- function(x) {
  x[
    order(
      match(as.character(x$chr), chr_order),
      as.integer(x$k)
    ),
    ,
    drop = FALSE
  ]
}

read_window_map <- function(results_dir, dataset) {
  wdir <- file.path(results_dir, dataset, "window")
  if (!dir.exists(wdir)) return(NULL)

  files <- list.files(
    wdir,
    pattern = "^window_chr.*\\.tsv$",
    full.names = TRUE
  )
  if (!length(files)) return(NULL)

  rows <- lapply(files, function(path) {
    w <- read_tsv(path)

    if (!all(c("k", "window_power") %in% names(w))) {
      return(NULL)
    }

    if (!"chr" %in% names(w)) {
      w$chr <- sub(
        "^window_chr(.*)\\.tsv$",
        "\\1",
        basename(path)
      )
    }

    w$window_rank <- rank(
      -w$window_power,
      ties.method = "min",
      na.last = "keep"
    )

    n_finite <- sum(is.finite(w$window_power))
    w$window_pct <- if (n_finite > 0) {
      100 * w$window_rank / n_finite
    } else {
      NA_real_
    }

    if (!"coverage" %in% names(w)) {
      w$coverage <- NA_real_
    }

    w[
      ,
      c(
        "chr", "k",
        "window_power", "window_rank",
        "window_pct", "coverage"
      ),
      drop = FALSE
    ]
  })

  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(NULL)

  do.call(rbind, rows)
}

discover_inputs <- function(opt) {
  ds <- opt$datasets

  if (is.null(ds)) {
    ds <- list.dirs(
      opt$results_dir,
      recursive = FALSE,
      full.names = FALSE
    )

    ds <- ds[
      dir.exists(
        file.path(opt$results_dir, ds, "consensus")
      )
    ]
  }

  if (!length(ds)) {
    abort(
      "No dataset consensus directories under ",
      opt$results_dir
    )
  }

  out <- list()

  for (id in ds) {
    cdir <- file.path(
      opt$results_dir,
      id,
      "consensus"
    )

    files <- list.files(
      cdir,
      pattern = "^consensus_spectrum_.*\\.tsv$",
      full.names = TRUE
    )

    for (path in files) {
      cond <- sub(
        "^consensus_spectrum_(.*)\\.tsv$",
        "\\1",
        basename(path)
      )

      d <- read_tsv(path)

      needed <- c(
        "chr", "N", "k", "period",
        "median_power_normalised",
        "prevalence_rank",
        "plv",
        "mean_phase"
      )

      miss <- setdiff(needed, names(d))
      if (length(miss)) {
        abort(
          path,
          " lacks: ",
          paste(miss, collapse = ", ")
        )
      }

      d$dataset <- id
      d$condition <- cond
      d$source_condition <- cond

      if (!all(
        c(
          "window_power",
          "window_rank",
          "window_pct"
        ) %in% names(d)
      )) {
        wm <- read_window_map(
          opt$results_dir,
          id
        )

        if (!is.null(wm)) {
          m <- match(
            paste(d$chr, d$k),
            paste(wm$chr, wm$k)
          )

          d$window_power <- wm$window_power[m]
          d$window_rank <- wm$window_rank[m]
          d$window_pct <- wm$window_pct[m]
          d$coverage <- wm$coverage[m]
        } else {
          d$window_power <- NA_real_
          d$window_rank <- NA_real_
          d$window_pct <- NA_real_
          d$coverage <- NA_real_
        }
      }

      sig_path <- file.path(
        cdir,
        paste0(
          "signature_",
          cond,
          ".tsv"
        )
      )

      d$selected_in_signature <- FALSE
      d$confirmed_in_cohort <- FALSE

      if (file.exists(sig_path)) {
        s <- read_tsv(sig_path)
        hit <- match(
          feature_key(d),
          feature_key(s)
        )

        d$selected_in_signature <- !is.na(hit)

        if ("signature_class" %in% names(s)) {
          d$confirmed_in_cohort <-
            !is.na(hit) &
            s$signature_class[hit] %in%
            "confirmed"
        }
      }

      out[[
        paste(
          id,
          cond,
          sep = "\r"
        )
      ]] <- d
    }
  }

  if (!length(out)) {
    abort(
      "No consensus spectra matched ",
      "the requested inputs"
    )
  }

  out
}

add_healthy_superclass <- function(inputs) {
  healthy_labels <- c(
    "Normal_histology",
    "Control_disease_cohort",
    "Control_external_study"
  )

  source_names <- names(inputs)[
    vapply(
      inputs,
      function(x) {
        x$condition[1] %in%
          healthy_labels
      },
      logical(1)
    )
  ]

  for (nm in source_names) {
    d <- inputs[[nm]]
    d$source_condition <- d$condition
    d$condition <- "Healthy"

    inputs[[
      paste0(
        nm,
        "\rHealthy"
      )
    ]] <- d
  }

  inputs
}

add_within_cohort_effect <- function(inputs) {
  by_dataset <- split(
    names(inputs),
    vapply(
      inputs,
      function(x) x$dataset[1],
      character(1)
    )
  )

  for (id in names(by_dataset)) {
    nms <- by_dataset[[id]]

    all_power <- lapply(
      nms,
      function(nm) {
        d <- inputs[[nm]]
        stats::setNames(
          d$median_power_normalised,
          feature_key(d)
        )
      }
    )

    names(all_power) <- vapply(
      nms,
      function(nm) {
        inputs[[nm]]$condition[1]
      },
      character(1)
    )

    for (nm in nms) {
      d <- inputs[[nm]]
      cond <- d$condition[1]

      others <- all_power[
        names(all_power) != cond
      ]

      key <- feature_key(d)

      if (!length(others)) {
        d$cohort_log2_enrichment <- NA_real_
      } else {
        bg <- vapply(
          key,
          function(k) {
            median_finite(
              vapply(
                others,
                function(v) {
                  unname(v[k])
                },
                numeric(1)
              )
            )
          },
          numeric(1)
        )

        eps <- 1e-12

        d$cohort_log2_enrichment <-
          log2(
            (
              d$median_power_normalised +
                eps
            ) /
              (
                bg +
                  eps
              )
          )
      }

      inputs[[nm]] <- d
    }
  }

  inputs
}

# apply_period_floor(), reachable_bh_q() and draws_for_bh() come from
# R/period_floor.R. They used to be defined here and nowhere else, which meant
# the consensus stage tested frequencies this script would later discard -- the
# multiplicity correction ran over an inflated family. One definition now, and
# consensus applies it before its null.
source(file.path(Sys.getenv("TSF_ROOT", unset = getwd()), "R", "period_floor.R"))


aggregate_condition <- function(
  inputs,
  condition,
  min_cohorts,
  window_cut,
  min_prevalence = 0.5,
  exclude_chromosomes = character(0),
  signature_q = 0.05,
  min_period = "off",
  period_margin = 2,
  margin_mode = "add",
  min_period_biological = 0
) {
  pieces <- inputs[
    vapply(
      inputs,
      function(x) {
        identical(
          x$condition[1],
          condition
        )
      },
      logical(1)
    )
  ]

  if (!length(pieces)) {
    return(NULL)
  }

  long <- bind_rows_fill(pieces)

  if (length(exclude_chromosomes)) {
    drop <- as.character(long$chr) %in%
      exclude_chromosomes

    if (any(drop)) {
      long <- long[
        !drop,
        ,
        drop = FALSE
      ]
    }
  }

  if (!nrow(long)) {
    return(NULL)
  }

  long <- apply_period_floor(
    long = long,
    label = condition,
    min_period = min_period,
    period_margin = period_margin,
    margin_mode = margin_mode,
    min_period_biological =
      min_period_biological
  )

  if (!nrow(long)) {
    return(NULL)
  }

  groups <- split(
    seq_len(nrow(long)),
    feature_key(long)
  )

  rows <- lapply(
    groups,
    function(i) {
      x <- long[
        i,
        ,
        drop = FALSE
      ]

      phases <- suppressWarnings(
        as.numeric(x$mean_phase)
      )
      phases <- phases[
        is.finite(phases)
      ]

      phase_vec <- if (length(phases)) {
        mean(
          exp(
            1i *
              phases
          )
        )
      } else {
        NA_complex_
      }

      effect <- median_finite(
        x$cohort_log2_enrichment
      )

      data.frame(
        condition = condition,
        chr = as.character(x$chr[1]),
        N = as.integer(x$N[1]),
        k = as.integer(x$k[1]),
        freq = x$k[1] / x$N[1],
        period = x$N[1] / x$k[1],

        n_cohorts = length(
          unique(x$dataset)
        ),

        cohorts = paste(
          sort(
            unique(x$dataset)
          ),
          collapse = ","
        ),

        source_conditions = paste(
          sort(
            unique(
              if (
                "source_condition" %in%
                  names(x)
              ) {
                as.character(
                  x$source_condition
                )
              } else {
                as.character(
                  x$condition
                )
              }
            )
          ),
          collapse = ","
        ),

        final_power = median_finite(
          x$median_power_normalised
        ),

        final_power_min = min_finite(
          x$median_power_normalised
        ),

        final_power_max = max_finite(
          x$median_power_normalised
        ),

        condition_log2_enrichment =
          effect,

        condition_specific_power =
          if (is.finite(effect)) {
            max(
              2^effect - 1,
              0
            )
          } else {
            NA_real_
          },

        median_prevalence =
          median_finite(
            x$prevalence_rank
          ),

        median_plv_within_cohort =
          median_finite(
            x$plv
          ),

        phase_coherence_between_cohorts =
          if (
            is.finite(
              Re(phase_vec)
            )
          ) {
            Mod(phase_vec)
          } else {
            NA_real_
          },

        mean_phase_between_cohorts =
          if (
            is.finite(
              Re(phase_vec)
            )
          ) {
            Arg(phase_vec)
          } else {
            NA_real_
          },

        n_signature_cohorts =
          sum(
            x$selected_in_signature,
            na.rm = TRUE
          ),

        n_confirmed_cohorts =
          sum(
            x$confirmed_in_cohort,
            na.rm = TRUE
          ),

        signature_replication_fraction =
          mean(
            x$selected_in_signature,
            na.rm = TRUE
          ),

        worst_window_pct =
          min_finite(
            x$window_pct
          ),

        median_window_pct =
          median_finite(
            x$window_pct
          ),

        stringsAsFactors = FALSE
      )
    }
  )

  out <- do.call(
    rbind,
    rows
  )

  out$window_suspect <-
    is.finite(
      out$worst_window_pct
    ) &
    out$worst_window_pct <=
      window_cut

  stouffer <- function(
    p,
    floor_p
  ) {
    p <- suppressWarnings(
      as.numeric(p)
    )
    p <- p[
      is.finite(p)
    ]

    if (!length(p)) {
      return(NA_real_)
    }

    p <- pmin(
      pmax(
        p,
        floor_p
      ),
      1 -
        floor_p
    )

    stats::pnorm(
      sum(
        stats::qnorm(
          1 - p
        )
      ) /
        sqrt(
          length(p)
        ),
      lower.tail = FALSE
    )
  }

  if ("p_null" %in% names(long)) {
    n_draws <- if (
      "n_null" %in% names(long)
    ) {
      max_finite(long$n_null)
    } else {
      NA_real_
    }

    floor_p <- if (
      is.finite(n_draws) &&
      n_draws > 0
    ) {
      1 /
        (
          n_draws +
            1
        )
    } else {
      1e-4
    }

    p_groups <- split(
      long$p_null,
      feature_key(long)
    )

    meta_p <- vapply(
      p_groups,
      stouffer,
      numeric(1),
      floor_p = floor_p
    )

    out$p_meta_null <-
      unname(
        meta_p[
          feature_key(out)
        ]
      )

    out$q_meta_null <-
      stats::p.adjust(
        out$p_meta_null,
        method = "BH"
      )

    out$n_null_draws <-
      n_draws

    eligible_k <- out$n_cohorts[
      out$n_cohorts >= min_cohorts
    ]

    if (length(eligible_k)) {
      k_min <- min(
        eligible_k,
        na.rm = TRUE
      )

      best_meta <- stouffer(
        rep(
          floor_p,
          k_min
        ),
        floor_p
      )

      conservative_singleton_q <- best_meta *
        nrow(out)

      # This is only the BH value one would obtain if a single feature alone
      # occupied the smallest p-value rank. It is NOT a hard lower bound on BH
      # q because multiple tied/extreme p-values receive larger ranks i and can
      # yield smaller p_(i) * m / i values. Therefore never use this diagnostic
      # to claim that significance is mathematically impossible.
      message(
        "  calibrated null: ",
        n_draws,
        " draws, ",
        k_min,
        " cohort(s), ",
        nrow(out),
        " frequencies; pointwise floor p ~= ",
        signif(floor_p, 3),
        "; conservative rank-1 BH diagnostic ~= ",
        signif(conservative_singleton_q, 3),
        ". Final significance is determined from the observed BH-adjusted ",
        "p-value distribution, including ties/ranks."
      )
    }
  } else {
    out$p_meta_null <- NA_real_
    out$q_meta_null <- NA_real_
    out$n_null_draws <- NA_real_
  }

  out$eligible <-
    out$n_cohorts >=
    min_cohorts

  k_cohorts <- max(
    out$n_cohorts,
    na.rm = TRUE
  )

  coherence_floor <- max(
    0.8,
    min(
      0.95,
      2 *
        sqrt(pi) /
        (
          2 *
            sqrt(k_cohorts)
        )
    )
  )

  out$coherence_floor <-
    coherence_floor

  out$meta_score <- with(
    out,
    condition_specific_power *
      median_prevalence *
      median_plv_within_cohort *
      phase_coherence_between_cohorts
  )

  out$meta_score[
    !is.finite(
      out$meta_score
    )
  ] <- 0

  out$signature_class <- "background"

  candidate <-
    out$eligible &
    is.finite(
      out$condition_log2_enrichment
    ) &
    out$condition_log2_enrichment > 0 &
    out$signature_replication_fraction >= 0.5 &
    !out$window_suspect &
    is.finite(
      out$median_prevalence
    ) &
    out$median_prevalence >=
      min_prevalence

  out$signature_class[
    candidate
  ] <- "candidate"

  # If a meta-null exists, only the calibrated q-value can promote a
  # candidate into the selected set. If the null is absent, candidates
  # remain exploratory and are written to the continuous spectrum but are
  # not silently upgraded to robust evidence.
  has_meta_null <- any(
    is.finite(
      out$q_meta_null
    )
  )

  if (has_meta_null) {
    out$signature_selected <-
      candidate &
      is.finite(
        out$q_meta_null
      ) &
      out$q_meta_null <=
        signature_q
  } else {
    out$signature_selected <-
      candidate
  }

  robust <-
    candidate &
    out$signature_selected &
    out$n_signature_cohorts >= 2L &
    is.finite(
      out$phase_coherence_between_cohorts
    ) &
    out$phase_coherence_between_cohorts >=
      coherence_floor &
    has_meta_null

  out$signature_class[
    robust
  ] <- "robust"

  out$signature_class[
    out$window_suspect &
      out$eligible
  ] <- "window_suspect"

  order_spectrum(out)
}

condition_evidence_status <- function(
  tab,
  min_cohorts
) {
  n_available <- max(
    tab$n_cohorts,
    na.rm = TRUE
  )

  if (!is.finite(n_available) ||
      n_available <
        min_cohorts) {
    return(
      "insufficient_cohorts_for_meta_analysis"
    )
  }

  if (!any(
    is.finite(
      tab$q_meta_null
    )
  )) {
    return(
      "missing_meta_null"
    )
  }

  "meta_analysis_available"
}

condition_evidence_message <- function(
  tab,
  condition,
  min_cohorts
) {
  status <- condition_evidence_status(
    tab,
    min_cohorts
  )

  n_available <- max(
    tab$n_cohorts,
    na.rm = TRUE
  )

  if (identical(
    status,
    "insufficient_cohorts_for_meta_analysis"
  )) {
    message(
      "  NOTE: ",
      condition,
      " is represented by only ",
      n_available,
      " cohort(s). A continuous spectrum ",
      "can be reported, but cross-cohort ",
      "meta-analytic replication cannot be ",
      "established. Increasing permutations ",
      "cannot create an independent cohort."
    )
    return(invisible(status))
  }

  if (identical(
    status,
    "missing_meta_null"
  )) {
    message(
      "  NOTE: no usable pointwise ",
      "permutation null is present for ",
      condition,
      ". Re-run ./tsf consensus with ",
      "an adequate --n-null."
    )
    return(invisible(status))
  }

  if (sum(
    tab$n_confirmed_cohorts,
    na.rm = TRUE
  ) == 0L) {
    n_draws_here <- max_finite(
      tab$n_null_draws
    )

    message(
      "  NOTE: no component retained after ",
      "chromosome/period filters was confirmed ",
      "by one cohort alone for ",
      condition,
      ". Cross-cohort evidence remains ",
      "estimable; ",
      if (is.finite(n_draws_here)) {
        paste0(
          n_draws_here,
          " null draws are present."
        )
      } else {
        "a pointwise meta-null is present."
      }
    )
  }

  invisible(status)
}

# -------------------------------------------------------------------------
# Invariant layer
# -------------------------------------------------------------------------

default_invariant_conditions <- function(
  inputs
) {
  conds <- unique(
    vapply(
      inputs,
      function(x) {
        as.character(
          x$condition[1]
        )
      },
      character(1)
    )
  )

  # Healthy is a synthetic superclass. Including it together with its source
  # classes would double-count the same observations.
  conds <- setdiff(
    conds,
    "Healthy"
  )

  preferred <- c(
    "Normal_histology",
    "Control_disease_cohort",
    "Control_external_study",
    "F0", "F1", "F2", "F3", "F4"
  )

  c(
    intersect(
      preferred,
      conds
    ),
    setdiff(
      sort(conds),
      preferred
    )
  )
}

build_invariant_spectrum <- function(
  inputs,
  conditions,
  min_cohorts = 2L,
  min_condition_fraction = 0.70,
  min_prevalence = 0.50,
  min_phase_coherence = 0.80,
  max_power_cv = 1.00,
  max_abs_log2_enrichment = 0.50,
  core_phase_coherence = 0.90,
  core_max_power_cv = 0.30,
  core_max_abs_log2_enrichment = 0.50,
  window_cut = 1,
  exclude_chromosomes = c("Y", "MT"),
  min_period = "off",
  period_margin = 2,
  margin_mode = "add",
  min_period_biological = 0
) {
  conditions <- unique(
    setdiff(
      conditions,
      "Healthy"
    )
  )

  if (!length(conditions)) {
    return(NULL)
  }

  pieces <- inputs[
    vapply(
      inputs,
      function(x) {
        x$condition[1] %in%
          conditions
      },
      logical(1)
    )
  ]

  if (!length(pieces)) {
    return(NULL)
  }

  long <- bind_rows_fill(
    pieces
  )

  # Protect against accidental synthetic duplication even if callers pass
  # inputs after add_healthy_superclass().
  long <- long[
    as.character(long$condition) !=
      "Healthy",
    ,
    drop = FALSE
  ]

  long <- long[
    as.character(long$condition) %in%
      conditions,
    ,
    drop = FALSE
  ]

  if (length(exclude_chromosomes)) {
    long <- long[
      !as.character(long$chr) %in%
        exclude_chromosomes,
      ,
      drop = FALSE
    ]
  }

  if (!nrow(long)) {
    return(NULL)
  }

  long <- apply_period_floor(
    long = long,
    label = "invariant layer",
    min_period = min_period,
    period_margin = period_margin,
    margin_mode = margin_mode,
    min_period_biological =
      min_period_biological
  )

  if (!nrow(long)) {
    return(NULL)
  }

  requested_conditions <- unique(
    conditions
  )
  n_conditions_total <- length(
    requested_conditions
  )

  groups <- split(
    seq_len(nrow(long)),
    feature_key(long)
  )

  rows <- lapply(
    groups,
    function(i) {
      x <- long[
        i,
        ,
        drop = FALSE
      ]

      power <- suppressWarnings(
        as.numeric(
          x$median_power_normalised
        )
      )

      valid_power <- is.finite(power)

      prevalence_raw <- suppressWarnings(
        as.numeric(
          x$prevalence_rank
        )
      )

      # Every final condition spectrum is continuous, so finite power alone
      # cannot define cross-condition "presence": it would make
      # condition_fraction almost always equal to one. A condition counts as
      # carrying this invariant component when its within-condition median
      # prevalence reaches the predefined invariant prevalence floor.
      condition_groups <- split(
        seq_len(nrow(x)),
        as.character(x$condition)
      )

      condition_present_flag <- vapply(
        condition_groups,
        function(j) {
          any(valid_power[j]) &&
            is.finite(
              median_finite(
                prevalence_raw[j]
              )
            ) &&
            median_finite(
              prevalence_raw[j]
            ) >=
              min_prevalence
        },
        logical(1)
      )

      cond_present <- names(
        condition_present_flag
      )[
        condition_present_flag
      ]

      cohort_present <- unique(
        as.character(
          x$dataset[
            valid_power
          ]
        )
      )

      pair_present <- unique(
        paste(
          x$dataset[
            valid_power
          ],
          x$condition[
            valid_power
          ],
          sep = "::"
        )
      )

      power_valid <- power[
        valid_power &
          power > 0
      ]

      phases <- suppressWarnings(
        as.numeric(
          x$mean_phase
        )
      )
      phases <- phases[
        is.finite(phases)
      ]

      phase_vec <- if (length(phases)) {
        mean(
          exp(
            1i *
              phases
          )
        )
      } else {
        NA_complex_
      }

      prevalence <- suppressWarnings(
        as.numeric(
          x$prevalence_rank
        )
      )

      plv <- suppressWarnings(
        as.numeric(
          x$plv
        )
      )

      enrichment <- if (
        "cohort_log2_enrichment" %in%
          names(x)
      ) {
        suppressWarnings(
          as.numeric(
            x$cohort_log2_enrichment
          )
        )
      } else {
        rep(
          NA_real_,
          nrow(x)
        )
      }

      enrichment <- enrichment[
        is.finite(enrichment)
      ]

      window_pct <- if (
        "window_pct" %in%
          names(x)
      ) {
        suppressWarnings(
          as.numeric(
            x$window_pct
          )
        )
      } else {
        rep(
          NA_real_,
          nrow(x)
        )
      }

      mean_power <- mean_finite(
        power_valid
      )

      sd_power <- sd_finite(
        power_valid
      )

      power_cv <- if (
        is.finite(mean_power) &&
        mean_power > 0 &&
        is.finite(sd_power)
      ) {
        sd_power /
          mean_power
      } else {
        NA_real_
      }

      data.frame(
        chr = as.character(
          x$chr[1]
        ),
        N = as.integer(
          x$N[1]
        ),
        k = as.integer(
          x$k[1]
        ),
        freq = x$k[1] /
          x$N[1],
        period = x$N[1] /
          x$k[1],

        n_conditions = length(
          cond_present
        ),
        n_conditions_total =
          n_conditions_total,
        condition_fraction =
          length(cond_present) /
          n_conditions_total,
        conditions = paste(
          sort(
            cond_present
          ),
          collapse = ","
        ),

        n_cohorts = length(
          cohort_present
        ),
        cohorts = paste(
          sort(
            cohort_present
          ),
          collapse = ","
        ),

        n_condition_cohort_pairs =
          length(
            pair_present
          ),

        median_power =
          median_finite(
            power_valid
          ),
        mean_power =
          mean_power,
        power_cv =
          power_cv,

        median_prevalence =
          median_finite(
            prevalence
          ),

        median_plv =
          median_finite(
            plv
          ),

        phase_coherence =
          if (
            is.finite(
              Re(phase_vec)
            )
          ) {
            Mod(phase_vec)
          } else {
            NA_real_
          },

        mean_phase =
          if (
            is.finite(
              Re(phase_vec)
            )
          ) {
            Arg(phase_vec)
          } else {
            NA_real_
          },

        median_abs_condition_log2_enrichment =
          if (
            length(enrichment)
          ) {
            stats::median(
              abs(enrichment)
            )
          } else {
            NA_real_
          },

        max_abs_condition_log2_enrichment =
          if (
            length(enrichment)
          ) {
            max(
              abs(enrichment)
            )
          } else {
            NA_real_
          },

        worst_window_pct =
          min_finite(
            window_pct
          ),

        stringsAsFactors = FALSE
      )
    }
  )

  out <- do.call(
    rbind,
    rows
  )

  out$window_suspect <-
    is.finite(
      out$worst_window_pct
    ) &
    out$worst_window_pct <=
      window_cut

  power_stability <- ifelse(
    is.finite(
      out$power_cv
    ),
    1 /
      (
        1 +
          out$power_cv
      ),
    0
  )

  enrichment_neutrality <- ifelse(
    is.finite(
      out$median_abs_condition_log2_enrichment
    ),
    1 /
      (
        1 +
          out$median_abs_condition_log2_enrichment
      ),
    0
  )

  cohort_replication <- pmin(
    out$n_cohorts /
      min_cohorts,
    1
  )

  prevalence_score <- pmin(
    pmax(
      out$median_prevalence,
      0
    ),
    1
  )

  prevalence_score[
    !is.finite(
      prevalence_score
    )
  ] <- 0

  phase_score <- pmin(
    pmax(
      out$phase_coherence,
      0
    ),
    1
  )

  phase_score[
    !is.finite(
      phase_score
    )
  ] <- 0

  out$invariant_score <-
    out$condition_fraction *
    cohort_replication *
    prevalence_score *
    phase_score *
    power_stability *
    enrichment_neutrality

  out$invariant_score[
    !is.finite(
      out$invariant_score
    )
  ] <- 0

  out$invariant_class <- "background"

  shared_candidate <-
    out$n_cohorts >=
      min_cohorts &
    out$condition_fraction >=
      min_condition_fraction &
    is.finite(
      out$median_prevalence
    ) &
    out$median_prevalence >=
      min_prevalence &
    !out$window_suspect

  out$invariant_class[
    shared_candidate
  ] <- "shared_candidate"

  shared_robust <-
    shared_candidate &
    is.finite(
      out$phase_coherence
    ) &
    out$phase_coherence >=
      min_phase_coherence &
    is.finite(
      out$power_cv
    ) &
    out$power_cv <=
      max_power_cv &
    is.finite(
      out$median_abs_condition_log2_enrichment
    ) &
    out$median_abs_condition_log2_enrichment <=
      max_abs_log2_enrichment

  out$invariant_class[
    shared_robust
  ] <- "shared_robust"

  # Core invariants are the strict tissue baseline: the component must be
  # prevalent in every requested biological condition, preserve phase, show
  # low power heterogeneity, and have no single condition with a large
  # enrichment/depletion. This is intentionally stricter than shared_robust.
  core_invariant <-
    shared_candidate &
    abs(out$condition_fraction - 1) < 1e-12 &
    is.finite(
      out$phase_coherence
    ) &
    out$phase_coherence >=
      core_phase_coherence &
    is.finite(
      out$power_cv
    ) &
    out$power_cv <=
      core_max_power_cv &
    is.finite(
      out$max_abs_condition_log2_enrichment
    ) &
    out$max_abs_condition_log2_enrichment <=
      core_max_abs_log2_enrichment

  out$invariant_class[
    core_invariant
  ] <- "core_invariant"

  out$shared_candidate <- shared_candidate
  out$shared_robust <- shared_robust
  out$core_invariant <- core_invariant

  order_spectrum(out)
}

select_invariant_components <- function(
  invariant_spectrum
) {
  out <- invariant_spectrum[
    invariant_spectrum$invariant_class %in%
      c(
        "shared_candidate",
        "shared_robust",
        "core_invariant"
      ),
    ,
    drop = FALSE
  ]

  if (!nrow(out)) {
    return(out)
  }

  class_rank <- match(
    out$invariant_class,
    c(
      "core_invariant",
      "shared_robust",
      "shared_candidate"
    )
  )

  out[
    order(
      class_rank,
      -out$invariant_score,
      -out$median_power
    ),
    ,
    drop = FALSE
  ]
}

open_plot_device <- function(
  path,
  width = 2400,
  height = 1800,
  res = 180
) {
  if (requireNamespace(
    "ragg",
    quietly = TRUE
  )) {
    ragg::agg_png(
      path,
      width = width,
      height = height,
      units = "px",
      res = res
    )
    return(path)
  }

  if (isTRUE(
    capabilities("cairo")
  )) {
    grDevices::png(
      path,
      width = width,
      height = height,
      res = res,
      type = "cairo"
    )
    return(path)
  }

  pdf_path <- sub(
    "\\.png$",
    ".pdf",
    path,
    ignore.case = TRUE
  )

  warning(
    "Neither ragg nor Cairo is ",
    "available; writing PDF instead: ",
    pdf_path,
    call. = FALSE
  )

  grDevices::pdf(
    pdf_path,
    width = width / res,
    height = height / res,
    onefile = TRUE
  )

  pdf_path
}

plot_signature_by_chromosome <- function(
  signature,
  path,
  condition,
  chr_order = c(
    as.character(1:22),
    "X",
    "Y"
  )
) {
  actual_path <- open_plot_device(
    path
  )
  on.exit(
    grDevices::dev.off(),
    add = TRUE
  )

  old <- graphics::par(
    no.readonly = TRUE
  )
  on.exit(
    graphics::par(old),
    add = TRUE
  )

  graphics::par(
    mfrow = c(6, 4),
    mar = c(
      2.6, 3.2,
      2.1, 0.8
    ),
    oma = c(
      4, 4, 4, 1
    )
  )

  amp <- if (nrow(signature)) {
    pmax(
      signature$final_power,
      1e-15
    )
  } else {
    numeric(0)
  }

  ylim <- if (length(amp)) {
    range(amp) *
      c(
        0.5,
        2
      )
  } else {
    c(
      1e-6,
      1
    )
  }

  xlim <- if (nrow(signature)) {
    range(
      signature$period
    ) *
      c(
        0.8,
        1.2
      )
  } else {
    c(
      1,
      100
    )
  }

  for (chr in chr_order) {
    s <- signature[
      as.character(signature$chr) ==
        chr,
      ,
      drop = FALSE
    ]

    graphics::plot(
      NA,
      xlim = xlim,
      ylim = ylim,
      log = "xy",
      xlab = "period (genes)",
      ylab = "power",
      main = paste0(
        "chr",
        chr,
        if (nrow(s)) {
          paste0(
            " (",
            nrow(s),
            ")"
          )
        } else {
          ""
        }
      ),
      col.main = if (nrow(s)) {
        "black"
      } else {
        "grey60"
      }
    )

    if (!nrow(s)) {
      next
    }

    col <- ifelse(
      s$signature_class ==
        "robust",
      "#B2182B",
      "#2166AC"
    )

    graphics::segments(
      s$period,
      ylim[1],
      s$period,
      pmax(
        s$final_power,
        1e-15
      ),
      col = col,
      lwd = 1.2
    )

    graphics::points(
      s$period,
      pmax(
        s$final_power,
        1e-15
      ),
      pch = 19,
      cex = 0.8,
      col = col
    )
  }

  graphics::mtext(
    paste(
      "Selected components:",
      condition
    ),
    outer = TRUE,
    side = 3,
    line = 1.5,
    cex = 1.35,
    font = 2
  )

  graphics::mtext(
    "period (genes per cycle, log scale)",
    outer = TRUE,
    side = 1,
    line = 1.4
  )

  graphics::mtext(
    "median normalised power across cohorts",
    outer = TRUE,
    side = 2,
    line = 2
  )

  graphics::mtext(
    "red: robust    blue: candidate    grey panel: no components selected",
    outer = TRUE,
    side = 1,
    line = 2.6,
    cex = 0.85
  )

  invisible(actual_path)
}

plot_genome_wide <- function(
  tab,
  signature,
  path,
  condition,
  chr_order = c(
    as.character(1:22),
    "X",
    "Y"
  )
) {
  actual_path <- open_plot_device(
    path,
    width = 3600,
    height = 1200
  )
  on.exit(
    grDevices::dev.off(),
    add = TRUE
  )

  old <- graphics::par(
    no.readonly = TRUE
  )
  on.exit(
    graphics::par(old),
    add = TRUE
  )

  graphics::par(
    mar = c(
      4.2, 4.5,
      3.2, 1
    )
  )

  present <- chr_order[
    chr_order %in%
      as.character(
        tab$chr
      )
  ]

  tab <- tab[
    order(
      match(
        as.character(tab$chr),
        present
      ),
      tab$k
    ),
    ,
    drop = FALSE
  ]

  tab$pos <- seq_len(
    nrow(tab)
  )

  counts <- table(
    factor(
      as.character(tab$chr),
      levels = present
    )
  )

  bounds <- cumsum(
    as.integer(counts)
  )

  mids <- c(
    0,
    utils::head(
      bounds,
      -1
    )
  ) +
    diff(
      c(
        0,
        bounds
      )
    ) /
    2

  y <- pmax(
    tab$final_power,
    1e-15
  )

  graphics::plot(
    tab$pos,
    y,
    type = "h",
    log = "y",
    lwd = 0.4,
    col = "grey78",
    xlab = "",
    ylab = "power (median across cohorts)",
    main = paste(
      "Full spectrum:",
      condition
    ),
    xaxt = "n"
  )

  if (length(bounds) > 1L) {
    graphics::abline(
      v = utils::head(
        bounds,
        -1
      ),
      col = "grey90",
      lwd = 0.6
    )
  }

  graphics::axis(
    1,
    at = mids,
    labels = present,
    tick = FALSE,
    las = 1,
    cex.axis = 0.8,
    line = -0.5
  )

  if (nrow(signature)) {
    m <- match(
      feature_key(signature),
      feature_key(tab)
    )
    ok <- !is.na(m)

    col <- ifelse(
      signature$signature_class[
        ok
      ] ==
        "robust",
      "#B2182B",
      "#2166AC"
    )

    graphics::points(
      tab$pos[
        m[ok]
      ],
      pmax(
        signature$final_power[
          ok
        ],
        1e-15
      ),
      pch = 19,
      cex = 0.7,
      col = col
    )
  }

  graphics::mtext(
    "chromosome (axis: running frequency index, not a genomic coordinate)",
    side = 1,
    line = 3,
    cex = 0.9
  )

  invisible(actual_path)
}

plot_condition <- function(
  tab,
  signature,
  path,
  condition
) {
  actual_path <- open_plot_device(
    path
  )
  on.exit(
    grDevices::dev.off(),
    add = TRUE
  )

  old <- graphics::par(
    no.readonly = TRUE
  )
  on.exit(
    graphics::par(old),
    add = TRUE
  )

  graphics::par(
    mfrow = c(6, 4),
    mar = c(
      2.5, 3.2,
      2.1, 0.8
    ),
    oma = c(
      3, 4, 4, 1
    )
  )

  for (chr in chr_order) {
    x <- tab[
      as.character(tab$chr) ==
        chr,
      ,
      drop = FALSE
    ]

    if (!nrow(x)) {
      graphics::plot.new()
      next
    }

    y <- pmax(
      x$final_power,
      1e-15
    )

    graphics::plot(
      x$k,
      y,
      type = "l",
      log = "y",
      lwd = 0.7,
      xlab = "k",
      ylab = "power",
      main = paste0(
        "chr",
        chr
      ),
      col = "grey25"
    )

    s <- signature[
      as.character(signature$chr) ==
        chr,
      ,
      drop = FALSE
    ]

    if (nrow(s)) {
      graphics::points(
        s$k,
        pmax(
          s$final_power,
          1e-15
        ),
        pch = 19,
        cex = 0.65,
        col = ifelse(
          s$signature_class ==
            "robust",
          "#B2182B",
          "#2166AC"
        )
      )
    }
  }

  graphics::mtext(
    paste(
      "Espectro final de potencia —",
      condition
    ),
    outer = TRUE,
    side = 3,
    line = 1.5,
    cex = 1.35,
    font = 2
  )

  graphics::mtext(
    "Frecuencia discreta dentro de cada cromosoma",
    outer = TRUE,
    side = 1,
    line = 1.2
  )

  graphics::mtext(
    "Mediana intercohorte de potencia normalizada (escala log)",
    outer = TRUE,
    side = 2,
    line = 2
  )

  invisible(actual_path)
}

plot_invariant_spectrum <- function(
  tab,
  path
) {
  actual_path <- open_plot_device(
    path
  )
  on.exit(
    grDevices::dev.off(),
    add = TRUE
  )

  old <- graphics::par(
    no.readonly = TRUE
  )
  on.exit(
    graphics::par(old),
    add = TRUE
  )

  graphics::par(
    mfrow = c(6, 4),
    mar = c(
      2.5, 3.2,
      2.1, 0.8
    ),
    oma = c(
      3, 4, 4, 1
    )
  )

  for (chr in chr_order) {
    x <- tab[
      as.character(tab$chr) ==
        chr,
      ,
      drop = FALSE
    ]

    if (!nrow(x)) {
      graphics::plot.new()
      next
    }

    y <- pmax(
      x$median_power,
      1e-15
    )

    graphics::plot(
      x$k,
      y,
      type = "l",
      log = "y",
      lwd = 0.7,
      xlab = "k",
      ylab = "power",
      main = paste0(
        "chr",
        chr
      ),
      col = "grey35"
    )

    selected <- x[
      x$invariant_class %in%
        c(
          "shared_candidate",
          "shared_robust",
          "core_invariant"
        ),
      ,
      drop = FALSE
    ]

    if (nrow(selected)) {
      graphics::points(
        selected$k,
        pmax(
          selected$median_power,
          1e-15
        ),
        pch = 19,
        cex = 0.65,
        col = ifelse(
          selected$invariant_class ==
            "core_invariant",
          "#762A83",
          ifelse(
            selected$invariant_class ==
              "shared_robust",
            "#B2182B",
            "#2166AC"
          )
        )
      )
    }
  }

  graphics::mtext(
    "Arquitectura espectral invariante",
    outer = TRUE,
    side = 3,
    line = 1.5,
    cex = 1.35,
    font = 2
  )

  graphics::mtext(
    "Frecuencia discreta dentro de cada cromosoma",
    outer = TRUE,
    side = 1,
    line = 1.2
  )

  graphics::mtext(
    "Mediana de potencia entre condición-cohorte (escala log)",
    outer = TRUE,
    side = 2,
    line = 2
  )

  graphics::mtext(
    "purple: core invariant    red: shared robust    blue: shared candidate",
    outer = TRUE,
    side = 1,
    line = 2.5,
    cex = 0.82
  )

  invisible(actual_path)
}

write_invariant_outputs <- function(
  invariant_spectrum,
  conditions,
  opt
) {
  invariant_dir <- file.path(
    opt$out_dir,
    "invariants"
  )

  dir.create(
    invariant_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  spectrum_path <- file.path(
    invariant_dir,
    "invariant_spectrum.tsv"
  )

  components_path <- file.path(
    invariant_dir,
    "invariant_components.tsv"
  )

  core_path <- file.path(
    invariant_dir,
    "invariant_core.tsv"
  )

  manifest_path <- file.path(
    invariant_dir,
    "invariant_manifest.tsv"
  )

  plot_path <- file.path(
    invariant_dir,
    "invariant_power_spectrum.png"
  )

  write_tsv(
    invariant_spectrum,
    spectrum_path
  )

  invariant_components <-
    select_invariant_components(
      invariant_spectrum
    )

  write_tsv(
    invariant_components,
    components_path
  )

  invariant_core <- invariant_spectrum[
    invariant_spectrum$core_invariant %in% TRUE,
    ,
    drop = FALSE
  ]

  if (nrow(invariant_core)) {
    invariant_core <- invariant_core[
      order(
        -invariant_core$invariant_score,
        -invariant_core$median_power
      ),
      ,
      drop = FALSE
    ]
  }

  write_tsv(
    invariant_core,
    core_path
  )

  actual_plot <- plot_invariant_spectrum(
    invariant_spectrum,
    plot_path
  )

  manifest <- data.frame(
    builder_version = TSF_BUILD_VERSION,
    n_conditions =
      length(conditions),
    conditions =
      paste(
        conditions,
        collapse = ","
      ),
    n_frequencies =
      nrow(
        invariant_spectrum
      ),
    n_shared_selected =
      sum(
        invariant_spectrum$shared_candidate,
        na.rm = TRUE
      ),
    n_shared_robust =
      sum(
        invariant_spectrum$shared_robust,
        na.rm = TRUE
      ),
    n_shared_candidate_only =
      sum(
        invariant_spectrum$invariant_class ==
          "shared_candidate",
        na.rm = TRUE
      ),
    n_core_invariants =
      sum(
        invariant_spectrum$core_invariant,
        na.rm = TRUE
      ),
    min_cohorts =
      opt$min_cohorts,
    min_condition_fraction =
      opt$invariant_min_condition_fraction,
    min_prevalence =
      opt$invariant_min_prevalence,
    min_phase_coherence =
      opt$invariant_min_phase_coherence,
    max_power_cv =
      opt$invariant_max_power_cv,
    max_abs_log2_enrichment =
      opt$invariant_max_abs_log2_enrichment,
    core_condition_fraction =
      1.0,
    core_phase_coherence =
      opt$invariant_core_phase_coherence,
    core_max_power_cv =
      opt$invariant_core_max_power_cv,
    core_max_abs_log2_enrichment =
      opt$invariant_core_max_abs_log2_enrichment,
    window_cut_pct =
      opt$window_cut,
    excluded_chromosomes =
      paste(
        opt$exclude_chromosomes,
        collapse = ","
      ),
    min_period =
      opt$min_period,
    period_margin =
      opt$period_margin,
    margin_mode =
      opt$margin_mode,
    min_period_biological =
      opt$min_period_biological,
    plot_file =
      basename(
        actual_plot
      ),
    stringsAsFactors = FALSE
  )

  write_tsv(
    manifest,
    manifest_path
  )

  messagef(
    "%-24s frequencies=%d shared_selected=%d shared_robust=%d core=%d",
    "INVARIANTS",
    nrow(
      invariant_spectrum
    ),
    sum(
      invariant_spectrum$shared_candidate,
      na.rm = TRUE
    ),
    sum(
      invariant_spectrum$shared_robust,
      na.rm = TRUE
    ),
    sum(
      invariant_spectrum$core_invariant,
      na.rm = TRUE
    )
  )

  invisible(
    list(
      spectrum = spectrum_path,
      components = components_path,
      core = core_path,
      manifest = manifest_path,
      plot = actual_plot
    )
  )
}

main <- function() {
  message("TissueSpectF final-spectrum builder version: ", TSF_BUILD_VERSION)

  opt <- parse_args(
    commandArgs(
      trailingOnly = TRUE
    )
  )

  dir.create(
    opt$out_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  # Compute within-cohort effects before adding the optional Healthy synthetic
  # superclass. This prevents Healthy duplicates from contaminating the
  # background used for the original conditions.
  inputs <- discover_inputs(opt)
  inputs <- add_within_cohort_effect(
    inputs
  )

  base_inputs <- inputs

  if (opt$healthy_superclass) {
    message(
      "NOTE: --healthy-superclass merges ",
      "three healthy classes whose equivalence ",
      "has not been demonstrated. ",
      "source_conditions records what went in."
    )

    inputs <- add_healthy_superclass(
      inputs
    )
  }

  conditions <- unique(
    vapply(
      inputs,
      function(x) {
        x$condition[1]
      },
      character(1)
    )
  )

  preferred <- c(
    "Healthy",
    "Normal_histology",
    "Control_disease_cohort",
    "Control_external_study",
    "F0", "F1", "F2", "F3", "F4"
  )

  conditions <- c(
    intersect(
      preferred,
      conditions
    ),
    setdiff(
      sort(conditions),
      preferred
    )
  )

  if (!is.null(opt$conditions)) {
    conditions <- opt$conditions[
      opt$conditions %in%
        conditions
    ]
  }

  manifest <- list()
  condition_tables <- list()

  # Aggregation is the expensive half: it reads every cohort's consensus table
  # for a condition and takes medians across them. It is also pure -- it reads
  # inputs and returns a data frame -- so it parallelises across conditions
  # cleanly. Every WRITE stays in the serial loop below: forked workers writing
  # to the same output directory is how a library ends up half-formed with no
  # error to show for it.
  #
  # Splitting on conditions caps the useful worker count at length(conditions),
  # typically 6. Asking for 32 cores will not make this faster than asking for
  # 6, and the message says so rather than letting the number look effective.
  n_cores <- max(1L, min(opt$cores, length(conditions)))
  if (opt$cores > n_cores) {
    messagef(
      "--cores %d requested; using %d, one per condition (aggregation splits on conditions only)",
      opt$cores, n_cores
    )
  }

  agg_one <- function(cond) {
    aggregate_condition(
      inputs = inputs,
      condition = cond,
      min_cohorts = opt$min_cohorts,
      window_cut = opt$window_cut,
      min_prevalence = opt$min_prevalence,
      exclude_chromosomes =
        opt$exclude_chromosomes,
      signature_q = opt$signature_q,
      min_period = opt$min_period,
      period_margin =
        opt$period_margin,
      margin_mode = opt$margin_mode,
      min_period_biological =
        opt$min_period_biological
    )
  }

  if (n_cores > 1L) {
    messagef("aggregating %d condition(s) on %d worker(s)",
             length(conditions), n_cores)
    tabs <- parallel::mclapply(conditions, agg_one, mc.cores = n_cores,
                               mc.set.seed = FALSE)
    # mclapply reports a worker failure as a try-error in the result rather
    # than by raising, so an unchecked run would silently drop a condition and
    # still write a manifest claiming success.
    failed <- vapply(tabs, function(x) inherits(x, "try-error"), logical(1))
    if (any(failed)) {
      abort("aggregation failed for condition(s): ",
            paste(conditions[failed], collapse = ", "), "\n  ",
            paste(unique(vapply(tabs[failed], conditionMessage, character(1))),
                  collapse = "\n  "))
    }
  } else {
    tabs <- lapply(conditions, agg_one)
  }
  names(tabs) <- conditions

  for (cond in conditions) {
    tab <- tabs[[cond]]

    if (is.null(tab)) {
      next
    }

    condition_tables[[cond]] <- tab

    full_path <- file.path(
      opt$out_dir,
      paste0(
        "condition_spectrum_",
        cond,
        ".tsv"
      )
    )

    write_tsv(
      tab,
      full_path
    )

    # Preserve both exploratory candidates and robust calls. The signature file
    # is therefore the selected non-background layer, while the complete
    # continuous spectrum remains in condition_spectrum_*.tsv.
    sig <- tab[
      tab$eligible %in%
        TRUE &
        tab$signature_class %in%
        c(
          "candidate",
          "robust"
        ),
      ,
      drop = FALSE
    ]

    sig <- sig[
      order(
        sig$signature_class !=
          "robust",
        -sig$meta_score,
        -sig$final_power
      ),
      ,
      drop = FALSE
    ]

    if (!is.na(opt$top) &&
        nrow(sig) >
          opt$top) {
      message(
        "  NOTE: ",
        nrow(sig),
        " components are candidate/robust for ",
        cond,
        "; writing the ",
        opt$top,
        " highest-scoring. Raise --top ",
        "to keep them all."
      )

      sig <- utils::head(
        sig,
        opt$top
      )
    }

    sig_path <- file.path(
      opt$out_dir,
      paste0(
        "condition_signature_",
        cond,
        ".tsv"
      )
    )

    write_tsv(
      sig,
      sig_path
    )

    png_path <- file.path(
      opt$out_dir,
      paste0(
        "power_spectrum_",
        cond,
        ".png"
      )
    )

    plot_path <- plot_condition(
      tab,
      sig,
      png_path,
      cond
    )

    sig_png <- plot_signature_by_chromosome(
      sig,
      file.path(
        opt$out_dir,
        paste0(
          "signature_peaks_",
          cond,
          ".png"
        )
      ),
      cond
    )

    genome_png <- plot_genome_wide(
      tab,
      sig,
      file.path(
        opt$out_dir,
        paste0(
          "genome_spectrum_",
          cond,
          ".png"
        )
      ),
      cond
    )

    evidence_status <-
      condition_evidence_status(
        tab,
        opt$min_cohorts
      )

    condition_evidence_message(
      tab,
      cond,
      opt$min_cohorts
    )

    manifest[[cond]] <- data.frame(
      builder_version = TSF_BUILD_VERSION,
      condition = cond,
      n_cohorts_max =
        max(
          tab$n_cohorts,
          na.rm = TRUE
        ),
      evidence_status =
        evidence_status,
      coherence_floor =
        tab$coherence_floor[1],
      min_prevalence =
        opt$min_prevalence,
      window_cut_pct =
        opt$window_cut,
      excluded_chromosomes =
        paste(
          opt$exclude_chromosomes,
          collapse = ","
        ),
      healthy_superclass =
        opt$healthy_superclass,
      n_frequencies =
        nrow(tab),
      n_robust =
        sum(
          tab$signature_class ==
            "robust",
          na.rm = TRUE
        ),
      n_candidates =
        sum(
          tab$signature_class ==
            "candidate",
          na.rm = TRUE
        ),
      n_signature_calibrated =
        sum(
          tab$signature_selected,
          na.rm = TRUE
        ),
      plot_signature =
        basename(sig_png),
      plot_genome =
        basename(genome_png),
      signature_q =
        opt$signature_q,
      min_period =
        opt$min_period,
      period_margin =
        opt$period_margin,
      margin_mode =
        opt$margin_mode,
      min_period_biological =
        opt$min_period_biological,
      n_window_suspect =
        sum(
          tab$signature_class ==
            "window_suspect",
          na.rm = TRUE
        ),
      plot_file =
        basename(plot_path),
      stringsAsFactors = FALSE
    )

    messagef(
      "%-24s cohorts=%d frequencies=%d robust=%d candidates=%d status=%s",
      cond,
      max(
        tab$n_cohorts,
        na.rm = TRUE
      ),
      nrow(tab),
      sum(
        tab$signature_class ==
          "robust",
        na.rm = TRUE
      ),
      sum(
        tab$signature_class ==
          "candidate",
        na.rm = TRUE
      ),
      evidence_status
    )
  }

  if (length(manifest)) {
    write_tsv(
      do.call(
        rbind,
        manifest
      ),
      file.path(
        opt$out_dir,
        "condition_library_manifest.tsv"
      )
    )
  }

  if (opt$build_invariants) {
    invariant_conditions <- if (
      is.null(
        opt$invariant_conditions
      )
    ) {
      default_invariant_conditions(
        base_inputs
      )
    } else {
      unique(
        setdiff(
          opt$invariant_conditions,
          "Healthy"
        )
      )
    }

    available_base_conditions <- unique(
      vapply(
        base_inputs,
        function(x) {
          x$condition[1]
        },
        character(1)
      )
    )

    missing_invariant_conditions <- setdiff(
      invariant_conditions,
      available_base_conditions
    )

    if (length(
      missing_invariant_conditions
    )) {
      message(
        "NOTE: invariant conditions not found ",
        "and therefore ignored: ",
        paste(
          missing_invariant_conditions,
          collapse = ", "
        )
      )

      invariant_conditions <- intersect(
        invariant_conditions,
        available_base_conditions
      )
    }

    if (length(
      invariant_conditions
    )) {
      invariant_spectrum <- build_invariant_spectrum(
        inputs = base_inputs,
        conditions = invariant_conditions,
        min_cohorts = opt$min_cohorts,
        min_condition_fraction =
          opt$invariant_min_condition_fraction,
        min_prevalence =
          opt$invariant_min_prevalence,
        min_phase_coherence =
          opt$invariant_min_phase_coherence,
        max_power_cv =
          opt$invariant_max_power_cv,
        max_abs_log2_enrichment =
          opt$invariant_max_abs_log2_enrichment,
        core_phase_coherence =
          opt$invariant_core_phase_coherence,
        core_max_power_cv =
          opt$invariant_core_max_power_cv,
        core_max_abs_log2_enrichment =
          opt$invariant_core_max_abs_log2_enrichment,
        window_cut =
          opt$window_cut,
        exclude_chromosomes =
          opt$exclude_chromosomes,
        min_period =
          opt$min_period,
        period_margin =
          opt$period_margin,
        margin_mode =
          opt$margin_mode,
        min_period_biological =
          opt$min_period_biological
      )

      if (!is.null(
        invariant_spectrum
      ) &&
          nrow(
            invariant_spectrum
          )) {
        write_invariant_outputs(
          invariant_spectrum,
          invariant_conditions,
          opt
        )
      } else {
        message(
          "NOTE: invariant layer produced no ",
          "frequencies after the configured filters."
        )
      }
    } else {
      message(
        "NOTE: no conditions are available ",
        "for the invariant layer."
      )
    }
  }

  message(
    "Wrote final condition spectra to: ",
    normalizePath(
      opt$out_dir
    )
  )
}

if (sys.nframe() == 0L) {
  main()
}
