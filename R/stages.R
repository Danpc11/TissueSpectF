# stages.R -- every pipeline stage as a callable function.
#
# The numbered scripts and the `tsf` CLI both call these, so there is exactly
# one implementation of each stage. Each takes (project, opt) where opt carries
# datasets, cond, branch and force, and each returns a short summary the caller
# can print or aggregate.

stage_names <- c("ingest", "spectra", "maxt", "condition", "consensus", "clean",
                 "stability", "peaks", "compare")

# Stages that can be skipped on a small machine. `maxt` is per sample and costs
# hours; `condition` answers the primary question at ~1/n of the cost, so a run
# with --from=condition --to=compare is a complete analysis minus the
# per-sample reproducibility figures.
optional_stages <- c("maxt")

#' Resolve which datasets a stage should run over.
stage_datasets <- function(opt) {
  if (length(opt$datasets)) opt$datasets
  else sub("\\.R$", "", list.files("config/datasets", pattern = "\\.R$"))
}

#' Worker budget for the local process manager.
#'
#' There is deliberately no scheduler detection here: this project runs on a
#' host whose resource manager exposes a local process budget.  An explicit
#' --cores value wins, followed by N_WORKERS, then a conservative automatic
#' default.  BLAS is kept at one thread so that each fork remains one CPU task.
local_workers <- function(opt, n_tasks = Inf, default_max = 8L) {
  valid_int <- function(x) {
    x <- suppressWarnings(as.integer(x))
    if (length(x) != 1L || is.na(x) || x < 1L) NA_integer_ else x
  }

  requested <- valid_int(opt$cores %||% NA_integer_)
  source <- "--cores"
  if (is.na(requested)) {
    requested <- valid_int(Sys.getenv("N_WORKERS", unset = NA_character_))
    source <- "N_WORKERS"
  }

  detected <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (length(detected) != 1L || is.na(detected) || detected < 1L) detected <- 1L
  if (is.na(requested)) {
    requested <- min(as.integer(default_max), max(1L, detected - 1L))
    source <- "automatic"
  }

  n_tasks <- suppressWarnings(as.integer(n_tasks))
  if (length(n_tasks) != 1L || is.na(n_tasks) || n_tasks < 1L) n_tasks <- detected
  workers <- max(1L, min(requested, detected, n_tasks))

  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
             MKL_NUM_THREADS = "1", BLIS_NUM_THREADS = "1",
             VECLIB_MAXIMUM_THREADS = "1")
  attr(workers, "source") <- source
  workers
}

# --- 01 ingest ---------------------------------------------------------------
stage_ingest <- function(project, opt) {
  ids <- stage_datasets(opt)
  present <- list()
  for (id in ids) {
    res <- ingest_dataset(id, project)
    present[[id]] <- res$audit$present
  }
  if (length(present) > 1) {
    common <- Reduce(intersect, present)
    all_conditions <- unique(unlist(present, use.names = FALSE))
    for (id in names(present)) {
      missing <- setdiff(all_conditions, present[[id]])
      if (length(missing)) {
        tsf_log(id, " does not contribute: ", paste(missing, collapse = ", "))
      }
    }
    if (length(common)) {
      tsf_log("Conditions shared by every selected dataset: ",
              paste(common, collapse = ", "))
    } else {
      tsf_warn("No condition is shared by every selected dataset. This is valid ",
               "for library construction: each cohort was ingested independently. ",
               "Cross-cohort synthesis must group cohorts by condition rather than ",
               "requiring one global intersection.")
    }
  }
  sprintf("%d dataset(s) ingested", length(ids))
}

# --- 02 spectra --------------------------------------------------------------
stage_spectra <- function(project, opt) {
  n <- 0L
  for (id in stage_datasets(opt)) {
    inp <- tsf_stage_inputs(project, id)
    tsf_log(id, ": spectra over ", length(inp$chrom_idx), " chromosome(s)")
    for (cond in tsf_conditions(inp$conditions, opt)) {
      spec <- compute_condition_spectra(inp$dataset, cond, inp$chrom_idx)
      if (is.null(spec)) next
      is_summary <- spec$sample %in% c("avg_signal", "median_signal")
      write_tsv_tsf(spec[is_summary, ], p_spectra_condition(inp$paths, cond))
      write_tsv_tsf(spec[!is_summary, ], p_spectra_samples(inp$paths, cond))
      tsf_log("  ", cond, ": ", sum(is_summary), " summary rows, ",
              sum(!is_summary), " sample rows")
      n <- n + 1L
    }
  }
  sprintf("%d condition spectra", n)
}

# --- 03 maxT -----------------------------------------------------------------
stage_maxt <- function(project, opt) {
  computed <- 0L; reused <- 0L
  for (id in stage_datasets(opt)) {
    inp <- tsf_stage_inputs(project, id)
    n_cores <- maxt_cores(inp$chrom_idx)
    tsf_log(id, ": maxT (B = ", project$maxt$B, ", ", n_cores, " core(s))")
    for (cond in tsf_conditions(inp$conditions, opt)) {
      out_path <- p_maxt(inp$paths, cond)
      if (!isTRUE(opt$force) && file.exists(out_path) && file.size(out_path) > 0) {
        tsf_log("  ", cond, ": reusing existing maxT (--force to recompute)")
        reused <- reused + 1L
        next
      }
      t0 <- Sys.time()
      res <- maxt_condition(inp$dataset, cond, inp$chrom_idx, project$maxt, n_cores)
      if (is.null(res)) { tsf_warn("  ", cond, ": no maxT result"); next }
      write_tsv_tsf(res, out_path)
      write_tsv_tsf(data.frame(condition = cond,
                               sample_id = attr(res, "expected_samples")),
                    file.path(inp$paths$maxt, sprintf("expected_samples_%s.tsv", cond)))
      tsf_log("  ", cond, ": ", nrow(res), " rows in ",
              round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
      computed <- computed + 1L
    }
  }
  sprintf("%d computed, %d reused", computed, reused)
}

# --- condition-level test ---------------------------------------------------
stage_condition <- function(project, opt) {
  n_sig <- 0L
  for (id in stage_datasets(opt)) {
    inp <- tsf_stage_inputs(project, id, need = "maxt")
    n_cores <- maxt_cores(inp$chrom_idx)
    tsf_log(id, ": condition-level test (B = ",
            project$maxt$condition_B %||% project$maxt$B, ", ", n_cores, " core(s))")
    for (cond in tsf_conditions(inp$conditions, opt)) {
      for (branch in opt$branches) {
        mi <- inp$maxt[[cond]]
        if (is.null(mi)) {
          tsf_log("  ", cond, ": no per-sample maxT; Stouffer will be omitted")
        }
        cs <- condition_significance(inp$dataset, cond, inp$chrom_idx, project$maxt,
                                     branch = branch, maxt_individual = mi,
                                     n_cores = n_cores)
        if (is.null(cs)) { tsf_warn("  ", cond, "/", branch, ": no result"); next }
        write_tsv_tsf(cs, file.path(inp$paths$base, "condition",
                                    sprintf("condition_significance_%s_%s.tsv",
                                            branch, cond)))
        hit <- sum(cs$q_condition <= 0.05, na.rm = TRUE)
        n_sig <- n_sig + hit
        tsf_log("  ", cond, "/", branch, ": ", hit, " peak(s) at q <= 0.05 of ",
                nrow(cs), " frequencies")
      }
    }
  }
  sprintf("%d condition-level peak(s) at q <= 0.05", n_sig)
}

# --- consensus spectrum ------------------------------------------------------
stage_consensus <- function(project, opt) {
  n_sig <- 0L
  for (id in stage_datasets(opt)) {
    inp <- tsf_stage_inputs(project, id, need = "maxt")
    detected <- suppressWarnings(parallel::detectCores(logical = TRUE))
    if (is.na(detected) || detected < 1L) detected <- 1L
    requested <- suppressWarnings(as.integer(opt$cores %||% NA_integer_))
    # Consensus workers each touch a large spectral table. Eight is a safer
    # default than blindly forking once per chromosome; --cores can raise or
    # lower it after checking available RAM.
    n_cores <- if (is.na(requested)) min(8L, max(1L, detected - 1L)) else
      max(1L, min(requested, detected))
    tsf_log(id, ": consensus spectra from per-sample spectra (", n_cores,
            " worker(s))")
    # The label-permuted null needs every sample of the dataset, not just the
    # condition's, so it is assembled once.
    pool <- do.call(rbind, lapply(inp$conditions, function(c)
      read_tsv_tsf(p_spectra_samples(inp$paths, c), required = FALSE)))
    # The null depends only on how many samples are drawn, so conditions of the
    # same size share it. Without this the same permutation set is recomputed
    # once per condition, which dominates the stage's cost on real data.
    null_cache <- list()
    # Samples are not always independent. When the metadata names a blocking
    # variable, the null draws whole blocks so its dependence structure matches
    # the data's.
    block_col <- project$consensus$permutation_block %||% NULL
    blocks <- NULL
    if (!is.null(block_col)) {
      if (block_col %in% colnames(inp$dataset$samples)) {
        blocks <- stats::setNames(as.character(inp$dataset$samples[[block_col]]),
                                  inp$dataset$samples$sample_id)
        tsf_log("  null blocked by '", block_col, "': ",
                length(unique(blocks)), " block(s) over ", length(blocks), " samples")
      } else {
        tsf_warn("  consensus$permutation_block names '", block_col,
                 "', which is not a column of samples.tsv; the null will treat ",
                 "samples as independent")
      }
    }
    for (cond in tsf_conditions(inp$conditions, opt)) {
      sp <- read_tsv_tsf(p_spectra_samples(inp$paths, cond), required = FALSE)
      if (is.null(sp) || !nrow(sp)) {
        tsf_warn("  ", cond, ": no per-sample spectra; run the spectra stage")
        next
      }
      cs <- consensus_spectrum(sp, maxt = inp$maxt[[cond]],
                               n_boot = project$consensus$n_boot %||% 500L,
                               alpha = project$maxt$alpha,
                               seed = project$maxt$seed,
                               quantile_cut = project$consensus$quantile_cut %||% 0.95,
                               n_cores = n_cores)
      if (is.null(cs)) { tsf_warn("  ", cond, ": no consensus"); next }

      n_here <- length(unique(sp$sample))
      key <- as.character(n_here)
      null_dist <- NULL
      if (!is.null(pool)) {
        if (is.null(null_cache[[key]])) {
          t0 <- Sys.time()
          null_cache[[key]] <- list(d = null_consensus_distribution(
            pool, n_here, n_null = project$consensus$n_null %||% 50L,
            seed = project$maxt$seed,
            quantile_cut = project$consensus$quantile_cut %||% 0.95,
            blocks = blocks, n_cores = n_cores))
          tsf_log("  null for n = ", n_here, ": global q95 = ",
                  signif(null_cache[[key]]$d$global %||% NA, 3), " over ",
                  null_cache[[key]]$d$n_null %||% 0, " draws (",
                  round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2),
                  " min)")
        }
        null_dist <- null_cache[[key]]$d
      }
      cs <- null_component_pvalues(cs, null_dist,
                                   null_q = project$consensus$null_q %||% 0.05)
      write_tsv_tsf(cs, file.path(inp$paths$base, "consensus",
                                  sprintf("consensus_spectrum_%s.tsv", cond)))
      sig <- consensus_signature(cs,
                                 null_q = project$consensus$null_q %||% 0.05,
                                 max_components = project$consensus$max_components %||% 50L,
                                 min_prevalence = project$consensus$min_prevalence %||% 0.5,
                                 plv_q = project$consensus$plv_q %||% 0.05)
      if (!is.null(sig)) {
        sig$condition <- cond
        write_tsv_tsf(sig, file.path(inp$paths$base, "consensus",
                                     sprintf("signature_%s.tsv", cond)))
        n_conf <- sum(sig$signature_class == "confirmed")
        n_sig <- n_sig + n_conf
        tsf_log("  ", cond, ": ", n_conf, " confirmed and ",
                nrow(sig) - n_conf, " exploratory component(s) of ", nrow(cs),
                " frequencies | top: chr", sig$chr[1], " k", sig$k[1],
                " period ", round(sig$period[1]), " genes, PLV ",
                round(sig$plv[1], 2), ", prevalence ", round(sig$prevalence[1], 2),
                " [", sig$signature_class[1], "]")
      } else {
        tsf_log("  ", cond, ": no component passes prevalence and phase alignment")
      }
    }
  }
  sprintf("%d confirmed signature component(s) across conditions", n_sig)
}

# --- CLEAN decomposition -----------------------------------------------------
stage_clean <- function(project, opt) {
  total <- 0L
  for (id in stage_datasets(opt)) {
    inp <- tsf_stage_inputs(project, id)
    n_cores <- maxt_cores(inp$chrom_idx)
    tsf_log(id, ": CLEAN decomposition (EBIC gamma = ",
            project$clean$ebic_gamma %||% 1, ")")
    for (cond in tsf_conditions(inp$conditions, opt)) {
      for (branch in opt$branches) {
        cp <- clean_condition(inp$dataset, cond, inp$chrom_idx, branch = branch,
                              per_sample = isTRUE(project$clean$per_sample),
                              clean_cfg = project$clean, n_cores = n_cores)
        if (is.null(cp)) { tsf_log("  ", cond, "/", branch, ": no component"); next }
        write_tsv_tsf(cp, file.path(inp$paths$base, "clean",
                                    sprintf("components_%s_%s.tsv", branch, cond)))
        tsf_log("  ", cond, "/", branch, ": ", nrow(cp), " component(s) over ",
                length(unique(cp$chr)), " chromosome(s), median period ",
                round(stats::median(cp$period)), " genes")
        total <- total + nrow(cp)
      }
    }
  }
  sprintf("%d component(s) extracted", total)
}

# --- 04 stability ------------------------------------------------------------
stage_stability <- function(project, opt) {
  total_stable <- 0L
  for (id in stage_datasets(opt)) {
    inp <- tsf_stage_inputs(project, id, need = c("maxt", "spectra"))
    tsf_log(id, ": stability (alpha = ", project$maxt$alpha,
            ", stable_frac = ", project$maxt$stable_frac, ")")
    for (cond in tsf_conditions(inp$conditions, opt)) {
      m <- inp$maxt[[cond]]
      if (is.null(m)) {
        # No per-sample maxT: fall back to the condition-level result alone.
        st <- stability_from_condition(inp$paths, cond, opt$branches[1])
        if (is.null(st)) {
          tsf_warn("  ", cond, ": neither maxT nor a condition test; run one of them")
          next
        }
        write_tsv_tsf(st, p_stability(inp$paths, cond))
        tsf_log("  ", cond, ": ", sum(st$is_stable), " selected (condition test only)")
        total_stable <- total_stable + sum(st$is_stable)
        for (branch in opt$branches) {
          pk <- condition_peak_table(st, inp$spectra[[cond]], branch)
          if (!is.null(pk)) write_tsv_tsf(pk, p_peaks(inp$paths, branch, cond))
        }
        next
      }
      expected <- read_tsv_tsf(file.path(inp$paths$maxt,
                                         sprintf("expected_samples_%s.tsv", cond)),
                               required = FALSE)
      expected_samples <- if (!is.null(expected)) expected$sample_id else unique(m$sample)
      st <- stable_peaks_maxt(m, expected_samples,
                              alpha = project$maxt$alpha,
                              stable_frac = project$maxt$stable_frac)
      st <- attach_condition_test(st, inp$paths, cond, opt$branches[1],
                                  project$stability_criterion %||% "condition")
      write_tsv_tsf(st, p_stability(inp$paths, cond))
      tsf_log("  ", cond, ": ", sum(st$is_stable), " stable of ", nrow(st), " peaks")
      total_stable <- total_stable + sum(st$is_stable)
      for (branch in opt$branches) {
        pk <- condition_peak_table(st, inp$spectra[[cond]], branch)
        if (!is.null(pk)) write_tsv_tsf(pk, p_peaks(inp$paths, branch, cond))
      }
    }
  }
  sprintf("%d stable peak(s) across conditions", total_stable)
}

# --- 05 peak genes -----------------------------------------------------------
stage_peaks <- function(project, opt) {
  written <- 0L
  for (id in stage_datasets(opt)) {
    inp <- tsf_stage_inputs(project, id)
    tsf_log(id, ": peak-gene reconstruction")
    for (cond in tsf_conditions(inp$conditions, opt)) {
      for (branch in opt$branches) {
        peaks <- read_tsv_tsf(p_peaks(inp$paths, branch, cond), required = FALSE)
        if (is.null(peaks) || !nrow(peaks)) {
          tsf_log("  ", cond, "/", branch,
                  ": no selected peak (the stability result may legitimately be zero)")
          next
        }
        out_dir <- p_peak_genes_dir(inp$paths, branch, cond)
        unlink(list.files(out_dir, pattern = "^pico_chr.*\\.tsv$", full.names = TRUE))
        n <- write_peak_gene_tables(peaks, inp$dataset$genes, inp$chrom_idx,
                                    out_dir, branch, cond)
        tsf_log("  ", cond, "/", branch, ": ", n, " peak file(s)")
        written <- written + n
      }
    }
  }
  sprintf("%d peak file(s)", written)
}

# --- 06 compare --------------------------------------------------------------
stage_compare <- function(project, opt) {
  ids <- stage_datasets(opt)
  if (length(ids) < 2) {
    tsf_warn("Comparison needs at least two datasets; skipped")
    return("skipped")
  }
  loaded <- lapply(ids, function(id) {
    inp <- tsf_stage_inputs(project, id, need = c("stability", "maxt"))
    inp$peaks <- stats::setNames(lapply(c("average", "median"), function(br)
      stats::setNames(lapply(inp$conditions, function(c)
        read_tsv_tsf(p_peaks(inp$paths, br, c), required = FALSE)), inp$conditions)),
      c("average", "median"))
    inp
  })
  names(loaded) <- ids

  # Datasets are only comparable within a tissue and a vocabulary. Crossing a
  # liver fibrosis stage against a lung disease grade would merge on a label
  # that means different things, so groups are formed first and each is compared
  # on its own.
  keys <- vapply(ids, function(id) {
    cfg <- load_dataset_config(id)
    paste(cfg$tissue %||% "unknown", cfg$vocabulary_spec$id, sep = "/")
  }, character(1))
  groups <- split(ids, keys)
  singles <- names(groups)[lengths(groups) < 2]
  if (length(singles)) {
    tsf_warn("Only one dataset in: ", paste(singles, collapse = ", "),
             " -- no cross-dataset replication possible there")
  }
  groups <- groups[lengths(groups) >= 2]
  if (!length(groups)) {
    tsf_warn("No tissue/vocabulary group has two datasets; nothing to compare")
    return("nothing to compare")
  }
  tsf_log("Comparison groups: ", paste(names(groups), collapse = " | "))

  # The expensive independent unit is a dataset.  Groups and branches stay
  # sequential; within a branch compare_one_group forks at most one worker per
  # cohort, then the parent performs the cross-dataset joins.  This gives one
  # global CPU budget and prevents nested parallelism.
  max_group_size <- max(lengths(groups))
  n_cores <- local_workers(opt, n_tasks = max_group_size)
  tsf_log("compare: ", n_cores, " local worker(s) [",
          attr(n_cores, "source"), "]; BLAS threads per worker = 1")

  # EVERY group is processed, each into its own directory. Taking only the
  # largest group would silently drop, say, the lung datasets whenever the liver
  # ones outnumbered them.
  results <- vapply(names(groups), function(gname) {
    compare_one_group(project, loaded[groups[[gname]]], groups[[gname]],
                      gname, opt, n_cores = n_cores)
  }, character(1))
  return(paste(results, collapse = "; "))
}

#' Compare the datasets of one tissue/vocabulary group.
compare_one_group <- function(project, loaded, ids, group_name, opt,
                              n_cores = 1L) {
  compare_dir <- file.path(project$results_dir, "comparison", group_name)
  ensure_dir(compare_dir)
  tsf_log("=== group ", group_name, ": ", paste(ids, collapse = ", "), " ===")
  n_replicated <- 0L

  voc <- load_dataset_config(ids[1])$vocabulary_spec
  declared <- tsf_levels(voc)
  available <- intersect(declared, unique(unlist(lapply(loaded, function(x)
    as.character(x$conditions)), use.names = FALSE)))
  contributors <- stats::setNames(vapply(available, function(cond)
    sum(vapply(loaded, function(x) cond %in% x$conditions, logical(1))), integer(1)),
    available)
  tsf_log("Condition coverage: ", paste0(names(contributors), "=",
          contributors, collapse = ", "))
  one_only <- names(contributors)[contributors < 2L]
  if (length(one_only)) {
    tsf_warn("No cross-dataset replication for condition(s) contributed by one ",
             "dataset only: ", paste(one_only, collapse = ", "))
  }
  ordered_voc <- isTRUE(voc$ordered)
  if (!ordered_voc) {
    tsf_log("Vocabulary '", voc$id, "' is unordered: reporting per-level ",
            "comparisons, no transitions")
  }
  n_replicated <- 0L

  for (branch in opt$branches) {
    tsf_log("branch: ", branch)
    workers_now <- max(1L, min(as.integer(n_cores), length(ids)))

    per_dataset <- parallel::mclapply(ids, function(id) {
      x <- loaded[[id]]
      dataset_conditions <- intersect(declared, as.character(x$conditions))
      dataset_progression <- intersect(tsf_progression(voc), dataset_conditions)
      sig <- constant_signature(x$stability, dataset_conditions)
      if (!is.null(sig)) {
        sig$branch <- branch
        write_tsv_tsf(sig, file.path(x$paths$base, "comparison",
                                     sprintf("constant_signature_%s.tsv", branch)))
        tsf_log("  ", id, ": constant signature = ", nrow(sig), " peak(s)")
      } else tsf_log("  ", id, ": constant signature is empty")
      conditions <- do.call(rbind, lapply(dataset_conditions, function(c) {
        p <- x$peaks[[branch]][[c]]
        if (is.null(p) || !nrow(p)) return(NULL)
        p$condition <- c
        p
      }))

      trans <- list()
      for (i in if (ordered_voc) seq_len(max(0L, length(dataset_progression) - 1L)) else integer(0)) {
        t <- transition_table(dataset_progression[i], dataset_progression[i + 1L], x$peaks[[branch]],
                              x$stability, x$maxt, branch)
        if (is.null(t)) {
          tsf_log("  ", id, ": ", dataset_progression[i], " -> ", dataset_progression[i + 1L],
                  ": no shared stable peak")
          next
        }
        tsf_log("  ", id, ": ", dataset_progression[i], " -> ", dataset_progression[i + 1L], ": ", nrow(t),
                " shared, ", sum(t$p_power_fdr <= 0.05, na.rm = TRUE),
                " with a significant power change")
        trans[[length(trans) + 1]] <- t
      }
      trans <- if (length(trans)) do.call(rbind, trans) else NULL
      if (!is.null(trans)) {
        write_tsv_tsf(trans, file.path(x$paths$base, "comparison",
                                       sprintf("transitions_%s.tsv", branch)))
      }
      list(signature = sig, conditions = conditions, transitions = trans)
    }, mc.cores = workers_now, mc.preschedule = TRUE, mc.set.seed = FALSE)
    names(per_dataset) <- ids

    failed <- vapply(per_dataset, inherits, logical(1), what = "try-error")
    if (any(failed)) {
      detail <- paste(vapply(per_dataset[failed], as.character, character(1)),
                      collapse = "; ")
      tsf_abort("compare failed for ", paste(ids[failed], collapse = ", "),
                ": ", detail)
    }

    signature_by_ds <- lapply(per_dataset, `[[`, "signature")
    conditions_by_ds <- lapply(per_dataset, `[[`, "conditions")
    transitions_by_ds <- lapply(per_dataset, `[[`, "transitions")

    sig_cross <- cross_datasets(signature_by_ds, c("chr", "N", "k"))
    if (!is.null(sig_cross)) {
      write_tsv_tsf(sig_cross, file.path(compare_dir,
                                         sprintf("constant_signature_shared_%s.tsv", branch)))
    }
    # A global intersection would let a healthy-only cohort erase every
    # fibrosis comparison, and a late-stage cohort erase F0--F2.  Cross each
    # biological condition independently, using only its contributors.
    cond_cross <- cross_by_level(conditions_by_ds, "condition", available,
                                 c("chr", "N", "k", "condition"))
    write_cross_level_tables(cond_cross, compare_dir,
                             paste0("conditions_shared_", branch))

    transition_levels <- unique(unlist(lapply(transitions_by_ds, function(x)
      if (is.null(x)) character(0) else as.character(x$transition)),
      use.names = FALSE))
    trans_cross <- cross_by_level(transitions_by_ds, "transition",
                                  transition_levels,
                                  c("chr", "N", "k", "transition"))
    if (length(trans_cross)) {
      for (nm in names(trans_cross)) {
        contributing <- ids[vapply(transitions_by_ds, function(x)
          !is.null(x) && any(x$transition == nm), logical(1))]
        trans_cross[[nm]] <- add_replication_flags(trans_cross[[nm]], contributing)
        if (!is.null(trans_cross[[nm]]$replicated)) {
          n_replicated <- n_replicated +
            sum(trans_cross[[nm]]$replicated, na.rm = TRUE)
        }
      }
      write_cross_level_tables(trans_cross, compare_dir,
                               paste0("transitions_shared_", branch))
    }
  }
  sprintf("%s: %d replicated", group_name, n_replicated)
}

#' Cross tables independently at each biological level.
cross_by_level <- function(tables_by_dataset, level_col, levels_now, key_cols) {
  out <- list()
  for (level in levels_now) {
    subset <- lapply(tables_by_dataset, function(x) {
      if (is.null(x) || !nrow(x) || !level_col %in% colnames(x)) return(NULL)
      y <- x[!is.na(x[[level_col]]) & x[[level_col]] == level, , drop = FALSE]
      if (nrow(y)) y else NULL
    })
    keep <- vapply(subset, function(x) !is.null(x) && nrow(x) > 0, logical(1))
    if (sum(keep) < 2L) {
      tsf_log("  ", level_col, " ", level, ": ", sum(keep),
              " non-empty dataset(s); no replication table")
      next
    }
    crossed <- cross_datasets(subset[keep], key_cols)
    if (!is.null(crossed) && nrow(crossed)) out[[level]] <- crossed
  }
  out
}

#' Write one unambiguous table per condition or transition.
write_cross_level_tables <- function(tables, out_dir, stem) {
  if (!length(tables)) return(invisible(0L))
  safe <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)
  for (nm in names(tables)) {
    write_tsv_tsf(tables[[nm]], file.path(out_dir,
      sprintf("%s_%s.tsv", stem, safe(nm))))
  }
  invisible(length(tables))
}

# --- reference library -------------------------------------------------------
stage_reference <- function(project, opt) {
  ids <- stage_datasets(opt)
  fps <- list(); kept_datasets <- list()
  for (id in ids) {
    inp <- tsf_stage_inputs(project, id)
    tsf_log(id, ": fingerprinting ", nrow(inp$dataset$samples), " samples")
    f <- fingerprint_dataset(inp$dataset, inp$chrom_idx,
                             k_max = project$fingerprint$k_max %||% 64L,
                             features = project$fingerprint$features %||% "amplitude")
    if (is.null(f)) { tsf_warn("  no fingerprint for ", id); next }
    # Kept so coverage calibration can re-fingerprint from masked grids: the
    # loss that matters is loss of genes, which only the expression matrix and
    # the grid can reproduce.
    kept_datasets[[id]] <- list(dataset = inp$dataset, chrom_idx = inp$chrom_idx)
    write_tsv_tsf(cbind(f$labels, as.data.frame(f$matrix)),
                  file.path(inp$paths$base, "fingerprints", "fingerprints.tsv"))
    fps[[id]] <- f
  }
  if (!length(fps)) tsf_abort("No fingerprints could be built")

  # The CANONICAL annotation grid travels with the reference -- not the genes
  # the first dataset happened to observe. Using genes.tsv would drop genes seen
  # only in the other cohorts, make query coverage look better than it is, and
  # make the reference depend on the order of `ids`.
  grids <- stats::setNames(lapply(names(fps), function(id) load_grid(id, project)),
                           names(fps))
  canonical_grid <- assert_compatible_grids(grids)
  prov <- grids[[1]]$provenance
  ref <- build_reference(fps, target = project$fingerprint$target %||% "condition",
                         n_features = project$fingerprint$n_features %||% 500L,
                         n_masks = project$fingerprint$n_masks %||% 10L,
                         datasets = kept_datasets,
                         max_queries_per_mask = project$fingerprint$max_queries_per_mask %||% 25L,
                         threshold_policy = project$fingerprint$threshold_policy %||% "pooled",
                         grid = canonical_grid,
                         params = list(
                           k_max = project$fingerprint$k_max %||% 64L,
                           features = project$fingerprint$features %||% "amplitude",
                           gene_universe = prov$gene_universe %||% "all",
                           annotation = prov$annotation_file %||% project$annotation_file,
                           species = prov$species %||% NA_character_,
                           genome_build = prov$genome_build %||% NA_character_,
                           annotation_release = prov$annotation_release %||% NA_character_,
                           grid_digest = prov$grid_digest %||% NA_character_,
                           expression_unit = prov$expression_unit %||% NA_character_,
                           chrom_levels = paste(project$chrom_levels, collapse = ","),
                           datasets = paste(names(fps), collapse = ",")))
  ensure_dir(file.path(project$results_dir, "reference"))
  path <- file.path(project$results_dir, "reference", "reference.rds")
  saveRDS(ref, path)

  if (!is.null(ref$validation)) {
    write_tsv_tsf(ref$validation$predictions,
                  file.path(project$results_dir, "reference",
                            "out_of_cohort_predictions.tsv"))
    write_tsv_tsf(as.data.frame(ref$validation$confusion),
                  file.path(project$results_dir, "reference", "confusion_matrix.tsv"))
    print(ref$validation$confusion)
  }
  tsf_log(reference_status(ref))
  sprintf("reference over %d sample(s), %d class(es)",
          nrow(ref$labels), length(ref$model$classes))
}

#' Spectral window per chromosome: what the gap pattern alone produces.
#'
#' This is the diagnostic for a peak near the Nyquist limit. If the sampling
#' pattern has power at the same frequency as an observed peak, that peak can be
#' an alias of structure elsewhere rather than structure at its own frequency.
#' Run it before interpreting any high-frequency result.
stage_window <- function(project, opt) {
  rows <- list()
  for (id in stage_datasets(opt)) {
    inp <- tsf_stage_inputs(project, id)
    out_dir <- file.path(inp$paths$base, "window")
    ensure_dir(out_dir)
    tsf_log(id, ": spectral window over ", length(inp$chrom_idx), " chromosome(s)")

    for (chr_now in names(inp$chrom_idx)) {
      ci <- inp$chrom_idx[[chr_now]]
      w <- spectral_window(ci$t, ci$N)
      w$chr <- chr_now
      write_tsv_tsf(w, file.path(out_dir, sprintf("window_chr%s.tsv", chr_now)))

      # Where do the stable peaks of this chromosome sit in the window?
      for (cond in tsf_conditions(inp$conditions, opt)) {
        pk <- read_tsv_tsf(p_peaks(inp$paths, opt$branches[1], cond), required = FALSE)
        if (is.null(pk) || !nrow(pk)) next
        pk <- pk[pk$chr == chr_now, , drop = FALSE]
        if (!nrow(pk)) next
        for (i in seq_len(nrow(pk))) {
          j <- match(pk$k[i], w$k)
          rows[[length(rows) + 1]] <- data.frame(
            dataset = id, condition = cond, chr = chr_now,
            N = ci$N, k = pk$k[i], period = pk$period[i],
            coverage_pct = round(100 * ci$coverage, 1),
            window_power = w$window_power[j],
            window_rank = rank(-w$window_power)[j],
            window_pct = round(100 * rank(-w$window_power)[j] / nrow(w), 2),
            stringsAsFactors = FALSE)
        }
      }
    }
  }
  if (!length(rows)) return("no stable peak to place in the window")
  tab <- do.call(rbind, rows)
  tab <- tab[order(tab$window_rank), ]
  write_tsv_tsf(tab, file.path(project$results_dir, "comparison",
                               "peaks_vs_spectral_window.tsv"))
  suspicious <- sum(tab$window_pct <= 1, na.rm = TRUE)
  if (suspicious) {
    tsf_warn(suspicious, " peak(s) sit in the top 1% of the spectral window -- ",
             "treat those as sampling artefacts until shown otherwise")
  }
  print(utils::head(tab, 10), row.names = FALSE)
  sprintf("%d peak(s) placed, %d in the window's top 1%%", nrow(tab), suspicious)
}

stage_functions <- list(
  fetch     = stage_fetch,
  ingest    = stage_ingest,
  spectra   = stage_spectra,
  maxt      = stage_maxt,
  condition = stage_condition,
  consensus = stage_consensus,
  clean     = stage_clean,
  stability = stage_stability,
  peaks     = stage_peaks,
  compare   = stage_compare,
  window    = stage_window,
  reference = stage_reference
)

#' What exists on disk for each dataset and stage.
pipeline_status <- function(project, opt) {
  rows <- list()
  for (id in stage_datasets(opt)) {
    p <- tsf_paths(project, id)
    interim <- file.path(project$interim_dir, id)
    samples <- read_tsv_tsf(file.path(interim, "samples.tsv"), required = FALSE)
    conds <- if (is.null(samples)) character(0) else unique(samples$condition)
    count <- function(dir, pattern) length(list.files(dir, pattern, recursive = TRUE))
    rows[[id]] <- data.frame(
      dataset = id,
      ingested = !is.null(samples),
      samples = if (is.null(samples)) 0L else nrow(samples),
      conditions = paste(conds, collapse = ","),
      spectra = count(p$spectra, "^spectra_condition_.*\\.tsv$"),
      maxt = count(p$maxt, "^maxt_individual_.*\\.tsv$"),
      stability = count(p$stability, "^maxt_stability_.*\\.tsv$"),
      peak_tables = count(p$peaks, "^peaks_.*\\.tsv$"),
      peak_gene_files = count(p$peak_genes, "^pico_chr.*\\.tsv$"),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}
