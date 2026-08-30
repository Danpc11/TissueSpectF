#!/usr/bin/env Rscript
# prepare_ae_data.R -- export per-sample spectra for TissueSpect-AE.
#
#   ./tsf ae-prepare --interim-dir interim --results-dir results_pc
#
# This is the entire contract between the R pipeline and the Python model. It is
# a directory of plain TSVs plus a manifest, so the exchange is inspectable with
# `head` and diffable in review; nothing about the model can silently depend on
# an R object it cannot see.
#
# WHAT IS AND IS NOT EXPORTED
#
# Every frequency of every sample is exported, not a peak list. The model has to
# see the whole spectrum, including the parts the statistical pipeline calls
# non-significant: a fingerprint does not need any component to be individually
# significant, and filtering first would remove exactly what makes classes
# separable.
#
# Coverage travels as data, and the observed mask travels with it. A frequency
# that could not be estimated is marked, never filled. The Python side is
# required to zero those positions only AFTER building the mask, purely as array
# padding, and to exclude them from every loss.

suppressPackageStartupMessages({
  source("R/utils_io.R")
  tsf_load_all("R")
})

prepare_ae_data <- function(project, opt) {
  out_dir <- opt$output_dir %||% file.path(project$results_dir, "autoencoder", "data")
  ensure_dir(out_dir)
  ids <- stage_datasets(opt)
  k_max <- opt$k_max %||% project$fingerprint$k_max %||% 64L

  spectra <- list(); samples <- list(); provenance <- list()

  for (id in ids) {
    inp <- tsf_stage_inputs(project, id)
    prov <- tryCatch(load_grid(id, project)$provenance, error = function(e) list())
    provenance[[id]] <- prov
    tsf_log(id, ": exporting per-sample spectra (k <= ", k_max, ")")

    for (cond in inp$conditions) {
      sp <- read_tsv_tsf(p_spectra_samples(inp$paths, cond), required = FALSE)
      if (is.null(sp) || !nrow(sp)) {
        tsf_warn("  ", cond, ": no per-sample spectra; run the spectra stage")
        next
      }
      sp <- sp[sp$k >= 1 & sp$k <= k_max, , drop = FALSE]
      if (!nrow(sp)) next
      lab <- inp$dataset$samples
      i <- match(sp$sample, lab$sample_id)
      sp$dataset_id <- id
      sp$tissue <- if ("tissue" %in% colnames(lab)) as.character(lab$tissue)[i] else NA_character_
      sp$state <- if ("state" %in% colnames(lab)) as.character(lab$state)[i] else NA_character_
      sp$condition <- as.character(cond)
      sp$class_id <- if ("class_id" %in% colnames(lab)) as.character(lab$class_id)[i] else NA_character_
      keep <- intersect(c("sample", "dataset_id", "tissue", "state", "condition",
                          "class_id", "chr", "N", "k", "freq", "period",
                          "amplitude", "power", "phase", "coverage"), colnames(sp))
      spectra[[length(spectra) + 1]] <- sp[, keep, drop = FALSE]
    }

    s <- inp$dataset$samples
    if (!"class_id" %in% colnames(s)) {
      tsf_abort(id, "'s samples.tsv has no class_id column, so it predates the ",
                "tissue::state::condition scheme. Re-run ingest before exporting: ",
                "a model trained on conditions without their tissue and state ",
                "cannot be extended to another tissue later without retraining.")
    }
    cols <- intersect(c("sample_id", "dataset_id", "tissue", "state", "condition",
                        "class_id", "cohort", "keep"), colnames(s))
    samples[[id]] <- s[, cols, drop = FALSE]
  }

  if (!length(spectra)) tsf_abort("Nothing to export: no per-sample spectra found")
  sp_all <- do.call(rbind, spectra)
  names(sp_all)[names(sp_all) == "sample"] <- "sample_id"
  sm_all <- do.call(rbind, samples)

  # --- alignment: refuse before exporting, not after training ----------------
  # Every cohort has to share the annotation, gene universe, grid and expression
  # unit. A model trained across two grids would be learning the grids.
  fields <- c("species", "genome_build", "annotation_release", "gene_universe",
              "grid_digest", "expression_unit")
  ref <- provenance[[ids[1]]]
  for (id in ids[-1]) {
    for (f in fields) {
      a <- ref[[f]] %||% NA_character_; b <- provenance[[id]][[f]] %||% NA_character_
      if (!identical(as.character(a), as.character(b))) {
        tsf_abort("Cohorts ", ids[1], " and ", id, " disagree on ", f, " (", a,
                  " vs ", b, "). They cannot train one model: a feature named ",
                  "chr7_k12 would mean different things in each.")
      }
    }
  }

  # grid_N must be one value per chromosome across the whole export.
  gn <- unique(sp_all[, c("chr", "N")])
  dup <- gn$chr[duplicated(gn$chr)]
  if (length(dup)) {
    tsf_abort("Chromosome(s) ", paste(unique(dup), collapse = ", "),
              " carry more than one grid_N across cohorts. Re-run ingest with ",
              "one annotation before exporting.")
  }
  chrom_order <- gn$chr[order(match(gn$chr, project$chrom_levels))]

  write_tsv_tsf(sp_all, file.path(out_dir, "spectra.tsv"))
  write_tsv_tsf(sm_all, file.path(out_dir, "samples.tsv"))
  write_tsv_tsf(gn[match(chrom_order, gn$chr), ], file.path(out_dir, "chromosomes.tsv"))

  manifest <- data.frame(
    key = c("annotation", "species", "genome_build", "annotation_release",
            "gene_universe", "grid_digest", "expression_unit", "k_max",
            "channels", "chromosome_order", "dataset_ids", "class_ids",
            "n_samples", "n_rows", "pipeline_version", "r_version",
            "created_at", "seed"),
    value = c(ref$annotation_file %||% project$annotation_file,
              ref$species %||% NA, ref$genome_build %||% NA,
              ref$annotation_release %||% NA, ref$gene_universe %||% NA,
              ref$grid_digest %||% NA, ref$expression_unit %||% NA,
              k_max,
              "log_power,amplitude,cos_phase,sin_phase,coverage,observed_mask",
              paste(chrom_order, collapse = ","),
              paste(ids, collapse = ","),
              paste(sort(unique(sp_all$class_id)), collapse = ","),
              length(unique(sp_all$sample_id)), nrow(sp_all),
              TSF_VERSION, paste(R.version$major, R.version$minor, sep = "."),
              format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
              opt$seed %||% project$maxt$seed),
    stringsAsFactors = FALSE)
  write_tsv_tsf(manifest, file.path(out_dir, "manifest.tsv"))

  tsf_log("Exported ", nrow(sp_all), " rows for ",
          length(unique(sp_all$sample_id)), " samples, ",
          length(unique(sp_all$class_id)), " classes, ",
          length(chrom_order), " chromosomes to ", out_dir)
  per_class <- table(sm_all$class_id[sm_all$keep %in% TRUE])
  for (cl in names(per_class)) tsf_log("  ", cl, ": ", per_class[[cl]])
  invisible(out_dir)
}

# No auto-run block here on purpose. This file is `source`d by scripts/tsf.R,
# and a trailing "if not interactive, run" guard fires during that source() --
# re-reading config/project.R from scratch and taking the CLI's own command name
# as a dataset id, so every --interim-dir and --results-dir override is lost.
# The entry point is `./tsf ae-prepare`.
