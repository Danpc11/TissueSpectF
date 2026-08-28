# selfcheck.R -- run the whole pipeline on synthetic data with a known answer.
#
# This is the forward-looking replacement for diffing against the old scripts.
# It builds a small GEO-shaped dataset in a temporary directory, injects one
# sinusoid on chromosome 1 whose amplitude grows monotonically with fibrosis
# stage, runs ingest -> compare, and asserts that:
#
#   1. labels resolve as declared (Control only where the config allows it)
#   2. the recovered peak is the injected (chr, N, k)
#   3. the recovered phase matches the injected phase
#   4. amplitude increases across every transition, in both datasets
#   5. the change replicates across datasets
#
# It exercises the same code the real run uses, so a regression in the FFT,
# the permutation test, the stability threshold or the comparison shows up here.
# It cannot tell you whether the biology is real -- only that the machinery
# recovers a signal it is known to contain.

SELFCHECK_K <- 6L
SELFCHECK_PHASE <- 0.5
SELFCHECK_GENES_PER_CHR <- 200L

write_synthetic_geo <- function(dir, amplitude_by_condition) {
  ensure_dir(dir)
  n_chr <- 3L
  n_genes <- n_chr * SELFCHECK_GENES_PER_CHR
  accs <- c("NC_000001.11", "NC_000002.12", "NC_000003.12")

  gz_write <- function(lines, path) {
    con <- gzfile(path, "wt"); writeLines(lines, con); close(con)
  }

  annot <- c("GeneID\tSymbol\tEnsemblGeneID\tChrAcc\tChrStart",
             vapply(seq_len(n_genes), function(i) {
               sprintf("%d\tSYM%d\tENSG%011d.3\t%s\t%d",
                       2000L + i, i, i,
                       accs[(i - 1) %/% SELFCHECK_GENES_PER_CHR + 1],
                       1000L * ((i - 1) %% SELFCHECK_GENES_PER_CHR + 1))
             }, character(1)))
  gz_write(annot, file.path(dir, "Human.GRCh38.p13.annot.tsv.gz"))

  make_dataset <- function(id, conds, titles, characteristics) {
    samples <- sprintf("GSM%s%03d", substr(id, 4, 6), seq_along(conds))
    quoted <- function(v) paste0('"', v, '"')
    lines <- c(
      paste0("!Sample_title\t", paste(quoted(titles), collapse = "\t")),
      paste0("!Sample_geo_accession\t", paste(quoted(samples), collapse = "\t")),
      vapply(characteristics, function(ch)
        paste0("!Sample_characteristics_ch1\t", paste(quoted(ch), collapse = "\t")),
        character(1)))
    gz_write(lines, file.path(dir, sprintf("%s_series_matrix.txt.gz", id)))

    set.seed(42L)
    rows <- character(n_genes + 1L)
    rows[1] <- paste0("GeneID\t", paste(samples, collapse = "\t"))
    for (i in seq_len(n_genes)) {
      chrom <- (i - 1) %/% SELFCHECK_GENES_PER_CHR
      pos <- (i - 1) %% SELFCHECK_GENES_PER_CHR
      vals <- vapply(conds, function(cnd) {
        periodic <- if (chrom == 0) {
          amplitude_by_condition[[cnd]] * 80 *
            cos(2 * pi * SELFCHECK_K * pos / SELFCHECK_GENES_PER_CHR + SELFCHECK_PHASE)
        } else 0
        max(0L, as.integer(200 + periodic + stats::rnorm(1, sd = 12)))
      }, numeric(1))
      rows[i + 1L] <- paste0(2000L + i, "\t", paste(vals, collapse = "\t"))
    }
    gz_write(rows, file.path(dir, sprintf("%s.tsv.gz", id)))
  }

  conds1 <- c(rep("Control", 6), rep("0", 5), rep("1", 5), rep("2", 5),
              rep("3", 5), rep("4", 5))
  make_dataset("GSE135251", conds1,
               sprintf("sample %d", seq_along(conds1)),
               list(paste0("disease: ", ifelse(conds1 == "Control", "Control", "NAFLD")),
                    paste0("fibrosis stage: ", ifelse(conds1 == "Control", "0", conds1))))

  conds2 <- c(rep("N", 6), rep("F0", 5), rep("F1", 5), rep("F2", 5),
              rep("F3", 5), rep("F4", 5))
  make_dataset("GSE162694", conds2,
               sprintf("Liver %s %d", conds2, seq_along(conds2)),
               list(paste0("fibrosis stage: ",
                           ifelse(conds2 == "N", "normal liver histology",
                                  substr(conds2, 2, 2)))))
  invisible(dir)
}

run_selfcheck <- function() {
  tmp <- file.path(tempdir(), paste0("tsf_selfcheck_", Sys.getpid()))
  geo <- file.path(tmp, "data")
  tsf_log("selfcheck workspace: ", tmp)

  amps <- list(Control = 0.2, `0` = 0.3, `1` = 0.6, `2` = 1.0, `3` = 1.4, `4` = 1.8,
               N = 0.3, F0 = 0.3, F1 = 0.65, F2 = 1.05, F3 = 1.35, F4 = 1.75)
  write_synthetic_geo(geo, amps)

  project <- load_project_config("config/project.R")
  project$geo_dir <- geo
  project$interim_dir <- file.path(tmp, "interim")
  project$results_dir <- file.path(tmp, "results")
  project$maxt$B <- 100L   # enough to separate an injected peak from noise

  opt <- list(datasets = c("GSE135251", "GSE162694"), cond = NULL,
              branch = NULL, branches = "average", force = TRUE)

  for (stage in stage_names) {
    tsf_log("selfcheck stage: ", stage)
    stage_functions[[stage]](project, opt)
  }

  failures <- 0L
  check <- function(label, ok) {
    ok <- isTRUE(ok)
    cat(if (ok) "  PASS  " else "  FAIL  ", label, "\n")
    if (!ok) failures <<- failures + 1L
  }

  s1 <- read_tsv_tsf(file.path(project$interim_dir, "GSE135251", "samples.tsv"))
  s2 <- read_tsv_tsf(file.path(project$interim_dir, "GSE162694", "samples.tsv"))
  check("GSE135251 has a Control group", sum(s1$condition == "Control") == 6)
  check("GSE135251 keeps F0 separate from Control", sum(s1$condition == "F0") == 5)
  check("GSE162694 has no Control", !any(s2$condition == "Control"))
  check("GSE162694 normal histology lands in F0", sum(s2$condition == "F0") == 11)

  sig <- read_tsv_tsf(file.path(project$results_dir, "GSE135251", "comparison",
                                "constant_signature_average.tsv"), required = FALSE)
  check("a constant signature was found", !is.null(sig) && nrow(sig) >= 1)
  if (!is.null(sig) && nrow(sig)) {
    check("the signature peak is the injected one",
          any(sig$chr == "1" & sig$N == SELFCHECK_GENES_PER_CHR & sig$k == SELFCHECK_K))
  }

  peaks <- read_tsv_tsf(file.path(project$results_dir, "GSE135251", "peaks",
                                  "average", "peaks_F4.tsv"), required = FALSE)
  if (!is.null(peaks)) {
    row <- peaks[peaks$chr == "1" & peaks$k == SELFCHECK_K, ]
    check("the recovered phase matches the injected phase",
          nrow(row) == 1 && abs(row$phase[1] - SELFCHECK_PHASE) < 0.2)
  } else check("the recovered phase matches the injected phase", FALSE)

  for (id in c("GSE135251", "GSE162694")) {
    tr <- read_tsv_tsf(file.path(project$results_dir, id, "comparison",
                                 "transitions_average.tsv"), required = FALSE)
    ok <- !is.null(tr) && nrow(tr) >= 4 && all(tr$log2_amplitude_ratio > 0, na.rm = TRUE)
    check(paste0(id, ": amplitude increases across every transition"), ok)
  }

  shared <- read_tsv_tsf(file.path(project$results_dir, "comparison",
                                   "transitions_shared_average.tsv"), required = FALSE)
  check("the increase replicates across datasets",
        !is.null(shared) && sum(shared$replicated, na.rm = TRUE) >= 3)

  unlink(tmp, recursive = TRUE)
  if (failures == 0L) {
    tsf_log("selfcheck passed: the pipeline recovers the injected signal end to end.")
    return(0L)
  }
  tsf_warn("selfcheck: ", failures, " check(s) failed.")
  1L
}
