# ingest.R -- turn any GEO dataset into the common on-disk format.
#
# Common format, written to <interim_dir>/<dataset_id>/ :
#
#   samples.tsv      sample_id, dataset_id, condition, fibrosis_stage, cohort, ...
#   genes.tsv        gene_id, gene_name, chr, start, gene_order
#   counts.tsv       gene_id + one integer column per sample
#   expression.tsv   gene_id + one asinh(TPM) column per sample
#   label_audit.tsv  every input sample, resolved or not, with the rule that fired
#   manifest.tsv     versions, config fingerprint, timestamp, counts per condition
#
# Everything downstream (spectra, maxT, comparison) reads only this format and
# never touches GEO again, so adding a third dataset means adding one config
# file, not copying 4,000 lines.

#' Read the phenotype table of a series matrix.
#'
#' Uses GEOquery when available, and falls back to parsing the !Sample_ lines
#' directly so the ingest layer stays testable without Bioconductor.
read_series_pheno <- function(path) {
  if (requireNamespace("GEOquery", quietly = TRUE)) {
    gse <- GEOquery::getGEO(filename = path, getGPL = FALSE)
    ph <- Biobase::pData(gse)
    ph[] <- lapply(ph, as.character)
    return(ph)
  }
  lfft_warn("GEOquery not installed; using the minimal series-matrix parser")
  parse_series_matrix_pheno(path)
}

parse_series_matrix_pheno <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  lines <- readLines(con, warn = FALSE)
  lines <- lines[grepl("^!Sample_", lines)]
  if (!length(lines)) lfft_abort("No !Sample_ lines in ", path)

  parsed <- lapply(lines, function(l) {
    parts <- strsplit(l, "\t", fixed = TRUE)[[1]]
    list(key = sub("^!Sample_", "", parts[1]),
         values = gsub('^"|"$', "", parts[-1]))
  })
  keys <- vapply(parsed, function(x) x$key, character(1))
  # GEOquery names characteristics columns after their "field: " prefix; do the
  # same so dataset configs can refer to "fibrosis stage" either way.
  is_char <- grepl("^characteristics", keys)
  keys[is_char] <- vapply(parsed[is_char], function(x) {
    fields <- unique(sub(":.*$", "", x$values[grepl(":", x$values)]))
    if (length(fields) == 1L) paste0(trimws(fields), ":ch1") else "characteristics_ch1"
  }, character(1))
  keys <- make.unique(keys, sep = "_")
  n <- max(vapply(parsed, function(x) length(x$values), integer(1)))
  df <- as.data.frame(
    lapply(parsed, function(x) c(x$values, rep(NA_character_, n - length(x$values)))),
    stringsAsFactors = FALSE, col.names = keys)
  colnames(df) <- keys
  df
}

#' Read a GEO count table and map its identifiers to Ensembl gene ids.
read_counts <- function(path, id_type = "ENTREZID") {
  dt <- read_tsv_lfft(path)
  id_col <- colnames(dt)[1]
  colnames(dt)[1] <- "source_id"
  lfft_log("Counts: ", nrow(dt), " rows x ", ncol(dt) - 1, " samples (id column: ", id_col, ")")
  attr(dt, "id_type") <- id_type
  dt
}

#' Gene coordinates from the NCBI annotation table shipped with the GEO counts.
read_gene_annotation <- function(path, chrom_levels) {
  annot <- read_tsv_lfft(path)
  ncbi_to_chr <- setNames(
    c(as.character(1:22), "X", "Y"),
    c("NC_000001.11","NC_000002.12","NC_000003.12","NC_000004.12","NC_000005.10",
      "NC_000006.12","NC_000007.14","NC_000008.11","NC_000009.12","NC_000010.11",
      "NC_000011.10","NC_000012.12","NC_000013.11","NC_000014.9","NC_000015.10",
      "NC_000016.10","NC_000017.11","NC_000018.10","NC_000019.10","NC_000020.11",
      "NC_000021.9","NC_000022.11","NC_000023.11","NC_000024.10"))

  gene_id <- sub("\\..*$", "", as.character(annot$EnsemblGeneID))
  chr <- unname(ncbi_to_chr[trimws(as.character(annot$ChrAcc))])
  start <- suppressWarnings(as.numeric(annot$ChrStart))

  out <- data.frame(gene_id = gene_id,
                    gene_name = as.character(annot$Symbol),
                    entrez_id = as.character(annot$GeneID),
                    chr = chr, start = start,
                    stringsAsFactors = FALSE)
  out <- out[!is.na(out$gene_id) & nzchar(out$gene_id) &
             out$chr %in% chrom_levels & !is.na(out$start), ]
  out <- out[!duplicated(out$gene_id), ]
  lfft_log("Annotation: ", nrow(out), " genes with chromosome and start")
  out
}

#' Gene order along each chromosome -- the FFT input axis.
#'
#' Computed once here instead of being re-derived with order(start) at every
#' call site, which is how the ordering silently diverged between stages.
add_gene_order <- function(genes, chrom_levels, min_genes_per_chr = 8L) {
  genes <- genes[genes$chr %in% chrom_levels, ]
  genes <- genes[order(match(genes$chr, chrom_levels), genes$start, genes$gene_id), ]
  genes$gene_order <- stats::ave(seq_len(nrow(genes)), genes$chr, FUN = seq_along)
  keep_chr <- names(which(table(genes$chr) >= min_genes_per_chr))
  dropped <- setdiff(unique(genes$chr), keep_chr)
  if (length(dropped)) {
    lfft_warn("Chromosomes with < ", min_genes_per_chr, " genes dropped: ",
              paste(dropped, collapse = ", "))
  }
  genes[genes$chr %in% keep_chr, ]
}

#' TPM followed by asinh, matching the original transform.
counts_to_expression <- function(count_mat, gene_length = NULL) {
  count_mat[is.na(count_mat)] <- 0
  if (is.null(gene_length)) {
    # No length available: CPM is the honest fallback, and it is recorded in the
    # manifest so nobody assumes TPM downstream.
    scaled <- t(t(count_mat) / pmax(colSums(count_mat), 1)) * 1e6
    attr(scaled, "unit") <- "CPM"
  } else {
    rpk <- count_mat / (gene_length / 1000)
    scaled <- t(t(rpk) / pmax(colSums(rpk, na.rm = TRUE), 1)) * 1e6
    attr(scaled, "unit") <- "TPM"
  }
  unit <- attr(scaled, "unit")
  out <- asinh(scaled)
  attr(out, "unit") <- paste0("asinh(", unit, ")")
  out
}

#' Keep genes expressed in a reasonable fraction of samples.
filter_expressed <- function(expr_mat, min_value, min_fraction) {
  keep <- rowMeans(expr_mat >= asinh(min_value)) >= min_fraction
  lfft_log("Expression filter: ", sum(keep), "/", length(keep), " genes kept")
  expr_mat[keep, , drop = FALSE]
}

#' Full ingest for one dataset.
ingest_dataset <- function(dataset_id, project, dataset_dir = "config/datasets") {
  cfg <- load_dataset_config(dataset_id, dataset_dir)
  lfft_log("=== Ingesting ", cfg$id, " (", cfg$description, ") ===")

  out_dir <- file.path(project$interim_dir, cfg$id)
  ensure_dir(out_dir)

  # ---- phenotype and labels -------------------------------------------------
  pheno <- read_series_pheno(file.path(project$geo_dir, cfg$series_matrix))
  labels <- harmonize_conditions(pheno, cfg)
  write_tsv_lfft(labels, file.path(out_dir, "label_audit.tsv"))
  audit <- audit_labels(labels, cfg)
  samples <- labels[labels$keep, , drop = FALSE]
  samples$condition <- factor(samples$condition, levels = LFFT_CONDITION_LEVELS)

  # ---- counts and annotation ------------------------------------------------
  counts <- read_counts(file.path(project$geo_dir, cfg$counts_file), cfg$count_id_type)
  genes <- read_gene_annotation(file.path(project$geo_dir, project$annotation_file),
                                project$chrom_levels)

  key <- if (identical(cfg$count_id_type, "ENTREZID")) "entrez_id" else "gene_id"
  merged <- merge(counts, genes, by.x = "source_id", by.y = key)
  if (!nrow(merged)) lfft_abort("No count row mapped to the annotation for ", cfg$id)
  merged <- merged[!duplicated(merged$gene_id), ]

  sample_cols <- intersect(colnames(counts)[-1], samples$sample_id)
  missing_cols <- setdiff(samples$sample_id, colnames(counts))
  if (length(missing_cols)) {
    lfft_warn(length(missing_cols), " labelled sample(s) absent from the count matrix: ",
              paste(utils::head(missing_cols, 5), collapse = ", "))
    samples <- samples[samples$sample_id %in% sample_cols, , drop = FALSE]
  }
  if (!length(sample_cols)) {
    lfft_abort("Sample ids in the series matrix do not match count columns for ", cfg$id,
               " -- check quoting/prefixes in sample_id_column.")
  }

  count_mat <- as.matrix(merged[, sample_cols, drop = FALSE])
  storage.mode(count_mat) <- "double"
  rownames(count_mat) <- merged$gene_id

  gene_len <- if ("Length" %in% colnames(merged)) suppressWarnings(as.numeric(merged$Length)) else NULL
  expr_mat <- counts_to_expression(count_mat, gene_len)
  unit <- attr(expr_mat, "unit")
  expr_mat <- filter_expressed(expr_mat, project$min_tpm, project$min_fraction)

  genes_out <- merged[merged$gene_id %in% rownames(expr_mat),
                      c("gene_id", "gene_name", "chr", "start")]
  genes_out <- add_gene_order(genes_out, project$chrom_levels, project$min_genes_per_chr)
  expr_mat <- expr_mat[genes_out$gene_id, , drop = FALSE]

  # ---- write the common format ---------------------------------------------
  write_tsv_lfft(samples, file.path(out_dir, "samples.tsv"))
  write_tsv_lfft(genes_out, file.path(out_dir, "genes.tsv"))
  write_tsv_lfft(data.frame(gene_id = rownames(count_mat), count_mat[, sample_cols, drop = FALSE],
                            check.names = FALSE),
                 file.path(out_dir, "counts.tsv"))
  write_tsv_lfft(data.frame(gene_id = rownames(expr_mat), expr_mat, check.names = FALSE),
                 file.path(out_dir, "expression.tsv"))
  write_run_manifest(file.path(out_dir, "manifest.tsv"), cfg$id, list(project, cfg),
                     extra = list(
                       n_samples = nrow(samples),
                       n_genes = nrow(genes_out),
                       expression_unit = unit,
                       has_control_cohort = cfg$has_control_cohort,
                       conditions_present = audit$present,
                       samples_per_condition = paste(names(audit$counts), audit$counts,
                                                     sep = "=", collapse = " ")))

  lfft_log("Wrote common format to ", out_dir)
  invisible(list(dir = out_dir, samples = samples, genes = genes_out,
                 expression = expr_mat, audit = audit))
}

#' Load a previously ingested dataset in the common format.
load_dataset <- function(dataset_id, project) {
  dir <- file.path(project$interim_dir, dataset_id)
  samples <- read_tsv_lfft(file.path(dir, "samples.tsv"))
  genes   <- read_tsv_lfft(file.path(dir, "genes.tsv"))
  expr    <- read_tsv_lfft(file.path(dir, "expression.tsv"))
  mat <- as.matrix(expr[, setdiff(colnames(expr), "gene_id"), drop = FALSE])
  rownames(mat) <- expr$gene_id
  samples$condition <- factor(samples$condition, levels = LFFT_CONDITION_LEVELS)
  list(id = dataset_id, samples = samples, genes = genes, expression = mat)
}
