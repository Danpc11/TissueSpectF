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
    ph <- tryCatch({
      gse <- GEOquery::getGEO(filename = path, getGPL = FALSE)
      p <- Biobase::pData(gse)
      p[] <- lapply(p, as.character)
      p
    }, error = function(e) {
      tsf_warn("GEOquery could not parse ", basename(path), " (", conditionMessage(e),
               "); falling back to the minimal parser")
      NULL
    })
    if (!is.null(ph) && nrow(ph)) return(ph)
  } else {
    tsf_warn("GEOquery not installed; using the minimal series-matrix parser")
  }
  parse_series_matrix_pheno(path)
}

parse_series_matrix_pheno <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  lines <- readLines(con, warn = FALSE)
  lines <- lines[grepl("^!Sample_", lines)]
  if (!length(lines)) tsf_abort("No !Sample_ lines in ", path)

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

#' Read a GEO count table.
#'
#' GEO publishes count matrices in whatever shape the submitters chose: tabs or
#' commas, Entrez or Ensembl or symbols in the first column, sometimes a second
#' header row carrying the real sample names above meaningless column labels.
#' Those are properties of a file, not of the analysis, so they are declared in
#' the dataset config rather than handled by a bespoke script per series.
#'
#' @param spec optional list with sep, id_column, symbol_column, sample_map_row,
#'   skip. `sample_map_row` names a second header row whose values are the
#'   sample identifiers the phenotype table uses; it is extracted, used to
#'   rename the columns, and removed BEFORE anything is coerced to numeric --
#'   left in place it turns every count column into character and the whole
#'   matrix silently becomes text.
read_counts <- function(path, id_type = "ENTREZID", spec = list()) {
  sep <- spec$sep %||% "\t"
  dt <- utils::read.delim(path, sep = sep, header = TRUE, check.names = FALSE,
                          stringsAsFactors = FALSE, quote = "\"",
                          comment.char = "", skip = spec$skip %||% 0L)
  sample_map <- NULL

  if (!is.null(spec$sample_map_row)) {
    row <- as.integer(spec$sample_map_row)
    if (row < 1L || row > nrow(dt)) {
      tsf_abort("sample_map_row = ", row, " is out of range for ", basename(path))
    }
    raw <- as.character(unlist(dt[row, ], use.names = FALSE))
    mapped <- raw

    # The descriptive row rarely matches the phenotype table verbatim. In
    # GSE142530 it reads "Control_Lille 389" and "Control_TPF 111484" while the
    # series matrix says "Control_389" and "Control_111484": a recruiting-site
    # token the accessions do not carry. Declared substitutions bring the two
    # into line, and every one of them is recorded next to the raw value in
    # count_column_map.tsv, so the correspondence can be checked rather than
    # trusted. A regular expression that silently rewrites sample names and
    # leaves no table behind is not a mapping, it is a guess.
    for (tr in spec$map_transform %||% list()) {
      mapped <- gsub(tr$pattern, tr$replacement %||% "", mapped, perl = TRUE)
    }
    mapped <- trimws(mapped)

    sample_map <- data.frame(count_column = colnames(dt), raw_label = raw,
                             mapped_id = mapped, stringsAsFactors = FALSE)
    dt <- dt[-row, , drop = FALSE]
    tsf_log("Counts: second header row consumed as the sample map",
            if (length(spec$map_transform %||% list()))
              paste0(" (", length(spec$map_transform), " substitution(s) applied)")
            else "")
  }

  id_col <- spec$id_column %||% 1L
  if (is.numeric(id_col)) id_col <- colnames(dt)[id_col]
  if (!id_col %in% colnames(dt)) {
    tsf_abort("counts id_column '", id_col, "' is not in ", basename(path),
              " (columns start: ", paste(utils::head(colnames(dt), 4),
                                         collapse = ", "), ")")
  }
  symbol_col <- spec$symbol_column
  if (is.numeric(symbol_col)) symbol_col <- colnames(dt)[symbol_col]
  value_cols <- setdiff(colnames(dt), c(id_col, symbol_col))

  # Columns a submitter left in the file but did not use. Naming them in the
  # config is better than letting them fail to match later, where the symptom is
  # an unexplained missing sample rather than a stated exclusion.
  # Matched against the column name AND the descriptive label, because a
  # submitter's "not.used" lives in the descriptive row while the column itself
  # is called RB_N3. Matching only the column name would leave it in, and the
  # symptom would be one sample that never maps to a phenotype.
  drop_cols <- spec$exclude_columns %||% character(0)
  if (length(drop_cols)) {
    pattern <- paste(drop_cols, collapse = "|")
    by_name <- grepl(pattern, value_cols, fixed = FALSE)
    by_label <- rep(FALSE, length(value_cols))
    if (!is.null(sample_map)) {
      lab <- sample_map$raw_label[match(value_cols, sample_map$count_column)]
      by_label <- !is.na(lab) & grepl(pattern, lab)
    }
    hit <- value_cols[by_name | by_label]
    if (length(hit)) {
      tsf_log("Counts: excluding declared column(s): ", paste(hit, collapse = ", "))
      value_cols <- setdiff(value_cols, hit)
      if (!is.null(sample_map)) sample_map <- sample_map[!sample_map$count_column %in% hit, ]
    }
  }
  for (cl in value_cols) dt[[cl]] <- suppressWarnings(as.numeric(dt[[cl]]))

  # The map describes SAMPLE columns. Leaving the id and symbol columns in it
  # puts two NA rows in count_column_map.tsv and makes any count of mapped
  # samples wrong by two.
  if (!is.null(sample_map)) {
    sample_map <- sample_map[sample_map$count_column %in% value_cols, , drop = FALSE]
  }

  out <- dt[, c(id_col, value_cols), drop = FALSE]
  colnames(out)[1] <- "source_id"
  out$source_id <- sub("\\..*$", "", trimws(as.character(out$source_id)))
  tsf_log("Counts: ", nrow(out), " rows x ", length(value_cols),
          " samples (id column: ", id_col, ", type: ", id_type, ")")
  attr(out, "id_type") <- id_type
  attr(out, "sample_map") <- sample_map
  out
}

#' Gene coordinates from the NCBI annotation table shipped with the GEO counts.
read_gene_annotation <- function(path, chrom_levels) {
  annot <- read_tsv_tsf(path)
  ncbi_to_chr <- setNames(
    c(as.character(1:22), "X", "Y"),
    c("NC_000001.11","NC_000002.12","NC_000003.12","NC_000004.12","NC_000005.10",
      "NC_000006.12","NC_000007.14","NC_000008.11","NC_000009.12","NC_000010.11",
      "NC_000011.10","NC_000012.12","NC_000013.11","NC_000014.9","NC_000015.10",
      "NC_000016.10","NC_000017.11","NC_000018.10","NC_000019.10","NC_000020.11",
      "NC_000021.9","NC_000022.11","NC_000023.11","NC_000024.10"))

  gene_id <- sub("\\..*$", "", as.character(annot$EnsemblGeneID))
  chr <- unname(ncbi_to_chr[trimws(as.character(annot$ChrAcc))])
  if (all(is.na(chr)) && "Chromosome" %in% colnames(annot)) {
    chr <- normalise_chrom_names(annot$Chromosome)
  }
  chr <- normalise_chrom_names(chr)
  start <- suppressWarnings(as.numeric(annot$ChrStart))

  out <- data.frame(gene_id = gene_id,
                    gene_name = as.character(annot$Symbol),
                    entrez_id = as.character(annot$GeneID),
                    chr = chr, start = start,
                    stringsAsFactors = FALSE)
  # Kept when the annotation provides them: transcript length enables TPM
  # instead of CPM, and gene biotype lets a run restrict the universe.
  if ("Length" %in% colnames(annot)) {
    out$gene_length <- suppressWarnings(as.numeric(annot$Length))
  }
  if ("GeneType" %in% colnames(annot)) {
    out$gene_type <- as.character(annot$GeneType)
  }
  out <- out[!is.na(out$gene_id) & nzchar(out$gene_id) &
             out$chr %in% chrom_levels & !is.na(out$start), ]
  out <- out[!duplicated(out$gene_id), ]
  tsf_log("Annotation: ", nrow(out), " genes with chromosome and start")
  out
}

#' DEPRECATED: rank among filtered genes. Kept only for reading old outputs.
#' The spectral axis is now the annotation grid; see build_reference_grid().
add_gene_order <- function(genes, chrom_levels, min_genes_per_chr = 8L) {
  genes <- genes[genes$chr %in% chrom_levels, ]
  genes <- genes[order(match(genes$chr, chrom_levels), genes$start, genes$gene_id), ]
  genes$gene_order <- stats::ave(seq_len(nrow(genes)), genes$chr, FUN = seq_along)
  keep_chr <- names(which(table(genes$chr) >= min_genes_per_chr))
  dropped <- setdiff(unique(genes$chr), keep_chr)
  if (length(dropped)) {
    tsf_warn("Chromosomes with < ", min_genes_per_chr, " genes dropped: ",
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
  bad <- rowSums(!is.finite(expr_mat)) > 0
  if (any(bad)) {
    tsf_warn(sum(bad), " gene(s) have non-finite expression and are dropped")
    expr_mat <- expr_mat[!bad, , drop = FALSE]
  }
  keep <- rowMeans(expr_mat >= asinh(min_value), na.rm = TRUE) >= min_fraction
  keep[is.na(keep)] <- FALSE
  tsf_log("Expression filter: ", sum(keep), "/", length(keep), " genes kept")
  expr_mat[keep, , drop = FALSE]
}

#' Labels straight from a metadata table, for source = "matrix".
#'
#' No rules to apply: the file already names the condition. Values outside the
#' declared vocabulary are still refused rather than silently accepted.
labels_from_metadata <- function(pheno, cfg) {
  id_col <- cfg$sample_id_column %||% "sample_id"
  cond_col <- cfg$condition_column %||% "condition"
  for (need in c(id_col, cond_col)) {
    if (!need %in% colnames(pheno)) {
      tsf_abort("metadata_file has no column '", need, "' (columns: ",
                paste(colnames(pheno), collapse = ", "), ")")
    }
  }
  levels_now <- tsf_levels(cfg$vocabulary_spec)
  cond <- clean_pheno_value(pheno[[cond_col]])
  bad <- !is.na(cond) & !(cond %in% levels_now)
  if (any(bad)) {
    tsf_warn(sum(bad), " sample(s) carry a condition outside vocabulary '",
             cfg$vocabulary_spec$id, "' (", paste(unique(cond[bad]), collapse = ", "),
             "); dropped")
    cond[bad] <- NA_character_
  }
  baseline <- cfg$vocabulary_spec$baseline
  data.frame(
    sample_id = clean_pheno_value(pheno[[id_col]]),
    dataset_id = cfg$id, tissue = cfg$tissue %||% NA_character_,
    vocabulary = cfg$vocabulary_spec$id,
    condition = cond,
    state = if (is.null(cfg$vocabulary_spec$states)) "healthy" else
      unname(cfg$vocabulary_spec$states[cond]),
    class_id = tsf_class_id(cfg$tissue, cond, cfg$vocabulary_spec),
    fibrosis_stage = NA_real_,
    fibrosis_stage_reported = NA_real_, label_rule = "metadata_column",
    label_mismatch = FALSE,
    cohort = tsf_cohort_role(cond, cfg$vocabulary_spec),
    label_resolved = !is.na(cond),
    in_scope = TRUE,
    condition_selected = TRUE,
    keep = !is.na(cond), stringsAsFactors = FALSE)
}

#' Full ingest for one dataset.
ingest_dataset <- function(dataset_id, project, dataset_dir = "config/datasets") {
  cfg <- load_dataset_config(dataset_id, dataset_dir)
  tsf_log("=== Ingesting ", cfg$id, " (", cfg$description, ") ===")

  out_dir <- file.path(project$interim_dir, cfg$id)
  ensure_dir(out_dir)

  # ---- phenotype and labels -------------------------------------------------
  # source = "geo"    : series matrix + condition rules
  # source = "matrix" : a plain metadata TSV with a column naming the condition,
  #                     for data that never went through GEO
  pheno <- if (identical(cfg$source, "matrix")) {
    read_tsv_tsf(file.path(project$geo_dir, cfg$metadata_file))
  } else {
    read_series_pheno(file.path(project$geo_dir, cfg$series_matrix))
  }
  labels <- if (identical(cfg$source, "matrix") && is.null(cfg$condition_rules)) {
    labels_from_metadata(pheno, cfg)
  } else {
    harmonize_conditions(pheno, cfg)
  }
  labels <- apply_sample_filter(labels, pheno, cfg)
  write_tsv_tsf(labels, file.path(out_dir, "label_audit.tsv"))
  audit <- audit_labels(labels, cfg)
  samples <- labels[labels$keep, , drop = FALSE]
  samples$condition <- factor(samples$condition, levels = tsf_levels(cfg$vocabulary_spec))

  # ---- counts and annotation ------------------------------------------------
  counts <- read_counts(file.path(project$geo_dir, cfg$counts_file),
                        cfg$count_id_type, cfg$counts_spec %||% list())
  sample_map <- attr(counts, "sample_map")
  if (!is.null(sample_map)) {
    keep_map <- sample_map$count_column %in% colnames(counts)
    colnames(counts)[match(sample_map$count_column[keep_map], colnames(counts))] <-
      sample_map$mapped_id[keep_map]
    # Written with the match recorded, so a failed mapping is visible in the
    # file rather than only as an empty intersection three steps later.
    sample_map$matched_phenotype <- sample_map$mapped_id %in% labels$sample_id
    write_tsv_tsf(sample_map, file.path(out_dir, "count_column_map.tsv"))
    n_matched <- sum(sample_map$matched_phenotype & keep_map)
    tsf_log("Count columns matched to phenotype: ", n_matched, "/",
            sum(keep_map))
    if (n_matched == 0L) {
      tsf_abort("No count column matches a phenotype sample id for ", cfg$id,
                ". The map is in ", file.path(out_dir, "count_column_map.tsv"),
                "; compare mapped_id there against samples.tsv and add a ",
                "counts_spec$map_transform entry.")
    }
  }
  genes <- read_annotation(project)

  key <- switch(cfg$count_id_type %||% "ENTREZID",
                ENTREZID = "entrez_id",
                ENSEMBL  = "gene_id",
                SYMBOL   = "gene_name",
                tsf_abort("count_id_type must be ENTREZID, ENSEMBL or SYMBOL; got '",
                          cfg$count_id_type, "'"))
  if (!key %in% colnames(genes)) {
    tsf_abort("The annotation has no ", key, " column, which count_id_type = ",
              cfg$count_id_type, " requires.")
  }
  merged <- merge(counts, genes, by.x = "source_id", by.y = key)
  if (!nrow(merged)) {
    tsf_abort("No count row mapped to the annotation for ", cfg$id,
              ". count_id_type is ", cfg$count_id_type, "; the file's first ",
              "identifiers are ", paste(utils::head(counts$source_id, 3),
                                        collapse = ", "), ".")
  }
  if (identical(cfg$count_id_type, "SYMBOL")) {
    dup <- sum(duplicated(merged$source_id))
    if (dup) tsf_log(dup, " duplicated symbol(s) collapsed to one gene each; ",
                     "symbols are not unique identifiers and this is the cost ",
                     "of a count table keyed on them")
  }

  # merge() names its output join column after by.x (source_id). When counts
  # are keyed directly on Ensembl, the annotation's gene_id is therefore
  # consumed by the join. Restore the canonical identifier BEFORE duplicate
  # removal: subsetting with duplicated(merged$gene_id) while gene_id is absent
  # produces a zero-row data frame and later looks like an expression filter
  # removed every gene.
  if (identical(cfg$count_id_type, "ENSEMBL") &&
      !"gene_id" %in% colnames(merged)) {
    merged$gene_id <- as.character(merged$source_id)
  }

  # The same merge behaviour applies to SYMBOL: gene_name is the right-hand
  # join key and is therefore consumed into source_id.  Restore it so every
  # common-format genes.tsv carries the annotation symbol, including cohorts
  # such as GSE276114 whose count matrix itself is keyed by gene symbol.
  if (identical(cfg$count_id_type, "SYMBOL") &&
      !"gene_name" %in% colnames(merged)) {
    merged$gene_name <- as.character(merged$source_id)
  }

  merged <- merged[!duplicated(merged$gene_id), ]
  # merge() consumed the join key into source_id; restore it under its own name
  # so a query keyed on Entrez ids can be matched without conversion.
  if (identical(cfg$count_id_type, "ENTREZID") && !"entrez_id" %in% colnames(merged)) {
    merged$entrez_id <- as.character(merged$source_id)
  }

  sample_cols <- intersect(colnames(counts)[-1], samples$sample_id)
  missing_cols <- setdiff(samples$sample_id, colnames(counts))
  if (length(missing_cols)) {
    tsf_warn(length(missing_cols), " labelled sample(s) absent from the count matrix: ",
              paste(utils::head(missing_cols, 5), collapse = ", "))
    samples <- samples[samples$sample_id %in% sample_cols, , drop = FALSE]
  }
  if (!length(sample_cols)) {
    tsf_abort("Sample ids in the series matrix do not match count columns for ", cfg$id,
               " -- check quoting/prefixes in sample_id_column.")
  }
  # A cohort that declares how many samples it should contribute gets checked.
  # Silently ingesting eleven where twelve were expected is how a count in a
  # manuscript stops matching the data behind it.
  if (!is.null(cfg$expected_n_samples)) {
    if (nrow(samples) != as.integer(cfg$expected_n_samples)) {
      tsf_abort(cfg$id, " yields ", nrow(samples), " samples but its config ",
                "declares ", cfg$expected_n_samples, ". Reconcile the two ",
                "before proceeding: see count_column_map.tsv and label_audit.tsv.")
    }
    tsf_log("Sample count matches the declared ", cfg$expected_n_samples)
  }

  count_mat <- as.matrix(merged[, sample_cols, drop = FALSE])
  storage.mode(count_mat) <- "double"
  rownames(count_mat) <- merged$gene_id

  gene_len <- if ("gene_length" %in% colnames(merged)) {
    suppressWarnings(as.numeric(merged$gene_length))
  } else {
    tsf_warn("No gene_length column in the annotation; expression will be CPM, not TPM")
    NULL
  }

  if (!is.null(gene_len)) {
    # A single NA or non-positive length would produce an all-NA row in the TPM
    # matrix, and rowMeans() downstream would then return NA for the filter.
    # Drop those genes outright: they keep their slot on the reference grid and
    # are simply unobserved, which is exactly how any other unmeasured gene is
    # treated. Never impute a length.
    valid <- is.finite(gene_len) & gene_len > 0
    if (!any(valid)) {
      tsf_warn("No usable gene length; expression will be CPM, not TPM")
      gene_len <- NULL
    } else if (any(!valid)) {
      tsf_log(sum(!valid), " gene(s) dropped for missing or non-positive length ",
              "(they stay on the grid as unobserved)")
      count_mat <- count_mat[valid, , drop = FALSE]
      gene_len <- gene_len[valid]
      merged <- merged[valid, , drop = FALSE]
    }
  }
  expr_mat <- counts_to_expression(count_mat, gene_len)
  unit <- attr(expr_mat, "unit")
  expr_mat <- filter_expressed(expr_mat, project$min_tpm, project$min_fraction)

  # ---- reference grid ------------------------------------------------------
  # The axis is every annotated gene of the allowed biotypes, NOT the genes that
  # survived the expression filter. Filtered-out genes keep their grid slot and
  # are simply unobserved; they are never imputed.
  grid <- build_reference_grid(genes, project$chrom_levels,
                               biotypes = project$gene_universe,
                               min_genes_per_chr = project$min_genes_per_chr)

  # entrez_id is kept because GEO count tables are commonly keyed on it, and a
  # query file has to be matchable without the user converting identifiers.
  keep_cols <- intersect(c("gene_id", "gene_name", "entrez_id", "chr", "start",
                           "gene_length", "gene_type"), colnames(merged))
  genes_out <- merged[merged$gene_id %in% rownames(expr_mat), keep_cols]
  genes_out <- merge(genes_out, grid[, c("gene_id", "grid_index", "grid_N")],
                     by = "gene_id")
  off_grid <- sum(!rownames(expr_mat) %in% genes_out$gene_id)
  if (off_grid) {
    tsf_log(off_grid, " expressed gene(s) are not on the reference grid ",
            "(biotype excluded or unannotated); dropped from the spectral axis")
  }
  genes_out <- genes_out[order(match(genes_out$chr, project$chrom_levels),
                               genes_out$grid_index), ]
  expr_mat <- expr_mat[genes_out$gene_id, , drop = FALSE]

  cov <- stats::aggregate(genes_out$gene_id, list(chr = genes_out$chr), length)
  cov$N <- grid$grid_N[match(cov$chr, grid$chr)]
  cov$coverage_pct <- round(100 * cov$x / cov$N, 1)
  tsf_log("Grid coverage: median ", stats::median(cov$coverage_pct), "%, range ",
          min(cov$coverage_pct), "-", max(cov$coverage_pct), "%")
  write_tsv_tsf(stats::setNames(cov, c("chr", "n_observed", "grid_N", "coverage_pct")),
                file.path(out_dir, "grid_coverage.tsv"))

  # ---- write the common format ---------------------------------------------
  # (expression unit is known by now: see counts_to_expression above)
  # The CANONICAL grid, i.e. the whole annotation universe, is written out
  # separately from genes.tsv. genes.tsv holds only the genes this dataset
  # observed; a reference built from it would silently inherit one cohort's
  # coverage as if it were the annotation.
  write_tsv_tsf(grid, file.path(out_dir, "grid.tsv"))
  write_tsv_tsf(data.frame(
    key = c("species", "genome_build", "annotation_release", "annotation_format",
            "gene_length_mode", "gene_universe", "annotation_file",
            "grid_genes", "grid_digest", "expression_unit"),
    value = c(project$species %||% NA_character_,
              project$genome_build %||% NA_character_,
              project$annotation_release %||% NA_character_,
              project$annotation_format %||% "ncbi",
              attr(genes, "length_mode") %||% NA_character_,
              project$gene_universe %||% "all",
              project$annotation_file,
              nrow(grid), grid_digest(grid), unit),
    stringsAsFactors = FALSE),
    file.path(out_dir, "grid_provenance.tsv"))

  write_tsv_tsf(samples, file.path(out_dir, "samples.tsv"))
  write_tsv_tsf(genes_out, file.path(out_dir, "genes.tsv"))
  write_tsv_tsf(data.frame(gene_id = rownames(count_mat), count_mat[, sample_cols, drop = FALSE],
                            check.names = FALSE),
                 file.path(out_dir, "counts.tsv"))
  write_tsv_tsf(data.frame(gene_id = rownames(expr_mat), expr_mat, check.names = FALSE),
                 file.path(out_dir, "expression.tsv"))
  write_run_manifest(file.path(out_dir, "manifest.tsv"), cfg$id, list(project, cfg),
                     extra = list(
                       n_samples = nrow(samples),
                       n_genes = nrow(genes_out),
                       grid_genes = nrow(grid),
                       gene_universe = project$gene_universe %||% "all",
                       expression_unit = unit,
                       has_control_cohort = cfg$has_control_cohort,
                       conditions_present = audit$present,
                       samples_per_condition = paste(names(audit$counts), audit$counts,
                                                     sep = "=", collapse = " ")))

  tsf_log("Wrote common format to ", out_dir)
  invisible(list(dir = out_dir, samples = samples, genes = genes_out,
                 expression = expr_mat, audit = audit))
}

#' Fingerprint of a grid: same genes, same order, same N.
#'
#' Two datasets may only be combined into one reference if this matches. It is
#' cheap and it is the thing that actually has to be identical -- comparing
#' annotation file names would pass for two different releases with the same
#' name and fail for the same release stored under different names.
grid_digest <- function(grid) {
  # Canonical serialisation: sorted by chromosome then grid index, so the digest
  # cannot depend on row order.
  g <- grid[order(as.character(grid$chr), grid$grid_index, grid$gene_id), ]
  key <- paste(g$gene_id, g$chr, g$grid_index, g$grid_N, sep = "|", collapse = "\n")

  # SHA-256 where available. The algorithm is recorded in the value, so two
  # digests are only ever compared when they were produced the same way; a
  # mismatch of algorithm falls back to comparing the grids themselves.
  if (requireNamespace("openssl", quietly = TRUE)) {
    return(paste0("sha256:", as.character(openssl::sha256(key))))
  }
  if (requireNamespace("digest", quietly = TRUE)) {
    return(paste0("sha256:", digest::digest(key, algo = "sha256", serialize = FALSE)))
  }
  tmp <- tempfile(); on.exit(unlink(tmp), add = TRUE)
  writeLines(key, tmp)
  sha <- suppressWarnings(tryCatch(
    system2("sha256sum", tmp, stdout = TRUE, stderr = FALSE), error = function(e) NULL))
  if (length(sha) && nzchar(sha[1])) {
    return(paste0("sha256:", sub("\\s.*$", "", sha[1])))
  }
  paste0("md5:", unname(tools::md5sum(tmp)))
}

#' Do two grids describe the same axis?
#'
#' Digests settle it when they were computed with the same algorithm; otherwise
#' the grids are compared directly rather than declared incompatible over a
#' missing package.
grids_identical <- function(ga, gb, digest_a = NULL, digest_b = NULL) {
  algo <- function(d) if (is.null(d) || is.na(d)) NA_character_ else sub(":.*$", "", d)
  if (!is.na(algo(digest_a)) && identical(algo(digest_a), algo(digest_b))) {
    return(identical(as.character(digest_a), as.character(digest_b)))
  }
  key <- function(g) {
    o <- g[order(as.character(g$chr), g$grid_index, g$gene_id), ]
    paste(o$gene_id, o$chr, o$grid_index, o$grid_N, sep = "|")
  }
  identical(key(ga), key(gb))
}

#' Read the canonical grid and its provenance for an ingested dataset.
load_grid <- function(dataset_id, project) {
  dir <- file.path(project$interim_dir, dataset_id)
  grid <- read_tsv_tsf(file.path(dir, "grid.tsv"), required = FALSE)
  if (is.null(grid)) {
    tsf_abort("No grid.tsv for ", dataset_id, " -- re-run the ingest stage. ",
              "The reference needs the annotation grid, not just the genes this ",
              "dataset observed.")
  }
  prov <- read_tsv_tsf(file.path(dir, "grid_provenance.tsv"), required = FALSE)
  p <- if (is.null(prov)) list() else
    stats::setNames(as.list(prov$value), prov$key)
  list(grid = grid, provenance = p)
}

#' Refuse to build one reference from incompatible grids.
assert_compatible_grids <- function(grids) {
  ids <- names(grids)
  fields <- c("species", "genome_build", "annotation_release", "gene_universe",
              "grid_digest")
  ref <- grids[[1]]$provenance
  for (id in ids[-1]) {
    p <- grids[[id]]$provenance
    if (!grids_identical(grids[[1]]$grid, grids[[id]]$grid,
                         ref$grid_digest, p$grid_digest)) {
      tsf_abort("Datasets ", ids[1], " and ", id, " do not share the same grid.")
    }
    for (f in setdiff(fields, "grid_digest")) {
      a <- ref[[f]] %||% NA_character_; b <- p[[f]] %||% NA_character_
      if (!identical(as.character(a), as.character(b))) {
        tsf_abort("Datasets ", ids[1], " and ", id, " disagree on ", f, " (",
                  a, " vs ", b, "). They cannot share a reference: a feature ",
                  "named chrX_k7 would mean different things in each.")
      }
    }
  }
  tsf_log("Grid compatibility: ", length(ids), " dataset(s) share ",
          nrow(grids[[1]]$grid), " grid genes (",
          ref$genome_build %||% "build unknown", ", ",
          ref$gene_universe %||% "all", ")")
  grids[[1]]$grid
}

#' Load a previously ingested dataset in the common format.
load_dataset <- function(dataset_id, project, dataset_dir = "config/datasets") {
  dir <- file.path(project$interim_dir, dataset_id)
  samples <- read_tsv_tsf(file.path(dir, "samples.tsv"))
  genes   <- read_tsv_tsf(file.path(dir, "genes.tsv"))
  expr    <- read_tsv_tsf(file.path(dir, "expression.tsv"))
  mat <- as.matrix(expr[, setdiff(colnames(expr), "gene_id"), drop = FALSE])
  rownames(mat) <- expr$gene_id
  # Condition order comes from the dataset's vocabulary, so transitions follow
  # the declared progression rather than alphabetical order.
  lv <- tryCatch(tsf_levels(load_dataset_config(dataset_id, dataset_dir)$vocabulary_spec),
                 error = function(e) sort(unique(samples$condition)))
  samples$condition <- factor(samples$condition, levels = lv[lv %in% samples$condition])
  list(id = dataset_id, samples = samples, genes = genes, expression = mat,
       tissue = if ("tissue" %in% colnames(samples)) samples$tissue[1] else NA_character_)
}
