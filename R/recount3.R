# recount3.R -- bring recount3 gene sums into the format `./tsf ingest` reads.
#
# WHY recount3
# ------------
# GTEx raw reads are controlled access, but recount3 publishes GTEx gene-level
# coverage sums openly, processed with the SAME pipeline (Monorail / STAR /
# GENCODE v26) it used for every SRA study it holds. A cohort taken from
# recount3 and a GTEx tissue taken from recount3 differ by biology and by
# library, not by aligner, annotation or quantifier. That is the only way a
# differential spectrum (sample minus tissue reference) is a statement about
# the sample rather than about two pipelines.
#
# WHAT IS WRITTEN
# ---------------
# For a project P (a GTEx tissue like LIVER, or an SRA study like SRP217231):
#   <geo_dir>/R3_<P>_reads.tsv.gz    Name + one column per sample, read counts
#   <geo_dir>/R3_<P>_pheno.tsv       sample_id, donor, tissue, plus the
#                                    exploded sample attributes (SRA) or the
#                                    GTEx tissue fields
#   config/datasets/R3_<P>.R         a dataset config, if none exists yet
# Then `./tsf ingest R3_<P>` runs the ordinary path: Ensembl ids matched to the
# project annotation, TPM from the annotation lengths, asinh, the grid.
#
# COUNTS FROM COVERAGE SUMS
# -------------------------
# recount3 stores, per gene, the sum of per-base coverage over the gene. A read
# of length L contributes L to that sum, so
#     reads ~= sum / average_input_read_length
# which is what recount3's own compute_read_counts() does. The read length is
# per sample in the metadata (recount_qc.star.average_input_read_length). Note
# that for paired-end data STAR reports the length of the PAIR, so the result
# is fragments, which is the right unit for TPM anyway.
#
# URLS
# ----
# base = https://duffel.rail.bio/recount3/human
#   GTEx  gene sums : {base}/data_sources/gtex/gene_sums/{last2}/{P}/gtex.gene_sums.{P}.G026.gz
#         metadata  : {base}/data_sources/gtex/metadata/{last2}/{P}/gtex.gtex.{P}.MD.gz
#                     {base}/data_sources/gtex/metadata/{last2}/{P}/gtex.recount_qc.{P}.MD.gz
#   SRA   gene sums : {base}/data_sources/sra/gene_sums/{last2}/{P}/sra.gene_sums.{P}.G026.gz
#         metadata  : {base}/data_sources/sra/metadata/{last2}/{P}/sra.sra.{P}.MD.gz
#                     {base}/data_sources/sra/metadata/{last2}/{P}/sra.recount_qc.{P}.MD.gz
# where {last2} is the last two characters of P. These are the paths the
# recount3 Bioconductor package builds; they are checked with a HEAD request
# before anything is downloaded, and a missing project is reported, not guessed.

RECOUNT3_BASE <- Sys.getenv("TSF_RECOUNT3_BASE", "https://duffel.rail.bio/recount3/human")

recount3_last2 <- function(p) substr(p, nchar(p) - 1L, nchar(p))

recount3_urls <- function(project, source = c("gtex", "sra")) {
  source <- match.arg(source)
  l2 <- recount3_last2(project)
  base <- RECOUNT3_BASE
  list(
    gene_sums = sprintf("%s/data_sources/%s/gene_sums/%s/%s/%s.gene_sums.%s.G026.gz",
                        base, source, l2, project, source, project),
    metadata  = sprintf("%s/data_sources/%s/metadata/%s/%s/%s.%s.%s.MD.gz",
                        base, source, l2, project, source, source, project),
    qc        = sprintf("%s/data_sources/%s/metadata/%s/%s/%s.recount_qc.%s.MD.gz",
                        base, source, l2, project, source, project)
  )
}

#' TRUE if the URL answers 200 to a HEAD request. Uses curl if present, else
#' url() -- both without downloading the body.
recount3_url_exists <- function(url) {
  if (nzchar(Sys.which("curl"))) {
    code <- suppressWarnings(system2("curl", c("-s", "-o", "/dev/null", "-I", "-L",
                                               "-w", "%{http_code}", shQuote(url)),
                                     stdout = TRUE))
    return(identical(code[length(code)], "200"))
  }
  ok <- tryCatch({ con <- url(url, "rb"); close(con); TRUE }, error = function(e) FALSE)
  ok
}

#' Which source holds this project, or NA.
recount3_locate <- function(project) {
  for (src in c("gtex", "sra")) {
    if (recount3_url_exists(recount3_urls(project, src)$gene_sums)) return(src)
  }
  NA_character_
}

recount3_download <- function(url, dest) {
  if (file.exists(dest) && file.size(dest) > 0) return(dest)
  ensure_dir(dirname(dest))
  tsf_log("  downloading ", basename(url))
  ok <- tryCatch(utils::download.file(url, dest, mode = "wb", quiet = TRUE) == 0L,
                 error = function(e) FALSE)
  if (!ok || !file.exists(dest) || file.size(dest) == 0) {
    unlink(dest)
    tsf_abort("download failed: ", url)
  }
  dest
}

read_recount3_md <- function(path) {
  utils::read.delim(path, sep = "\t", header = TRUE, check.names = FALSE,
                    stringsAsFactors = FALSE, quote = "", comment.char = "",
                    na.strings = c("", "NA"))
}

#' Explode recount3's `sra.sample_attributes` ("key;;value|key;;value") into
#' columns, so a dataset config can use column_match / fibrosis_stage rules
#' on them exactly as it does on a GEO series matrix.
explode_sample_attributes <- function(x) {
  x[is.na(x)] <- ""
  parsed <- lapply(strsplit(x, "|", fixed = TRUE), function(pairs) {
    kv <- strsplit(pairs, ";;", fixed = TRUE)
    kv <- kv[vapply(kv, length, integer(1)) == 2L]
    stats::setNames(vapply(kv, `[`, "", 2L), tolower(trimws(vapply(kv, `[`, "", 1L))))
  })
  keys <- unique(unlist(lapply(parsed, names)))
  keys <- keys[nzchar(keys)]
  if (!length(keys)) return(NULL)
  out <- as.data.frame(lapply(keys, function(k) {
    vapply(parsed, function(p) if (k %in% names(p)) p[[k]] else NA_character_, "")
  }), stringsAsFactors = FALSE, check.names = FALSE)
  colnames(out) <- keys
  out
}

#' Fetch one recount3 project and write the ingest inputs.
#'
#' @param project  "LIVER" (GTEx) or "SRP217231" (SRA); case matters.
#' @param project_dir  where geo_dir and config/datasets live (TSF_ROOT).
#' @param tissue  tissue label for the config (e.g. "liver").
#' @param vocabulary  vocabulary id for the config.
#' @param max_samples  optional cap (GTEx LIVER has ~250 samples; fine as is).
fetch_recount3_project <- function(project, geo_dir, cache_dir = file.path(geo_dir, "recount3_cache"),
                                   tissue = NA_character_, vocabulary = "case_control",
                                   max_samples = NULL, write_config = TRUE,
                                   config_dir = "config/datasets") {
  src <- recount3_locate(project)
  if (is.na(src)) {
    tsf_abort("recount3 does not hold a project named '", project, "' under gtex or sra. ",
              "For a GEO series, resolve its SRP with scripts/resolve_recount3.R; ",
              "if the SRP is not in recount3 the cohort must come from GEO counts ",
              "(the ordinary ./tsf fetch path) and the pipeline difference must be ",
              "declared in the write-up.")
  }
  u <- recount3_urls(project, src)
  id <- paste0("R3_", project)
  tsf_log("=== recount3 ", src, "/", project, " -> ", id, " ===")

  f_sums <- recount3_download(u$gene_sums, file.path(cache_dir, basename(u$gene_sums)))
  f_md   <- recount3_download(u$metadata,  file.path(cache_dir, basename(u$metadata)))
  f_qc   <- recount3_download(u$qc,        file.path(cache_dir, basename(u$qc)))

  # gene sums: header lines start with "##", then a tab table gene_id x sample
  sums <- utils::read.delim(f_sums, sep = "\t", header = TRUE, check.names = FALSE,
                            stringsAsFactors = FALSE, comment.char = "#", quote = "")
  gene_col <- colnames(sums)[1]
  gene_ids <- sub("\\..*$", "", as.character(sums[[gene_col]]))
  mat <- as.matrix(sums[, -1, drop = FALSE])
  rownames(mat) <- gene_ids
  # PAR_Y copies and versioned duplicates: keep the first occurrence.
  mat <- mat[!duplicated(rownames(mat)), , drop = FALSE]

  md <- read_recount3_md(f_md)
  qc <- read_recount3_md(f_qc)
  key <- if ("external_id" %in% colnames(md)) "external_id" else colnames(md)[1]
  md <- md[match(colnames(mat), md[[key]]), , drop = FALSE]
  qc <- qc[match(colnames(mat), qc[[if ("external_id" %in% colnames(qc)) "external_id" else colnames(qc)[1]]]), , drop = FALSE]

  rl_col <- grep("average_input_read_length$", colnames(qc), value = TRUE)[1]
  if (is.na(rl_col)) tsf_abort("recount_qc metadata lacks average_input_read_length")
  read_len <- suppressWarnings(as.numeric(qc[[rl_col]]))
  bad <- !is.finite(read_len) | read_len <= 0
  if (any(bad)) {
    tsf_warn(sum(bad), " sample(s) without a read length; dropped")
    mat <- mat[, !bad, drop = FALSE]; md <- md[!bad, , drop = FALSE]; read_len <- read_len[!bad]
  }
  counts <- round(sweep(mat, 2, read_len, "/"))
  storage.mode(counts) <- "integer"

  # phenotype table
  if (src == "gtex") {
    pheno <- data.frame(
      sample_id = colnames(counts),
      donor     = md[["gtex.subjid"]] %||% sub("^(GTEX-[^-]+).*", "\\1", colnames(counts)),
      tissue    = md[["gtex.smtsd"]] %||% project,
      tissue_group = md[["gtex.smts"]] %||% NA_character_,
      sex       = md[["gtex.sex"]] %||% NA,
      age       = md[["gtex.age"]] %||% NA,
      rin       = md[["gtex.smrin"]] %||% NA,
      condition = "Control_external_study",
      stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    attrs <- explode_sample_attributes(md[["sra.sample_attributes"]])
    pheno <- data.frame(
      sample_id  = colnames(counts),
      run        = md[["sra.run_acc"]] %||% colnames(counts),
      experiment = md[["sra.experiment_acc"]] %||% NA_character_,
      biosample  = md[["sra.sample_acc.x"]] %||% md[["sra.sample_acc"]] %||% NA_character_,
      sample_title = md[["sra.sample_title"]] %||% NA_character_,
      donor      = md[["sra.sample_acc.x"]] %||% colnames(counts),
      stringsAsFactors = FALSE, check.names = FALSE)
    if (!is.null(attrs)) pheno <- cbind(pheno, attrs)
    if (!"condition" %in% colnames(pheno)) pheno$condition <- NA_character_
  }

  if (!is.null(max_samples) && ncol(counts) > max_samples) {
    set.seed(1L)
    # donor-disjoint by construction here: one sample per donor in GTEx tissues
    keep <- sort(sample.int(ncol(counts), max_samples))
    counts <- counts[, keep, drop = FALSE]; pheno <- pheno[keep, , drop = FALSE]
  }

  ensure_dir(geo_dir)
  reads_file <- paste0(id, "_reads.tsv.gz")
  pheno_file <- paste0(id, "_pheno.tsv")
  con <- gzfile(file.path(geo_dir, reads_file), "w")
  utils::write.table(data.frame(Name = rownames(counts), counts, check.names = FALSE),
                     con, sep = "\t", quote = FALSE, row.names = FALSE)
  close(con)
  utils::write.table(pheno, file.path(geo_dir, pheno_file), sep = "\t",
                     quote = FALSE, row.names = FALSE, na = "")
  tsf_log("  ", ncol(counts), " sample(s), ", nrow(counts), " gene(s) -> ", reads_file)

  if (write_config) {
    cfg_path <- file.path(config_dir, paste0(id, ".R"))
    if (!file.exists(cfg_path)) {
      write_recount3_config(cfg_path, id, project, src, reads_file, pheno_file,
                            tissue, vocabulary, pheno)
      tsf_log("  wrote ", cfg_path, if (src == "sra")
        " -- EDIT its condition_rules before ingesting: the sample attributes are exploded into columns of the pheno file" else "")
    } else {
      tsf_log("  config exists, left untouched: ", cfg_path)
      old <- tryCatch(source(cfg_path, local = TRUE)$value, error = function(e) NULL)
      if (is.null(old) || !identical(old$source, "matrix") || is.null(old$metadata_file) ||
          !identical(old$count_id_type, "ENSEMBL")) {
        tsf_warn("  ", cfg_path, " was written by an earlier generator and ingest will not read it ",
                 "correctly. Run: Rscript scripts/validate_recount3_configs.R --migrate")
      }
    }
  }
  invisible(list(id = id, source = src, reads = file.path(geo_dir, reads_file),
                 pheno = file.path(geo_dir, pheno_file), n_samples = ncol(counts)))
}

write_recount3_config <- function(path, id, project, src, reads_file, pheno_file,
                                  tissue, vocabulary, pheno) {
  ensure_dir(dirname(path))
  tissue_r <- if (is.na(tissue)) "NA_character_" else sprintf('"%s"', tissue)
  attr_cols <- setdiff(colnames(pheno), c("sample_id", "run", "experiment", "biosample",
                                          "sample_title", "donor", "condition"))
  rules <- if (src == "gtex") {
'  condition_rules = list(
    list(id = "gtex_healthy", type = "column_match", column = "^condition$",
         values = c("Control_external_study"), assign = "Control_external_study")
  ),'
  } else {
sprintf('  # recount3 exploded the SRA sample attributes into these columns:
  #   %s
  # Write rules on them as for a GEO series matrix (column_match, fibrosis_stage).
  # Until then every sample is unlabelled and ingest keeps none.
  condition_rules = list(
    # list(id = "biopsy_fibrosis_stage", type = "fibrosis_stage",
    #      column = "fibrosis stage", normal_terms = c("normal"))
  ),', paste(attr_cols, collapse = ", "))
  }
  txt <- sprintf('# Generated by R/recount3.R from recount3 %s/%s. Same pipeline as every other
# recount3 project (Monorail, GENCODE v26 gene sums -> reads via read length).
list(
  id            = "%s",
  tissue        = %s,
  vocabulary    = "%s",
  description   = "recount3 %s/%s",
  counts_file   = "%s",
  source        = "matrix",
  metadata_file = "%s",
  sample_id_column = "sample_id",
  counts_spec   = list(sep = "\\t", id_column = "Name"),
  count_id_type = "ENSEMBL",
  expression_unit = "counts",
  has_control_cohort = %s,
  donor_column  = "donor",
%s
  notes = "%s"
)
', src, project, id, tissue_r, vocabulary, src, project, reads_file, pheno_file,
    if (src == "gtex") "TRUE" else "FALSE", rules,
    if (src == "gtex") "Postmortem tissue: the tissue reference, never a clinical control group."
    else "Check the pheno columns and complete condition_rules before ingest.")
  writeLines(txt, path)
}
