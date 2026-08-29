# config.R -- configuration is data, not code scattered through the pipeline.
#
# Configs are plain R files that return a list, so there is no yaml/jsonlite
# dependency on the cluster and they can be validated and unit-tested.

TSF_CHROM_LEVELS <- c(as.character(1:22), "X", "Y")

load_project_config <- function(path = "config/project.R") {
  cfg <- source(path, local = TRUE)$value
  required <- c("geo_dir", "interim_dir", "results_dir")
  missing <- setdiff(required, names(cfg))
  if (length(missing)) tsf_abort("project config missing: ", paste(missing, collapse = ", "))
  cfg$chrom_levels <- cfg$chrom_levels %||% TSF_CHROM_LEVELS
  cfg
}

load_dataset_config <- function(dataset_id, dir = "config/datasets") {
  path <- file.path(dir, paste0(dataset_id, ".R"))
  if (!file.exists(path)) tsf_abort("No config for dataset ", dataset_id, " at ", path)
  cfg <- source(path, local = TRUE)$value
  cfg <- attach_vocabulary(cfg)
  validate_dataset_config(cfg, dataset_id)
}

#' Load the condition vocabulary a dataset declares.
#'
#' A dataset may name a vocabulary file, override its levels, or (for
#' tissue-atlas style designs) supply condition_levels directly.
load_vocabulary <- function(name, dir = "config/vocabularies") {
  path <- file.path(dir, paste0(name, ".R"))
  if (!file.exists(path)) {
    tsf_abort("No vocabulary '", name, "' at ", path, ". Available: ",
              paste(sub("\\.R$", "", list.files(dir, pattern = "\\.R$")),
                    collapse = ", "))
  }
  source(path, local = TRUE)$value
}

attach_vocabulary <- function(cfg) {
  voc <- load_vocabulary(cfg$vocabulary %||% TSF_DEFAULT_VOCABULARY)
  if (!is.null(cfg$condition_levels)) voc$levels <- cfg$condition_levels
  if (!is.null(cfg$tissue)) voc$tissue <- cfg$tissue
  if (is.null(voc$levels) || !length(voc$levels)) {
    tsf_abort("Vocabulary '", voc$id, "' declares no levels and dataset ", cfg$id,
              " supplies no condition_levels")
  }
  cfg$vocabulary_spec <- voc
  cfg$tissue <- cfg$tissue %||% voc$tissue
  cfg
}

validate_dataset_config <- function(cfg, dataset_id) {
  # Required fields depend on where the data comes from and on whether the
  # vocabulary even has a baseline. A tissue atlas has no notion of control, so
  # demanding has_control_cohort from it would be asking a question that does
  # not apply.
  cfg$source <- cfg$source %||% "geo"
  required <- switch(cfg$source,
    geo = c("id", "tissue", "counts_file", "series_matrix", "condition_rules"),
    matrix = c("id", "tissue", "counts_file", "metadata_file"),
    tsf_abort("Unknown source '", cfg$source, "' in ", dataset_id,
              ". Supported: geo, matrix"))
  missing <- setdiff(required, names(cfg))
  if (length(missing)) {
    tsf_abort("Dataset config ", dataset_id, " (source: ", cfg$source,
              ") missing: ", paste(missing, collapse = ", "))
  }
  has_baseline <- !is.null(cfg$vocabulary_spec$baseline)
  if (has_baseline && is.null(cfg$has_control_cohort)) {
    tsf_abort("Vocabulary '", cfg$vocabulary_spec$id, "' has a baseline (",
              cfg$vocabulary_spec$baseline, "), so ", dataset_id,
              " must declare has_control_cohort. That single field is what keeps ",
              "the baseline label from meaning two different things across cohorts.")
  }
  if (!has_baseline) cfg$has_control_cohort <- cfg$has_control_cohort %||% FALSE
  if (!identical(cfg$id, dataset_id)) {
    tsf_abort("Config id (", cfg$id, ") does not match filename (", dataset_id, ")")
  }
  if (!is.logical(cfg$has_control_cohort) || is.na(cfg$has_control_cohort)) {
    tsf_abort("has_control_cohort must be TRUE or FALSE in ", dataset_id,
               " -- this is the field that keeps Control from meaning two things.")
  }
  ids <- vapply(cfg$condition_rules, function(r) r$id %||% "", character(1))
  if (any(!nzchar(ids)) || anyDuplicated(ids)) {
    tsf_abort("condition_rules need unique, non-empty ids in ", dataset_id)
  }
  cfg$sample_id_column <- cfg$sample_id_column %||% "geo_accession"
  cfg
}

tsf_source_all <- function(dir = "R") {
  files <- sort(list.files(dir, pattern = "\\.R$", full.names = TRUE))
  # utils_io first: everything else uses its helpers.
  files <- c(files[grepl("utils_io\\.R$", files)], files[!grepl("utils_io\\.R$", files)])
  invisible(lapply(files, source))
}
