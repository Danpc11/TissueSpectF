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
  validate_dataset_config(cfg, dataset_id)
}

validate_dataset_config <- function(cfg, dataset_id) {
  required <- c("id", "counts_file", "series_matrix", "condition_rules",
                "has_control_cohort")
  missing <- setdiff(required, names(cfg))
  if (length(missing)) {
    tsf_abort("Dataset config ", dataset_id, " missing: ", paste(missing, collapse = ", "))
  }
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
