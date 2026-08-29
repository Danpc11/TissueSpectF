# utils_io.R -- I/O primitives shared by every stage of the pipeline.
#
# Design note: the original scripts used save_tsv(overwrite = FALSE), which
# silently kept files from previous runs. Here overwrite defaults to TRUE and a
# run manifest records what was written, so two runs can never be interleaved
# without leaving a trace.

tsf_log <- function(..., level = "INFO") {
  cat(sprintf("[%s] %-5s %s\n", format(Sys.time(), "%H:%M:%S"), level,
              paste0(..., collapse = "")))
  invisible(NULL)
}

tsf_warn <- function(...) tsf_log(..., level = "WARN")

tsf_abort <- function(...) {
  stop(paste0(..., collapse = ""), call. = FALSE)
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

#' Write a data.frame as TSV.
#'
#' @param overwrite TRUE by default. Set FALSE only for artefacts that are
#'   genuinely expensive and immutable, and say so at the call site.
write_tsv_tsf <- function(df, path, overwrite = TRUE) {
  ensure_dir(dirname(path))
  if (!overwrite && file.exists(path) && file.size(path) > 0L) {
    tsf_warn("skip (already exists): ", path)
    return(invisible(path))
  }
  utils::write.table(df, path, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = TRUE, na = "NA")
  invisible(path)
}

read_tsv_tsf <- function(path, required = TRUE) {
  if (!file.exists(path)) {
    if (required) tsf_abort("Missing required file: ", path)
    return(NULL)
  }
  utils::read.delim(path, sep = "\t", header = TRUE, check.names = FALSE,
                    stringsAsFactors = FALSE, quote = "", comment.char = "")
}

#' Strip GEO characteristics prefixes and stray quoting.
#'
#' GEO series matrices store values as "fibrosis stage: 3" and sometimes wrap
#' accessions in quotes. Both scripts needed this; only one of them had it.
clean_pheno_value <- function(x) {
  v <- as.character(x)
  v <- sub("^[^:]+:\\s*", "", v)
  v <- gsub("^['\"]+|['\"]+$", "", trimws(v))
  v[v %in% c("", "NA", "na", "null", "NULL", "-")] <- NA_character_
  v
}

#' Record what a run produced, so outputs are traceable to code + config.
write_run_manifest <- function(path, dataset_id, config, extra = list()) {
  entries <- c(
    list(
      dataset_id      = dataset_id,
      run_timestamp   = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      pipeline_version = TSF_VERSION,
      r_version       = paste(R.version$major, R.version$minor, sep = "."),
      config_digest   = digest_config(config)
    ),
    lapply(extra, function(x) paste(x, collapse = ", "))
  )
  df <- data.frame(key = names(entries),
                   value = vapply(entries, function(x) paste(as.character(x), collapse = ", "),
                                  character(1)),
                   stringsAsFactors = FALSE)
  write_tsv_tsf(df, path, overwrite = TRUE)
}

#' Cheap, dependency-free config fingerprint (no digest package required).
digest_config <- function(config) {
  txt <- paste(utils::capture.output(str(config, max.level = 4)), collapse = "\n")
  sum(utf8ToInt(txt) * seq_along(utf8ToInt(txt))) %% .Machine$integer.max
}

#' Load every module in R/, dependencies first.
#'
#' Sourcing a hand-written list of file names meant that adding a module and
#' forgetting to add it here produced "could not find function" at run time,
#' while the unit tests kept passing because they source what they need
#' explicitly. Globbing the directory removes that failure mode entirely: a file
#' that exists is loaded.
tsf_module_order <- function(dir = "R") {
  first <- c("utils_io", "config", "labels", "paths")
  files <- sort(sub("\\.R$", "", basename(list.files(dir, pattern = "\\.R$"))))
  ordered <- c(intersect(first, files), setdiff(files, first))
  file.path(dir, paste0(ordered, ".R"))
}

tsf_load_all <- function(dir = "R") {
  invisible(lapply(tsf_module_order(dir), source))
}

TSF_VERSION <- "TissueSpectF 0.1.0"
