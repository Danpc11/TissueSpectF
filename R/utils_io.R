# utils_io.R -- I/O primitives shared by every stage of the pipeline.
#
# Design note: the original scripts used save_tsv(overwrite = FALSE), which
# silently kept files from previous runs. Here overwrite defaults to TRUE and a
# run manifest records what was written, so two runs can never be interleaved
# without leaving a trace.

#' Default for NULL. Defined here so any module works when sourced alone.
`%||%` <- function(a, b) if (is.null(a)) b else a

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
#' Is this file a module, or a script that happens to live in R/?
#'
#' tsf_load_all() sources every .R file in R/. That is fine for a module, which
#' only defines things, and destructive for a script, which DOES things at the
#' top level: a copy of scripts/clean_results.R placed in R/ deleted the results
#' tree on every single `./tsf` invocation, before argument parsing had even
#' run, so the deletion appeared in the log above an unrelated usage error.
#'
#' The signatures of a script are a shebang and a top-level commandArgs() call.
#' Neither belongs in a module, so seeing either is enough to refuse -- loudly,
#' naming the file and where it should go, rather than sourcing it and finding
#' out what it does.
tsf_script_not_module <- function(path) {
  head <- tryCatch(readLines(path, n = 40, warn = FALSE),
                   error = function(e) character(0))
  if (!length(head)) return(NA_character_)
  if (grepl("^#!", head[1])) return("it starts with a shebang")
  body <- head[!grepl("^\\s*#", head)]
  if (any(grepl("commandArgs\\s*\\(", body))) {
    return("it reads commandArgs() at the top level")
  }
  NA_character_
}

tsf_module_order <- function(dir = "R") {
  # stages.R calls into contrast.R and differential.R, and alphabetical order
  # would put "stages" before them only by luck. Named explicitly so a new
  # module cannot silently break the load order.
  first <- c("utils_io", "config", "labels", "paths", "grid", "period_floor",
             "fingerprint", "contrast", "differential")
  paths <- list.files(dir, pattern = "\\.R$", full.names = TRUE)

  # Refuse before sourcing, not after. A script in here is a bug in the tree
  # layout, and the cost of guessing wrong is whatever that script does.
  for (p in paths) {
    why <- tsf_script_not_module(p)
    if (!is.na(why)) {
      stop(p, " looks like a script, not a module: ", why, ".\n",
           "  Everything in ", dir, "/ is sourced on every run, so a script ",
           "here executes before any argument is parsed.\n",
           "  Move it to scripts/ and call it explicitly.",
           call. = FALSE)
    }
  }

  files <- sort(sub("\\.R$", "", basename(paths)))
  ordered <- c(intersect(first, files), setdiff(files, first))
  file.path(dir, paste0(ordered, ".R"))
}

# Worker-count policy, in one place. It lived in stages.R, but it is a utility
# rather than a stage: three chromosome-parallel stages, consensus and compare
# all need it, and having it inside stages.R meant nothing could ask about the
# policy without loading every stage function in the pipeline.
#
# Resolution: --cores, then N_WORKERS, then automatic. `n_tasks` is the real
# ceiling -- a worker with no task is not a faster run -- and `default_max`
# caps only the automatic path.
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

#' Order-independent fingerprint of a set of frequency keys.
#'
#' Used to decide whether a cached permutation null belongs to the family being
#' tested now. It needs to be stable and to change when the family changes; it
#' does not need to be cryptographic, since it is only ever compared against a
#' cache this same installation wrote.
#'
#' No external dependency: the pipeline runs on base R, and adding `digest` for
#' one cache key would be the first package in the tree.
digest_keys <- function(keys) {
  keys <- sort(unique(as.character(keys)))
  if (!length(keys)) return("empty")
  # Position-weighted sum of character codes, taken modulo a large prime so the
  # value stays in double range for any family size.
  acc <- 0
  for (i in seq_along(keys)) {
    codes <- utf8ToInt(keys[i])
    acc <- (acc * 131 + sum(codes * seq_along(codes)) + i) %% 2147483647
  }
  sprintf("n%d-h%010.0f", length(keys), acc)
}

tsf_load_all <- function(dir = "R") {
  invisible(lapply(tsf_module_order(dir), source))
}

TSF_VERSION <- "TissueSpectF 0.1.0"
