# paths.R -- one place that knows where every stage writes.
#
# Stages communicate only through these files, never through session state.
# That is what makes any stage re-runnable on its own: the old
# regenerar_picos_genes_condicion.R existed because the gene reconstruction
# could only be rebuilt from objects living in a finished R session.

tsf_paths <- function(project, dataset_id) {
  base <- file.path(project$results_dir, dataset_id)
  list(
    base       = base,
    spectra    = file.path(base, "spectra"),
    maxt       = file.path(base, "maxt"),
    stability  = file.path(base, "stability"),
    peaks      = file.path(base, "peaks"),
    peak_genes = file.path(base, "peaks_genes"),
    compare    = file.path(project$results_dir, "comparison")
  )
}

p_spectra_samples   <- function(p, cond) file.path(p$spectra, sprintf("spectra_samples_%s.tsv", cond))
p_spectra_condition <- function(p, cond) file.path(p$spectra, sprintf("spectra_condition_%s.tsv", cond))
p_maxt              <- function(p, cond) file.path(p$maxt, sprintf("maxt_individual_%s.tsv", cond))
p_stability         <- function(p, cond) file.path(p$stability, sprintf("maxt_stability_%s.tsv", cond))
p_peaks             <- function(p, branch, cond) file.path(p$peaks, branch, sprintf("peaks_%s.tsv", cond))
p_peak_genes_dir    <- function(p, branch, cond) file.path(p$peak_genes, branch, cond)

#' Load the ingested dataset plus whatever a stage needs from earlier stages.
tsf_stage_inputs <- function(project, dataset_id, need = character(0)) {
  ds <- load_dataset(dataset_id, project)
  p <- tsf_paths(project, dataset_id)
  conds <- levels(droplevels(ds$samples$condition))
  out <- list(dataset = ds, paths = p, conditions = conds,
              chrom_idx = chromosome_index(ds$genes, project$chrom_levels,
                                           project$min_genes_per_chr))
  if ("spectra" %in% need) {
    out$spectra <- stats::setNames(lapply(conds, function(c)
      read_tsv_tsf(p_spectra_condition(p, c), required = FALSE)), conds)
  }
  if ("maxt" %in% need) {
    out$maxt <- stats::setNames(lapply(conds, function(c)
      read_tsv_tsf(p_maxt(p, c), required = FALSE)), conds)
  }
  if ("stability" %in% need) {
    out$stability <- stats::setNames(lapply(conds, function(c)
      read_tsv_tsf(p_stability(p, c), required = FALSE)), conds)
  }
  out
}

#' Shared CLI parsing: dataset ids, --cond, --branch, --force.
tsf_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  flag <- function(name) {
    hit <- grep(paste0("^--", name, "="), args, value = TRUE)
    if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else NULL
  }
  list(
    datasets = args[!grepl("^--", args)],
    cond     = flag("cond"),
    branch   = flag("branch"),
    force    = any(args == "--force")
  )
}

#' Datasets to process: those named on the command line, or all configured.
tsf_dataset_ids <- function(cli) {
  if (length(cli$datasets)) cli$datasets
  else sub("\\.R$", "", list.files("config/datasets", pattern = "\\.R$"))
}

#' Conditions to process, honouring --cond=F2.
tsf_conditions <- function(available, cli) {
  if (is.null(cli$cond)) return(available)
  want <- strsplit(cli$cond, ",")[[1]]
  missing <- setdiff(want, available)
  if (length(missing)) tsf_abort("Condition not present in this dataset: ",
                                 paste(missing, collapse = ", "))
  want
}

tsf_source_pipeline <- function() {
  for (f in c("utils_io", "config", "labels", "fetch", "ingest", "paths", "grid",
              "spectrum", "maxt", "condition_test", "stability", "peaks_genes", "compare")) {
    source(file.path("R", paste0(f, ".R")))
  }
}
