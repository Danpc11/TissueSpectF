#!/usr/bin/env Rscript
# Dependency-free tests for the label layer. Run: Rscript tests/test_labels.R
source("R/utils_io.R"); source("R/config.R"); source("R/labels.R")
source("R/grid.R"); source("R/ingest.R"); source("R/fingerprint.R")
source("R/reference.R"); source("R/bundle.R")

failures <- 0L
check <- function(label, expr) {
  ok <- isTRUE(tryCatch(expr, error = function(e) { cat("   error: ", conditionMessage(e), "\n"); FALSE }))
  cat(if (ok) "  PASS  " else "  FAIL  ", label, "\n")
  if (!ok) failures <<- failures + 1L
}

cfg_135 <- load_dataset_config("GSE135251")
cfg_162 <- load_dataset_config("GSE162694")

# --- GSE135251: Control cohort and F0 must stay separate --------------------
pheno_135 <- data.frame(
  geo_accession = c('"GSM1"', "GSM2", "GSM3", "GSM4"),
  title = c("control 1", "NAFLD 2", "NAFLD 3", "NAFLD 4"),
  `disease state` = c("disease: Control", "disease: NAFLD", "disease: NAFLD", "disease: NAFLD"),
  `fibrosis stage` = c("fibrosis stage: 0", "fibrosis stage: 0",
                       "fibrosis stage: 3", "fibrosis stage: not available"),
  check.names = FALSE, stringsAsFactors = FALSE)

lab135 <- harmonize_conditions(pheno_135, cfg_135)

check("quoted accessions are cleaned", identical(lab135$sample_id[1], "GSM1"))
check("non-disease cohort -> Control", identical(lab135$condition[1], "Control"))
check("NAFLD with stage 0 -> F0, not Control", identical(lab135$condition[2], "F0"))
check("stage 3 -> F3", identical(lab135$condition[3], "F3"))
check("unparseable stage -> NA, never 'FNA'",
      is.na(lab135$condition[4]) && !any(grepl("FNA", lab135$condition)))
check("keep flag marks the unresolved sample", identical(lab135$keep, c(TRUE, TRUE, TRUE, FALSE)))
check("cohort column", identical(lab135$cohort[1:2], c("control", "disease")))

# --- GSE162694: no control cohort; normal histology is F0 -------------------
pheno_162 <- data.frame(
  geo_accession = c("GSM10", "GSM11", "GSM12", "GSM13"),
  title = c("Liver N 01", "Liver F2 02", "Liver F4 03", "Liver 04"),
  `fibrosis stage` = c("fibrosis stage: normal liver histology", "fibrosis stage: 2",
                       "fibrosis stage: 4", "fibrosis stage: 1"),
  check.names = FALSE, stringsAsFactors = FALSE)

lab162 <- harmonize_conditions(pheno_162, cfg_162)

check("title token N -> F0 (not Control)", identical(lab162$condition[1], "F0"))
check("no Control label is ever produced", !any(lab162$condition %in% "Control"))
check("title token F2", identical(lab162$condition[2], "F2"))
check("fallback to fibrosis column when the title has no token",
      identical(lab162$condition[4], "F1") && identical(lab162$label_rule[4], "biopsy_fibrosis_stage"))
check("normal liver histology parses as stage 0",
      identical(parse_fibrosis_stage("fibrosis stage: normal liver histology",
                                     c("normal liver histology")), 0))
check("the same sample gets the same label from either route",
      identical(lab162$condition[1], "F0") &&
        identical(paste0("F", parse_fibrosis_stage(pheno_162$`fibrosis stage`[1],
                                                   c("normal liver histology"))), "F0"))

# --- guardrails --------------------------------------------------------------
check("config with has_control_cohort=FALSE rejects a Control rule",
      inherits(tryCatch({
        bad <- cfg_162
        bad$condition_rules <- list(list(id = "x", type = "column_match",
                                         column = "title", values = c("Liver N 01"),
                                         assign = "Control"))
        harmonize_conditions(pheno_162, bad); FALSE
      }, error = function(e) e), "error"))

check("ambiguous title with two stage tokens is unresolved by the title rule",
      is.na(parse_title_stage("sample F1 vs F3")))

check("stage out of range is NA", is.na(parse_fibrosis_stage("7")))

check("comparable_conditions is the ordered intersection",
      identical(comparable_conditions(list(a = c("F0","F1","Control"), b = c("F1","F0"))),
                c("F0", "F1")))

# --- input schemas and vocabularies without a baseline -----------------------
tmp_v <- file.path(tempdir(), "voc"); tmp_d <- file.path(tempdir(), "ds")
dir.create(tmp_v, showWarnings = FALSE); dir.create(tmp_d, showWarnings = FALSE)
invisible(file.copy(list.files("config/vocabularies", full.names = TRUE), tmp_v, overwrite = TRUE))
local_cfg <- function(txt, name) {
  writeLines(txt, file.path(tmp_d, paste0(name, ".R")))
  old <- load_vocabulary
  assign("load_vocabulary",
         function(n, dir = tmp_v) source(file.path(dir, paste0(n, ".R")), local = TRUE)$value,
         envir = globalenv())
  on.exit(assign("load_vocabulary", old, envir = globalenv()), add = TRUE)
  tryCatch(load_dataset_config(name, tmp_d), error = function(e) e)
}

check("a tissue atlas needs no has_control_cohort", {
  cfg <- local_cfg('list(id="ATLAS", tissue="multi", source="matrix",
    vocabulary="tissue_atlas", condition_levels=c("liver","lung"),
    counts_file="c.tsv", metadata_file="m.tsv")', "ATLAS")
  !inherits(cfg, "error") && identical(cfg$has_control_cohort, FALSE) })

check("source = matrix does not require a series matrix", {
  cfg <- local_cfg('list(id="MAT", tissue="lung", source="matrix",
    vocabulary="case_control", has_control_cohort=TRUE,
    counts_file="c.tsv", metadata_file="m.tsv")', "MAT")
  !inherits(cfg, "error") && identical(cfg$source, "matrix") })

check("a vocabulary with a baseline still demands has_control_cohort", {
  cfg <- local_cfg('list(id="LIV", tissue="liver", vocabulary="liver_fibrosis",
    counts_file="c.tsv", series_matrix="s.txt",
    condition_rules=list(list(id="r", type="fibrosis_stage", column="fib")))', "LIV")
  inherits(cfg, "error") })

check("an unknown source is refused", {
  cfg <- local_cfg('list(id="ODD", tissue="x", source="salmon",
    vocabulary="case_control", has_control_cohort=TRUE)', "ODD")
  inherits(cfg, "error") })

check("incompatible grids cannot share a reference", {
  g <- data.frame(gene_id = paste0("g", 1:10), chr = "1",
                  grid_index = 1:10, grid_N = 10L, stringsAsFactors = FALSE)
  a <- list(grid = g, provenance = list(species = "Homo sapiens",
            genome_build = "GRCh38", annotation_release = "NCBI",
            gene_universe = "^protein-coding$", grid_digest = grid_digest(g)))
  b <- a; b$provenance$genome_build <- "GRCh37"
  inherits(tryCatch(assert_compatible_grids(list(A = a, B = b)),
                    error = function(e) e), "error") })

check("identical grids are accepted", {
  g <- data.frame(gene_id = paste0("g", 1:10), chr = "1",
                  grid_index = 1:10, grid_N = 10L, stringsAsFactors = FALSE)
  a <- list(grid = g, provenance = list(species = "Homo sapiens",
            genome_build = "GRCh38", annotation_release = "NCBI",
            gene_universe = "^protein-coding$", grid_digest = grid_digest(g)))
  !inherits(tryCatch(assert_compatible_grids(list(A = a, B = a)),
                     error = function(e) e), "error") })

check("a changed grid changes its digest", {
  g <- data.frame(gene_id = paste0("g", 1:10), chr = "1",
                  grid_index = 1:10, grid_N = 10L, stringsAsFactors = FALSE)
  h <- g; h$grid_index[3] <- 99L
  grid_digest(g) != grid_digest(h) })

check("chromosome names are normalised", {
  identical(normalise_chrom_names(c("chr1", "1", "CHR1", "chrM", "MT")),
            c("1", "1", "1", "MT", "MT")) })

# --- integration: the CLI's own load path resolves every stage ---------------
check("every stage function exists after the CLI load path", {
  # Loading exactly as scripts/tsf.R does, then checking that each declared
  # stage resolves. A module added to R/ but forgotten in a hand-written source
  # list used to fail only at run time, in the middle of a long pipeline.
  env <- new.env()
  sys.source("R/utils_io.R", envir = env)
  invisible(lapply(get("tsf_module_order", envir = env)(),
                   function(f) sys.source(f, envir = env)))
  fns <- get("stage_functions", envir = env)
  stages <- get("stage_names", envir = env)
  all(stages %in% names(fns)) &&
    all(vapply(fns[stages], is.function, logical(1))) })

check("the consensus entry points are reachable from the load path", {
  env <- new.env()
  sys.source("R/utils_io.R", envir = env)
  invisible(lapply(get("tsf_module_order", envir = env)(),
                   function(f) sys.source(f, envir = env)))
  all(vapply(c("consensus_spectrum", "consensus_signature", "phase_locking",
               "prevalence_from_rank", "signature_features"),
             function(f) is.function(get0(f, envir = env)), logical(1))) })

check("no module is missing from the load order", {
  env <- new.env(); sys.source("R/utils_io.R", envir = env)
  files <- basename(get("tsf_module_order", envir = env)())
  setequal(files, list.files("R", pattern = "\\.R$")) })

# --- the distributable bundle ------------------------------------------------
check("the bundle carries every module the app sources", {
  # Whatever app.R sources must be in BUNDLE_MODULES, or the bundle breaks on
  # someone else's machine and nowhere else.
  src <- readLines("app/app.R", warn = FALSE)
  line <- grep('^for \\(f in c\\("utils_io"', src, value = TRUE)
  needed <- regmatches(line, gregexpr('"[a-z_]+"', line))[[1]]
  needed <- gsub('"', "", needed)
  all(needed %in% BUNDLE_MODULES) })

check("bundled modules do not reach outside the bundle", {
  # A module that sources ../R/ or reads config/project.R would work in the
  # repository and fail in the bundle.
  bad <- unlist(lapply(BUNDLE_MODULES, function(m) {
    txt <- readLines(file.path("R", paste0(m, ".R")), warn = FALSE)
    grep('config/project\\.R|\\.\\./R|load_project_config|tsf_stage_inputs',
         txt, value = TRUE)
  }))
  length(bad) == 0 })

check("the launcher and provenance are generated", {
  ref <- list(built = "x", version = "y", target = "condition",
              model = list(classes = c("A", "B")),
              labels = data.frame(a = 1:3),
              grid = data.frame(gene_id = 1:10),
              params = list(k_max = 64L, features = "amplitude"),
              validation = NULL)
  length(bundle_launcher_sh()) > 5 && length(bundle_launcher_bat()) > 5 &&
    any(grepl("NOT validated", bundle_readme(ref))) &&
    any(grepl("grid digest", bundle_manifest(ref, "reference.rds"))) })

cat("\n", if (failures == 0L) "All tests passed." else paste(failures, "test(s) failed."), "\n")
quit(status = if (failures == 0L) 0L else 1L)
