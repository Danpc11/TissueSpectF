#!/usr/bin/env Rscript
# Dependency-free tests for the label layer. Run: Rscript tests/test_labels.R
source("R/utils_io.R"); source("R/config.R"); source("R/labels.R")
source("R/grid.R"); source("R/ingest.R"); source("R/annotation.R")
source("R/fingerprint.R")
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
check("non-disease cohort -> Control_disease_cohort",
      identical(lab135$condition[1], "Control_disease_cohort"))
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

check("title token N -> Normal_histology, not F0 and not Control",
      identical(lab162$condition[1], "Normal_histology"))
check("no Control label is ever produced", !any(lab162$condition %in% "Control"))
check("the three healthy states are distinct classes under one state", {
  # Merging them would assume what is worth testing: that a healthy liver looks
  # the same whichever study recruited it.
  ids <- vapply(c("Control_disease_cohort", "Control_external_study",
                  "Normal_histology"),
                function(c) tsf_class_id("liver", c, cfg_135$vocabulary_spec),
                character(1))
  all(startsWith(ids, "liver::healthy::")) && length(unique(ids)) == 3L })
check("a stage becomes a disease class", {
  identical(tsf_class_id("liver", "F2", cfg_135$vocabulary_spec),
            "liver::disease::NAFLD_fibrosis_F2") })
check("the progression excludes the healthy states", {
  identical(tsf_progression(cfg_135$vocabulary_spec), paste0("F", 0:4)) })
check("title token F2", identical(lab162$condition[2], "F2"))
check("fallback to fibrosis column when the title has no token",
      identical(lab162$condition[4], "F1") && identical(lab162$label_rule[4], "biopsy_fibrosis_stage"))
check("normal liver histology parses as stage 0 only when a config asks it to",
      identical(parse_fibrosis_stage("fibrosis stage: normal liver histology",
                                     c("normal liver histology")), 0) &&
        is.na(parse_fibrosis_stage("fibrosis stage: normal liver histology")))
check("the phenotype field and the title token agree on the healthy class", {
  # Both routes must reach Normal_histology; before this change one said F0 and
  # the other said stage 0, which is how the two groups got merged.
  identical(lab162$condition[1], "Normal_histology") &&
    identical(lab162$label_rule[1], "normal_histology_field") })

# --- guardrails --------------------------------------------------------------
check("config with has_control_cohort=FALSE rejects a Control rule",
      inherits(tryCatch({
        bad <- cfg_162
        bad$condition_rules <- list(list(id = "x", type = "column_match",
                                         column = "title", values = c("Liver N 01"),
                                         assign = "Control_disease_cohort"))
        harmonize_conditions(pheno_162, bad); FALSE
      }, error = function(e) e), "error"))

check("ambiguous title with two stage tokens is unresolved by the title rule",
      is.na(parse_title_stage("sample F1 vs F3")))

check("stage out of range is NA", is.na(parse_fibrosis_stage("7")))

check("comparable_conditions is the ordered intersection",
      identical(comparable_conditions(list(a = c("F0","F1","Control"), b = c("F1","F0"))),
                c("F0", "F1")))

check("a compound rule needs every sub-condition to hold", {
  ph <- data.frame(
    geo_accession = paste0("GSM", 1:5),
    `steatosis grade` = c("steatosis grade: 0", "steatosis grade: 0",
                          "steatosis grade: 1", "steatosis grade: 0",
                          "steatosis grade: 0"),
    `cytological ballooning grade` = paste0("cytological ballooning grade: ",
                                            c(0, 1, 0, 0, 0)),
    `fibrosis stage` = paste0("fibrosis stage: ", c(0, 0, 0, 2, 0)),
    check.names = FALSE, stringsAsFactors = FALSE)
  rule <- list(id = "nh", type = "compound", assign = "Normal_histology",
               all_of = list(list(column = "steatosis grade", values = "0"),
                             list(column = "cytological ballooning grade", values = "0"),
                             list(column = "fibrosis stage", values = "0")))
  got <- apply_condition_rule(rule, ph)
  # Only samples 1 and 5 satisfy all three: 2 has ballooning, 3 has steatosis,
  # 4 has fibrosis. A liver with fibrosis is not histologically normal whatever
  # its steatosis score.
  identical(got, c("Normal_histology", NA, NA, NA, "Normal_histology")) })

check("a compound rule with a missing column does not fire", {
  ph <- data.frame(geo_accession = "GSM1", `steatosis grade` = "0",
                   check.names = FALSE, stringsAsFactors = FALSE)
  rule <- list(id = "nh", type = "compound", assign = "Normal_histology",
               all_of = list(list(column = "steatosis grade", values = "0"),
                             list(column = "absent column", values = "0")))
  all(is.na(apply_condition_rule(rule, ph))) })

check("the compound rule is ordered before the stage rule", {
  cfg <- load_dataset_config("GSE130970")
  ids <- vapply(cfg$condition_rules, function(r) r$id, character(1))
  identical(ids[1], "normal_histology_scores") &&
    identical(ids[2], "biopsy_fibrosis_stage") })

# --- sample filters and per-dataset class restriction ------------------------
check("a sample filter keeps only what matches", {
  ph <- data.frame(geo_accession = paste0("G", 1:4),
                   etiology = c("etiology: MASLD", "etiology: HBV",
                                "etiology: MASH", "etiology: alcohol"),
                   check.names = FALSE, stringsAsFactors = FALSE)
  lab <- data.frame(condition = "F4", keep = rep(TRUE, 4), stringsAsFactors = FALSE)
  cfg <- list(id = "X", sample_filter = list(
    list(column = "etiolog", pattern = "masld|mash")))
  identical(apply_sample_filter(lab, ph, cfg)$keep, c(TRUE, FALSE, TRUE, FALSE)) })

check("an unlabelled sample is dropped, not admitted", {
  ph <- data.frame(geo_accession = c("G1", "G2"),
                   etiology = c("etiology: MASLD", "etiology: "),
                   check.names = FALSE, stringsAsFactors = FALSE)
  lab <- data.frame(condition = "F4", keep = c(TRUE, TRUE), stringsAsFactors = FALSE)
  cfg <- list(id = "X", sample_filter = list(
    list(column = "etiolog", pattern = "masld")))
  identical(apply_sample_filter(lab, ph, cfg)$keep, c(TRUE, FALSE)) })

check("a missing filter column aborts instead of admitting everything", {
  ph <- data.frame(geo_accession = "G1", other = "x", stringsAsFactors = FALSE)
  lab <- data.frame(condition = "F4", keep = TRUE, stringsAsFactors = FALSE)
  cfg <- list(id = "X", sample_filter = list(list(column = "etiolog", pattern = "masld")))
  inherits(tryCatch(apply_sample_filter(lab, ph, cfg), error = function(e) e), "error") })

check("keep_conditions restricts what a dataset contributes", {
  ph <- data.frame(geo_accession = paste0("G", 1:4), stringsAsFactors = FALSE)
  lab <- data.frame(condition = c("F0", "F3", "F4", "F1"), keep = rep(TRUE, 4),
                    stringsAsFactors = FALSE)
  cfg <- list(id = "X", keep_conditions = c("F3", "F4"))
  identical(apply_sample_filter(lab, ph, cfg)$keep, c(FALSE, TRUE, TRUE, FALSE)) })

check("the two new cohort configs load and declare their restrictions", {
  a <- load_dataset_config("GSE276114"); b <- load_dataset_config("GSE142530")
  identical(a$keep_conditions, c("F3", "F4")) && isFALSE(a$has_control_cohort) &&
    identical(b$keep_conditions, "Control_external_study") &&
    isTRUE(b$has_control_cohort) &&
    # GSE276114 needs two filters: etiology, and dropping the F0-2 bin, which
    # spans three classes and resolves to none of them.
    length(a$sample_filter) == 2 && length(b$sample_filter) == 1 })

check("an etiology filter is not fooled by a similarly named column", {
  # "^disease" would match "disease group:ch1" before "disease:ch1" and filter
  # on the stage bin while believing it filtered on etiology. The anchor with
  # the colon is what prevents it.
  ph <- data.frame(geo_accession = c("G1", "G2"),
                   `disease group:ch1` = c("F4", "F4"),
                   `disease:ch1` = c("MASLD", "CVH"),
                   check.names = FALSE, stringsAsFactors = FALSE)
  cfg <- list(id = "X", sample_filter = list(
    list(column = "^disease:", values = "MASLD")))
  lab <- data.frame(condition = "F4", keep = c(TRUE, TRUE), stringsAsFactors = FALSE)
  identical(apply_sample_filter(lab, ph, cfg)$keep, c(TRUE, FALSE)) })

check("a declared substitution reconciles labels and leaves a trail", {
  # The count file says "Control_Lille 389"; the series matrix says
  # "Control_389". The substitution is declared, and the raw label survives in
  # the map so the correspondence can be checked rather than trusted.
  f <- tempfile(fileext = ".csv")
  writeLines(c("gene_id,symbol,RB_N1,RB_N2,RB_N3",
               "NA,NA,Control_Lille 389,Control_TPF 111484,not.used",
               "ENSG00000000003.1,TSPAN6,10,20,5"), f)
  d <- read_counts(f, "ENSEMBL",
                   list(sep = ",", id_column = 1L, symbol_column = 2L,
                        sample_map_row = 1L,
                        map_transform = list(list(pattern = "_(Lille|TPF)[[:space:]]+",
                                                  replacement = "_")),
                        exclude_columns = c("not.used")))
  m <- attr(d, "sample_map")
  setequal(m$mapped_id, c("Control_389", "Control_111484")) &&
    "Control_Lille 389" %in% m$raw_label &&
    !any(grepl("not.used", colnames(d))) })

check("the Ensembl version suffix is stripped", {
  f <- tempfile(fileext = ".csv")
  writeLines(c("gene_id,symbol,S1", "ENSG00000000003.14,TSPAN6,10"), f)
  identical(read_counts(f, "ENSEMBL",
                        list(sep = ",", id_column = 1L, symbol_column = 2L))$source_id,
            "ENSG00000000003") })

check("cohort_roles keep a second control class out of the disease group", {
  # Comparing against the baseline alone called Control_external_study disease,
  # because it is not the baseline. That column feeds reports and balancing.
  v <- load_vocabulary("liver_fibrosis")
  got <- tsf_cohort_role(c("Control_disease_cohort", "Control_external_study",
                           "Normal_histology", "F2"), v)
  identical(got, c("control", "control", "within_disease_normal", "disease")) })

check("a vocabulary without cohort_roles falls back to the baseline", {
  v <- list(baseline = "Control", cohort_roles = NULL)
  identical(tsf_cohort_role(c("Control", "Case"), v), c("control", "disease")) })

check("an undeclared class is recorded as unknown, not guessed", {
  v <- load_vocabulary("liver_fibrosis")
  identical(tsf_cohort_role("Made_up", v), "unknown") })

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

# --- GTF / GENCODE annotation ------------------------------------------------
mini_gtf <- function() {
  f <- tempfile(fileext = ".gtf")
  a <- function(g, n, ty) sprintf('gene_id "%s"; gene_type "%s"; gene_name "%s";', g, ty, n)
  writeLines(c(
    "##description: synthetic",
    paste("chr1\tH\tgene\t1000\t5000\t.\t+\t.", a("ENSG1.5", "AAA", "protein_coding"), sep = "\t"),
    paste("chr1\tH\texon\t1000\t1500\t.\t+\t.", a("ENSG1.5", "AAA", "protein_coding"), sep = "\t"),
    paste("chr1\tH\texon\t1400\t2000\t.\t+\t.", a("ENSG1.5", "AAA", "protein_coding"), sep = "\t"),
    paste("chr1\tH\texon\t4000\t4200\t.\t+\t.", a("ENSG1.5", "AAA", "protein_coding"), sep = "\t"),
    paste("chr1\tH\tgene\t8000\t9000\t.\t+\t.", a("ENSG2.1", "BBB", "lncRNA"), sep = "\t"),
    paste("chr1\tH\texon\t8000\t9000\t.\t+\t.", a("ENSG2.1", "BBB", "lncRNA"), sep = "\t"),
    paste("chrY\tH\tgene\t100\t200\t.\t+\t.", a("ENSG1.5_PAR_Y", "AAA", "protein_coding"), sep = "\t"),
    paste("chrY\tH\texon\t100\t200\t.\t+\t.", a("ENSG1.5_PAR_Y", "AAA", "protein_coding"), sep = "\t"),
    paste("GL000009.2\tH\tgene\t1\t100\t.\t+\t.", a("ENSG4.1", "SCAF", "protein_coding"), sep = "\t"),
    paste("GL000009.2\tH\texon\t1\t100\t.\t+\t.", a("ENSG4.1", "SCAF", "protein_coding"), sep = "\t")), f)
  f
}
CHR_LEVELS <- c(as.character(1:22), "X", "Y")

check("a GTF is read into the common annotation shape", {
  a <- read_gtf_annotation(mini_gtf(), CHR_LEVELS)
  all(c("gene_id", "gene_name", "chr", "start", "end", "strand",
        "gene_length", "gene_type") %in% colnames(a)) &&
    identical(sort(a$gene_id), c("ENSG1", "ENSG2")) })

check("chromosome names and gene versions are normalised", {
  a <- read_gtf_annotation(mini_gtf(), CHR_LEVELS)
  all(a$chr == "1") && !any(grepl("\\.", a$gene_id)) })

check("exonic length is the union of exons, not their sum", {
  # 1000-1500 and 1400-2000 overlap: 1001 together, not 1102. Plus 4000-4200.
  a <- read_gtf_annotation(mini_gtf(), CHR_LEVELS, length_mode = "exonic")
  a$gene_length[a$gene_name == "AAA"] == 1202 })

check("the genomic span is not a transcript length", {
  a <- read_gtf_annotation(mini_gtf(), CHR_LEVELS, length_mode = "span")
  a$gene_length[a$gene_name == "AAA"] == 4001 })

check("PAR_Y features are dropped from every feature type", {
  # Stripping the version makes ENSG1.5_PAR_Y collide with its X copy, so a
  # PAR_Y exon left behind is summed into the X gene's length -- silently, and
  # only in the length, which nothing downstream would flag.
  a <- read_gtf_annotation(mini_gtf(), CHR_LEVELS)
  nrow(a) == 2 && a$gene_length[a$gene_name == "AAA"] == 1202 })

check("scaffolds are dropped", {
  a <- read_gtf_annotation(mini_gtf(), CHR_LEVELS)
  !("SCAF" %in% a$gene_name) })

check("an identifier map adds Entrez ids a GTF does not carry", {
  # GENCODE has Ensembl ids and symbols, never Entrez, and three cohorts publish
  # counts keyed on Entrez. Without the map the join is empty and ingest stops.
  m <- tempfile(fileext = ".tsv")
  writeLines(c("GeneID\tSymbol\tEnsemblGeneID",
               "111\tAAA\tENSG1.9", "222\tBBB\tENSG2.3"), m)
  p <- list(geo_dir = dirname(mini_gtf()), annotation_format = "gtf",
            chrom_levels = CHR_LEVELS, gene_universe = "^protein_coding$",
            id_map = list(file = basename(m), ensembl_column = "EnsemblGeneID",
                          entrez_column = "GeneID"))
  p$annotation_file <- basename(mini_gtf())
  p$geo_dir <- dirname(m)
  file.copy(mini_gtf(), file.path(dirname(m), p$annotation_file), overwrite = TRUE)
  a <- read_annotation(p)
  identical(sort(a$entrez_id[!is.na(a$entrez_id)]), c("111", "222")) })

check("a map that matches nothing aborts instead of ingesting an empty join", {
  m <- tempfile(fileext = ".tsv")
  writeLines(c("GeneID\tEnsemblGeneID", "111\tSOMETHING_ELSE"), m)
  p <- list(geo_dir = dirname(m), annotation_format = "gtf",
            chrom_levels = CHR_LEVELS, gene_universe = "^protein_coding$",
            id_map = list(file = basename(m), ensembl_column = "EnsemblGeneID",
                          entrez_column = "GeneID"))
  p$annotation_file <- basename(mini_gtf())
  file.copy(mini_gtf(), file.path(dirname(m), p$annotation_file), overwrite = TRUE)
  inherits(tryCatch(read_annotation(p), error = function(e) e), "error") })

check("each format gets its own protein-coding pattern", {
  identical(default_gene_universe("ncbi"), "^protein-coding$") &&
    identical(default_gene_universe("gtf"), "^protein_coding$") })

check("a universe pattern from the wrong format aborts with the right one", {
  # GENCODE writes protein_coding, NCBI writes protein-coding. A pattern for one
  # matches nothing in the other, and "nothing" would build an empty grid.
  p <- list(geo_dir = dirname(mini_gtf()), annotation_format = "gtf",
            chrom_levels = CHR_LEVELS, gene_universe = "^protein-coding$")
  p$annotation_file <- basename(mini_gtf())
  e <- tryCatch(read_annotation(p), error = function(e) conditionMessage(e))
  is.character(e) && grepl("\\^protein_coding\\$", e) })

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
