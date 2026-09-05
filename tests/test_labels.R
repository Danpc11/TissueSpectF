#!/usr/bin/env Rscript
# Dependency-free tests for the label layer. Run: Rscript tests/test_labels.R
source("R/utils_io.R"); source("R/config.R"); source("R/labels.R")
source("R/grid.R"); source("R/ingest.R"); source("R/annotation.R")
source("R/fingerprint.R")
source("R/reference.R"); source("R/bundle.R")
source("R/maxt.R")   # maxt_cores(); local_workers() comes from utils_io.R

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
check("the three healthy groups map to one Controles class", {
  # This asserted the opposite until the source publications settled it: all
  # three are liver tissue without disease. Two of them also existed in a
  # single cohort each, so leave-one-cohort-out could never learn them --
  # GSE142530 was skipped as a hold-out and both classes won zero predictions.
  ids <- vapply(c("Control_disease_cohort", "Control_external_study",
                  "Normal_histology"),
                function(c) tsf_class_id("liver", c, cfg_135$vocabulary_spec),
                character(1))
  all(ids == "liver::healthy::Controles") })

check("F0 is not folded into Controles", {
  # F0 is a histological stage in a NAFLD patient. Merging it with the controls
  # would assume what the Controles vs F0 contrast exists to test.
  f0 <- tsf_class_id("liver", "F0", cfg_135$vocabulary_spec)
  f0 != tsf_class_id("liver", "Normal_histology", cfg_135$vocabulary_spec) &&
    startsWith(f0, "liver::disease::") })

check("the raw labels survive the merge", {
  # The merge is in `conditions` only. Ingest, the audit trail and cohort_roles
  # must still be able to say which cohort a control came from.
  v <- cfg_135$vocabulary_spec
  all(c("Control_disease_cohort", "Control_external_study",
        "Normal_histology") %in% names(v$cohort_roles)) &&
    identical(unname(v$cohort_roles[["Normal_histology"]]),
              "within_disease_normal") })

check("the baseline names a class that exists after the merge", {
  v <- cfg_135$vocabulary_spec
  v$baseline %in% unname(v$conditions) })
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
check("config with has_control_cohort=FALSE rejects a Control rule", {
  # The flag is set explicitly here rather than borrowed from a real dataset
  # config: this test used cfg_162 as its FALSE example, and broke the moment
  # that dataset legitimately became a control contributor. A guardrail test
  # should fail when the guardrail breaks, not when a config changes.
  bad <- cfg_162
  bad$has_control_cohort <- FALSE
  bad$condition_rules <- list(list(id = "x", type = "column_match",
                                   column = "title", values = c("Liver N 01"),
                                   assign = "Control_disease_cohort"))
  inherits(tryCatch({ harmonize_conditions(pheno_162, bad); FALSE },
                    error = function(e) e), "error") })

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

# --- reproducibility of the shipped config -----------------------------------
#
# config/project.R used to default to one author's scratch directory, so a fresh
# clone failed at `./tsf check` with paths belonging to a machine the user has
# never seen. These assertions keep that from coming back: a default path is a
# choice about layout, never a location.

# --- R/ holds modules, not scripts -------------------------------------------
#
# tsf_load_all() sources every .R file in R/ on every invocation. A module only
# defines things, so that is safe. A script DOES things, and a copy of
# scripts/clean_results.R placed in R/ deleted the results tree on every run,
# before argument parsing -- which is why the deletion appeared in the log
# above an unrelated usage error. The boundary is now enforced at load time and
# asserted here.

check("no file in R/ is a script", {
  bad <- vapply(list.files("R", pattern = "[.]R$", full.names = TRUE),
                function(p) !is.na(tsf_script_not_module(p)), logical(1))
  !any(bad) })

# The probe file is removed inside the expression, not with on.exit(): on.exit
# registers on check()'s frame and does not fire where it is needed, and a
# probe left behind in R/ makes every later assertion fail for the wrong
# reason -- which is exactly what happened the first time this was written.
probe_refused <- function(name, lines) {
  p <- file.path("R", name)
  writeLines(lines, p)
  refused <- inherits(tryCatch(tsf_module_order("R"), error = function(e) e),
                      "error")
  unlink(p)
  refused && !file.exists(p) }

check("a shebang in R/ is refused before it is sourced", {
  probe_refused("zz_probe_shebang.R",
                c("#!/usr/bin/env Rscript", "cat(\"side effect\")")) })

check("a top-level commandArgs() in R/ is refused", {
  probe_refused("zz_probe_args.R",
                c("a <- commandArgs(trailingOnly = TRUE)", "unlink(\"x\")")) })

check("the guard does not reject a real module", {
  # utils_io.R is the most module-like file there is; if the guard trips on it
  # the pattern is too broad and the pipeline will not load at all.
  is.na(tsf_script_not_module("R/utils_io.R")) &&
    is.na(tsf_script_not_module("R/grid.R")) })

check("clean_results.R lives in scripts/, not R/", {
  file.exists("scripts/clean_results.R") &&
    !file.exists("R/clean_results.R") })

check("every option the CLI accepts is wired to something", {
  # --seed was in OPTION_ALIASES and set() by nothing: accepted without
  # complaint, then discarded. Every permutation and bootstrap derives from
  # maxt$seed, so a run passing --seed got the default 42 and identical
  # results -- a failure that looks like the seed not mattering.
  src <- readLines("scripts/tsf.R", warn = FALSE)
  alias_block <- src[seq(grep("^OPTION_ALIASES", src)[1],
                         grep("^\\)", src)[grep("^\\)", src) >
                                grep("^OPTION_ALIASES", src)[1]][1])]
  targets <- unique(gsub('.*"([a-z_0-9]+)".*', "\\1",
                         regmatches(alias_block,
                                    regexpr('= "[a-z_0-9]+"', alias_block))))
  body <- paste(src, collapse = " ")
  # Paths and scope options are consumed directly rather than through set().
  # Options a stage reads straight off `opt` rather than through set(): they
  # are still wired, just not into the project config. `grepl("opt$<name>")`
  # below catches them; the list is for the ones consumed before that point.
  consumed_directly <- c("geo_dir", "interim_dir", "results_dir", "config",
                         "from", "to", "cond", "branch", "force", "dry_run",
                         "log", "help", "query", "reference", "input_unit",
                         "out", "cores", "datasets")
  body <- paste(body, paste(readLines("R/stages.R", warn = FALSE),
                            collapse = " "))
  unwired <- Filter(function(t) {
    if (t %in% consumed_directly) return(FALSE)
    !grepl(paste0('"', t, '"'), body, fixed = TRUE) ||
      !grepl(paste0("opt$", t), body, fixed = TRUE)
  }, targets)
  length(unwired) == 0 })

check("--seed reaches maxt$seed", {
  out <- paste(suppressWarnings(system2("./tsf",
    c("consensus", "--seed", "99", "--config", "config/project.R",
      "--geo-dir", "d", "--interim-dir", "i", "--results-dir", "r",
      "--dry-run"), stdout = TRUE, stderr = TRUE)), collapse = " ")
  grepl("maxt$seed = 99", out, fixed = TRUE) })

check("path requirements are per command, not in bulk", {
  # An earlier version demanded all three paths from every command. That broke
  # `./tsf fetch --geo-dir data` -- fetch downloads and has no results tree --
  # and the CI step that dry-runs each stage. Demanding a path a command never
  # touches is a false requirement, and it teaches people to set variables to
  # satisfy a check rather than to say something true.
  fetch_ok <- suppressWarnings(system2("./tsf",
    c("fetch", "--geo-dir", "data", "--dry-run"),
    stdout = TRUE, stderr = TRUE))
  run_out <- suppressWarnings(system2("./tsf",
    c("run", "--to", "ingest", "--dry-run"),
    stdout = TRUE, stderr = TRUE))
  is.null(attr(fetch_ok, "status")) &&
    !is.null(attr(run_out, "status")) })

check("a missing-path error names every flag that is missing", {
  out <- paste(suppressWarnings(system2("./tsf",
    c("run", "--to", "ingest", "--dry-run"),
    stdout = TRUE, stderr = TRUE)), collapse = " ")
  all(vapply(c("--geo-dir", "--interim-dir", "--results-dir"),
             function(f) grepl(f, out, fixed = TRUE), logical(1))) })

check("CI supplies the paths its stage loop needs", {
  # The workflow lagged the policy change and started failing at the first
  # stage. That is what CI is for, but the fix belongs in the workflow.
  wf <- paste(readLines(".github/workflows/tests.yml", warn = FALSE),
              collapse = " ")
  grepl("TSF_GEO_DIR", wf, fixed = TRUE) &&
    grepl("TSF_INTERIM_DIR", wf, fixed = TRUE) &&
    grepl("TSF_RESULTS_DIR", wf, fixed = TRUE) })

check("config/project.R names no paths at all", {
  # Not "no absolute paths" but no paths: geo_dir, interim_dir and results_dir
  # have no default, so a run cannot inherit an output location from a file
  # nobody read. The invocation is the record.
  p <- load_project_config("config/project.R")
  is.null(p$geo_dir) && is.null(p$interim_dir) && is.null(p$results_dir) })

check("status asks for the trees it reads and not for the ones it does not", {
  # status reports what is on disk under interim and results. It does not read
  # the GEO downloads, so demanding --geo-dir from it would be the bulk check
  # this design replaced. This assertion is the per-command policy stated in
  # the one place it is easy to regress.
  out <- paste(suppressWarnings(system2("./tsf", "status",
                                        stdout = TRUE, stderr = TRUE)),
               collapse = " ")
  grepl("not set", out, fixed = TRUE) &&
    grepl("--interim-dir", out, fixed = TRUE) &&
    grepl("--results-dir", out, fixed = TRUE) &&
    !grepl("--geo-dir", out, fixed = TRUE) })

check("the environment supplies a path when no flag does", {
  old <- Sys.getenv("TSF_RESULTS_DIR", unset = NA)
  on.exit(if (is.na(old)) Sys.unsetenv("TSF_RESULTS_DIR") else
    Sys.setenv(TSF_RESULTS_DIR = old), add = TRUE)
  Sys.setenv(TSF_RESULTS_DIR = "run_probe")
  identical(load_project_config("config/project.R")$results_dir, "run_probe") })

check("an empty environment variable supplies nothing", {
  old <- Sys.getenv("TSF_RESULTS_DIR", unset = NA)
  on.exit(if (is.na(old)) Sys.unsetenv("TSF_RESULTS_DIR") else
    Sys.setenv(TSF_RESULTS_DIR = old), add = TRUE)
  Sys.setenv(TSF_RESULTS_DIR = "")
  is.null(load_project_config("config/project.R")$results_dir) })

check("the clean guard refuses the repository itself", {
  out <- suppressWarnings(system2("Rscript",
    c("scripts/clean_results.R", "--results-dir", "."),
    stdout = TRUE, stderr = TRUE))
  !is.null(attr(out, "status")) && attr(out, "status") != 0L })

check("clean will not delete a non-empty tree unconfirmed", {
  # results_dir points at the active run tree, so the default target of
  # `make clean` is real work. An unanswered prompt must leave it alone.
  tmp <- file.path(tempdir(), "tsf_clean_probe")
  dir.create(file.path(tmp, "sub"), recursive = TRUE, showWarnings = FALSE)
  writeLines("data", file.path(tmp, "sub", "keep.tsv"))
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  out <- suppressWarnings(system2("Rscript",
    c("scripts/clean_results.R", "--results-dir", tmp),
    stdin = "/dev/null", stdout = TRUE, stderr = TRUE))
  file.exists(file.path(tmp, "sub", "keep.tsv")) &&
    any(grepl("Nothing was deleted", paste(out, collapse = " "),
              fixed = TRUE)) })

check("clean with --force needs no prompt", {
  tmp <- file.path(tempdir(), "tsf_force_probe")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  writeLines("data", file.path(tmp, "gone.tsv"))
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  suppressWarnings(system2("Rscript",
    c("scripts/clean_results.R", "--results-dir", tmp, "--force"),
    stdin = "/dev/null", stdout = TRUE, stderr = TRUE))
  !file.exists(file.path(tmp, "gone.tsv")) })

check("--cores reaches the chromosome-parallel stages", {
  # maxt_cores() used to read detectCores() directly and ignore both --cores
  # and N_WORKERS, so the documented flag silently did nothing on maxT,
  # condition and CLEAN -- three stages, including the slow one. The source
  # attribute is what distinguishes an honoured flag from a fallback.
  chrom <- stats::setNames(vector("list", 24), paste0("c", 1:24))
  auto <- maxt_cores(chrom)
  flagged <- maxt_cores(chrom, list(cores = 2))
  identical(attr(auto, "source"), "automatic") &&
    identical(attr(flagged, "source"), "--cores") })

check("maxt_cores agrees with local_workers on the same request", {
  # One worker-count policy, not two. If these ever disagree, some stages are
  # obeying a different rule than the flag documents.
  chrom <- stats::setNames(vector("list", 12), paste0("c", 1:12))
  opt <- list(cores = 3)
  identical(as.integer(maxt_cores(chrom, opt)),
            as.integer(local_workers(opt, n_tasks = 12L, default_max = 12L))) })

check("the chromosome count caps the workers", {
  # These stages split on chromosomes, so the ceiling is structural: with 4
  # chromosomes the 5th core has nothing to do, whatever the flag says.
  chrom <- stats::setNames(vector("list", 4), paste0("c", 1:4))
  as.integer(maxt_cores(chrom, list(cores = 64))) <= 4L })

check("--config reads a different file and a flag still overrides it", {
  # A derived config keeps one run tree's settings under one name instead of a
  # row of flags retyped per command. It must inherit the defaults, and the
  # command line must still win over it, or the precedence documented in the
  # usage is a fiction.
  tmp <- tempfile(fileext = ".R")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(c(
    'base <- local({',
    '  .tsf_root <- Sys.getenv("TSF_ROOT", unset = getwd())',
    '  source(file.path(.tsf_root, "config", "project.R"), local = TRUE)$value',
    '})',
    'modifyList(base, list(results_dir = "results_test_tree"))'
  ), tmp)
  cfg <- load_project_config(tmp)
  base <- load_project_config("config/project.R")
  identical(cfg$results_dir, "results_test_tree") &&
    identical(cfg$geo_dir, base$geo_dir) &&
    length(cfg) == length(base) })

check("a missing --config file fails by name", {
  out <- suppressWarnings(system2("./tsf",
    c("status", "--config", "config/does_not_exist.R"),
    stdout = TRUE, stderr = TRUE))
  any(grepl("does_not_exist.R", paste(out, collapse = " "), fixed = TRUE)) })

check("every run logs the results tree it resolved", {
  # Writing to the wrong tree is silent by nature, so the resolved path is
  # announced before any work happens.
  out <- suppressWarnings(system2("./tsf",
    c("window", "--results-dir", "results_probe_tree", "--dry-run"),
    stdout = TRUE, stderr = TRUE))
  any(grepl("results_probe_tree", paste(out, collapse = " "), fixed = TRUE)) })

check("the help text lists every stage that exists", {
  # The stage list in the usage string used to be typed out by hand, and had
  # drifted: it omitted `consensus`, so the help told the user about a pipeline
  # one stage shorter than the real one. It is now built from stage_names, and
  # this asserts they cannot come apart again.
  # stage_names lives in R/stages.R, which this suite does not source; read it
  # from there so the assertion compares the help against the real definition
  # rather than against a second hand-written list.
  src <- paste(readLines("R/stages.R", warn = FALSE), collapse = " ")
  decl <- regmatches(src, regexpr('stage_names <- c\\([^)]*\\)', src))
  stages <- regmatches(decl, gregexpr('"[a-z]+"', decl))[[1]]
  stages <- gsub('"', "", stages)
  usage <- paste(system2("./tsf", "--help", stdout = TRUE, stderr = TRUE),
                 collapse = " ")
  length(stages) >= 9L &&
    all(vapply(stages, function(x) grepl(x, usage, fixed = TRUE), logical(1))) })

check("the help says datasets are positional", {
  # `--datasets GSE135251` is the natural guess and is wrong; the error for it
  # should say so rather than only listing valid options.
  out <- suppressWarnings(system2("./tsf", c("window", "--datasets", "GSE1"),
                                  stdout = TRUE, stderr = TRUE))
  any(grepl("positional", paste(out, collapse = " "), fixed = TRUE)) })

# --- exclusion by accession ---------------------------------------------------
#
# sample_filter evaluates one column across all samples, which is right for
# "this series mixes MASLD with HBV" and wrong for "these two samples are
# mislabelled": excluding fibrosis stages 1 and 2 to drop two miscoded controls
# would also delete every real F1 and F2 in the cohort.

check("exclude_samples drops exactly the named samples", {
  lab <- data.frame(sample_id = paste0("GSM", 1:5),
                    condition = c("Control", "Control", "F1", "F2", "F3"),
                    in_scope = TRUE, filtered_out = FALSE,
                    label_resolved = TRUE, condition_selected = TRUE,
                    stringsAsFactors = FALSE)
  ph <- data.frame(geo_accession = lab$sample_id, stringsAsFactors = FALSE)
  out <- apply_sample_filter(lab, ph,
    list(id = "T", exclude_samples = c("GSM1", "GSM2")))
  identical(out$in_scope, c(FALSE, FALSE, TRUE, TRUE, TRUE)) &&
    all(out$filtered_out[1:2]) })

check("the real F1 and F2 survive a control-only exclusion", {
  # The failure this replaces: two miscoded controls removed, and 100 genuine
  # mid-stage samples removed with them.
  lab <- data.frame(sample_id = paste0("GSM", 1:5),
                    condition = c("Control", "Control", "F1", "F2", "F3"),
                    in_scope = TRUE, filtered_out = FALSE,
                    label_resolved = TRUE, condition_selected = TRUE,
                    stringsAsFactors = FALSE)
  ph <- data.frame(geo_accession = lab$sample_id, stringsAsFactors = FALSE)
  out <- apply_sample_filter(lab, ph,
    list(id = "T", exclude_samples = c("GSM1", "GSM2")))
  all(out$in_scope[out$condition %in% c("F1", "F2", "F3")]) })

check("a misspelled accession aborts instead of excluding nothing", {
  lab <- data.frame(sample_id = paste0("GSM", 1:3), condition = "Control",
                    in_scope = TRUE, filtered_out = FALSE,
                    label_resolved = TRUE, condition_selected = TRUE,
                    stringsAsFactors = FALSE)
  ph <- data.frame(geo_accession = lab$sample_id, stringsAsFactors = FALSE)
  e <- tryCatch(apply_sample_filter(lab, ph,
    list(id = "T", exclude_samples = c("GSM1", "GSM_typo"))),
    error = function(e) e)
  inherits(e, "error") && grepl("GSM_typo", conditionMessage(e), fixed = TRUE) })

check("the two miscoded GSE135251 controls are excluded by the config", {
  cfg <- source("config/datasets/GSE135251.R")$value
  identical(sort(cfg$exclude_samples), c("GSM3998224", "GSM3998341")) })

check("the control guardrail compares raw labels against the baseline CLASS", {
  # The guardrail compared a raw label with a class name. Once the three
  # healthy groups merged into "Controles", nothing matched and it stopped
  # firing -- a dataset declaring has_control_cohort = FALSE could assign
  # controls unnoticed.
  src <- paste(readLines("R/labels.R", warn = FALSE), collapse = " ")
  grepl("control_labels", src, fixed = TRUE) })

# --- licensing ---------------------------------------------------------------
#
# One licence, CC BY-NC 4.0, everywhere. This repository briefly declared MIT in
# README.md and R/bundle.R while LICENSE said CC BY-NC, which is worse than
# either: the bundle told recipients they could use commercially what LICENSE
# forbade. The assertions below exist because that drift was invisible to
# review and visible only to whoever read the wrong file.

LICENCE_ID <- "CC-BY-NC-4.0"
LICENCE_PROSE <- "CC BY-NC 4.0"

check("LICENSE declares CC BY-NC 4.0 by SPDX id and by name", {
  txt <- paste(readLines("LICENSE", warn = FALSE), collapse = " ")
  grepl(LICENCE_ID, txt, fixed = TRUE) &&
    grepl(LICENCE_PROSE, txt, fixed = TRUE) })

check("no file still claims MIT", {
  files <- c("LICENSE", "README.md", "R/bundle.R",
             "scripts/sonify_tissuespectf.py", "Makefile")
  files <- files[file.exists(files)]
  claims <- unlist(lapply(files, function(f) {
    hits <- grep("MIT", readLines(f, warn = FALSE), value = TRUE)
    # A line explaining what the licence is NOT may name MIT; a line asserting
    # the licence may not.
    hits[!grepl("not|rather than|instead of|Apache|GPL", hits)]
  }))
  length(claims) == 0 })

check("every SPDX header agrees with LICENSE", {
  files <- c(list.files("R", pattern = "[.]R$", full.names = TRUE),
             list.files("scripts", pattern = "[.](R|py)$", full.names = TRUE),
             list.files("ml", pattern = "[.]py$", full.names = TRUE),
             "app/app.R", "tsf")
  files <- files[file.exists(files)]
  bad <- unlist(lapply(files, function(f) {
    spdx <- grep("SPDX-License-Identifier", readLines(f, n = 15, warn = FALSE),
                 value = TRUE)
    spdx[!grepl(LICENCE_ID, spdx, fixed = TRUE)]
  }))
  length(bad) == 0 })

check("the bundle does not promise commercial use", {
  # The bundle is a release asset a collaborator downloads and reads. Under NC
  # it must not tell them the opposite of what LICENSE says.
  src <- paste(readLines("R/bundle.R", warn = FALSE), collapse = " ")
  grepl(LICENCE_PROSE, src, fixed = TRUE) &&
    !grepl("commercially or not", src, fixed = TRUE) })

check("the bundle ships the licence", {
  src <- readLines("R/bundle.R", warn = FALSE)
  any(grepl('file.copy("LICENSE"', src, fixed = TRUE)) })

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
