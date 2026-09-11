# test_differential.R -- reference profile, differential signal, recount3 parsing.
suppressWarnings({ source("R/utils_io.R"); tsf_load_all("R") })
fails <- 0L
check <- function(name, expr) {
  ok <- tryCatch(isTRUE(expr), error = function(e) { message("  ERROR ", conditionMessage(e)); FALSE })
  cat(if (ok) "  PASS   " else "  FAIL   ", name, "\n"); if (!ok) fails <<- fails + 1L
}

# --- explode_sample_attributes ------------------------------------------------
x <- c("tissue;;liver|fibrosis stage;;3|Sex;;male", "tissue;;liver|fibrosis stage;;0", NA)
a <- explode_sample_attributes(x)
check("sample attributes explode into lower-case columns",
      identical(sort(colnames(a)), c("fibrosis stage", "sex", "tissue")) &&
        identical(a[["fibrosis stage"]], c("3", "0", NA)) && is.na(a$sex[2]))

# --- recount3 URL shape -------------------------------------------------------
u <- recount3_urls("LIVER", "gtex")
check("gtex url uses last two characters as the shard",
      grepl("/gene_sums/ER/LIVER/gtex.gene_sums.LIVER.G026.gz$", u$gene_sums))
u2 <- recount3_urls("SRP217231", "sra")
check("sra url uses last two characters as the shard",
      grepl("/gene_sums/31/SRP217231/sra.gene_sums.SRP217231.G026.gz$", u2$gene_sums) &&
        grepl("sra.sra.SRP217231.MD.gz$", u2$metadata))

# --- differential signal ------------------------------------------------------
set.seed(1)
tmp <- tempfile(); dir.create(tmp)
project <- list(interim_dir = tmp)
mk <- function(id, n, shift = 0) {
  m <- matrix(rnorm(100 * n, 5, 1), 100, n,
              dimnames = list(paste0("ENSG", 1:100, ".", 1:100 %% 3), paste0(id, "_s", 1:n)))
  m[1:10, ] <- m[1:10, ] + shift
  d <- file.path(tmp, id); dir.create(d)
  write_tsv_tsf(data.frame(gene_id = rownames(m), m, check.names = FALSE), file.path(d, "expression.tsv"))
  write_tsv_tsf(data.frame(gene_id = rownames(m), chr = "1", grid_index = 1:100, grid_N = 100L),
                file.path(d, "genes.tsv"))
  write_tsv_tsf(data.frame(sample_id = colnames(m), condition = "A", dataset_id = id),
                file.path(d, "samples.tsv"))
  invisible(m)
}
mk("REF", 30); mc <- mk("COH", 8, shift = 2)
# load_dataset needs a config; use a local stand-in with the same fields
load_dataset <- function(id, project, ...) {
  d <- file.path(project$interim_dir, id)
  e <- read_tsv_tsf(file.path(d, "expression.tsv"))
  m <- as.matrix(e[, -1, drop = FALSE]); rownames(m) <- e$gene_id
  list(id = id, samples = read_tsv_tsf(file.path(d, "samples.tsv")),
       genes = read_tsv_tsf(file.path(d, "genes.tsv")), expression = m)
}
ref <- build_reference_profile("REF", project, min_samples = 20L)
check("the profile is a median per gene with unversioned ids",
      nrow(ref) == 100 && !any(grepl("\\.", ref$gene_id)) &&
        abs(ref$ref_median[ref$gene_id == "ENSG50"] - stats::median(load_dataset("REF", project)$expression[50, ])) < 1e-9)
check("fewer reference samples than min_samples aborts",
      inherits(try(build_reference_profile("COH", project, min_samples = 20L), silent = TRUE), "try-error"))

d <- subtract_reference_profile(mc, ref)
check("the shifted genes deviate by the shift and the others by ~0",
      abs(mean(d[1:10, ]) - 2) < 0.3 && abs(mean(d[11:100, ])) < 0.2)
ref_part <- ref[1:90, ]
d2 <- subtract_reference_profile(mc, ref_part)
check("genes absent from the profile become NA, never zero",
      all(is.na(d2[91:100, ])) && attr(d2, "n_unmatched") == 10)

apply_reference_profile("COH", project, ref_part, "test")
e2 <- read_tsv_tsf(file.path(tmp, "COH", "expression.tsv"))
g2 <- read_tsv_tsf(file.path(tmp, "COH", "genes.tsv"))
check("apply rewrites expression.tsv and genes.tsv consistently and keeps the raw copy",
      nrow(e2) == 90 && nrow(g2) == 90 && file.exists(file.path(tmp, "COH", "expression_raw.tsv")) &&
        abs(mean(as.matrix(e2[1:10, -1])) - 2) < 0.3)
apply_reference_profile("COH", project, ref, "test2")
e3 <- read_tsv_tsf(file.path(tmp, "COH", "expression.tsv"))
check("re-applying starts from the raw copy, so subtractions do not stack",
      nrow(e3) == 100 && abs(mean(as.matrix(e3[1:10, -1])) - 2) < 0.3)
restore_raw_expression("COH", project)
e4 <- read_tsv_tsf(file.path(tmp, "COH", "expression.tsv"))
check("restore brings back the ingest output", abs(mean(as.matrix(e4[1:10, -1])) - 7) < 0.3)

# --- query path picks the profile up from the environment --------------------
pf <- file.path(tmp, "profile.tsv"); write_tsv_tsf(ref, pf)
Sys.setenv(TSF_REFERENCE_PROFILE = pf)
q <- stats::setNames(rep(1000, 100), rownames(mc))
qs <- query_signal(q, unit = "tpm")
check("query_signal subtracts the active profile", isTRUE(attr(qs, "reference_profile")) &&
        abs(mean(qs - (asinh(1000) - ref$ref_median))) < 1e-9)
Sys.unsetenv("TSF_REFERENCE_PROFILE")
check("without the variable query_signal is unchanged", is.null(attr(query_signal(q, "tpm"), "reference_profile")))

if (fails) { cat(" ", fails, "check(s) failed\n"); quit(status = 1) } else cat(" All tests passed. \n")
