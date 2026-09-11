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

# --- provenance survives the TSV round trip (finding 5) ------------------------
pf <- file.path(tmp, "profile.tsv"); write_reference_profile(ref, pf, list(grid_axis = "gene"))
ref2 <- read_reference_profile(pf)
man <- read_tsv_tsf(sub("\\.tsv$", "_manifest.tsv", pf))
check("n_ref_samples and source datasets are columns, and a manifest with md5 is written",
      identical(ref2$n_ref_samples[1], 30L) && identical(ref2$source_datasets[1], "REF") &&
        identical(man$value[man$key == "md5"], unname(tools::md5sum(pf))))
apply_reference_profile("COH", project, ref2, "profile.tsv", profile_path = pf)
ap <- read_reference_profile_applied(file.path(tmp, "COH"))
check("the applied marker records path, digest and n_ref_samples",
      identical(ap$profile_digest, unname(tools::md5sum(pf))) && ap$n_ref_samples == 30L &&
        identical(ap$profile_path, pf))

# --- a re-ingest after the correction is detected (finding 6) -----------------
Sys.sleep(1.2)
fresh <- mk("COH", 8, shift = 4)           # "ingest --force": expression.tsv is raw again, newer
check("a rewritten expression.tsv makes the dataset raw again",
      is.null(suppressWarnings(read_reference_profile_applied(file.path(tmp, "COH")))))
apply_reference_profile("COH", project, ref2, "profile.tsv", profile_path = pf)
e5 <- read_tsv_tsf(file.path(tmp, "COH", "expression.tsv"))
check("re-applying after a re-ingest uses the NEW ingest output, not the stale backup",
      abs(mean(as.matrix(e5[1:10, -1])) - 4) < 0.3)

# --- the query is corrected from the reference object, on both axes (3, 4) ----
grid <- data.frame(gene_id = sub("\\..*$", "", rownames(mc)), chr = "1", grid_index = 1:100, grid_N = 100L)
ref_obj <- list(reference_profile = list(profile = ref2[, c("gene_id", "ref_median")],
                                         digest = ap$profile_digest, grid_axis = "gene"),
                params = list(differential = TRUE))
y <- asinh(rep(1000, 100)); ids <- rownames(mc)
Sys.unsetenv("TSF_REFERENCE_PROFILE")
yc <- apply_query_reference_profile(y, ids, ref_obj)
check("a query against a differential library is corrected from the stored profile",
      identical(attr(yc, "reference_profile"), ap$profile_digest) &&
        abs(mean(yc - (asinh(1000) - ref2$ref_median))) < 1e-9)
raw_obj <- list(reference_profile = NULL, params = list())
check("a raw library leaves the query alone", identical(apply_query_reference_profile(y, ids, raw_obj), y))
Sys.setenv(TSF_REFERENCE_PROFILE = pf)
check("TSF_REFERENCE_PROFILE set against a RAW library is refused",
      inherits(try(apply_query_reference_profile(y, ids, raw_obj), silent = TRUE), "try-error"))
other <- file.path(tmp, "other.tsv"); write_tsv_tsf(ref2[1:50, ], other)
Sys.setenv(TSF_REFERENCE_PROFILE = other)
check("a variable pointing at a DIFFERENT profile than the library's is refused",
      inherits(try(apply_query_reference_profile(y, ids, ref_obj), silent = TRUE), "try-error"))
Sys.unsetenv("TSF_REFERENCE_PROFILE")
check("a query mostly off the profile's positions is refused",
      inherits(try(apply_query_reference_profile(y, paste0("ENSX", 1:100), ref_obj), silent = TRUE), "try-error"))
bin_obj <- ref_obj; bin_obj$reference_profile$profile$gene_id <- paste0("bin_1_", 1:100)
check("on the bp axis the profile is keyed by bin id and the same function applies",
      abs(mean(apply_query_reference_profile(y, paste0("bin_1_", 1:100), bin_obj) -
                 (asinh(1000) - ref2$ref_median))) < 1e-9)

if (fails) { cat(" ", fails, "check(s) failed\n"); quit(status = 1) } else cat(" All tests passed. \n")
