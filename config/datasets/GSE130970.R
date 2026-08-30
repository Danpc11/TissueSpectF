# GSE130970 -- NAFLD liver biopsies (Hoang et al.), 78 samples.
#
# NO CONTROL COHORT, and no explicit histology label either. There is no
# disease/control column: the phenotype table carries only tissue, sex, age,
# lobular inflammation, ballooning, steatosis, NAS and fibrosis stage. So a
# histologically normal biopsy can only be recognised from the scores.
#
# WHICH SAMPLES ARE Normal_histology, AND WHY
#
# NAFLD is defined by steatosis. A biopsy with no steatosis and no
# ballooning is not fatty liver disease however it is filed, so those samples
# do not belong in F0 alongside patients who have NAFLD but have not yet
# fibrosed. Among the 25 samples at fibrosis stage 0:
#
#   steatosis 0 AND ballooning 0   ->  6   Normal_histology
#   the rest                       ->  19  F0
#
# Two of those six carry NAS = 1 from grade-1 lobular inflammation, which is
# nonspecific and common in otherwise normal liver. The other four have NAS = 0.
# GSE162694's Normal_histology group is NAS = 0 throughout, so the two classes
# are close but not identically defined -- worth remembering when their spectra
# are compared. Restricting to `nas = "0"` instead would give 4 samples and an
# exact match to that cohort's definition; it is one line below, and which one
# applies is a declared choice rather than an accident.
#
# One further sample has steatosis 0 but fibrosis above 0. It stays in its
# fibrosis class: the compound rule requires stage 0 as well, because a liver
# with fibrosis is not histologically normal whatever its steatosis score.
#
# Sample sizes: Normal_histology = 6, F0 = 19, F1 = 28, F2 = 9, F3 = 14, F4 = 2.
# F4 with two samples cannot support a consensus spectrum (phase alignment
# needs roughly log(n_frequencies/q) samples); it will be reported as present
# but not testable rather than silently averaged in.
list(
  id            = "GSE130970",
  description   = "NAFLD liver biopsies, single cohort, no non-disease controls",
  tissue        = "liver",
  vocabulary    = "liver_fibrosis",
  counts_file   = "GSE130970.tsv.gz",
  series_matrix = "GSE130970_series_matrix.txt.gz",
  count_id_type = "ENTREZID",

  has_control_cohort = FALSE,

  sample_id_column = "geo_accession",
  fibrosis_column  = "fibrosis stage",
  normal_histology_terms = c("normal liver histology", "normal"),

  # Covariates worth carrying for the confounder diagnostics: this is the only
  # one of the three cohorts that publishes sex and age per sample.
  covariate_columns = c(sex = "^Sex", age = "age at biopsy",
                        nas = "nafld activity score",
                        steatosis = "steatosis grade",
                        ballooning = "cytological ballooning grade",
                        inflammation = "lobular inflammation grade"),

  condition_rules = list(
    # Runs FIRST: the general stage rule would otherwise swallow these six into
    # F0 before anything looked at the histology scores.
    list(
      id     = "normal_histology_scores",
      type   = "compound",
      assign = "Normal_histology",
      all_of = list(
        list(column = "steatosis grade", values = "0"),
        list(column = "cytological ballooning grade", values = "0"),
        list(column = "fibrosis stage", values = "0")
        # Add list(column = "nafld activity score", values = "0") to use the
        # stricter NAS = 0 definition (4 samples) matching GSE162694 exactly.
      )
    ),
    list(
      id           = "biopsy_fibrosis_stage",
      type         = "fibrosis_stage",
      column       = "fibrosis stage",
      normal_terms = c("normal liver histology", "normal")
    )
  )
)
