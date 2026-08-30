# GSE130970 -- NAFLD liver biopsies (Hoang et al.), 78 samples.
#
# NO CONTROL COHORT. Every sample is a liver biopsy scored on the NAFLD system:
# there is no disease/control field, and the phenotype table carries only
# tissue, sex, age, lobular inflammation, ballooning, steatosis, NAS and
# fibrosis stage. Some samples score 0 on activity and steatosis, but they come
# from the same biopsy series as the rest, so by the rule this project applies
# everywhere they are F0, not Control: Control is a statement about the cohort a
# subject belongs to, never about how normal a biopsy looks.
#
# Fibrosis stage is already numeric 0-4 with no textual "normal liver
# histology" values, so the stage rule applies directly and no title parsing is
# needed.
#
# Sample sizes per stage: F0 = 25, F1 = 28, F2 = 9, F3 = 14, F4 = 2.
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
    list(
      id           = "biopsy_fibrosis_stage",
      type         = "fibrosis_stage",
      column       = "fibrosis stage",
      normal_terms = c("normal liver histology", "normal")
    )
  )
)
