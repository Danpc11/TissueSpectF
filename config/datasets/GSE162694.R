# GSE162694 -- single NAFLD biopsy cohort (Pantano et al.). There is NO
# non-disease control group: the samples labelled "N" in the title are biopsies
# with normal liver histology, i.e. fibrosis stage 0 within the same cohort.
#
# This is the fix for the main incongruence in the original scripts: they mapped
# "N" to Healthy via the title, while the fallback mapped the same biology
# ("normal liver histology") to F0 -- so one sample's label depended on which
# branch fired, and "Healthy" meant something different from GSE135251's Healthy.
list(
  id            = "GSE162694",
  tissue      = "liver",
  vocabulary  = "liver_fibrosis",
  description   = "NAFLD liver biopsies, single cohort, no non-disease controls",
  counts_file   = "GSE162694.tsv.gz",
  series_matrix = "GSE162694_series_matrix.txt.gz",
  count_id_type = "ENTREZID",

  has_control_cohort = FALSE,   # no non-NAFLD cohort; the N group is Normal_histology

  sample_id_column = "geo_accession",
  fibrosis_column  = "fibrosis stage",
  # Deliberately no normal_histology_terms: "normal liver histology" must NOT
  # be parsed as stage 0 here. It is its own class.
  covariate_columns = c(sex = "^Sex", age = "^age", nas = "nas score"),

  condition_rules = list(
    # The phenotype field is decisive, so it goes first; the title token is the
    # fallback for anything it does not cover.
    list(
      id     = "normal_histology_field",
      type   = "column_match",
      column = "fibrosis stage",
      values = c("normal liver histology"),
      assign = "Normal_histology"
    ),
    list(
      id             = "title_stage_token",
      type           = "title_token",
      column         = "title",
      pattern        = "^F[0-4]$",
      normal_tokens  = c("N"),
      normal_assign  = "Normal_histology"
    ),
    list(
      id           = "biopsy_fibrosis_stage",
      type         = "fibrosis_stage",
      column       = "fibrosis stage"
    )
  ),

  notes = paste(
    "Expected: Normal_histology 31, F0 35, F1 30, F2 27, F3 8, F4 12.",
    "Titles carrying more than one F-token are treated as ambiguous and fall",
    "through to the fibrosis-stage rule instead of taking the first match."
  )
)
