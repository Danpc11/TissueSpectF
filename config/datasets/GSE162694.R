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
  description   = "NAFLD liver biopsies, single cohort, no non-disease controls",
  counts_file   = "GSE162694.tsv.gz",
  series_matrix = "GSE162694_series_matrix.txt.gz",
  count_id_type = "ENTREZID",

  has_control_cohort = FALSE,

  sample_id_column = "geo_accession",
  fibrosis_column  = "fibrosis stage",
  normal_histology_terms = c("normal liver histology", "normal"),

  condition_rules = list(
    list(
      id             = "title_stage_token",
      type           = "title_token",
      column         = "title",
      pattern        = "^F[0-4]$",
      normal_tokens  = c("N"),
      normal_assign  = "F0"      # NOT Control: same cohort, stage 0
    ),
    list(
      id           = "biopsy_fibrosis_stage",
      type         = "fibrosis_stage",
      column       = "fibrosis stage",
      normal_terms = c("normal liver histology", "normal")
    )
  ),

  notes = paste(
    "Titles carrying more than one F-token are treated as ambiguous and fall",
    "through to the fibrosis-stage rule instead of taking the first match."
  )
)
