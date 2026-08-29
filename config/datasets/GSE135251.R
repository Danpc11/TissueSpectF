# GSE135251 -- NAFLD biopsy cohort (Govaere et al.) WITH a separate group of
# non-NAFLD control livers. Control and F0 are therefore different populations
# and must stay separate: F0 means NAFLD without fibrosis.
list(
  id            = "GSE135251",
  tissue      = "liver",
  vocabulary  = "liver_fibrosis",
  description   = "NAFLD liver biopsies + non-NAFLD controls",
  counts_file   = "GSE135251.tsv.gz",
  series_matrix = "GSE135251_series_matrix.txt.gz",
  count_id_type = "ENTREZID",

  has_control_cohort = TRUE,

  sample_id_column = "geo_accession",
  fibrosis_column  = "fibrosis stage",
  normal_histology_terms = c("normal liver histology", "normal"),

  # Ordered: the first rule that resolves a sample wins.
  condition_rules = list(
    list(
      id     = "non_disease_cohort",
      type   = "column_match",
      column = "^disease",
      values = c("Control", "control", "normal"),
      assign = "Control"
    ),
    list(
      id            = "biopsy_fibrosis_stage",
      type          = "fibrosis_stage",
      column        = "fibrosis stage",
      normal_terms  = c("normal liver histology", "normal")
    )
  ),

  notes = paste(
    "Samples with a non-numeric fibrosis stage are left unresolved and dropped",
    "with a report. The previous implementation produced the literal label",
    "'FNA' for these."
  )
)
