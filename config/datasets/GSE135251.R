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

  # Two of the ten non-NAFLD controls carry incidental fibrosis. The `disease`
  # column calls them Control and so does `group in paper`; only
  # `fibrosis stage` disagrees:
  #
  #   GSM3998224   disease: Control   fibrosis stage: 1
  #   GSM3998341   disease: Control   fibrosis stage: 2
  #
  # `non_disease_cohort` matches on `disease` and is first, so it wins before
  # the stage is read, and both would enter Controles as healthy liver.
  #
  # Excluded by accession, not by a rule on `fibrosis stage`: sample_filter
  # evaluates ONE column per entry across ALL samples, so excluding stages 1
  # and 2 would also delete the 47 F1 and 53 F2 of the NAFLD cohort.
  #
  # This departs from the publication, which groups all ten as controls. The
  # control class here is histologically clean, not the paper's recruitment
  # group, and that is a claim to state rather than to assume.
  exclude_samples = c("GSM3998224", "GSM3998341"),

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
      assign = "Control_disease_cohort"
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
