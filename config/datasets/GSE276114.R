# GSE276114 -- mixed-etiology liver fibrosis, 177 liver samples.
#
# Verified against the series matrix. The phenotype table carries exactly two
# informative fields, and both matter:
#
#   disease        ARLD 14   CVH 82   MASLD 81
#   disease group  F0-2 39   F3 24    F4 114
#
#   disease x disease group:
#     ARLD  F4 14
#     CVH   F0-2 13   F3 11   F4 58
#     MASLD F0-2 26   F3 13   F4 42
#
# ONLY MASLD, AND ONLY F3/F4
#
# Alcohol-related and chronic viral hepatitis fibrosis are different diseases.
# Their fibrosis is not the same object as MASLD fibrosis -- different injury,
# different cell populations, different tissue context -- so pooling them into
# NAFLD_fibrosis_F4 would put one class label on a mixture of three diseases.
# `sample_filter` admits MASLD only, which leaves 81.
#
# THE STAGE FIELD IS A BIN, NOT A STAGE
#
# `disease group` reads F0-2, F3, F4 -- not 0,1,2,3,4. "F0-2" spans three of
# this project's classes and cannot be resolved to any one of them; assigning it
# to F0, F1 or F2 would be inventing a stage the data does not record. Those 26
# MASLD samples are therefore filtered out explicitly rather than left to fail
# label resolution, so the count of unresolved samples stays a signal that
# something is wrong rather than a routine 32%.
#
# That leaves F3 = 13 and F4 = 42, which is what this cohort contributes.
#
# CAUTION FOR INTERPRETATION
#
# 114 of the 177 samples are cirrhosis. At that scale the tissue is usually
# explants and resections rather than needle biopsies, so this cohort's F4 may
# differ systematically in tissue source from the other three -- and the
# difference is aligned with the class. A component appearing in F4 only after
# this cohort joins must be checked against the dataset-prediction diagnostic
# before it is read as biology.
list(
  id            = "GSE276114",
  description   = "Mixed-etiology liver fibrosis; MASLD F3/F4 only",
  tissue        = "liver",
  vocabulary    = "liver_fibrosis",
  counts_file   = "GSE276114.tsv.gz",          # VERIFY the name GEO actually uses
  series_matrix = "GSE276114_series_matrix.txt.gz",
  count_id_type = "ENTREZID",                  # VERIFY against the counts header

  has_control_cohort = FALSE,

  sample_id_column = "geo_accession",

  sample_filter = list(
    # "^disease:" and not "^disease": the parser names characteristics columns
    # <field>:ch1, and a bare "^disease" would match "disease group:ch1" first
    # -- filtering on the stage bin while believing it filtered on etiology.
    list(column = "^disease:", values = c("MASLD")),
    # F0-2 is a bin spanning three classes; it contributes nothing resolvable.
    list(column = "disease group", values = c("F3", "F4"))
  ),

  keep_conditions = c("F3", "F4"),

  # The field holds the label directly, so no numeric parsing is involved.
  condition_rules = list(
    list(
      id     = "disease_group_stage",
      type   = "column_match",
      column = "disease group",
      values = c("F3"),
      assign = "F3"
    ),
    list(
      id     = "disease_group_stage_f4",
      type   = "column_match",
      column = "disease group",
      values = c("F4"),
      assign = "F4"
    )
  ),

  notes = "Verified: MASLD F3 = 13, F4 = 42. Stop if it differs."
)
