# GSE142530 -- alcohol-related liver disease series; only its controls are used.
#
# Verified against the series matrix. 28 liver biopsies:
#
#   disease state   Normal 12   Alcoholic Hepatitis 10   Alcohol-related Cirrhosis 6
#
# The twelve are titled Control_<id> while the others are Severe AH_<id> and
# Cirrhosis without AH_<id>. `disease state` is a cohort statement here -- these
# are subjects without alcohol-related liver disease, not biopsies from the
# diseased series that happened to look normal -- so they map to Control, the
# same class as GSE135251's ten. That gives liver::healthy::Control a second
# cohort and lifts it out of single_cohort/provisional.
#
# THE ARGUMENT AGAINST, WHICH IS NOT NOTHING
#
# These are controls for an ALCOHOL study, not for MASLD. The class is defined
# at the level of tissue and state -- a healthy liver is a healthy liver -- and
# that is exactly what makes a tissue-level class useful. But the controls of a
# different disease study are recruited differently, sampled differently and
# processed in a different laboratory, and all of that is aligned with the
# cohort. If Control's spectrum turns out to depend on which cohort it came
# from, that is a finding about the method, not about liver.
#
# The alcohol-related samples are NOT used. ARLD fibrosis is a different disease
# from MASLD fibrosis, and this project does not have an alcohol vocabulary yet.
# When it does, they become liver::disease::ALD_* without any change here.
#
# Twelve samples is also near the floor for phase alignment (about nine are
# needed to make a Rayleigh test reachable), so expect this class to sit close
# to exploratory.
list(
  id            = "GSE142530",
  description   = "Alcohol-related liver disease series; normal controls only",
  tissue        = "liver",
  vocabulary    = "liver_fibrosis",
  counts_file   = "GSE142530.tsv.gz",          # VERIFY the name GEO actually uses
  series_matrix = "GSE142530_series_matrix.txt.gz",
  count_id_type = "ENTREZID",                  # VERIFY against the counts header

  has_control_cohort = TRUE,

  sample_id_column = "geo_accession",

  sample_filter = list(
    list(column = "disease state", values = c("Normal"))
  ),

  keep_conditions = c("Control"),

  condition_rules = list(
    list(
      id     = "normal_disease_state",
      type   = "column_match",
      column = "disease state",
      values = c("Normal"),
      assign = "Control"
    )
  ),

  notes = "Verified: Control = 12. The 16 alcohol-related samples are excluded."
)
