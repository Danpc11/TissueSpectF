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
# NOT THE SAME CLASS AS GSE135251's CONTROLS
#
# They map to Control_external_study, not to Control_disease_cohort. Merging the
# two would assume what is worth testing -- that a healthy liver looks the same
# whichever study recruited it -- and would do it by choosing a label. If the two
# spectra agree, a Healthy_consensus class can be defended afterwards.
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
# Eleven samples is also near the floor for phase alignment (about nine are
# needed to make a Rayleigh test reachable), so expect this class to sit close
# to exploratory.
list(
  id            = "GSE142530",
  description   = "Alcohol-related liver disease series; normal controls only",
  tissue        = "liver",
  vocabulary    = "liver_fibrosis",
  # A comma-separated file with Ensembl ids in the first column, symbols in the
  # second, and a DESCRIPTIVE SECOND ROW carrying the real sample names
  # (Control_Lille 389 and so on) above columns named RB_N1, RB_N2, ...
  #
  # That second row is why a rename is not enough: left in place it makes every
  # count column character, and the whole matrix silently becomes text. The
  # reader extracts it as the sample map, renames the columns with it, and
  # removes it before anything is coerced to numeric. The map is written to
  # count_column_map.tsv so the correspondence is auditable rather than implied.
  counts_file   = "GSE142530_Annoted-RNAseq-with-SampleIDs.csv.gz",
  series_matrix = "GSE142530_series_matrix.txt.gz",
  count_id_type = "ENSEMBL",

  counts_spec = list(
    sep = ",", id_column = 1L, symbol_column = 2L, sample_map_row = 1L,
    # The descriptive row carries a recruiting-site token the series matrix does
    # not: "Control_Lille 389" there against "Control_389" here. The two
    # substitutions below reconcile them, and count_column_map.tsv records the
    # raw label, the mapped id and whether it matched, so the correspondence is
    # auditable rather than asserted.
    map_transform = list(
      list(pattern = "_(Lille|TPF)\\s+", replacement = "_")
    ),
    # One column reads not.used. It is excluded by name rather than left to fail
    # matching, where the symptom would be an unexplained missing sample.
    exclude_columns = c("not.used")
  ),

  # Eleven, not twelve: the file carries twelve control labels but one column is
  # not.used. Ingest aborts if the count changes, so a number quoted anywhere
  # else cannot drift away from the data.
  expected_n_samples = 11L,

  has_control_cohort = TRUE,

  # After the rename the columns carry the descriptive names, which match the
  # series matrix titles rather than the accessions.
  sample_id_column = "title",

  sample_filter = list(
    list(column = "disease state", values = c("Normal"))
  ),

  keep_conditions = c("Control_external_study"),

  condition_rules = list(
    list(
      id     = "normal_disease_state",
      type   = "column_match",
      column = "disease state",
      values = c("Normal"),
      assign = "Control_external_study"
    )
  ),

  notes = paste(
    "Series lists 12 Normal samples; the count matrix has one not.used column,",
    "so 11 are usable. The 16 alcohol-related samples are excluded."
  )
)
