# Fibrosis staging of liver biopsies (NAFLD/NASH cohorts).
#
# Control is a COHORT statement: liver from a subject outside the disease
# cohort. It is never inferred from a fibrosis stage of 0 or from a description
# of normal histology. A dataset with no such cohort declares
# has_control_cohort = FALSE and its stage-0 samples are F0.
list(
  id       = "liver_fibrosis",
  tissue   = "liver",
  levels   = c("Control", paste0("F", 0:4)),
  ordered  = TRUE,       # transitions are meaningful between adjacent levels
  baseline = "Control",
  description = "METAVIR-style fibrosis stage plus a non-disease control group"
)
