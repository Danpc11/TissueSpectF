# The minimal vocabulary: any tissue, any disease, two groups.
#
# Use this for a dataset that only distinguishes affected from unaffected. It is
# unordered, so no transition tables are produced -- there is only one contrast
# and `compare` reports it as such.
list(
  id       = "case_control",
  tissue   = NA_character_,     # set per dataset
  levels   = c("Control", "Case"),
  ordered  = FALSE,
  baseline = "Control",
  description = "Two-group design with no intermediate states"
)
