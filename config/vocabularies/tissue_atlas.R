# One "condition" per tissue, for building a cross-tissue reference.
#
# Unordered by construction: liver is not between lung and kidney. With this
# vocabulary the pipeline answers "what is the characteristic spectrum of each
# tissue" rather than "how does it change with stage", and `compare` reports
# which components are shared across tissues and which are specific to one.
#
# Levels are declared per dataset, since which tissues a study covers is a
# property of that study. GTEx-style data would list them all here.
list(
  id       = "tissue_atlas",
  tissue   = NA_character_,
  levels   = NULL,          # dataset config supplies condition_levels
  ordered  = FALSE,
  baseline = NULL,
  description = "One level per tissue; no ordering, no transitions"
)
