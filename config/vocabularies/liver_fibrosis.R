# Liver fibrosis staging, with the healthy states kept apart.
#
# THREE LEVELS OF IDENTITY: class_id = tissue::state::condition
#
# `state` separates healthy tissue from diseased tissue; `condition` names the
# specific class within that state. This is what lets two different kinds of
# healthy liver coexist without one absorbing the other:
#
#   Control            a subject OUTSIDE the disease cohort. A cohort statement,
#                      not a histology one -- in GSE135251 two of the ten
#                      controls carry incidental fibrosis of stage 1 and 2.
#   Normal_histology   a biopsy with normal histology and no NAFLD activity
#                      (NAS = 0) taken WITHIN a biopsy series. In GSE162694 the
#                      phenotype table distinguishes these 31 samples from the
#                      35 with fibrosis stage 0, whose NAS runs 1-5.
#
# Folding Normal_histology into F0 would merge histologically normal livers with
# NAFLD patients who simply have no fibrosis yet. Calling it Control would merge
# two different definitions of healthy. Keeping both under state = healthy makes
# the comparison between them an answerable question instead of an assumption.
#
# `progression` is the ordered subset. Transitions and ordinal trends run over
# it alone: F0 -> F1 is a step, Control -> Normal_histology is not.
list(
  id       = "liver_fibrosis",
  tissue   = "liver",
  levels   = c("Control", "Normal_histology", paste0("F", 0:4)),
  states   = c(Control = "healthy", Normal_histology = "healthy",
               F0 = "disease", F1 = "disease", F2 = "disease",
               F3 = "disease", F4 = "disease"),
  conditions = c(Control = "Control", Normal_histology = "Normal_histology",
                 F0 = "NAFLD_fibrosis_F0", F1 = "NAFLD_fibrosis_F1",
                 F2 = "NAFLD_fibrosis_F2", F3 = "NAFLD_fibrosis_F3",
                 F4 = "NAFLD_fibrosis_F4"),
  progression = paste0("F", 0:4),
  ordered  = TRUE,
  baseline = "Control",
  description = "METAVIR-style fibrosis stage, plus two distinct healthy states"
)
