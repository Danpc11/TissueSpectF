# Liver fibrosis staging, with the healthy states kept apart.
#
# THREE LEVELS OF IDENTITY: class_id = tissue::state::condition
#
# `state` separates healthy tissue from diseased tissue; `condition` names the
# specific class within that state. This is what lets two different kinds of
# healthy liver coexist without one absorbing the other:
#
# THREE HEALTHY CLASSES, NOT ONE
#
#   Normal_histology         a biopsy with normal histology and no NAFLD
#                            activity (NAS = 0) taken WITHIN a biopsy series.
#                            GSE162694 (31) and GSE130970 (6).
#   Control_disease_cohort   a subject outside the disease cohort of a NAFLD
#                            study. A cohort statement, not a histology one:
#                            two of GSE135251's ten carry incidental fibrosis
#                            of stage 1 and 2, so this class is not a strictly
#                            healthy extreme and must not be used as one
#                            without auditing those two out.
#   Control_external_study   the control group of a DIFFERENT disease study.
#                            GSE142530's twelve are controls for an alcohol
#                            cohort, recruited, sampled and processed under that
#                            study's protocol.
#
# All three sit under state = healthy and none absorbs the others. Merging them
# into one Control would assume exactly what is worth testing: that a healthy
# liver looks the same whichever study recruited it. If their spectra turn out
# to agree, a Healthy_consensus class can be created and defended. Equivalence
# is demonstrated first; it is not asserted by sharing a label.
#
# Folding Normal_histology into F0 would likewise merge histologically normal
# livers with NAFLD patients who simply have no fibrosis yet.
#
# `progression` is the ordered subset. Transitions and ordinal trends run over
# it alone: F0 -> F1 is a step, Control -> Normal_histology is not.
list(
  id       = "liver_fibrosis",
  tissue   = "liver",
  levels   = c("Control_disease_cohort", "Control_external_study",
               "Normal_histology", paste0("F", 0:4)),
  states   = c(Control_disease_cohort = "healthy",
               Control_external_study = "healthy",
               Normal_histology = "healthy",
               F0 = "disease", F1 = "disease", F2 = "disease",
               F3 = "disease", F4 = "disease"),
  conditions = c(Control_disease_cohort = "Control_disease_cohort",
                 Control_external_study = "Control_external_study",
                 Normal_histology = "Normal_histology",
                 F0 = "NAFLD_fibrosis_F0", F1 = "NAFLD_fibrosis_F1",
                 F2 = "NAFLD_fibrosis_F2", F3 = "NAFLD_fibrosis_F3",
                 F4 = "NAFLD_fibrosis_F4"),
  # What role a class plays, independent of which one is the baseline. Deriving
  # `cohort` by comparing against the baseline alone marked
  # Control_external_study as disease, because it is not the baseline -- a
  # conclusion drawn from a naming convention rather than from the biology.
  cohort_roles = c(Control_disease_cohort = "control",
                   Control_external_study = "control",
                   Normal_histology       = "within_disease_normal",
                   F0 = "disease", F1 = "disease", F2 = "disease",
                   F3 = "disease", F4 = "disease"),

  progression = paste0("F", 0:4),
  ordered  = TRUE,
  baseline = "Control_disease_cohort",
  description = "METAVIR-style fibrosis stage, plus two distinct healthy states"
)
