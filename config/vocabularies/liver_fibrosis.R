# Liver fibrosis staging, with the three healthy groups as ONE control class.
#
# THREE LEVELS OF IDENTITY: class_id = tissue::state::condition
#
# `state` separates healthy tissue from diseased tissue; `condition` names the
# specific class within that state. This is what lets two different kinds of
# healthy liver coexist without one absorbing the other:
#
# THREE HEALTHY GROUPS, ONE CLASS
#
# The raw labels stay separate -- ingest, the audit trail and `cohort_roles`
# still record which cohort a control came from -- but they map to a single
# condition, `Controles`. The source publications place all three in the same
# group: liver tissue without disease.
#
# An earlier version kept them as three classes, on the argument that a healthy
# liver might not look the same whichever study recruited it. Two things
# decided against it, neither of which depends on any accuracy figure:
#
#   Control_disease_cohort exists only in GSE135251 and Control_external_study
#   only in GSE142530, so leave-one-cohort-out has no second cohort to learn
#   either from. In a real run GSE142530 was skipped as a hold-out entirely
#   ("fewer than two shared classes"), 10 samples were dropped from validation,
#   and both classes sat among the centroids competing for every prediction
#   while winning zero.
#
#   `states` already marked all three healthy. The split lived in the label.
#
# Merged, Controles has 56 samples across FOUR cohorts and GSE142530 becomes a
# usable hold-out.
#
# WHAT IT COSTS: the class is heterogeneous in clinical context.
# Normal_histology is a NAFLD-cohort patient with no fibrosis yet, not a donor
# liver; Control_external_study is the control arm of an alcohol study under a
# different protocol. And two of GSE135251's ten controls carried incidental
# fibrosis of stage 1 and 2 -- those are excluded by accession in that dataset's
# config, which is why the count is 56 and not 58.
#
# The three groups, and where they come from:
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
# F0 is NOT folded in. F0 is a histological stage in a NAFLD patient -- a
# measurement -- while these three are the absence of the disease context.
# Merging them would assume what the Controles vs F0 contrast exists to test.
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
  # The merge lives here and nowhere else: three raw labels, one condition.
  conditions = c(Control_disease_cohort = "Controles",
                 Control_external_study = "Controles",
                 Normal_histology       = "Controles",
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
  # Must name a class that exists after the merge.
  baseline = "Controles",
  description = "METAVIR-style fibrosis stage, with the three healthy groups merged into one Controles class"
)
