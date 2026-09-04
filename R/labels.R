# labels.R -- the single place where a sample gets its condition label.
#
# ---------------------------------------------------------------------------
# CONTROLLED VOCABULARY
# ---------------------------------------------------------------------------
# Control : liver from a subject that does NOT belong to the disease cohort
#           (healthy donor / non-NAFLD surgical control). It is a *cohort*
#           statement, not a histology statement.
# F0..F4  : fibrosis stage within the disease cohort, by biopsy.
#
# The rule that resolves the incongruence found in the original scripts:
#
#   A sample is NEVER labelled Control because its fibrosis stage is 0, or
#   because its histology is described as normal. Control requires explicit
#   evidence of a separate, non-diseased cohort in the phenotype table.
#
# Consequence: GSE135251 has Control (disease == "Control") *and* F0 (NAFLD
# without fibrosis) as distinct groups. GSE162694 is a single biopsy cohort:
# its "N" / "normal liver histology" samples are F0, not Control, and the
# dataset declares has_control_cohort = FALSE. Any cross-dataset comparison
# then runs only over conditions present on both sides (see compare layer).
# ---------------------------------------------------------------------------

# Default vocabulary, used when a dataset config names none. Nothing in the
# pipeline may assume these particular levels: the vocabulary is data, loaded
# from config/vocabularies/, so a new tissue or design needs a config file and
# no code change. See load_vocabulary() in config.R.
TSF_DEFAULT_VOCABULARY <- "liver_fibrosis"

#' The ordered subset over which transitions and trends are defined.
#'
#' Not every level is a step of a progression. Control -> Normal_histology is
#' not a transition, and treating it as one would report a "change" between two
#' unrelated healthy states. A vocabulary declares `progression` explicitly;
#' when it does not, the ordered levels are the progression.
tsf_progression <- function(voc) {
  if (!isTRUE(voc$ordered)) return(character(0))
  voc$progression %||% tsf_levels(voc)
}

#' class_id = tissue::state::condition
#'
#' Three levels because two are not enough. `state` separates healthy from
#' diseased tissue; `condition` names the class within it. With two levels,
#' healthy liver from a non-disease cohort and healthy liver from a biopsy
#' series would have to share one label or invent unrelated ones. With three,
#' both sit under liver::healthy and stay distinguishable, and TCGA's
#' adjacent-normal will fit later without renaming anything.
tsf_class_id <- function(tissue, condition, voc) {
  state <- if (is.null(voc$states)) "healthy" else
    unname(voc$states[condition])
  state[is.na(state)] <- "unknown"
  cond_name <- if (is.null(voc$conditions)) condition else {
    nm <- unname(voc$conditions[condition]); nm[is.na(nm)] <- condition[is.na(nm)]; nm
  }
  ifelse(is.na(condition), NA_character_,
         paste(tissue %||% "unknown", state, cond_name, sep = "::"))
}

#' The role a class plays: control, disease, or something in between.
#'
#' Derived from the vocabulary's `cohort_roles` when it declares them. The old
#' rule -- equal to the baseline means control, anything else means disease --
#' broke as soon as there was more than one healthy class: Control_external_study
#' is not the baseline, so it was labelled disease. That column feeds reports and
#' any future balancing, so the error would have propagated quietly.
#'
#' Without cohort_roles the baseline comparison is kept, which is right for a
#' two-group vocabulary and is what those declare.
tsf_cohort_role <- function(condition, voc) {
  roles <- voc$cohort_roles
  if (!is.null(roles)) {
    out <- unname(roles[condition])
    out[is.na(condition)] <- NA_character_
    unknown <- !is.na(condition) & is.na(out)
    if (any(unknown)) {
      tsf_warn("No cohort_role declared for: ",
               paste(unique(condition[unknown]), collapse = ", "),
               "; recorded as unknown rather than guessed")
      out[unknown] <- "unknown"
    }
    return(out)
  }
  baseline <- voc$baseline %||% "Control"
  ifelse(is.na(condition), NA_character_,
         ifelse(condition == baseline, "control", "disease"))
}

#' Levels of a dataset's vocabulary, in order.
tsf_levels <- function(x) {
  if (is.character(x)) return(x)
  if (!is.null(x$levels)) return(x$levels)
  if (!is.null(x$condition_levels)) return(x$condition_levels)
  tsf_abort("No condition levels declared")
}

#' Numeric fibrosis stage from a raw phenotype value.
#'
#' Accepts "3", "F3", "stage 3", "fibrosis stage: 3" and the textual
#' descriptions listed in `normal_terms` (mapped to 0). Anything else is NA --
#' never silently pasted into a label. The original GSE135251 script produced
#' the literal condition "FNA" here.
parse_fibrosis_stage <- function(x, normal_terms = character(0)) {
  v <- clean_pheno_value(x)
  low <- tolower(v)
  v[!is.na(low) & low %in% tolower(normal_terms)] <- "0"
  v <- sub("^[Ff](?=[0-4]$)", "", v, perl = TRUE)
  v <- sub("^stage\\s*", "", v, ignore.case = TRUE)
  num <- suppressWarnings(as.numeric(v))
  num[!is.na(num) & (num < 0 | num > 4 | num != round(num))] <- NA_real_
  num
}

#' Extract a stage token from a free-text sample title.
#'
#' Tokenises on non-alphanumerics and looks for F0..F4. `normal_tokens` (e.g.
#' "N") are mapped to `normal_assign`, which is a *stage* label, not Control.
parse_title_stage <- function(title, pattern = "^F[0-4]$",
                              normal_tokens = character(0),
                              normal_assign = "F0") {
  toks <- strsplit(as.character(title), "[^A-Za-z0-9]+")
  vapply(toks, function(tk) {
    tk_up <- toupper(tk[nzchar(tk)])
    hit <- tk_up[grepl(pattern, tk_up)]
    if (length(hit) == 1L) return(hit)
    if (length(hit) > 1L) return(NA_character_)   # ambiguous title -> unresolved
    if (length(normal_tokens) && any(tk_up %in% toupper(normal_tokens))) {
      return(normal_assign)
    }
    NA_character_
  }, character(1))
}

#' Apply one declarative rule from a dataset config.
#'
#' Rule types:
#'   column_match  : exact (case-insensitive) match on a phenotype column
#'   title_token   : stage token parsed from the title
#'   fibrosis_stage: numeric stage from a phenotype column
apply_condition_rule <- function(rule, pheno) {
  n <- nrow(pheno)
  out <- rep(NA_character_, n)

  if (identical(rule$type, "compound")) {
    # Every sub-condition must hold. Needed when a class is defined by a
    # combination of scores rather than by one field: GSE130970 publishes no
    # disease/control column, so a histologically normal biopsy can only be
    # recognised by its scores being zero together.
    hit <- rep(TRUE, n)
    for (sub in rule$all_of) {
      col <- find_pheno_column(pheno, sub$column)
      if (is.null(col)) {
        tsf_warn("compound rule '", rule$id, "' needs column matching '",
                 sub$column, "', which is absent; the rule cannot fire")
        return(out)
      }
      v <- tolower(clean_pheno_value(pheno[[col]]))
      hit <- hit & !is.na(v) & v %in% tolower(as.character(sub$values))
    }
    out[hit] <- rule$assign
    return(out)
  }

  if (identical(rule$type, "column_match")) {
    col <- find_pheno_column(pheno, rule$column)
    if (is.null(col)) return(out)
    v <- tolower(clean_pheno_value(pheno[[col]]))
    hit <- !is.na(v) & v %in% tolower(rule$values)
    out[hit] <- rule$assign

  } else if (identical(rule$type, "title_token")) {
    col <- find_pheno_column(pheno, rule$column %||% "title")
    if (is.null(col)) return(out)
    out <- parse_title_stage(pheno[[col]],
                             pattern       = rule$pattern %||% "^F[0-4]$",
                             normal_tokens = rule$normal_tokens %||% character(0),
                             normal_assign = rule$normal_assign %||% "F0")

  } else if (identical(rule$type, "fibrosis_stage")) {
    col <- find_pheno_column(pheno, rule$column)
    if (is.null(col)) return(out)
    stage <- parse_fibrosis_stage(pheno[[col]], rule$normal_terms %||% character(0))
    out[!is.na(stage)] <- paste0("F", stage[!is.na(stage)])

  } else {
    tsf_abort("Unknown condition rule type: ", rule$type)
  }
  out
}

find_pheno_column <- function(pheno, pattern) {
  if (is.null(pattern)) return(NULL)
  hits <- colnames(pheno)[grepl(pattern, colnames(pheno), ignore.case = TRUE)]
  if (length(hits)) hits[1] else NULL
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Resolve the condition of every sample, in declared rule order.
#'
#' Returns one row per input sample with the label, the numeric stage, the rule
#' that fired, and an explicit `keep` flag. Nothing is guessed and nothing is
#' dropped silently: unresolved samples come back with condition = NA and a
#' reason, and the caller decides.
harmonize_conditions <- function(pheno, dataset_config) {
  rules <- dataset_config$condition_rules
  if (!length(rules)) tsf_abort("Dataset config declares no condition_rules")

  n <- nrow(pheno)
  condition <- rep(NA_character_, n)
  rule_used <- rep(NA_character_, n)

  for (rule in rules) {
    pending <- is.na(condition)
    if (!any(pending)) break
    got <- apply_condition_rule(rule, pheno)
    take <- pending & !is.na(got)
    condition[take] <- got[take]
    rule_used[take] <- rule$id
  }

  levels_now <- tsf_levels(dataset_config$vocabulary_spec %||%
                             dataset_config$condition_levels)
  bad <- !is.na(condition) & !(condition %in% levels_now)
  if (any(bad)) {
    tsf_warn(sum(bad), " sample(s) resolved to a label outside the vocabulary (",
              paste(unique(condition[bad]), collapse = ", "), "); dropped.")
    condition[bad] <- NA_character_
    rule_used[bad] <- NA_character_
  }

  # `condition` here is the RAW label a rule assigned; `baseline` is a CLASS
  # name. Comparing them directly worked only while every healthy group was its
  # own class. Once the three merged into "Controles", no raw label equalled
  # the baseline and this guardrail silently stopped firing -- a dataset
  # declaring has_control_cohort = FALSE could assign controls unnoticed.
  #
  # Map the raw labels to their classes first, so the comparison happens in one
  # space. Every raw label whose class is the baseline counts.
  vocab <- dataset_config$vocabulary_spec %||% list()
  baseline <- vocab$baseline %||% "Control"
  control_labels <- if (length(vocab$conditions)) {
    names(vocab$conditions)[unname(vocab$conditions) %in% baseline]
  } else {
    baseline
  }
  if (!length(control_labels)) control_labels <- baseline

  if (isFALSE(dataset_config$has_control_cohort) &&
      any(condition %in% control_labels, na.rm = TRUE)) {
    offending <- unique(condition[condition %in% control_labels])
    tsf_abort("Dataset ", dataset_config$id, " declares has_control_cohort = FALSE ",
              "but a rule assigned ", paste(offending, collapse = ", "),
              ", which maps to the baseline class ", baseline,
              ". Fix the config, not the data.")
  }

  stage_col <- find_pheno_column(pheno, dataset_config$fibrosis_column)
  fibrosis_stage <- if (!is.null(stage_col)) {
    parse_fibrosis_stage(pheno[[stage_col]], dataset_config$normal_histology_terms %||% character(0))
  } else rep(NA_real_, n)

  # For F-labelled samples the stage is implied by the label; keep both and
  # report disagreement rather than picking one silently.
  implied <- suppressWarnings(as.numeric(sub("^F", "", condition)))
  mismatch <- !is.na(implied) & !is.na(fibrosis_stage) & implied != fibrosis_stage
  if (any(mismatch)) {
    tsf_warn(sum(mismatch), " sample(s) disagree between label and fibrosis column; ",
              "label kept, see label_audit.tsv")
  }

  data.frame(
    sample_id      = clean_pheno_value(pheno[[dataset_config$sample_id_column %||% "geo_accession"]]),
    dataset_id     = dataset_config$id,
    tissue         = dataset_config$tissue %||% NA_character_,
    vocabulary     = (dataset_config$vocabulary_spec %||% list())$id %||% NA_character_,
    state          = if (is.null(dataset_config$vocabulary_spec$states)) "healthy"
                     else unname(dataset_config$vocabulary_spec$states[condition]),
    class_id       = tsf_class_id(dataset_config$tissue, condition,
                                  dataset_config$vocabulary_spec %||% list()),
    condition      = condition,
    fibrosis_stage = ifelse(is.na(implied), fibrosis_stage, implied),
    fibrosis_stage_reported = fibrosis_stage,
    label_rule     = rule_used,
    label_mismatch = mismatch,
    cohort         = tsf_cohort_role(condition,
                                     dataset_config$vocabulary_spec %||% list()),
    label_resolved = !is.na(condition),
    in_scope       = TRUE,
    condition_selected = TRUE,
    keep           = !is.na(condition),
    stringsAsFactors = FALSE
  )
}

#' Drop samples a dataset should not contribute at all.
#'
#' Two mechanisms, both explicit, because a cohort can be useful for part of
#' what it contains and wrong for the rest:
#'
#'   sample_filter    keep only samples whose fields match (or exclude those
#'                    that match). This is how a mixed-etiology series
#'                    contributes its MASLD samples and not its HBV, ALD or
#'                    cholestatic ones -- those are different diseases whose
#'                    fibrosis is not the same object.
#'   keep_conditions  keep only certain classes. A series can be sound for
#'                    advanced fibrosis and too thin or too selected for the
#'                    early stages.
#'
#' A sample whose filter field is missing or unreadable is DROPPED, never kept
#' on the assumption that it belongs. Silently admitting a hepatitis B cirrhosis
#' into a MASLD class would be the most expensive error available here.
apply_sample_filter <- function(labels, pheno, dataset_config) {
  if (!"label_resolved" %in% colnames(labels)) {
    labels$label_resolved <- !is.na(labels$condition)
  }
  if (!"in_scope" %in% colnames(labels)) labels$in_scope <- TRUE
  if (!"condition_selected" %in% colnames(labels)) {
    labels$condition_selected <- TRUE
  }
  if (!"filtered_out" %in% colnames(labels)) labels$filtered_out <- FALSE

  # Named exclusions, applied before the column filters.
  #
  # sample_filter evaluates ONE column per entry, across ALL samples. That is
  # right for "this series mixes MASLD with HBV, keep the MASLD" and wrong for
  # "these two specific samples are mislabelled": a filter excluding fibrosis
  # stages 1 and 2 to drop two miscoded controls would also delete every real
  # F1 and F2 in the cohort. Naming the accessions is narrower and auditable --
  # the config says which samples left and the log says whether they were found.
  #
  # Refuses to proceed if a named sample is absent. A typo in an accession would
  # otherwise exclude nothing and leave the run looking clean.
  excl <- dataset_config$exclude_samples
  if (!is.null(excl) && length(excl)) {
    excl <- as.character(excl)
    hit <- labels$sample_id %in% excl
    missing <- setdiff(excl, labels$sample_id)

    # PARTIAL match means a typo: some accessions resolved, so this is the right
    # series and the rest should have resolved too. Abort -- excluding one of
    # two mislabelled samples is worse than excluding neither, because the run
    # looks corrected.
    #
    # NO match means a different sample universe: a synthetic self-check, a
    # subset, or a config applied to another series. Warn and carry on; there is
    # nothing here to exclude and nothing to be wrong about.
    if (length(missing) && any(hit)) {
      tsf_abort("exclude_samples for ", dataset_config$id, " names ",
                length(missing), " sample(s) not in this series: ",
                paste(missing, collapse = ", "),
                ", while ", sum(hit), " other(s) matched. A partial match is a ",
                "typo, not a different series: check the accessions.")
    }
    if (length(missing) && !any(hit)) {
      tsf_warn("exclude_samples for ", dataset_config$id, ": none of the ",
               length(excl), " named sample(s) are in this series (",
               paste(excl, collapse = ", "), "). Nothing was excluded. This is ",
               "expected for synthetic or subset data, and a mistake anywhere ",
               "else.")
    }
    if (any(hit)) {
      tsf_log("  exclude_samples: dropped ", sum(hit), " named sample(s) (",
              paste(excl, collapse = ", "), ")")
      labels$in_scope[hit] <- FALSE
      labels$filtered_out[hit] <- TRUE
    }
  }

  filt <- dataset_config$sample_filter
  if (!is.null(filt)) {
    for (f in filt) {
      col <- find_pheno_column(pheno, f$column)
      if (is.null(col)) {
        tsf_abort("sample_filter for ", dataset_config$id, " needs a column ",
                  "matching '", f$column, "', which is absent. Refusing to ",
                  "proceed: without it every sample would be admitted, ",
                  "including the ones the filter exists to exclude.")
      }
      v <- tolower(clean_pheno_value(pheno[[col]]))
      ok <- if (!is.null(f$values)) {
        !is.na(v) & v %in% tolower(as.character(f$values))
      } else if (!is.null(f$pattern)) {
        !is.na(v) & grepl(f$pattern, v, ignore.case = TRUE)
      } else if (!is.null(f$exclude_values)) {
        !is.na(v) & !(v %in% tolower(as.character(f$exclude_values)))
      } else {
        tsf_abort("sample_filter entry needs values, pattern or exclude_values")
      }
      previously_in_scope <- labels$in_scope
      dropped <- sum(previously_in_scope & !ok)
      if (dropped) {
        tsf_log("  filter '", f$column, "': dropped ", dropped, " sample(s) (",
                paste(utils::head(sort(unique(v[previously_in_scope & !ok])), 6),
                      collapse = ", "), ")")
      }
      labels$in_scope <- labels$in_scope & ok
      labels$filtered_out <- labels$filtered_out | !ok
    }
  }

  keep_cond <- dataset_config$keep_conditions
  if (!is.null(keep_cond)) {
    out <- labels$in_scope & labels$label_resolved &
      !(labels$condition %in% keep_cond)
    if (any(out)) {
      tsf_log("  keep_conditions: dropped ", sum(out), " sample(s) in ",
              paste(sort(unique(labels$condition[out])), collapse = ", "),
              " (this dataset contributes only ",
              paste(keep_cond, collapse = ", "), ")")
    }
    labels$condition_selected[out] <- FALSE
  }

  labels$keep <- labels$label_resolved & labels$in_scope &
    labels$condition_selected
  labels
}

#' Summarise labelling and refuse to continue on a bad parse.
audit_labels <- function(labels, dataset_config, max_unresolved_frac = 0.05,
                         min_n_per_condition = 5L) {
  levels_now <- tsf_levels(dataset_config$vocabulary_spec %||%
                             dataset_config$condition_levels)
  if (!"label_resolved" %in% colnames(labels)) {
    labels$label_resolved <- !is.na(labels$condition)
  }
  if (!"in_scope" %in% colnames(labels)) labels$in_scope <- TRUE
  if (!"condition_selected" %in% colnames(labels)) {
    labels$condition_selected <- TRUE
  }

  scope <- labels$in_scope & labels$condition_selected
  n_scope <- sum(scope)
  resolved_scope <- scope & labels$label_resolved
  unresolved <- sum(scope & !labels$label_resolved)
  filtered <- sum(!labels$in_scope)
  condition_excluded <- sum(labels$in_scope & !labels$condition_selected)

  tsf_log("Labels for ", dataset_config$id, ": ",
          sum(resolved_scope), "/", n_scope,
          " resolved within analysis scope; ", filtered, " filtered out",
          if (condition_excluded > 0L)
            paste0("; ", condition_excluded,
                   " resolved sample(s) excluded by keep_conditions")
          else "")

  tab <- table(factor(labels$condition[resolved_scope], levels = levels_now))
  for (cnd in names(tab)) tsf_log("  ", cnd, ": ", tab[[cnd]])

  if (n_scope == 0L) {
    tsf_abort("No samples remain within the analysis scope for ",
              dataset_config$id, ". Check sample_filter and keep_conditions.")
  }
  if (unresolved / n_scope > max_unresolved_frac) {
    tsf_abort(unresolved, "/", n_scope, " in-scope samples unresolved in ",
               dataset_config$id,
               " (> ", max_unresolved_frac * 100, "%). Check condition_rules.")
  }
  low <- names(tab)[tab > 0 & tab < min_n_per_condition]
  if (length(low)) {
    tsf_warn("Low sample size (n < ", min_n_per_condition, ") in: ",
              paste(low, collapse = ", "),
              " -- underpowered for per-condition spectra.")
  }
  empty <- names(tab)[tab == 0]
  if (length(empty)) {
    tsf_log("Conditions absent from ", dataset_config$id, ": ",
             paste(empty, collapse = ", "),
             " (recorded; cross-dataset comparison will skip them)")
  }
  invisible(list(counts = tab, unresolved = unresolved,
                 filtered_out = filtered,
                 condition_excluded = condition_excluded,
                 n_in_scope = n_scope,
                 present = names(tab)[tab > 0]))
}

#' Conditions usable for a cross-dataset comparison.
#'
#' Replaces the silent Reduce(intersect) failure: if a condition is missing on
#' one side, that is reported here rather than showing up as an empty signature.
comparable_conditions <- function(present_by_dataset, levels_now = NULL) {
  common <- Reduce(intersect, present_by_dataset)
  if (!is.null(levels_now)) common <- levels_now[levels_now %in% common]
  for (nm in names(present_by_dataset)) {
    missing <- setdiff(unlist(present_by_dataset), present_by_dataset[[nm]])
    if (length(missing)) {
      tsf_warn(nm, " lacks: ", paste(missing, collapse = ", "))
    }
  }
  if (!length(common)) tsf_abort("No condition is present in every dataset.")
  tsf_log("Comparable conditions: ", paste(common, collapse = ", "))
  common
}
