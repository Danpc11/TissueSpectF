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

LFFT_CONDITION_LEVELS <- c("Control", paste0("F", 0:4))

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
    lfft_abort("Unknown condition rule type: ", rule$type)
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
  if (!length(rules)) lfft_abort("Dataset config declares no condition_rules")

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

  bad <- !is.na(condition) & !(condition %in% LFFT_CONDITION_LEVELS)
  if (any(bad)) {
    lfft_warn(sum(bad), " sample(s) resolved to a label outside the vocabulary (",
              paste(unique(condition[bad]), collapse = ", "), "); dropped.")
    condition[bad] <- NA_character_
    rule_used[bad] <- NA_character_
  }

  if (isFALSE(dataset_config$has_control_cohort) && any(condition %in% "Control", na.rm = TRUE)) {
    lfft_abort("Dataset ", dataset_config$id, " declares has_control_cohort = FALSE ",
               "but a rule assigned Control. Fix the config, not the data.")
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
    lfft_warn(sum(mismatch), " sample(s) disagree between label and fibrosis column; ",
              "label kept, see label_audit.tsv")
  }

  data.frame(
    sample_id      = clean_pheno_value(pheno[[dataset_config$sample_id_column %||% "geo_accession"]]),
    dataset_id     = dataset_config$id,
    condition      = condition,
    fibrosis_stage = ifelse(is.na(implied), fibrosis_stage, implied),
    fibrosis_stage_reported = fibrosis_stage,
    label_rule     = rule_used,
    label_mismatch = mismatch,
    cohort         = ifelse(is.na(condition), NA_character_,
                            ifelse(condition == "Control", "control", "disease")),
    keep           = !is.na(condition),
    stringsAsFactors = FALSE
  )
}

#' Summarise labelling and refuse to continue on a bad parse.
audit_labels <- function(labels, dataset_config, max_unresolved_frac = 0.05,
                         min_n_per_condition = 5L) {
  n <- nrow(labels)
  unresolved <- sum(!labels$keep)
  lfft_log("Labels for ", dataset_config$id, ": ", n - unresolved, "/", n, " resolved")
  tab <- table(factor(labels$condition[labels$keep], levels = LFFT_CONDITION_LEVELS))
  for (cnd in names(tab)) lfft_log("  ", cnd, ": ", tab[[cnd]])

  if (unresolved / max(n, 1L) > max_unresolved_frac) {
    lfft_abort(unresolved, "/", n, " samples unresolved in ", dataset_config$id,
               " (> ", max_unresolved_frac * 100, "%). Check condition_rules.")
  }
  low <- names(tab)[tab > 0 & tab < min_n_per_condition]
  if (length(low)) {
    lfft_warn("Low sample size (n < ", min_n_per_condition, ") in: ",
              paste(low, collapse = ", "),
              " -- underpowered for per-condition spectra.")
  }
  empty <- names(tab)[tab == 0]
  if (length(empty)) {
    lfft_log("Conditions absent from ", dataset_config$id, ": ",
             paste(empty, collapse = ", "),
             " (recorded; cross-dataset comparison will skip them)")
  }
  invisible(list(counts = tab, unresolved = unresolved,
                 present = names(tab)[tab > 0]))
}

#' Conditions usable for a cross-dataset comparison.
#'
#' Replaces the silent Reduce(intersect) failure: if a condition is missing on
#' one side, that is reported here rather than showing up as an empty signature.
comparable_conditions <- function(present_by_dataset) {
  common <- Reduce(intersect, present_by_dataset)
  common <- LFFT_CONDITION_LEVELS[LFFT_CONDITION_LEVELS %in% common]
  for (nm in names(present_by_dataset)) {
    missing <- setdiff(unlist(present_by_dataset), present_by_dataset[[nm]])
    if (length(missing)) {
      lfft_warn(nm, " lacks: ", paste(missing, collapse = ", "))
    }
  }
  if (!length(common)) lfft_abort("No condition is present in every dataset.")
  lfft_log("Comparable conditions: ", paste(common, collapse = ", "))
  common
}
