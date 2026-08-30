# app.R -- local desktop app for TissueSpectF.
#
#   ./tsf app          (or: Rscript -e 'shiny::runApp("app")')
#
# Runs entirely on the machine it is started on: no upload leaves the computer,
# no server is contacted. Shiny is the only dependency beyond base R, and it is
# needed only for this file -- the CLI and the whole pipeline work without it.
#
# DESIGN CONSTRAINT: a matcher always returns a best match. This interface must
# never present one without the evidence that it means anything. The reference's
# out-of-cohort accuracy is shown first, above the result, in language that
# survives being read in a hurry.

library(shiny)

# Only what a query needs. No ingest, no stages: the reference is
# self-contained, so the app never rebuilds a grid from local config.
for (f in c("utils_io", "grid", "fingerprint", "reference")) {
  source(file.path("..", "R", paste0(f, ".R")))
}

# The reference path comes from the caller (`./tsf app --results-dir ...` sets
# TSF_APP_REFERENCE), never from re-reading config/project.R. Re-reading it was
# how the app ended up looking for a reference on the cluster's default path
# while the user had pointed the CLI somewhere else entirely. It can also be
# changed from the interface.
initial_ref_path <- Sys.getenv("TSF_APP_REFERENCE", "")

ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-family: -apple-system, Segoe UI, Roboto, sans-serif;
           max-width: 900px; margin: 0 auto; padding: 24px; }
    .banner { padding: 14px 18px; border-radius: 6px; margin: 16px 0;
              font-weight: 500; line-height: 1.5; }
    .ok    { background: #e6f4ea; border-left: 5px solid #34a853; }
    .weak  { background: #fef7e0; border-left: 5px solid #f9ab00; }
    .bad   { background: #fce8e6; border-left: 5px solid #ea4335; }
    .muted { color: #5f6368; font-size: 0.92em; }
    table  { width: 100%; border-collapse: collapse; margin-top: 10px; }
    th, td { text-align: left; padding: 7px 10px; border-bottom: 1px solid #e8eaed; }
    .top   { font-weight: 600; }
  "))),

  h2("TissueSpectF"),
  p(class = "muted",
    "Identify an expression profile by the shape of its chromosome-ordered ",
    "spectrum. Everything runs locally; nothing is uploaded."),

  textInput("ref_path", "Reference library (.rds)", value = initial_ref_path,
            width = "100%"),

  uiOutput("reference_banner"),

  hr(),
  fileInput("query", "Counts file (TSV: gene id column, then one column per sample)",
            accept = c(".tsv", ".txt", ".csv", ".gz")),
  selectInput("unit", "What the values are",
              choices = c("Raw counts" = "counts", "TPM" = "tpm", "CPM" = "cpm",
                          "Already log/asinh transformed" = "logged"),
              selected = "counts"),
  checkboxInput("show_all", "Show every class, not just the top five", FALSE),

  uiOutput("result"),

  hr(),
  p(class = "muted",
    "Gene ids must be Ensembl or Entrez, matching the annotation the reference ",
    "was built on. Raw counts are expected; normalisation happens here. ",
    "A match is a similarity to a class centroid, not a diagnosis.")
)

server <- function(input, output, session) {

  ref <- reactive({
    p <- input$ref_path
    if (is.null(p) || !nzchar(p) || !file.exists(p)) return(NULL)
    r <- readRDS(p)
    if (is.null(r$grid)) return(structure(list(), class = "stale_reference"))
    r
  })

  output$reference_banner <- renderUI({
    r <- ref()
    if (is.null(r)) {
      return(div(class = "banner bad",
                 strong("No reference library at that path."), br(),
                 "Build one with ", code("./tsf run --to spectra"), " then ",
                 code("./tsf reference"),
                 ", or point the box above at an existing reference.rds."))
    }
    if (inherits(r, "stale_reference")) {
      return(div(class = "banner bad",
                 strong("Reference too old."), br(),
                 "It carries no grid, so a query cannot be scored on the grid it ",
                 "was built with. Rebuild it: ", code("./tsf reference"), "."))
    }
    v <- r$validation
    if (is.null(v)) {
      return(div(class = "banner bad",
                 strong("Uncalibrated reference."), br(),
                 "Built from a single cohort, so there is no out-of-cohort ",
                 "validation. Any match below is a suggestion, not an ",
                 "identification."))
    }
    lift <- v$accuracy - v$baseline
    cls <- if (lift <= 0.02) "bad" else if (v$accuracy < 0.5) "weak" else "ok"
    head <- if (lift <= 0.02) {
      "This reference does not beat guessing."
    } else if (v$accuracy < 0.5) {
      "Weak reference: better than chance, wrong more often than right."
    } else {
      "Reference validated out of cohort."
    }
    calib <- v$calibration
    reject_note <- if (!is.null(calib) && !is.na(calib$separability_auc)) {
      sprintf(" A query scoring below the per-class %.0fth percentile of correct
               held-out matches is reported UNKNOWN; correct and incorrect
               matches separate with AUC %.2f. That threshold bounds how often a
               true member is rejected, not how often an out-of-domain sample is
               accepted -- no out-of-domain sample was in the validation.",
              100 * calib$quantile, calib$separability_auc)
    } else ""
    excl <- if (!is.null(v$excluded) && v$excluded$samples > 0) {
      sprintf(" %d sample(s) in %d class(es) were excluded from validation for
               not being shared across cohorts.", v$excluded$samples, v$excluded$classes)
    } else ""
    div(class = paste("banner", cls),
        strong(head), br(),
        sprintf("Trained and tested across %d independent cohorts. Accuracy on
                 held-out cohorts: %.1f%%, against %.1f%% for the training-fold
                 majority class. Target: %s.%s%s",
                v$n_datasets, 100 * v$accuracy, 100 * v$baseline, v$target,
                excl, reject_note))
  })

  query_fp <- reactive({
    req(input$query)
    r <- req(ref())
    # Unit compatibility is checked BEFORE any fingerprint is built, exactly as
    # the CLI does. A TPM query against a CPM reference is not a warning case:
    # length normalisation changes the spectrum, so the comparison is void.
    unit_error <- tryCatch({ assert_unit_compatible(r, input$unit); NULL },
                          error = function(e) conditionMessage(e))
    if (!is.null(unit_error)) return(list(unit_error = unit_error))

    q <- read_tsv_tsf(input$query$datapath)
    ids <- sub("\\..*$", "", as.character(q[[1]]))
    value_cols <- colnames(q)[-1]

    out <- lapply(value_cols, function(col) {
      # fingerprint_query builds the observed positions from the genes this file
      # actually contains. Genes it lacks are absent, never zero.
      fq <- fingerprint_query(q[[col]], ids, r, unit = input$unit)
      if (is.null(fq)) return(NULL)
      proj <- project_to_reference(fq$vector, r)
      # Coverage is checked BEFORE scoring. Computing a classification and then
      # warning about it invites the reader to take the number anyway.
      if (fq$coverage < 0.2 || proj$feature_coverage < 0.5) {
        return(list(name = col, result = NULL, coverage = fq$coverage,
                    feature_coverage = proj$feature_coverage,
                    id_type = fq$id_type, collapsed = fq$collapsed))
      }
      res <- match_query(r$model, proj$vector, proj$available)
      if (is.null(res)) return(NULL)
      res <- apply_rejection(res, r$validation$calibration,
                             coverage = proj$feature_coverage)
      list(name = col, result = res, coverage = fq$coverage,
           feature_coverage = proj$feature_coverage, id_type = fq$id_type,
           collapsed = fq$collapsed, n_features_used = res$n_features_used)
    })
    out[!vapply(out, is.null, logical(1))]
  })

  output$result <- renderUI({
    if (is.null(input$query)) return(NULL)
    r <- req(ref())
    matches <- query_fp()

    if (!is.null(matches$unit_error)) {
      return(div(class = "banner bad",
                 strong("Incompatible expression unit."), br(), matches$unit_error))
    }
    if (!length(matches)) {
      return(div(class = "banner bad",
                 strong("Nothing scorable in this file."), br(),
                 "No gene identifier matched the reference grid. Ensembl or ",
                 "Entrez ids are expected, matching the annotation the ",
                 "reference was built on."))
    }

    blocks <- lapply(matches, function(m) {
      if (is.null(m$result)) {
        return(tagList(h4(m$name),
          div(class = "banner bad",
              strong("LOW COVERAGE — not scored."), br(),
              sprintf("%.1f%% of the reference grid and %.1f%% of the features
                       the model uses are present in this file. Below 20%% and
                       50%% respectively, a similarity would be computed over so
                       little of the spectrum that it would not mean anything.",
                      100 * m$coverage, 100 * m$feature_coverage))))
      }
      res <- m$result
      s <- res$scores
      if (!isTRUE(input$show_all)) s <- utils::head(s, 5)

      # Decisions that mean "not classified" suppress the class table entirely,
      # as the CLI does. Showing the top class next to a warning invites the
      # reader to take the number anyway.
      not_classified <- res$decision %in% c("UNCALIBRATED_COVERAGE", "LOW_COVERAGE")
      if (identical(res$decision, "UNCALIBRATED_COVERAGE")) {
        return(tagList(h4(m$name),
          div(class = "banner bad",
              strong("NOT CLASSIFIED — no threshold for this coverage."), br(),
              sprintf("This query falls in the %s coverage band, and no
                       rejection threshold was calibrated for it. The
                       full-coverage threshold does not apply: the similarity
                       distribution shifts as coverage falls, so reusing it
                       would be the wrong number rather than a rough one.",
                      res$coverage_band %||% "unknown"))))
      }
      if (identical(res$decision, "LOW_COVERAGE")) {
        return(tagList(h4(m$name),
          div(class = "banner bad",
              strong("NOT CLASSIFIED — coverage below 50%."), br(),
              sprintf("Gene coverage %.1f%%. Below half the grid the similarity
                       is computed over too little of the spectrum to mean
                       anything.", 100 * m$coverage))))
      }

      verdict <- NULL
      if (identical(res$decision, "UNKNOWN")) {
        verdict <- div(class = "banner bad",
                       strong("UNKNOWN — outside the domain of this reference."),
                       br(),
                       sprintf("Closest class %s at similarity %.3f, below the
                                %.3f calibrated for it.", res$best, res$similarity,
                               res$threshold %||% NA_real_))
      } else if (identical(res$decision, "UNCALIBRATED")) {
        verdict <- div(class = "banner weak",
                       "No rejection threshold could be calibrated, so nothing ",
                       "rules out that this sample belongs to no class at all.")
      } else if (res$p_shuffle > 0.05) {
        verdict <- div(class = "banner bad",
                       "A randomly shuffled copy of this profile scores as well. ",
                       "There is no usable spectral shape here.")
      } else if (is.na(res$margin) || res$margin < 0.02) {
        verdict <- div(class = "banner weak",
                       "The top two classes are within 0.02. This call is not ",
                       "separable.")
      }

      tagList(
        h4(m$name), verdict,
        tags$table(
          tags$tr(tags$th("Class"), tags$th("Similarity")),
          lapply(seq_len(nrow(s)), function(i) {
            tags$tr(class = if (i == 1 && !identical(res$decision, "UNKNOWN")) "top" else "",
                    tags$td(s$class[i]), tags$td(sprintf("%.3f", s$similarity[i])))
          })
        ),
        p(class = "muted",
          sprintf("Margin %.3f | p against a shuffled query %.3f | grid coverage
                   %.1f%% (%s ids) | %d model features used (%.0f%%)%s",
                  res$margin, res$p_shuffle, 100 * m$coverage, m$id_type,
                  m$n_features_used, 100 * m$feature_coverage,
                  if (m$collapsed > 0)
                    sprintf(" | %d duplicate row(s) summed", m$collapsed) else ""))
      )
    })
    tagList(blocks)
  })
}

shinyApp(ui, server)
