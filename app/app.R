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

for (f in c("utils_io", "config", "labels", "ingest", "paths", "grid",
            "spectrum", "fingerprint", "reference")) {
  source(file.path("..", "R", paste0(f, ".R")))
}

project <- load_project_config(file.path("..", "config", "project.R"))
ref_path <- file.path(project$results_dir, "reference", "reference.rds")

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

  uiOutput("reference_banner"),

  hr(),
  fileInput("query", "Counts file (TSV: gene id column, then one column per sample)",
            accept = c(".tsv", ".txt", ".csv", ".gz")),
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
    if (!file.exists(ref_path)) return(NULL)
    readRDS(ref_path)
  })

  grid_ctx <- reactive({
    ids <- sub("\\.R$", "", list.files(file.path("..", "config", "datasets"),
                                       pattern = "\\.R$"))
    inp <- tsf_stage_inputs(project, ids[1])
    list(chrom_idx = inp$chrom_idx, genes = inp$dataset$genes,
         terms = fingerprint_terms(inp$chrom_idx))
  })

  output$reference_banner <- renderUI({
    r <- ref()
    if (is.null(r)) {
      return(div(class = "banner bad",
                 strong("No reference library."), br(),
                 "Build one first: ", code("./tsf run --to=spectra"), " then ",
                 code("./tsf reference"), "."))
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
    div(class = paste("banner", cls),
        strong(head), br(),
        sprintf("Trained and tested across %d independent cohorts. Accuracy on
                 held-out cohorts: %.1f%%, against %.1f%% for always guessing the
                 commonest class. Target: %s.",
                v$n_datasets, 100 * v$accuracy, 100 * v$baseline, v$target))
  })

  query_fp <- reactive({
    req(input$query)
    r <- req(ref())
    ctx <- grid_ctx()
    q <- read_tsv_tsf(input$query$datapath)
    ids <- sub("\\..*$", "", as.character(q[[1]]))
    hit <- match(ctx$genes$gene_id, ids)
    if (sum(!is.na(hit)) < 0.2 * nrow(ctx$genes) &&
        "entrez_id" %in% colnames(ctx$genes)) {
      hit <- match(ctx$genes$entrez_id, ids)
    }
    coverage <- mean(!is.na(hit))
    value_cols <- colnames(q)[-1]

    out <- lapply(value_cols, function(col) {
      counts <- suppressWarnings(as.numeric(q[[col]]))[hit]
      counts[!is.finite(counts)] <- 0
      y <- asinh(counts / max(sum(counts, na.rm = TRUE), 1) * 1e6)
      fp <- fingerprint_vector(y, ctx$chrom_idx, ctx$terms,
                               k_max = project$fingerprint$k_max %||% 64L,
                               features = project$fingerprint$features %||% "amplitude")
      if (is.null(fp)) return(NULL)
      fp <- (fp - mean(fp)) / max(stats::sd(fp), .Machine$double.eps)
      vv <- stats::setNames(rep(0, length(r$feature_space)), r$feature_space)
      shared <- intersect(names(fp), r$feature_space)
      vv[shared] <- fp[shared]
      list(name = col, result = match_query(r$model, vv), n_shared = length(shared))
    })
    list(coverage = coverage, matches = out[!vapply(out, is.null, logical(1))],
         n_features = length(r$feature_space))
  })

  output$result <- renderUI({
    if (is.null(input$query)) return(NULL)
    r <- req(ref())
    qf <- query_fp()

    if (qf$coverage < 0.2) {
      return(div(class = "banner bad",
                 strong("Identifiers not recognised."), br(),
                 sprintf("Only %.1f%% of the reference grid was found in this
                          file. Check that the gene ids are Ensembl or Entrez.",
                         100 * qf$coverage)))
    }

    blocks <- lapply(qf$matches, function(m) {
      s <- m$result$scores
      if (!isTRUE(input$show_all)) s <- utils::head(s, 5)
      warn <- NULL
      if (is.na(m$result$margin) || m$result$margin < 0.02) {
        warn <- div(class = "banner weak",
                    "The top two classes are within 0.02 of each other. ",
                    "This call is not separable.")
      }
      if (m$result$p_shuffle > 0.05) {
        warn <- div(class = "banner bad",
                    "A randomly shuffled version of this same profile scores as ",
                    "well. There is no usable spectral shape here.")
      }
      tagList(
        h4(m$name),
        tags$table(
          tags$tr(tags$th("Class"), tags$th("Similarity")),
          lapply(seq_len(nrow(s)), function(i) {
            tags$tr(class = if (i == 1) "top" else "",
                    tags$td(s$class[i]), tags$td(sprintf("%.3f", s$similarity[i])))
          })
        ),
        p(class = "muted",
          sprintf("Margin over runner-up: %.3f | p against a shuffled query: %.3f
                   | %d of %d reference features present",
                  m$result$margin, m$result$p_shuffle, m$n_shared, qf$n_features)),
        warn
      )
    })

    tagList(
      p(class = "muted",
        sprintf("Grid coverage of this file: %.1f%%", 100 * qf$coverage)),
      blocks
    )
  })
}

shinyApp(ui, server)
