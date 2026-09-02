# bundle.R -- package the app as a folder anyone can download and run.
#
#   ./tsf bundle --out TissueSpectF-app
#
# The result depends on nothing else: not on this repository, not on
# config/project.R, not on the interim or results directories, not on the GEO
# downloads. It carries the reference, the four R modules a query needs, a
# launcher, and a README. That is the whole point -- a collaborator should be
# able to receive a zip and identify a sample without installing a pipeline.
#
# Only `shiny` is fetched at first run, and only if it is missing.

BUNDLE_MODULES <- c("utils_io", "grid", "fingerprint", "reference")

bundle_app <- function(project, out_dir, reference_path = NULL) {
  ref_path <- reference_path %||%
    file.path(project$results_dir, "reference", "reference.rds")
  if (!file.exists(ref_path)) {
    tsf_abort("No reference at ", ref_path,
              ". Build one first: ./tsf reference")
  }
  ref <- readRDS(ref_path)
  if (is.null(ref$grid)) {
    tsf_abort("That reference carries no grid, so it cannot score a query on ",
              "its own. Rebuild it: ./tsf reference")
  }

  ensure_dir(out_dir)
  ensure_dir(file.path(out_dir, "R"))
  for (m in BUNDLE_MODULES) {
    file.copy(file.path("R", paste0(m, ".R")),
              file.path(out_dir, "R", paste0(m, ".R")), overwrite = TRUE)
  }
  file.copy(ref_path, file.path(out_dir, "reference.rds"), overwrite = TRUE)
  file.copy("app/app.R", file.path(out_dir, "app.R"), overwrite = TRUE)

  # A distributable with no statement of rights leaves the recipient guessing.
  # Everything in the bundle is MIT -- the four modules, the app and the
  # reference -- so the MIT text travels with it. The sonification is not in
  # the bundle, which is why LICENSE-ART.md is not needed here.
  if (file.exists("LICENSE")) {
    file.copy("LICENSE", file.path(out_dir, "LICENSE"), overwrite = TRUE)
  } else {
    tsf_warn("No LICENSE file to include in the bundle; recipients will have ",
             "no statement of their rights.")
  }

  writeLines(bundle_launcher_sh(), file.path(out_dir, "run.sh"))
  Sys.chmod(file.path(out_dir, "run.sh"), "755")
  writeLines(bundle_launcher_bat(), file.path(out_dir, "run.bat"))
  writeLines(bundle_readme(ref), file.path(out_dir, "README.txt"))
  writeLines(bundle_manifest(ref, ref_path), file.path(out_dir, "REFERENCE.txt"))

  zip_path <- paste0(out_dir, ".zip")
  ok <- tryCatch({
    utils::zip(zip_path, out_dir, flags = "-rq")
    file.exists(zip_path)
  }, error = function(e) FALSE, warning = function(w) FALSE)

  tsf_log("Bundle written to ", out_dir,
          if (ok) paste0(" and ", zip_path) else
            " (no zip: the `zip` command was not available)")
  tsf_log("It needs nothing from this repository. Send the folder, run run.sh.")
  invisible(list(dir = out_dir, zip = if (ok) zip_path else NA_character_))
}

bundle_launcher_sh <- function() c(
  "#!/usr/bin/env bash",
  "# Start the TissueSpectF app. Everything runs on this machine.",
  "set -euo pipefail",
  'cd "$(dirname "$0")"',
  "",
  "if ! command -v Rscript >/dev/null 2>&1; then",
  '  echo "R is not installed. Get it from https://cran.r-project.org and run this again." >&2',
  "  exit 1",
  "fi",
  "",
  "Rscript -e 'if (!requireNamespace(\"shiny\", quietly = TRUE)) {",
  "  message(\"Installing shiny (once)...\");",
  "  install.packages(\"shiny\", repos = \"https://cloud.r-project.org\")",
  "}'",
  "",
  'echo "Opening in your browser. Press Ctrl-C here to stop."',
  "Rscript -e 'shiny::runApp(\".\", launch.browser = TRUE)'")

bundle_launcher_bat <- function() c(
  "@echo off",
  "REM Start the TissueSpectF app. Everything runs on this machine.",
  "cd /d \"%~dp0\"",
  "where Rscript >nul 2>nul",
  "if errorlevel 1 (",
  "  echo R is not installed. Get it from https://cran.r-project.org",
  "  pause",
  "  exit /b 1",
  ")",
  "Rscript -e \"if (!requireNamespace('shiny', quietly=TRUE)) install.packages('shiny', repos='https://cloud.r-project.org')\"",
  "Rscript -e \"shiny::runApp('.', launch.browser=TRUE)\"",
  "pause")

bundle_readme <- function(ref) {
  v <- ref$validation
  status <- reference_status(ref)
  c("TissueSpectF - spectral tissue matcher",
    "======================================",
    "",
    "Run:",
    "  macOS / Linux:  ./run.sh",
    "  Windows:        double-click run.bat",
    "",
    "It opens in your browser. Nothing is uploaded; the whole thing runs on",
    "this computer. R is required (https://cran.r-project.org); shiny is",
    "installed automatically the first time.",
    "",
    "Licence: MIT. See LICENSE. Everything in this folder is free to use,",
    "modify and redistribute, commercially or not.",
    "",
    "What to feed it",
    "---------------",
    "A TSV: first column gene identifiers (Ensembl or Entrez), then one column",
    "per sample. Raw counts by default - tell the app if the values are TPM,",
    "CPM or already log-transformed.",
    "",
    "How much to trust a match",
    "-------------------------",
    status,
    "",
    "The app refuses to classify a sample whose gene coverage is below 50%,",
    "and reports UNKNOWN when the best similarity falls under the threshold",
    "calibrated for that coverage. A match is a similarity to a class centroid.",
    "It is not a diagnosis.",
    "",
    if (is.null(v)) "This reference was NOT validated across cohorts." else
      sprintf("Validated across %d independent cohorts on the target '%s'.",
              v$n_datasets, v$target),
    "",
    "See REFERENCE.txt for exactly what this reference was built from.")
}

bundle_manifest <- function(ref, ref_path) {
  p <- ref$params
  v <- ref$validation
  c("Reference provenance",
    "====================",
    paste("built           ", ref$built),
    paste("pipeline        ", ref$version),
    paste("source file     ", basename(ref_path)),
    paste("datasets        ", p$datasets %||% NA),
    paste("target          ", ref$target),
    paste("classes         ", paste(ref$model$classes, collapse = ", ")),
    paste("samples         ", nrow(ref$labels)),
    "",
    "Grid and scale",
    "--------------",
    paste("species         ", p$species %||% NA),
    paste("genome build    ", p$genome_build %||% NA),
    paste("annotation      ", p$annotation %||% NA, "/", p$annotation_release %||% NA),
    paste("gene universe   ", p$gene_universe %||% NA),
    paste("grid genes      ", nrow(ref$grid)),
    paste("grid digest     ", p$grid_digest %||% NA),
    paste("expression unit ", p$expression_unit %||% NA),
    paste("k_max           ", p$k_max),
    paste("features        ", p$features),
    "",
    "Validation",
    "----------",
    if (is.null(v)) "none - single cohort" else c(
      sprintf("out-of-cohort accuracy   %.1f%%", 100 * v$accuracy),
      sprintf("majority-class baseline  %.1f%%", 100 * v$baseline),
      if (!is.null(v$calibration$bands)) c(
        "", "per coverage band:",
        apply(v$calibration$bands, 1, function(r)
          sprintf("  %-9s accuracy %5s  threshold %6s  rejects %5s of members",
                  r[["band"]],
                  round(as.numeric(r[["accuracy"]]), 3),
                  round(as.numeric(r[["threshold_applied"]]), 3),
                  round(as.numeric(r[["expected_rejection_of_members"]]), 3))))
      else "no coverage-band calibration"))
}
