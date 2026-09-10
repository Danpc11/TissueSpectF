#!/usr/bin/env Rscript
# shared_gene_mask.R -- la interseccion de genes que TODAS las cohortes midieron.
#
#   Rscript scripts/shared_gene_mask.R --interim-dir interim_bp250 \
#     --datasets GSE130970,GSE135251,GSE162694,GSE276114 \
#     --out interim_bp250/shared_gene_mask.tsv
#
# POR QUE HACE FALTA
# ------------------
# Con --grid-axis bp, ingest agrega varios genes por bin. Pero decide QUE genes
# entran con filter_expressed(), que usa rowMeans sobre las muestras de esa
# cohorte: cada una retiene un conjunto distinto, asi que el bin
# chr1:10.0-10.1 Mb puede calcularse con tres genes en una cohorte y con uno en
# otra. La posicion es comun y el valor no.
#
# Y una consulta individual no puede reproducir ese filtro: no conoce la
# distribucion de expresion de la cohorte. Aplicar la interseccion SOLO a la
# consulta es peor que no aplicarla, porque entonces hay tres
# representaciones: A, B, y A n B.
#
# EL FLUJO DE DOS PASADAS
# -----------------------
#   1. ./tsf ingest ... --grid-axis bp            (sin mascara)
#      escribe retained_genes.tsv por cohorte, con los ids REALES
#   2. Rscript scripts/shared_gene_mask.R ...     (este script)
#   3. ./tsf ingest ... --gene-mask <archivo>     (re-ingesta con la mascara)
#      y desde ahi cohortes y consultas comparten el conjunto
#
# El paso 3 es una re-ingesta completa, no un parche: los bins hay que volver a
# agregarlos, y con ellos los espectros, la validacion y los centroides.
#
# POR QUE genes.tsv NO SIRVE
# --------------------------
# Con eje bp, genes.tsv guarda `bin_<chr>_<index>` como gene_id, porque el bin
# es la posicion del eje. Leer la mascara de ahi da ids de BIN, y compararlos
# contra ids de gen no interseca nunca -- 0 de 15 en el caso medido, y la
# consulta abortaba en la primera referencia real de dos cohortes.

suppressWarnings({
  source("R/utils_io.R")
  source("R/config.R")
})

args <- commandArgs(trailingOnly = TRUE)
flag <- function(n, d = NULL) {
  h <- grep(paste0("^", n, "="), args, value = TRUE)
  if (length(h)) return(sub(paste0("^", n, "="), "", h[1]))
  i <- match(n, args)
  if (!is.na(i) && length(args) > i) args[i + 1] else d
}

interim <- flag("--interim-dir", Sys.getenv("TSF_INTERIM_DIR", ""))
datasets <- strsplit(flag("--datasets", ""), ",")[[1]]
out <- flag("--out", "")
mode <- flag("--mode", "intersect")

if (!nzchar(interim)) tsf_abort("Pasa --interim-dir <dir>.")
if (!length(datasets)) tsf_abort("Pasa --datasets A,B,C.")
if (!nzchar(out)) out <- file.path(interim, "shared_gene_mask.tsv")
if (!mode %in% c("intersect", "union")) {
  tsf_abort("--mode debe ser 'intersect' o 'union', no '", mode, "'")
}

sets <- list()
for (ds in datasets) {
  f <- file.path(interim, ds, "retained_genes.tsv")
  if (!file.exists(f)) {
    tsf_abort("No hay ", f, ". Corre primero `./tsf ingest ... --grid-axis bp` ",
              "sin mascara: ese paso lo escribe. Y no uses genes.tsv, que con ",
              "eje bp guarda ids de bin.")
  }
  t <- read_tsv_tsf(f, required = TRUE)
  if (!"gene_id" %in% names(t)) {
    tsf_abort(f, " no tiene columna gene_id.")
  }
  sets[[ds]] <- unique(as.character(t$gene_id))
  tsf_log(ds, ": ", length(sets[[ds]]), " gen(es) retenido(s)")
}

if (length(sets) < 2L) {
  tsf_warn("Una sola cohorte: la mascara es su propio conjunto y no resuelve ",
           "nada. El problema aparece al comparar cohortes.")
}

shared <- if (identical(mode, "intersect")) {
  Reduce(intersect, sets)
} else {
  Reduce(union, sets)
}
un <- length(Reduce(union, sets))

tsf_log("")
tsf_log(mode, ": ", length(shared), " de ", un, " gen(es) en la union (",
        round(100 * length(shared) / max(un, 1)), "%)")
if (identical(mode, "intersect")) {
  # Lo que cuesta, dicho: cada cohorte pierde los genes que otra no midio.
  for (ds in names(sets)) {
    lost <- length(setdiff(sets[[ds]], shared))
    tsf_log("  ", ds, " pierde ", lost, " de ", length(sets[[ds]]),
            " (", round(100 * lost / max(length(sets[[ds]]), 1)), "%)")
  }
  if (length(shared) < 0.5 * un) {
    tsf_warn("La interseccion es menos de la mitad de la union. Con cohortes ",
             "de plataformas muy distintas eso descarta mucha senal; ",
             "considera --mode union, que conserva los genes y deja que ",
             "bin_min_coverage decida por muestra cuales bins son usables.")
  }
}

ref <- read_tsv_tsf(file.path(interim, datasets[1], "retained_genes.tsv"))
res <- ref[match(shared, as.character(ref$gene_id)), , drop = FALSE]
res$dataset_id <- NULL
write_tsv_tsf(res, out)
tsf_log("Escrito ", out)
tsf_log("")
tsf_log("Siguiente paso -- RE-INGESTA, no un parche: los bins hay que volver a ")
tsf_log("agregarlos, y con ellos los espectros, la validacion y los centroides.")
tsf_log("  ./tsf ingest <datasets> --grid-axis bp --gene-mask ", out, " --force")
