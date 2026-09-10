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
# --mode se acepta pero solo "intersect" es coherente. Ver la nota de abajo.
mode <- flag("--mode", "intersect")

if (!nzchar(interim)) tsf_abort("Pasa --interim-dir <dir>.")
if (!length(datasets)) tsf_abort("Pasa --datasets A,B,C.")
if (!nzchar(out)) out <- file.path(interim, "shared_gene_mask.tsv")
if (!identical(mode, "intersect")) {
  # --mode union NO FUNCIONA, y recomendarlo era un error.
  #
  # La mascara solo puede QUITAR genes: en la segunda ingesta cada cohorte
  # conserva los de la union que realmente tiene, asi que los conjuntos siguen
  # siendo distintos y stage_reference() aborta igual al comparar interseccion
  # contra union. La union no aporta nada que la cohorte no tuviera ya.
  #
  # Si la interseccion descarta demasiado, las salidas reales son otras: subir
  # min_tpm o min_fraction para que el filtro por cohorte sea menos dependiente
  # de la profundidad, ampliar el bin para que cada uno tenga mas genes
  # anotados y sobreviva a perder algunos, o quitar la cohorte que mas
  # discrepa. Ninguna se resuelve con la union.
  tsf_abort("--mode '", mode, "' no esta soportado. Solo 'intersect' produce ",
            "conjuntos identicos entre cohortes, que es lo que el eje bp ",
            "necesita: la mascara solo puede quitar genes, asi que una union ",
            "deja a cada cohorte con los que ya tenia y los conjuntos siguen ",
            "difiriendo.")
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

shared <- Reduce(intersect, sets)
un <- length(Reduce(union, sets))

tsf_log("")
tsf_log("interseccion: ", length(shared), " de ", un, " gen(es) en la union (",
        round(100 * length(shared) / max(un, 1)), "%)")
# Lo que cuesta, dicho: cada cohorte pierde los genes que otra no midio.
for (ds in names(sets)) {
  lost <- length(setdiff(sets[[ds]], shared))
  tsf_log("  ", ds, " pierde ", lost, " de ", length(sets[[ds]]),
          " (", round(100 * lost / max(length(sets[[ds]]), 1)), "%)")
}
if (length(shared) < 0.5 * un) {
  tsf_warn("La interseccion es menos de la mitad de la union: se descarta ",
           "mucha senal. La union NO es la salida --la mascara solo puede ",
           "quitar genes, asi que cada cohorte se quedaria con los que ya ",
           "tenia y los conjuntos seguirian difiriendo. Las salidas reales: ",
           "bajar min_tpm o min_fraction para que el filtro por cohorte ",
           "dependa menos de la profundidad, ampliar --bin-size para que cada ",
           "bin tenga mas genes anotados y sobreviva a perder algunos, o ",
           "quitar la cohorte que mas discrepa.")
}
if (length(shared) == 0L) {
  tsf_abort("La interseccion esta vacia: ninguna posicion la miden todas las ",
            "cohortes. Con eje bp no hay malla comun posible; revisa los ",
            "filtros de expresion y la compatibilidad de plataformas.")
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
