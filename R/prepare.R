# prepare.R -- la preparación gen → eje, compartida por ingest y por consultas.
#
# EL PROBLEMA
# -----------
# Las dos rutas hacían las mismas operaciones EN ORDEN DISTINTO:
#
#   ingest:   counts -> TPM -> asinh -> filtrar genes -> agregar por bin
#   consulta: counts -> agregar por bin -> CPM -> asinh
#
# Y no conmutan:
#
#   mean(asinh(TPM_i))  !=  asinh(TPM_del_bin)
#
# asinh es concavo para valores positivos, asi que promediar despues de
# transformar da menos que transformar despues de promediar, y la diferencia
# crece con la dispersion de los genes del bin. Un bin con un gen a 1000 TPM y
# otro a 1 da 3.8 por una ruta y 6.9 por la otra.
#
# Habia un segundo desfase, mas silencioso: ingest usa TPM con longitud de gen
# y la consulta usaba CPM sin ella. Dos genes con el mismo conteo y longitudes
# distintas reciben valores distintos en la referencia e iguales en la
# consulta.
#
# El resultado: aunque consulta y referencia usaran las mismas posiciones --lo
# que costo el arreglo anterior-- el VALOR de cada posicion se calculaba
# distinto, y la similitud comparaba dos cosas que no son la misma cantidad.
#
# LA SOLUCION
# -----------
# Una sola funcion, `prepare_axis_values()`, con el orden fijado en un solo
# sitio. Las dos rutas la llaman. Si el orden cambia, cambia para las dos o no
# cambia para ninguna.

#' Normalizar, transformar, enmascarar y agregar, en ese orden.
#'
#' @param counts matriz genes x muestras, o un vector para una consulta
#' @param gene_length longitudes en pb, o NULL para CPM. Pasarlo o no cambia la
#'   unidad, y la unidad viaja en el atributo para que nadie suponga TPM cuando
#'   fue CPM.
#' @param unit "counts" normaliza; "cpm" y "tpm" ya vienen normalizados y solo
#'   se transforman; "logged" no se toca.
#' @param bin_key clave de bin por gen, o NULL para no agregar. La agregacion
#'   va AL FINAL, despues de la transformacion, porque es el orden de ingest y
#'   el que hay que reproducir.
#' @param bin_aggregate "mean" o "sum"
#' @param bin_annotated genes ANOTADOS por bin, con nombres de bin. El
#'   denominador de la cobertura del bin: sin el la cobertura saldria 1 siempre
#'   porque solo se ven los genes que llegaron.
#' @param bin_min_coverage un bin por debajo queda NA en esa muestra
#' @param mask logico por gen, TRUE = usable. Se aplica ANTES de agregar, que
#'   es donde ingest lo aplica.
#'
#' @return matriz o vector con el atributo "unit", y cuando hay agregacion los
#'   atributos "bin_coverage" y "n_dropped"
prepare_axis_values <- function(counts, gene_length = NULL, unit = "counts",
                                bin_key = NULL, bin_aggregate = "mean",
                                bin_annotated = NULL, bin_min_coverage = 0.5,
                                mask = NULL) {
  vec_in <- !is.matrix(counts)
  m <- if (vec_in) matrix(counts, ncol = 1L,
                          dimnames = list(names(counts), "query")) else counts

  # 1. NORMALIZAR. La unidad se decide aqui y viaja en el atributo.
  # QUE ES NO MEDIDO SE DECIDE ANTES DE NORMALIZAR.
  #
  # `m[is.na(m)] <- 0` convertia los NA en ceros, y despues la cobertura del
  # bin los contaba como medidos: con 3 o 5 genes anotados y uno solo presente
  # la cobertura salia 1 y no se descartaba nada. Medido: 0 de 16 bins
  # descartados con umbral 0.9.
  #
  # La mascara de ausencia se guarda aqui y se aplica despues del asinh, para
  # que un cero de normalizacion --que es un conteo real de cero-- no se
  # confunda con un gen que la muestra no midio.
  absent <- !is.finite(m)

  if (identical(unit, "counts")) {
    m[is.na(m)] <- 0
    if (is.null(gene_length)) {
      scaled <- t(t(m) / pmax(colSums(m), 1)) * 1e6
      u <- "CPM"
    } else {
      if (length(gene_length) != nrow(m)) {
        tsf_abort("prepare_axis_values: ", length(gene_length),
                  " longitud(es) para ", nrow(m), " gen(es)")
      }
      rpk <- m / (gene_length / 1000)
      scaled <- t(t(rpk) / pmax(colSums(rpk, na.rm = TRUE), 1)) * 1e6
      u <- "TPM"
    }
  } else if (unit %in% c("cpm", "tpm")) {
    scaled <- m
    u <- toupper(unit)
  } else if (identical(unit, "logged")) {
    scaled <- m
    u <- "logged"
  } else {
    tsf_abort("prepare_axis_values: unidad '", unit,
              "' desconocida. Use counts, cpm, tpm o logged.")
  }

  # 2. TRANSFORMAR. Antes de agregar, porque es el orden de ingest: asinh es
  #    concavo y promediar despues de transformar no es transformar despues de
  #    promediar.
  y <- if (identical(u, "logged")) scaled else asinh(scaled)

  # 3. ENMASCARAR. Primero lo ausente, luego la mascara del llamador.
  y[absent] <- NA_real_
  if (!is.null(mask)) {
    if (length(mask) != nrow(y)) {
      tsf_abort("prepare_axis_values: mascara de ", length(mask),
                " para ", nrow(y), " gen(es)")
    }
    y[!mask, ] <- NA_real_
  }

  if (is.null(bin_key)) {
    out <- if (vec_in) as.numeric(y[, 1]) else y
    attr(out, "unit") <- if (identical(u, "logged")) u else paste0("asinh(", u, ")")
    return(out)
  }

  # 4. AGREGAR, al final.
  if (length(bin_key) != nrow(y)) {
    tsf_abort("prepare_axis_values: ", length(bin_key), " clave(s) de bin ",
              "para ", nrow(y), " gen(es)")
  }
  if (!bin_aggregate %in% c("mean", "sum")) {
    tsf_abort("prepare_axis_values: bin_aggregate debe ser 'mean' o 'sum', ",
              "no '", bin_aggregate, "'")
  }
  idx <- split(seq_len(nrow(y)), bin_key)
  bins <- names(idx)

  n_annot <- if (is.null(bin_annotated)) {
    vapply(idx, length, integer(1))
  } else {
    v <- as.integer(bin_annotated[bins])
    v[is.na(v)] <- vapply(idx, length, integer(1))[is.na(v)]
    v
  }

  ns <- ncol(y)
  agg <- matrix(NA_real_, length(idx), ns, dimnames = list(bins, colnames(y)))
  meas <- matrix(0L, length(idx), ns, dimnames = dimnames(agg))
  for (i in seq_along(idx)) {
    blk <- y[idx[[i]], , drop = FALSE]
    meas[i, ] <- colSums(is.finite(blk))
    v <- if (identical(bin_aggregate, "sum")) colSums(blk, na.rm = TRUE)
         else colMeans(blk, na.rm = TRUE)
    v[meas[i, ] == 0L] <- NA_real_
    agg[i, ] <- v
  }

  cov <- sweep(meas, 1L, pmax(n_annot, 1L), "/")
  drop <- !is.finite(cov) | cov < bin_min_coverage
  agg[drop] <- NA_real_

  out <- if (vec_in) stats::setNames(as.numeric(agg[, 1]), bins) else agg
  attr(out, "unit") <- if (identical(u, "logged")) u else paste0("asinh(", u, ")")
  attr(out, "bin_coverage") <- cov
  attr(out, "bin_measured") <- meas
  attr(out, "bin_annotated") <- n_annot
  attr(out, "n_dropped") <- sum(drop, na.rm = TRUE)
  out
}

#' LA definicion de cobertura, unica para calibracion y consulta.
#'
#' Habia TRES denominadores distintos en el repo:
#'
#'   calibracion  n_observed / (grid_size %||% n_observed)   -> 1 por defecto,
#'                porque grid_size nunca se pasaba: el log decia "covers X%" y
#'                calculaba 100% por construccion
#'   consulta v1  n_observed / sum(grid_N)                    -> los bins
#'                fisicos del cromosoma
#'   consulta v2  n_observed / nrow(grid)                     -> las posiciones
#'                de la referencia
#'
#' El umbral de rechazo se elegia con una banda calibrada con un denominador y
#' se aplicaba a una consulta medida con otro. Una consulta con TODAS las
#' posiciones de la referencia podia salir al 6% y rechazarse por <50%.
#'
#' La definicion es: **fraccion de las posiciones de la referencia que estan
#' observadas**. Es la cantidad que el umbral necesita --"dada esta fraccion de
#' la referencia, cuanta similitud hace falta"-- y la unica que significa lo
#' mismo en los dos ejes.
#'
#' La fraccion del CROMOSOMA es otra cantidad, honesta para decir cuanto del
#' genoma cubren los datos, y `genomic_coverage()` la calcula aparte. No debe
#' alimentar umbrales.
#'
#' @param n_observed posiciones observadas
#' @param grid la malla de la referencia
reference_coverage <- function(n_observed, grid) {
  n <- nrow(grid)
  if (!is.finite(n) || n <= 0) return(NA_real_)
  min(max(n_observed / n, 0), 1)
}

#' Fraccion del cromosoma observada. Se reporta, no elige umbrales.
genomic_coverage <- function(n_observed, grid) {
  if (!"grid_N" %in% names(grid)) return(NA_real_)
  Ns <- vapply(split(grid$grid_N, as.character(grid$chr)),
               function(v) as.numeric(v[1]), numeric(1))
  tot <- sum(Ns, na.rm = TRUE)
  if (!is.finite(tot) || tot <= 0) return(NA_real_)
  min(max(n_observed / tot, 0), 1)
}

#' Re-normalizar tras enmascarar, como haria una consulta real.
#'
#' La calibracion de cobertura quitaba posiciones DESPUES de normalizar; una
#' consulta real normaliza con los conteos que tiene. No es la misma
#' perturbacion, y la diferencia es grande: al 26% de cobertura la media en
#' escala asinh difiere en 1.53, que es log(1/0.26) = 1.35 como predice la
#' teoria --quitar una fraccion f de genes baja el total y escala los que
#' quedan por 1/f, y asinh(x/f) - asinh(x) ~ log(1/f) para x grande.
#'
#' El umbral se calibraba con una perturbacion y se aplicaba a otra.
#'
#' No hacen falta los conteos: asinh es invertible. sinh() devuelve el TPM (o
#' CPM), se enmascara, se re-normaliza a 1e6 sobre lo que queda, y se vuelve a
#' aplicar asinh. Eso reproduce exactamente lo que hace una consulta parcial.
#'
#' @param y valores en asinh(TPM) o asinh(CPM)
#' @param keep logico, TRUE = posicion observada por la consulta simulada
#' @param unit "logged" no se toca: sin normalizacion no hay nada que rehacer
renormalise_after_mask <- function(y, keep, unit = "asinh(TPM)") {
  if (identical(unit, "logged")) return(ifelse(keep, y, NA_real_))
  if (length(keep) != length(y)) {
    tsf_abort("renormalise_after_mask: ", length(keep), " indicador(es) para ",
              length(y), " valor(es)")
  }
  lin <- sinh(y)                       # de vuelta a TPM/CPM
  lin[!keep] <- NA_real_
  tot <- sum(lin, na.rm = TRUE)
  if (!is.finite(tot) || tot <= 0) return(rep(NA_real_, length(y)))
  asinh(lin / tot * 1e6)
}
