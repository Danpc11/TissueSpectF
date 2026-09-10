# background.R -- fondo 1/f por muestra, y picos probados contra él.
#
# EL PROBLEMA
# -----------
# La expresión a lo largo de un cromosoma no es ruido blanco: los dominios de
# coexpresión le dan estructura tipo 1/f, con más potencia en los períodos
# largos. Cualquier prueba que suponga un fondo plano cuenta esa pendiente como
# señal.
#
# Medido: la g de Fisher, cuyo nulo analítico supone ruido blanco gaussiano en
# rejilla completa, declara periodicidad significativa en el 100% de las
# réplicas de puro ruido 1/f. No es un problema de calibración, es que el
# fondo del test no es el fondo de los datos.
#
# La permutación completa tiene la misma debilidad en otra forma: barajar las
# posiciones destruye TODA la autocorrelación, incluida la pendiente 1/f, así
# que el nulo sale más plano que la realidad y los picos parecen más fuertes de
# lo que son. Los esquemas de bloque lo mitigan, y por eso `--primary-scheme
# all` importa, pero no lo eliminan.
#
# LO QUE HACE ESTA RUTA
# ---------------------
# Estima el fondo DE LA PROPIA MUESTRA como una función suave del período y
# prueba cada frecuencia contra él. Nada viene de fuera, así que no hay varianza
# externa importada -- que es el defecto de restar un fondo estimado de otros
# datos: var(residuo) = var(señal) + var(estimación), y tratar el fondo como
# exacto es anticonservador.
#
# LO QUE NO HACE
# --------------
# No corrige la ventana espectral. La ventana es un invariante técnico por
# dataset --depende de qué posiciones se observaron-- y un fondo suave en
# período no puede representar sus lóbulos, que son picudos. La etapa `window`
# sigue siendo el control para eso.

#' Fondo espectral suave, ajustado en log-log.
#'
#' Un espectro de ruido rojo va como potencia ~ período^alpha, que en log-log es
#' una recta. Se ajusta con regresión robusta iterativa para que un pico real no
#' arrastre el fondo hacia arriba y se enmascare a sí mismo.
#'
#' @param period períodos, en las unidades del eje
#' @param power potencia por frecuencia
#' @param span 0 para una recta global; >0 para un ajuste local de ese ancho
#'   en décadas de log-período, que sigue una pendiente que cambia con la escala
#' @param iters iteraciones de la reponderación robusta
#' @return vector del fondo ajustado, en las unidades de `power`
spectral_background <- function(period, power, span = 0.5, iters = 4L) {
  ok <- is.finite(period) & is.finite(power) & period > 0 & power > 0
  if (sum(ok) < 12L) return(rep(NA_real_, length(power)))
  lp <- log(period[ok]); ly <- log(power[ok])

  fit_at <- function(w) {
    if (span <= 0) {
      # Una recta global: el caso 1/f puro.
      m <- stats::lm.wfit(cbind(1, lp), ly, w)
      as.numeric(cbind(1, lp) %*% m$coefficients)
    } else {
      # Local: para cada punto, una recta ponderada por distancia en
      # log-período. Sigue una pendiente que cambia entre escalas, que es lo
      # que un espectro real hace.
      vapply(seq_along(lp), function(i) {
        d <- abs(lp - lp[i]) / span
        k <- pmax(0, 1 - d^3)^3          # tricúbica
        ww <- w * k
        if (sum(ww > 0) < 4L) return(ly[i])
        m <- stats::lm.wfit(cbind(1, lp), ly, ww)
        sum(m$coefficients * c(1, lp[i]))
      }, numeric(1))
    }
  }

  w <- rep(1, length(lp))
  bg <- fit_at(w)
  for (it in seq_len(iters)) {
    r <- ly - bg
    # Reponderación de un lado solamente: un punto MUY POR ENCIMA del fondo es
    # un pico candidato y debe dejar de tirar del ajuste, pero uno muy por
    # debajo es un hueco del espectro y sí informa sobre el fondo. Bajar el
    # peso de los dos por igual sesgaría el fondo hacia arriba justo donde hay
    # señal, que es el error clásico de restar una línea base.
    s <- stats::mad(r, center = 0)
    if (!is.finite(s) || s <= 0) break
    w <- ifelse(r > 0, 1 / (1 + (r / (3 * s))^2), 1)
    bg <- fit_at(w)
  }
  out <- rep(NA_real_, length(power))
  out[ok] <- exp(bg)
  out
}

#' Exceso sobre el fondo, y su p-valor.
#'
#' @param sp salida de tsf_spectrum(): necesita `period` y `power`
#' @param span pasa a spectral_background()
#' @param dof grados de libertad del estimador en cada frecuencia. Un
#'   periodograma tiene 2; el multitaper con K tapers tiene aproximadamente
#'   2K, y de ahí sale su menor varianza. Poner 2 con multitaper daría
#'   p-valores demasiado grandes: el test sería válido pero ciego.
#' @return `sp` con background, excess, p_background y q_background
peaks_over_background <- function(sp, span = 0.5, dof = 2L) {
  sp$background <- spectral_background(sp$period, sp$power, span = span)
  sp$excess <- sp$power / sp$background
  # Bajo el fondo, 2*dof*I/S sigue una chi cuadrada con 2*dof grados de
  # libertad: el resultado estándar para un periodograma contra un espectro
  # conocido. El fondo aquí está ESTIMADO, no conocido, así que el p es
  # ligeramente anticonservador -- de los ~1000 puntos que lo ajustan, la
  # incertidumbre en cada frecuencia es pequeña pero no cero, y eso va dicho.
  ok <- is.finite(sp$excess) & sp$excess > 0

  # NULO EMPÍRICO, NO ANALÍTICO.
  #
  # La version anterior usaba chi cuadrada con 2*dof grados de libertad, el
  # resultado estandar de un periodograma contra un espectro conocido. Sobre
  # ruido 1/f al 26% de cobertura daba un 96-99% de falsos positivos: apenas
  # mejor que el 100% de la g de Fisher, y por las mismas dos razones -- las
  # ordenadas de una rejilla incompleta no son chi cuadradas independientes, y
  # el fondo se estima de los mismos datos que se prueban.
  #
  # El fondo sigue sirviendo, pero para BLANQUEAR: tras dividir por el, el
  # exceso deberia distribuirse igual en todas las frecuencias, y ese es el
  # punto. Asi que el nulo se toma de los propios excesos -- una frecuencia es
  # un pico si su exceso es extremo respecto al de las demas frecuencias de su
  # cromosoma. No supone forma alguna y se calibra solo.
  #
  # `dof` ya no entra en el p-valor: se conserva en la salida porque describe
  # el estimador y sirve para interpretar la magnitud del exceso.
  sp$p_background <- NA_real_
  if (sum(ok) >= 24L) {
    e <- sp$excess[ok]
    # (1 + #{e_j >= e_i}) / (n + 1): la misma forma que un p de permutacion, y
    # con el mismo piso de 1/(n+1). Con ~250 frecuencias el piso es 0.004, asi
    # que BH sobre ellas alcanza q = 0.05 solo si varias empatan en el piso --
    # exactamente la restriccion que el resto del pipeline ya reporta.
    r <- rank(-e, ties.method = "min")
    sp$p_background[ok] <- r / (length(e) + 1)
  }
  # NO SE APLICA BH SOBRE p_background.
  #
  # El p empírico se calcula por RANGO entre las m frecuencias del cromosoma,
  # así que ya es una afirmación familiar: "esta es la mas extrema de las m".
  # Corregirlo otra vez con BH lo multiplica por m, y como su piso es 1/(m+1)
  # el resultado en rango 1 es q ~ 1 para cualquier señal, por fuerte que sea.
  # Medido: con una amplitud de 8 sobre ruido 1/f el pico tenia exceso 37,
  # rango 1 de 255, y q = 0.99. La potencia del test era 0% -- no conservador,
  # ciego.
  #
  # `p_background` se reporta tal cual y es el valor a usar. `q_background` se
  # conserva sólo para no romper lo que ya lo lee, y vale lo mismo.
  sp$q_background <- sp$p_background
  sp$dof <- dof
  sp
}

#' Grados de libertad implícitos en un estimador.
#'
#' Separado para que peaks_over_background() no tenga que adivinar: pasar 2 con
#' un espectro multitaper convierte un test válido en uno ciego.
estimator_dof <- function(estimator) {
  est <- (estimator %||% list())$estimator %||% "periodogram"
  if (identical(est, "multitaper")) {
    as.integer((estimator %||% list())$mt_k %||% 5L)
  } else {
    1L
  }
}

#' Nulo del exceso, agrupado entre muestras.
#'
#' POR QUÉ HACE FALTA
#' ------------------
#' El nulo empírico dentro de un solo espectro es válido y a la vez impotente:
#' su p mínimo es 1/(m+1), así que con m = 250 frecuencias el piso es 0.004 y BH
#' en rango 1 exige m*p_min = 1.0. Medido: 0% de falsos positivos y 0% de
#' detección de una señal real. Un test que nunca rechaza está calibrado y no
#' sirve.
#'
#' Agrupando el exceso de S muestras el piso pasa a 1/(S*m+1) -- con 30 muestras
#' y 250 frecuencias, 1.3e-4 -- y BH vuelve a ser alcanzable.
#'
#' LO QUE ASUME, Y ES FUERTE
#' -------------------------
#' Que el exceso blanqueado es intercambiable entre muestras: que tras dividir
#' por su propio fondo, dos muestras de la misma condición tienen la misma
#' distribución de excesos. Si una muestra tiene mucho más ruido que las demás,
#' domina el nulo agrupado y encoge los p-valores de todas.
#'
#' Es la misma suposición que hace `p_pointwise` en R/maxt.R, y ahí ya costó un
#' error: con un solo sorteo la desviación estándar era indefinida y el nulo
#' agrupado fabricó 4.950 descubrimientos. De ahí el mínimo de muestras.
#'
#' @param specs lista de salidas de peaks_over_background(), una por muestra
#' @param min_samples por debajo de esto devuelve NA en vez de un p apretado
#'   con un nulo que no puede sostenerlo
pooled_background_null <- function(specs, min_samples = 8L) {
  specs <- Filter(function(s) !is.null(s) && "excess" %in% names(s), specs)
  if (length(specs) < min_samples) {
    tsf_warn("pooled_background_null: ", length(specs), " muestra(s), hacen ",
             "falta ", min_samples, ". Se devuelve NA en vez de un p-valor ",
             "que el nulo no puede sostener.")
    return(lapply(specs, function(s) {
      s$p_pooled <- NA_real_; s$q_pooled <- NA_real_; s
    }))
  }
  pool <- unlist(lapply(specs, function(s) s$excess[is.finite(s$excess)]),
                 use.names = FALSE)
  n <- length(pool)
  pool <- sort(pool)
  lapply(specs, function(s) {
    e <- s$excess
    ok <- is.finite(e)
    s$p_pooled <- NA_real_
    # findInterval sobre el pool ordenado: cuántos del nulo agrupado igualan o
    # superan el exceso observado.
    ge <- n - findInterval(e[ok] - .Machine$double.eps, pool)
    s$p_pooled[ok] <- (1 + ge) / (n + 1)
    s$q_pooled <- NA_real_
    if (any(ok)) s$q_pooled[ok] <- stats::p.adjust(s$p_pooled[ok], "BH")
    s$n_pooled <- n
    s
  })
}

#' Marcar sólo las frecuencias que están por encima de su fondo.
#'
#' Añade `over_background`, lógico, que es la columna a filtrar. Nada se
#' elimina: una frecuencia descartada sigue en la tabla con su exceso y su p,
#' porque saber que una escala NO destaca es tan informativo como lo contrario
#' y borrarla haría imposible reconstruir el espectro completo.
#'
#' `alpha` se aplica a `p_background`, que ya es familiar por rango entre las m
#' frecuencias del cromosoma. NO se le vuelve a aplicar BH: el p se calcula
#' sobre las mismas m que BH corregiría, así que en rango 1 el resultado seria
#' q ~ 1 para cualquier señal -- medido, con amplitud 8 el pico tenia exceso 37
#' y q = 0.99, y la potencia del test era 0%.
#'
#' `min_excess` es un piso de tamaño de efecto además del p. Con ~250
#' frecuencias por cromosoma el p mas pequeno alcanzable es 1/251, asi que en
#' un cromosoma sin senal el rango 1 siempre tiene p = 0.004 y cruza cualquier
#' alpha razonable. El exceso dice CUANTO destaca, y sin un piso el marcado
#' devolveria una frecuencia por cromosoma por construccion.
mark_over_background <- function(sp, alpha = 0.05, min_excess = NULL,
                                 excess_null = NULL, verbose = FALSE) {
  if (!"p_background" %in% names(sp)) {
    tsf_abort("mark_over_background: falta p_background; corre ",
              "peaks_over_background() primero")
  }

  # EL RANGO SOLO NO SIRVE, MEDIDO.
  #
  # p_background es un rango entre las m frecuencias del cromosoma, asi que
  # `p <= alpha` marca el alpha superior POR CONSTRUCCION: en ruido 1/f puro
  # marca 25 de 511 (4.9%) y con una senal real marca las mismas 25. El rango
  # no distingue porque siempre hay un 5% superior.
  #
  # Lo que distingue es el EXCESO. Sobre 1/f al 26% de cobertura su
  # distribucion es mediana 1.3, q95 5.2, q99 8.0, maximo 11.6; el pico real
  # de amplitud 6 tenia exceso 37. Se separan en el extremo, no en el rango.
  #
  # `excess_null` son los excesos de un nulo --de otras muestras, de
  # permutaciones, o de pooled_background_null()-- y de ahi sale el umbral. Sin
  # el se usa `min_excess`, y si tampoco se da se aborta en vez de marcar el
  # 5% superior de cualquier cosa.
  # EL UMBRAL ES FAMILIAR, NO POR FRECUENCIA.
  #
  # Un cuantil sobre TODOS los excesos del nulo no sirve, y por una razon
  # circular: p_background ES el rango del exceso, asi que `p <= alpha` y
  # `exceso >= q95(todos)` seleccionan el mismo conjunto -- el alpha superior,
  # con senal o sin ella. Medido: 25 de 511 en ruido puro y 25 de 511 con una
  # senal de amplitud 6.
  #
  # `excess_null` debe ser una LISTA, un vector de excesos por realizacion del
  # nulo, y el umbral es el cuantil del MAXIMO de cada una. Esa es la logica de
  # maxT: en ruido puro el 5% de los cromosomas tendra alguna frecuencia por
  # encima, no el 5% de las frecuencias.
  if (is.null(min_excess) && !is.null(excess_null)) {
    if (!is.list(excess_null)) {
      tsf_abort("mark_over_background: excess_null tiene que ser una LISTA de ",
                "vectores, uno por realizacion del nulo. Un vector plano da un ",
                "umbral por frecuencia, que es circular: p_background ya es el ",
                "rango del exceso, asi que seleccionaria el alpha superior de ",
                "cualquier espectro.")
    }
    mx <- vapply(excess_null, function(v) {
      v <- v[is.finite(v)]
      if (!length(v)) NA_real_ else max(v)
    }, numeric(1))
    mx <- mx[is.finite(mx)]
    if (length(mx) < 20L) {
      tsf_abort("mark_over_background: ", length(mx), " realizacion(es) del ",
                "nulo; hacen falta al menos 20 para el cuantil ", 1 - alpha, ".")
    }
    min_excess <- as.numeric(stats::quantile(mx, 1 - alpha))
    # `verbose` y no un tsf_log incondicional: la funcion se llama una vez por
    # (muestra, cromosoma), y en una corrida real eso son miles de lineas
    # identicas que entierran cualquier mensaje que si informe.
    if (isTRUE(verbose)) {
      tsf_log("  umbral familiar de exceso: ", signif(min_excess, 3),
              " (cuantil ", 1 - alpha, " del maximo de ", length(mx),
              " realizacion(es))")
    }
  }
  if (is.null(min_excess)) {
    tsf_abort("mark_over_background: se necesita `excess_null` para calibrar ",
              "el umbral, o `min_excess` explicito. Marcar solo por ",
              "p_background <= alpha selecciona el alpha superior de cualquier ",
              "espectro, con senal o sin ella.")
  }

  # Solo el exceso: anadir `p_background <= alpha` no filtra nada mas, porque
  # el p es el rango del mismo exceso.
  sp$over_background <- is.finite(sp$excess) & sp$excess >= min_excess
  attr(sp, "min_excess") <- min_excess
  sp
}

#' Cuántas frecuencias sobreviven, por cromosoma y por muestra.
#'
#' Para el log: una fraccion alta es la senal de alarma que ya nos ahorro dos
#' falsos hallazgos en este proyecto. Si el 29% de las frecuencias "destaca",
#' lo que destaca es el criterio.
background_summary <- function(sp) {
  if (!"over_background" %in% names(sp)) return(NULL)
  key <- if ("chr" %in% names(sp)) as.character(sp$chr) else "all"
  n <- tapply(sp$over_background, key, function(v) sum(v, na.rm = TRUE))
  m <- tapply(sp$over_background, key, length)
  data.frame(chr = names(n), n_over = as.integer(n), n_freq = as.integer(m),
             frac = round(as.numeric(n) / as.numeric(m), 4),
             stringsAsFactors = FALSE)
}
