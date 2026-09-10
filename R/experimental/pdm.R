# pdm.R -- minimización de dispersión de fase, y GLS ponderado por varianza.
#
# DOS SUPUESTOS QUE EL PIPELINE HACE Y NO HABÍA VERIFICADO
# --------------------------------------------------------
# LA FORMA. Lomb-Scargle es óptimo para AJUSTAR UN SINUSOIDE, que no es lo
# mismo que ser óptimo para encontrar el periodo de una señal genérica. Si la
# organización cromosómica produce dominios con bordes abruptos --bloques de
# expresión alta y baja, que es lo que un compartimento A/B o un dominio de
# replicación parece-- entonces el estimador está buscando la forma equivocada
# y reparte la potencia del bloque entre su fundamental y sus armónicos.
#
# El PDM no supone forma: dobla las posiciones al periodo de prueba, agrupa por
# fase y compara la dispersión dentro de los grupos contra la total. Un periodo
# correcto agrupa valores parecidos sea cual sea la forma de la onda.
#
# EL PESO. Los datos son conteos, así que cada gen tiene una varianza
# estimable, y el GLS los trata todos igual. Un gen con 5 lecturas y otro con
# 5000 pesan lo mismo en el ajuste, cuando el segundo lleva mucha más
# información. Ponderar por 1/varianza es estándar en el Lomb-Scargle de
# astronomía, donde cada medición trae su barra de error.

#' Estadístico theta de minimización de dispersión de fase.
#'
#' theta = (varianza dentro de los grupos de fase) / (varianza total).
#' Cerca de 1 no hay estructura a ese periodo; muy por debajo de 1 la hay.
#' Es un MÍNIMO y no un máximo, al contrario que el periodograma.
#'
#' @param t posiciones observadas
#' @param y valores en esas posiciones
#' @param periods periodos a probar, en las mismas unidades que t
#' @param n_bins grupos de fase. Pocos pierden resolución de forma, muchos
#'   dejan grupos casi vacíos cuya varianza no se puede estimar. 10 es el
#'   compromiso habitual.
#' @param min_per_bin un grupo con menos de esto no aporta; su varianza sería
#'   ruido. Con cobertura del 26% y 10 grupos hacen falta unas 30 posiciones
#'   por cromosoma para que ninguno quede vacío.
pdm_theta <- function(t, y, periods, n_bins = 10L, min_per_bin = 3L) {
  t <- as.numeric(t); y <- as.numeric(y)
  ok <- is.finite(t) & is.finite(y)
  t <- t[ok]; y <- y[ok]
  n <- length(y)
  if (n < n_bins * min_per_bin) {
    return(rep(NA_real_, length(periods)))
  }
  # Varianza total con el mismo denominador que la de los grupos, para que
  # theta sea 1 cuando no hay estructura y no un valor sesgado por n.
  vtot <- sum((y - mean(y))^2) / (n - 1)
  if (!is.finite(vtot) || vtot <= 0) return(rep(NA_real_, length(periods)))

  vapply(periods, function(P) {
    if (!is.finite(P) || P <= 0) return(NA_real_)
    ph <- (t %% P) / P
    b <- pmin(floor(ph * n_bins) + 1L, n_bins)
    # Sumas por grupo en una pasada, sin split(): con 24 cromosomas y miles de
    # periodos el coste de crear listas domina.
    cnt <- tabulate(b, nbins = n_bins)
    if (any(cnt < min_per_bin)) return(NA_real_)
    s1 <- vapply(seq_len(n_bins), function(j) sum(y[b == j]), numeric(1))
    s2 <- vapply(seq_len(n_bins), function(j) sum(y[b == j]^2), numeric(1))
    # Suma de cuadrados dentro de cada grupo, agregada
    ssw <- sum(s2 - s1^2 / cnt)
    dfw <- n - n_bins
    if (dfw <= 0) return(NA_real_)
    (ssw / dfw) / vtot
  }, numeric(1))
}

#' Espectro PDM sobre la misma familia de periodos que el GLS.
#'
#' Devuelve las mismas columnas de identificación que gls_spectrum() más
#' `theta`, para que las dos rutas se puedan cruzar por (k, N).
#'
#' NO devuelve un p-valor. theta no tiene distribución nula conocida sobre una
#' rejilla incompleta con estructura 1/f, y las aproximaciones publicadas
#' suponen ruido blanco -- el mismo supuesto que hace fallar a la g de Fisher
#' con el 100% de falsos positivos en estos datos. El nulo tiene que venir de
#' permutación, y por eso esta función sólo entrega el estadístico.
pdm_spectrum <- function(y, terms, n_bins = 10L, min_per_bin = 3L) {
  y <- as.numeric(y)
  if (anyNA(y) || any(!is.finite(y))) {
    tsf_abort("pdm_spectrum() recibió ", sum(!is.finite(y)), " valor(es) no ",
              "finito(s). Pasá la señal por gls_observed() primero.")
  }
  N <- terms$N
  per <- N / terms$k
  th <- pdm_theta(terms$t_index, y, per, n_bins = n_bins,
                  min_per_bin = min_per_bin)
  m <- length(terms$k)
  structure(list(k = terms$k, freq = terms$k / N, period = per,
                 N = rep(N, m), n_observed = rep(terms$n, m),
                 theta = th,
                 # 1 - theta para que, como la potencia, más grande sea más
                 # estructura. Facilita reusar el código de picos y de nulos,
                 # que en todo el pipeline busca máximos.
                 pdm_strength = 1 - th),
            class = "data.frame", row.names = seq_len(m))
}

#' GLS ponderado por varianza.
#'
#' @param w pesos, uno por posición observada; típicamente 1/varianza. Se
#'   normalizan para que sumen n, así que la escala de la potencia no cambia y
#'   los resultados siguen siendo comparables con los del GLS sin pesos.
#'
#' Los datos son conteos y cada gen tiene una varianza estimable, pero el GLS
#' los trata todos igual: un gen con 5 lecturas pesa lo mismo que uno con 5000.
#' En el Lomb-Scargle de astronomía cada medición trae su barra de error y el
#' ajuste la usa; esto es lo mismo.
#'
#' No usa la FFT. La FFT exige pesos uniformes, y con pesos hay que evaluar las
#' sumas trigonométricas directamente: O(n*m) en vez de O(N log N). Con n de
#' unos cientos y m de unos miles es del orden de un segundo por espectro, así
#' que sirve para un análisis dirigido y no para las mil permutaciones de maxT.
gls_weighted <- function(y, terms, w) {
  y <- as.numeric(y); w <- as.numeric(w)
  if (length(w) != length(y)) {
    tsf_abort("gls_weighted: ", length(w), " peso(s) para ", length(y),
              " valor(es)")
  }
  if (anyNA(y) || any(!is.finite(y)) || anyNA(w) || any(!is.finite(w))) {
    tsf_abort("gls_weighted: valores o pesos no finitos. Pasá la señal por ",
              "gls_observed() y filtrá los pesos con las mismas posiciones.")
  }
  if (any(w < 0)) tsf_abort("gls_weighted: pesos negativos")
  if (sum(w) <= 0) tsf_abort("gls_weighted: todos los pesos son cero")

  n <- length(y)
  w <- w * n / sum(w)          # suman n: la potencia queda en la misma escala
  t <- terms$t_index
  N <- terms$N
  k <- terms$k

  W <- sum(w)
  ybar <- sum(w * y) / W
  yc <- y - ybar
  YY <- sum(w * yc^2)

  out <- vapply(k, function(kk) {
    ang <- 2 * pi * kk * t / N
    C <- cos(ang); S <- sin(ang)
    # Media flotante ponderada: los mismos productos que el GLS sin pesos,
    # con w dentro de cada suma.
    Cw <- sum(w * C) / W; Sw <- sum(w * S) / W
    Ch <- C - Cw; Sh <- S - Sw
    CC <- sum(w * Ch^2); SS <- sum(w * Sh^2); CS <- sum(w * Ch * Sh)
    YC <- sum(w * yc * Ch); YS <- sum(w * yc * Sh)
    D <- CC * SS - CS^2
    if (!is.finite(D) || abs(D) < .Machine$double.eps^0.5) return(c(0, 0))
    a <- (YC * SS - YS * CS) / D
    b <- (YS * CC - YC * CS) / D
    c(a, b)
  }, numeric(2))

  a <- out[1, ]; b <- out[2, ]
  amp <- sqrt(a^2 + b^2)
  m <- length(k)
  structure(list(k = k, freq = k / N, period = N / k,
                 N = rep(N, m), n_observed = rep(n, m),
                 amplitude = amp,
                 phase = atan2(-b, a),
                 power = amp^2 / 2,
                 power_normalised = if (YY > 0)
                   pmin(pmax(amp^2 * n / (2 * YY), 0), 1) else rep(0, m),
                 weight_ess = rep(sum(w)^2 / sum(w^2), m)),
            class = "data.frame", row.names = seq_len(m))
}

#' Pesos a partir de conteos, con el modelo de varianza de un conteo.
#'
#' @param counts conteos crudos por posición
#' @param unit "log" si la señal se modeló en escala log --que es el caso tras
#'   asinh(TPM)-- o "raw"
#'
#' Para un conteo con media mu, Var(mu) ~ mu bajo Poisson, así que en escala
#' log la varianza delta-método es ~1/mu y el peso es mu. Con sobredispersión
#' binomial negativa, Var = mu + mu^2/theta y el peso baja: el termino mu^2
#' domina en los genes muy expresados, asi que ponderar por mu puro les da
#' demasiada importancia.
#'
#' `theta` por defecto es 10, un valor tipico de RNA-seq de tejido. No se
#' estima de los datos aqui a proposito: estimarlo por gen con estas n haria el
#' peso mas ruidoso que el dato que corrige.
counts_to_weights <- function(counts, unit = c("log", "raw"), theta = 10) {
  unit <- match.arg(unit)
  mu <- pmax(as.numeric(counts), 0.5)   # piso: un conteo de 0 no da peso cero
  v_raw <- mu + mu^2 / theta            # binomial negativa
  w <- if (identical(unit, "log")) mu^2 / v_raw else 1 / v_raw
  w[!is.finite(w) | w <= 0] <- min(w[is.finite(w) & w > 0], na.rm = TRUE)
  w
}
