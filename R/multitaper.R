# multitaper.R -- multitaper Lomb-Scargle on an incomplete grid.
#
# WHY
# ---
# The periodogram is a poor estimator of a spectrum in three specific ways, and
# all three matter here.
#
#   INCONSISTENT. Its variance does not fall as data accumulate: with one
#   realisation per (sample, chromosome) the estimate at each frequency has
#   roughly 100% relative error however many genes the chromosome holds. That is
#   why nothing survives a permutation null -- the observed value is as noisy as
#   the null draws.
#
#   BIASED in finite samples.
#
#   LEAKY. Power from one frequency appears at others through the sidelobes of
#   the spectral window, which on a 26%-covered grid are large.
#
# Thomson's multitaper fixes all three by averaging periodograms computed under
# K orthogonal tapers: the variance falls roughly as 1/K while the bias and the
# leakage are controlled by the time-bandwidth product. It assumes uniform
# sampling, which this grid is not, so the tapers are evaluated at the OBSERVED
# positions and the Lomb-Scargle (here GLS) machinery does the rest.
#
# WHAT IT COSTS
# -------------
# Resolution. Averaging K tapers with time-bandwidth NW smooths the spectrum
# over a band of about 2*NW/N in frequency, so two components closer than that
# merge into one. With NW = 3 and K = 5 on N = 2000 that is a band of 3
# frequency bins -- acceptable when the question is "is there power near this
# scale", wrong when it is "is the period 30.1 or 30.4 genes".
#
# The trade is deliberate: this project has never had a resolution problem, it
# has had a variance problem.

#' Discrete prolate spheroidal sequences (Slepian tapers).
#'
#' The DPSS of length N and time-bandwidth NW are the eigenvectors of a
#' symmetric tridiagonal matrix, so base R's eigen() is enough and no package
#' is needed. Solved as the full symmetric problem rather than by a specialised
#' tridiagonal routine: N here is a few thousand at most, and correctness on
#' the first try is worth more than the speed.
#'
#' @return N x K matrix, columns ordered by decreasing concentration
# Los tapers dependen solo de (N, NW, K), no de los datos, y calcularlos exige
# la descomposicion propia de una matriz N x N: 7.15 s con N = 2000. Sin cache
# eso se repite en cada permutacion y maxT pasaria de 50 minutos a 118 dias.
# Medido, no estimado.
.tsf_mt_warned <- new.env(parent = emptyenv())

.dpss_cache <- new.env(parent = emptyenv())

dpss_tapers <- function(N, NW = 3, K = NULL) {
  N <- as.integer(N)
  ck <- paste(N, NW, K %||% "auto", sep = "|")
  hit <- .dpss_cache[[ck]]
  if (!is.null(hit)) return(hit)
  if (N < 8L) tsf_abort("dpss_tapers: N must be at least 8, got ", N)
  if (is.null(K)) K <- max(1L, floor(2 * NW) - 1L)
  K <- min(as.integer(K), N - 1L)
  W <- NW / N

  t <- seq_len(N) - 1L
  diagv <- ((N - 1 - 2 * t) / 2)^2 * cos(2 * pi * W)
  offv <- t[-1] * (N - t[-1]) / 2

  A <- matrix(0, N, N)
  diag(A) <- diagv
  idx <- cbind(seq_len(N - 1L), 2:N)
  A[idx] <- offv
  A[idx[, c(2, 1)]] <- offv

  e <- eigen(A, symmetric = TRUE)
  V <- e$vectors[, seq_len(K), drop = FALSE]
  # Sign convention: even-order tapers positive at the centre, odd-order with a
  # positive first half. Without it the sign is arbitrary per call and phases
  # computed from different tapers cannot be compared.
  for (j in seq_len(K)) {
    v <- V[, j]
    flip <- if (j %% 2L == 1L) sum(v) < 0 else v[which.max(abs(v))] < 0
    if (flip) V[, j] <- -v
  }
  assign(ck, V, envir = .dpss_cache)
  V
}

#' Vaciar el cache de tapers.
#'
#' Cada entrada es una matriz N x K de doubles: con N = 2000 y K = 5 son
#' 80 KB, y 24 cromosomas caben en 2 MB. No hace falta vaciarlo en una corrida
#' normal; existe para los tests, que verifican el cache y necesitan poder
#' medir el caso frio.
dpss_cache_clear <- function() {
  rm(list = ls(envir = .dpss_cache), envir = .dpss_cache)
  invisible(NULL)
}

#' Multitaper GLS spectrum.
#'
#' @param y values at the observed positions
#' @param terms output of gls_terms() for those positions
#' @param NW time-bandwidth product; 3 to 4 is the usual range
#' @param K tapers to average; defaults to 2*NW - 1, the number that stay
#'   well concentrated
#' @param adaptive weight each taper by its own concentration. The higher-order
#'   tapers leak more, so weighting by eigenvalue keeps them from dominating a
#'   frequency where the true power is small.
#'
#' @return the same columns as gls_spectrum(), plus:
#'   n_tapers   how many tapers contributed
#'   power_sd   spread of power across tapers, which is the estimator's own
#'              uncertainty at that frequency -- the quantity a single
#'              periodogram cannot report at all
gls_multitaper <- function(y, terms, NW = 3, K = NULL, adaptive = TRUE) {
  y <- as.numeric(y)
  if (anyNA(y) || any(!is.finite(y))) {
    tsf_abort("gls_multitaper() received ", sum(!is.finite(y)), " non-finite ",
              "value(s). Route the signal through gls_observed() first.")
  }
  N <- terms$N

  # EL ANCHO DE BANDA TIENE QUE CABER EN LA MALLA.
  #
  # Promediar K tapers con banda NW suaviza sobre 2*NW bins de frecuencia. En
  # una malla corta eso se traga el pico: medido con N = 200 y un componente en
  # k = 10, NW = 3 cubre de k = 7 a k = 13 y el pico cae al rango 2 mientras el
  # periodograma lo deja en el 1. Con NW <= 2.5 vuelve al rango 1.
  #
  # La regla: el ancho no debe pasar de la quinta parte de las frecuencias
  # disponibles, o el suavizado deja de ser suavizado y pasa a ser borrado.
  # Lo que importa NO es el ancho frente al total de frecuencias, sino frente
  # a la frecuencia de cada componente: un pico en k = 10 con ancho 6 queda
  # cubierto de k = 7 a k = 13, el 60% de su propia frecuencia, mientras uno en
  # k = 100 solo pierde el 6%. Las frecuencias BAJAS --los periodos largos, que
  # es donde vive la estructura cromosomica-- son las vulnerables.
  #
  # Se reporta por frecuencia en `bandwidth_fraction` en vez de avisar una vez,
  # porque el mismo espectro tiene frecuencias seguras y frecuencias borradas.
  # Once per (N, NW) rather than once per call: this function runs inside the
  # permutation loop, so the same warning was printed thousands of times per
  # stage and buried everything else in the log.
  n_low <- sum(terms$k < 2 * NW * 2)
  warn_key <- paste("mt", terms$N, NW)
  if (n_low > 0 && !isTRUE(.tsf_mt_warned[[warn_key]])) {
    .tsf_mt_warned[[warn_key]] <- TRUE
    tsf_warn("multitaper: con NW = ", NW, " el suavizado cubre ", 2 * NW,
             " bins, asi que las ", n_low, " frecuencia(s) con k < ",
             4 * NW, " quedan dentro de su propia banda y su potencia se ",
             "mezcla con la de sus vecinas. Medido: un pico en k = 10 con ",
             "NW = 3 cae del rango 1 al 2; con NW <= 2.5 vuelve al 1. Ver ",
             "`bandwidth_fraction` en la salida.")
  }

  V <- dpss_tapers(N, NW = NW, K = K)
  K <- ncol(V)

  # The tapers live on the full grid 1..N; the signal lives at t_index. Taking
  # V at those rows is what carries multitaper onto an incomplete grid: each
  # taper becomes a weight per observed position.
  Vobs <- V[terms$t_index, , drop = FALSE]

  # Concentration of each taper ON THE OBSERVED SET, not on the full grid.
  # A taper concentrated in a region the sampling never covers contributes
  # almost nothing and should not be weighted as though it did.
  wt <- colSums(Vobs^2)
  wt[!is.finite(wt) | wt <= 0] <- 0
  if (!any(wt > 0)) tsf_abort("gls_multitaper: no taper has support on the observed positions")
  if (!adaptive) wt[wt > 0] <- 1
  wt <- wt / sum(wt)

  # UNA FFT DE MATRIZ, no K llamadas.
  #
  # gls_spectrum() ya usa FFT, asi que el costo por taper es una FFT de
  # longitud N mas unas operaciones vectoriales. Llamarla K veces repite la
  # sobrecarga de R --validacion, asignacion de vectores, indexado-- K veces,
  # y eso es lo que costaba 21.6x cuando el trabajo real es 5x.
  #
  # Aqui las K senales tapereadas se colocan en una matriz N x K de una vez y
  # mvfft() hace las K transformadas en una sola llamada a la FFT compilada.
  # Todo lo que no depende de y --Cw, Sw, CCh, SSh, CSh, D-- se reusa tal cual.
  Z <- matrix(0, nrow = N, ncol = K)
  keep <- which(wt > 0)
  for (j in keep) {
    yj <- y * Vobs[, j]
    # Centrar por taper: un taper no preserva la media, asi que la senal
    # tapereada arrastra un termino DC que la media flotante absorberia en las
    # frecuencias bajas.
    Z[terms$t_index, j] <- yj - mean(yj)
  }
  if (!length(keep)) tsf_abort("gls_multitaper: todo taper degenero")

  F <- stats::mvfft(Z)
  kk <- terms$k + 1L
  YC <- Re(F[kk, keep, drop = FALSE])
  YS <- -Im(F[kk, keep, drop = FALSE])
  # La suma y la varianza de cada senal tapereada y centrada: la suma es cero
  # por construccion, asi que los terminos de media flotante se anulan y YCh,
  # YSh son YC, YS.
  D <- terms$D
  bad <- !is.finite(D) | abs(D) < .Machine$double.eps^0.5
  denom <- ifelse(bad, 1, D)
  A <- (YC * terms$SSh - YS * terms$CSh) / denom
  B <- (YS * terms$CCh - YC * terms$CSh) / denom
  A[bad, ] <- 0
  B[bad, ] <- 0

  AMP <- sqrt(A^2 + B^2)
  PWR <- AMP^2 / 2
  PHS <- exp(1i * atan2(-B, A))

  w <- wt[keep] / sum(wt[keep])
  am <- as.numeric(AMP %*% w)
  pw <- as.numeric(PWR %*% w)
  ph <- Arg(as.complex(PHS %*% w))
  psd <- if (length(keep) > 1L) apply(PWR, 1L, stats::sd) else rep(NA_real_, nrow(PWR))

  # Las columnas que no dependen de y se construyen directamente. Llamar a
  # gls_spectrum() solo para copiarlas costaba un 8% del total, y era la parte
  # que dejaba el factor en 9x cuando el trabajo real es 5x.
  # structure() y no data.frame(): data.frame() valida cada columna y genera
  # nombres de fila, y eso costaba MAS que la llamada a gls_spectrum que
  # reemplaza -- el factor subio de 9.0x a 11.7x al usarlo. Medido.
  m <- length(terms$k)
  cols <- list(k = terms$k, freq = terms$k / N, period = N / terms$k,
               N = rep(N, m), n_observed = rep(terms$n, m))
  if (!is.null(terms$window_power)) cols$window_power <- terms$window_power

  # power_normalised es un cociente a la varianza de la propia senal y no es
  # lineal en las senales tapereadas, asi que se recalcula de la amplitud
  # promediada en vez de promediarse entre tapers.
  yy <- sum((y - mean(y))^2)
  pn <- if (yy > 0) pmin(pmax(am^2 * terms$n / (2 * yy), 0), 1) else rep(0, length(am))

  cols$amplitude <- am
  cols$phase <- ph
  cols$power <- pw
  cols$power_normalised <- pn
  cols$n_tapers <- rep(length(keep), m)
  # Que fraccion de su propia frecuencia cubre el suavizado. Por encima de ~0.5
  # la potencia de esa frecuencia esta mezclada con la de sus vecinas y su
  # rango no significa lo que parece.
  cols$bandwidth_fraction <- pmin(2 * NW / pmax(terms$k, 1), 1)
  cols$power_sd <- psd
  structure(cols, class = "data.frame", row.names = seq_len(m))
}

#' Single dispatch point for the spectral estimator.
#'
#' Every stage calls this instead of gls_spectrum() directly, so the choice of
#' estimator is made once from the config and cannot diverge between stages --
#' a spectrum computed by one estimator and a null computed by the other would
#' not be a permutation test of anything.
#'
#' Defaults to "periodogram" so an existing results tree still means what it
#' said.
tsf_spectrum <- function(y, terms, project = NULL) {
  est <- (project %||% list())$estimator %||% "periodogram"
  if (identical(est, "multitaper")) {
    gls_multitaper(y, terms,
                   NW = (project %||% list())$mt_nw %||% 3,
                   K  = (project %||% list())$mt_k  %||% 5)
  } else if (identical(est, "periodogram")) {
    gls_spectrum(y, terms)
  } else {
    tsf_abort("estimator must be 'periodogram' or 'multitaper', got '", est, "'")
  }
}

#' Estimator settings from a project config, as a small list.
#'
#' Passed by value to every stage that computes a spectrum, so the choice
#' cannot diverge between them.
estimator_spec <- function(project) {
  list(estimator = (project %||% list())$estimator %||% "periodogram",
       mt_nw = (project %||% list())$mt_nw %||% 3,
       mt_k = (project %||% list())$mt_k %||% 5)
}
