# grid.R -- reference gene grid and least-squares (Lomb-Scargle) spectrum.
#
# WHY
# ---
# The FFT axis used to be the rank of a gene among those that survived the
# expression filter. That axis is dataset-specific (14,495 vs 16,269 genes), so
# N differed per dataset and "period in genes" meant different things on each
# side. Worse, dropping a gene made its neighbours adjacent, compressing the
# axis wherever coverage was poor.
#
# Here the axis is fixed by the annotation: every gene of the allowed biotypes
# on a chromosome, ordered by start, indexed 1..N. N is a property of the
# annotation, identical for every dataset and every condition. Genes without an
# expression measurement are simply absent from the observed positions -- they
# are never imputed and never zero-filled as data.
#
# HOW THE SPECTRUM IS COMPUTED
# ----------------------------
# For each frequency we fit, by least squares,
#     y(t) ~ mu + a*cos(2*pi*k*t/N) + b*sin(2*pi*k*t/N)
# over the OBSERVED positions only (generalised Lomb-Scargle, Zechmeister &
# Kuerster 2009, floating mean). This is still Fourier analysis: on a complete
# grid the fit reproduces the DFT coefficients exactly, term by term. What it
# adds is a correct treatment of gaps.
#
# The sums are evaluated with two FFTs. Writing z for the signal zero-filled at
# missing positions and w for the 0/1 presence indicator,
#     sum(y*cos) = Re(FFT(z))_k        sum(cos)   = Re(FFT(w))_k
#     sum(y*sin) = -Im(FFT(z))_k       sum(sin)   = -Im(FFT(w))_k
#     sum(cos^2) = (n + Re(FFT(w))_2k)/2      etc.
# The zero-filling here is a computational device, not imputation: the missing
# positions enter the normalisation through w, which is exactly what naive
# zero-filling fails to do. The result is identical to fitting only the
# observed points, at O(N log N) instead of O(N * n_frequencies).
#
# WHAT GAPS COST
# --------------
# With gaps the sinusoids are no longer orthogonal over the observed positions,
# so frequency bins are not independent, Parseval no longer holds exactly, and
# the Nyquist limit is not uniquely defined. The first is why the permutation
# null is mandatory rather than optional; the third is why every peak is
# reported alongside the spectral window (see spectral_window), which shows how
# much of the apparent structure the sampling pattern alone can produce.

#' Canonical chromosome names.
#'
#' Annotations disagree on "chr1" vs "1" vs "CHR1", and MT/M/chrM. Everything
#' downstream keys on the chromosome name, so a single spelling has to be
#' imposed at the point the annotation is read -- not patched later.
normalise_chrom_names <- function(x) {
  v <- toupper(trimws(as.character(x)))
  v <- sub("^CHR", "", v)
  v[v %in% c("M", "MT", "MITO")] <- "MT"
  v
}

#' Reference grid: every annotated gene of the allowed biotypes, per chromosome.
#'
#' @param annot annotation table with gene_id, chr, start and optionally gene_type
#' @param biotypes regular expression matched against gene_type; NULL keeps all
#' @param axis "gene" places genes at consecutive ranks; "bp" places them in
#'   fixed base-pair bins of `bin_size`.
#'
#'   THE GENE-RANK AXIS IS CONFOUNDED WITH THE BIOLOGY. Genes are not evenly
#'   spaced: gene density varies more than tenfold along a chromosome, and it
#'   varies WITH the things a spectral analysis would want to explain -- dense
#'   regions are open chromatin, early-replicating, GC-rich; sparse regions are
#'   lamina-associated, late-replicating gene deserts. So the rank-to-megabase
#'   factor is not a constant: a period of 30 genes can be 1 Mb in a dense
#'   region and 15 Mb in a desert, within one chromosome.
#'
#'   The consequence is that "period = 30 genes" does not name a physical
#'   scale, and no mechanism can be assigned to it: TADs are 0.1-1 Mb,
#'   replication timing domains 0.4-0.8 Mb, lamina-associated domains
#'   0.1-10 Mb. Those are base-pair statements and the rank axis cannot reach
#'   them.
#'
#'   Coverage does not suffer -- it improves, and becomes honest. On the rank
#'   axis "coverage" is expressed genes over annotated genes, and the axis is
#'   DEFINED by gene positions, so it is high by construction. On a bp axis it
#'   is bins-with-data over bins-in-the-chromosome, which is how much of the
#'   chromosome was actually observed: about 83% at 100 kb on chr1 against
#'   65-77% on the rank axis. The gaps are the regime the floating-mean GLS
#'   exists for.
#'
#' @param bin_size bin width in base pairs, for axis = "bp". 100 kb gives TAD
#'   resolution; 250 kb is near-complete coverage with less resolution. Run
#'   both: a band that appears at one width and not the other is a property of
#'   the binning.
#' Longitudes de los cromosomas de GRCh38, en pares de base.
#'
#' Hacen falta porque `max(start)` de los genes anotados NO es la longitud del
#' cromosoma: termina en el último gen del universo seleccionado. Usar eso como
#' N desplaza frecuencias, periodos, fase y cobertura, y el sesgo es peor en los
#' acrocéntricos --13, 14, 15, 21, 22-- donde el brazo corto es heterocromatina
#' sin genes codificantes.
#'
#' Y cambia con el universo génico: con `^protein-coding$` el último gen cae en
#' un sitio y con lncRNA incluidos en otro, así que dos corridas del mismo
#' cromosoma tendrían N distintos y sus espectros no serían comparables.
GRCH38_CHROM_LENGTHS <- c(
  "1" = 248956422, "2" = 242193529, "3" = 198295559, "4" = 190214555,
  "5" = 181538259, "6" = 170805979, "7" = 159345973, "8" = 145138636,
  "9" = 138394717, "10" = 133797422, "11" = 135086622, "12" = 133275309,
  "13" = 114364328, "14" = 107043718, "15" = 101991189, "16" = 90338345,
  "17" = 83257441, "18" = 80373285, "19" = 58617616, "20" = 64444167,
  "21" = 46709983, "22" = 50818468, "X" = 156040895, "Y" = 57227415,
  "MT" = 16569, "M" = 16569)

#' @param chrom_lengths longitudes en pb, con nombres de cromosoma. Por defecto
#'   GRCh38. Un cromosoma ausente de la tabla cae al último gen anotado con una
#'   advertencia, en vez de fallar -- un ensamblado no humano debe poder correr.
build_reference_grid <- function(annot, chrom_levels, biotypes = NULL,
                                 min_genes_per_chr = 8L,
                                 axis = c("gene", "bp"),
                                 bin_size = 100000L,
                                 chrom_lengths = GRCH38_CHROM_LENGTHS) {
  axis <- match.arg(axis)
  g <- annot[annot$chr %in% chrom_levels & !is.na(annot$start), , drop = FALSE]

  if (!is.null(biotypes) && "gene_type" %in% colnames(g)) {
    observed <- sort(unique(g$gene_type))
    keep <- grepl(biotypes, g$gene_type, ignore.case = TRUE)
    tsf_log("Gene biotypes present: ", paste(observed, collapse = ", "))
    tsf_log("Grid keeps ", sum(keep), "/", nrow(g), " annotated genes (pattern: ",
            biotypes, ")")
    if (!any(keep)) tsf_abort("No gene matches gene_universe = '", biotypes,
                              "'. Observed biotypes: ", paste(observed, collapse = ", "))
    g <- g[keep, , drop = FALSE]
  } else if (!is.null(biotypes)) {
    tsf_warn("gene_universe was set but the annotation has no gene_type column; ",
             "the grid keeps every annotated gene")
  }

  g <- g[order(match(g$chr, chrom_levels), g$start, g$gene_id), ]
  g <- g[!duplicated(g$gene_id), ]

  if (identical(axis, "gene")) {
    # ADVERTENCIA INEQUIVOCA, no una nota en la documentacion.
    #
    # El default sigue siendo "gene" para que un arbol de resultados existente
    # siga significando lo que decia, pero eso deja un hueco entre la facilidad
    # de uso y la recomendacion cientifica: una corrida estandar produce
    # periodos que alguien puede interpretar despues como distancias.
    #
    # La densidad genica varia mas de diez veces a lo largo de un cromosoma Y
    # varia CON el estado de cromatina que un resultado espectral querria
    # explicar, asi que el factor rango->Mb no es constante y un periodo en
    # genes no nombra una distancia.
    tsf_warn("Eje de RANGO DE GEN: un periodo en genes NO es una distancia ",
             "genomica y no se puede comparar con TADs, LADs, dominios de ",
             "replicacion ni citobandas. La densidad genica varia mas de 10x ",
             "por cromosoma y covaria con el estado de cromatina, asi que el ",
             "factor rango->Mb no es constante. Para interpretacion biologica: ",
             "--grid-axis bp --bin-size 100000")
    g$grid_index <- stats::ave(seq_len(nrow(g)), g$chr, FUN = seq_along)
    n_by_chr <- table(g$chr)
    g$grid_N <- as.integer(n_by_chr[g$chr])
  } else {
    bin_size <- as.numeric(bin_size)
    if (!is.finite(bin_size) || bin_size < 1000) {
      tsf_abort("bin_size must be at least 1000 bp; got ", bin_size)
    }
    # Bin 1 is the first bin of the chromosome, not the first bin that contains
    # a gene: N must be the chromosome's length in bins, or the axis would
    # shrink to the span of the observed genes and a gene desert at either end
    # would silently disappear from the grid.
    g$grid_index <- as.integer(floor(g$start / bin_size)) + 1L

    # N es la LONGITUD DEL CROMOSOMA en bins, no el último gen anotado.
    # max(start) termina el eje donde acaba el universo génico seleccionado, y
    # eso desplaza frecuencias, periodos y fase -- y hace que dos universos
    # distintos den N distintos para el mismo cromosoma.
    chrs <- unique(as.character(g$chr))
    have <- chrs[chrs %in% names(chrom_lengths)]
    miss <- setdiff(chrs, have)
    L <- stats::setNames(rep(NA_real_, length(chrs)), chrs)
    L[have] <- as.numeric(chrom_lengths[have])
    if (length(miss)) {
      fallback <- tapply(g$start, as.character(g$chr), max)
      L[miss] <- as.numeric(fallback[miss])
      tsf_warn("Sin longitud declarada para: ", paste(miss, collapse = ", "),
               ". Se usa el último gen anotado, que subestima el cromosoma y ",
               "desplaza los periodos. Pasá `chrom_lengths` para el ensamblado ",
               "correcto.")
    }
    g$grid_N <- as.integer(ceiling(L[as.character(g$chr)] / bin_size))
    g$chrom_length <- as.numeric(L[as.character(g$chr)])

    per_chr <- tapply(g$grid_index, g$chr, function(x) length(unique(x)))
    cov <- per_chr / tapply(g$grid_N, g$chr, function(x) x[1])
    tsf_log("Grid axis: base pairs, ", format(bin_size, big.mark = ","),
            " bp bins. Occupied bins per chromosome: median ",
            round(100 * stats::median(cov), 1), "%, range ",
            round(100 * min(cov), 1), "-", round(100 * max(cov), 1), "%")
    tsf_log("  ", nrow(g), " gene(s) into ",
            sum(unlist(per_chr)), " occupied bin(s) of ",
            sum(vapply(split(g$grid_N, g$chr), function(x) x[1], numeric(1))),
            " total; a bin with several genes is aggregated at ingest, an ",
            "empty bin stays UNOBSERVED and is never zero-filled")
  }

  n_by_chr <- table(g$chr)
  short <- names(n_by_chr)[n_by_chr < min_genes_per_chr]
  if (length(short)) {
    tsf_warn("Chromosomes below ", min_genes_per_chr, " grid genes, dropped: ",
             paste(short, collapse = ", "))
    g <- g[!g$chr %in% short, ]
  }
  tsf_log("Reference grid: ", nrow(g), " genes over ", length(unique(g$chr)),
          " chromosomes")
  g
}

#' Observed grid positions per chromosome for a dataset.
#'
#' Returns, per chromosome, the row indices into `genes` and their positions on
#' the reference grid, plus N and the coverage fraction.
grid_index <- function(genes, chrom_levels, min_observed = 8L,
                       min_coverage = 0.05) {
  out <- list()
  for (chr_now in chrom_levels) {
    idx <- which(genes$chr == chr_now)
    if (!length(idx)) next
    N <- genes$grid_N[idx[1]]
    ord <- idx[order(genes$grid_index[idx])]
    t <- genes$grid_index[ord]
    coverage <- length(t) / N
    if (length(t) < min_observed || coverage < min_coverage) {
      tsf_warn("chr", chr_now, ": ", length(t), "/", N, " observed (",
               round(100 * coverage), "%), skipped")
      next
    }
    out[[chr_now]] <- list(rows = ord, t = as.integer(t), N = as.integer(N),
                           coverage = coverage)
  }
  out
}

# --- the spectrum -------------------------------------------------------------

#' Sums needed by the GLS fit, from the presence indicator alone.
#'
#' Depends only on which positions are observed, so it is computed once per
#' chromosome and reused across samples and permutations.
gls_window_terms <- function(t, N) {
  k <- seq_len(floor((N - 1) / 2))
  w <- numeric(N); w[t] <- 1
  W <- stats::fft(w)
  n <- length(t)

  idx  <- k + 1L
  idx2 <- (2L * k) %% N + 1L          # wraps for k > N/4, which is correct

  Cw <- Re(W[idx]);  Sw <- -Im(W[idx])
  C2 <- Re(W[idx2]); S2 <- -Im(W[idx2])

  CC <- (n + C2) / 2
  SS <- (n - C2) / 2
  CS <- S2 / 2

  # Floating-mean (hatted) versions
  CCh <- CC - Cw^2 / n
  SSh <- SS - Sw^2 / n
  CSh <- CS - Cw * Sw / n
  D <- CCh * SSh - CSh^2

  list(k = k, n = n, N = N, Cw = Cw, Sw = Sw,
       CCh = CCh, SSh = SSh, CSh = CSh, D = D,
       window_power = (Cw^2 + Sw^2) / n^2)
}

#' Generalised Lomb-Scargle spectrum of one signal on a gappy grid.
#'
#' @param y values at the observed positions, in the order of `t`
#' @param terms output of gls_window_terms for those positions
#' @return data.frame with k, freq, period, amplitude, phase, power,
#'   power_normalised (the classic 0-1 GLS statistic) and window_power
gls_spectrum <- function(y, terms) {
  y <- as.numeric(y)
  # Callers are expected to have removed unmeasured positions with
  # gls_observed(); a non-finite value reaching here would be counted as
  # observed by w and fitted as a zero. Refuse rather than quietly imputing.
  if (anyNA(y) || any(!is.finite(y))) {
    tsf_abort("gls_spectrum() received ", sum(!is.finite(y)), " non-finite ",
              "value(s) at observed positions. Route the signal through ",
              "gls_observed() so unmeasured genes are dropped from the fit ",
              "instead of being zero-filled.")
  }
  n <- terms$n; N <- terms$N; k <- terms$k

  z <- numeric(N); z[terms$t_index] <- y
  Z <- stats::fft(z)
  YC <- Re(Z[k + 1L]); YS <- -Im(Z[k + 1L])
  Ysum <- sum(y)
  YY <- sum((y - Ysum / n)^2)

  YCh <- YC - Ysum * terms$Cw / n
  YSh <- YS - Ysum * terms$Sw / n

  D <- terms$D
  bad <- !is.finite(D) | abs(D) < .Machine$double.eps^0.5
  a <- (YCh * terms$SSh - YSh * terms$CSh) / D
  b <- (YSh * terms$CCh - YCh * terms$CSh) / D
  a[bad] <- 0; b[bad] <- 0

  amp <- sqrt(a^2 + b^2)
  phase <- atan2(-b, a)                       # y ~ amp * cos(2*pi*k*t/N + phase)
  power <- amp^2 / 2                          # same scale as the old one-sided power
  pnorm_gls <- (terms$SSh * YCh^2 + terms$CCh * YSh^2 - 2 * terms$CSh * YCh * YSh) /
    (YY * D)
  pnorm_gls[bad | !is.finite(pnorm_gls)] <- 0
  pnorm_gls <- pmin(pmax(pnorm_gls, 0), 1)

  data.frame(k = k, freq = k / N, period = N / k,
             amplitude = amp, phase = phase, power = power,
             power_normalised = pnorm_gls,
             window_power = terms$window_power,
             N = N, n_observed = n,
             stringsAsFactors = FALSE)
}

#' Prepare the reusable terms for a chromosome (positions + window).
#'
#' PHASE CONVENTION: positions are taken as t0 = grid_index - 1, so the fitted
#' model is amp * cos(2*pi*k*t0/N + phase) with t0 = 0..N-1. This matches
#' reconstruct_signal(), which evaluates the sinusoid over 0:(N-1), and matches
#' the DFT convention of R's fft(), where bin k has exponent -2i*pi*k*(j-1)/N.
#' Comparing against a design matrix built on 1-based positions produces a
#' constant phase offset of 2*pi*k/N that leaves the amplitude untouched.
gls_prepare <- function(t, N) {
  terms <- gls_window_terms(t, N)
  terms$t_index <- as.integer(t)
  terms
}

#' Restrict a signal to the positions that actually carry a measurement.
#'
#' WHY THIS EXISTS
#' ---------------
#' gls_spectrum() replaces a non-finite value with zero before the FFT. For the
#' positions that are absent from `t` that is the documented computational
#' device: they never enter the fit, because the normalisation goes through the
#' presence indicator w. For a position that IS in `t` but whose value is NA it
#' is something else entirely -- w counts it as observed, so a zero lands at a
#' fixed grid position and is treated as a measurement of zero expression. That
#' is exactly the failure the reference grid was built to avoid, and it would
#' arrive silently.
#'
#' Ingest drops genes with no usable value, so in a normal run every observed
#' position is finite and this function returns the shared terms untouched at no
#' cost. It matters in the case ingest cannot catch: a gene finite in most
#' samples but NA in one (a duplicate-identifier group that collapses to all-NA
#' for that sample alone) survives the expression filter and reaches here with a
#' hole in one column.
#'
#' When that happens the position is removed from the observed set and the
#' window terms are rebuilt for the reduced set, so the missingness of this
#' signal is described by its own w. The cost is one extra pair of FFTs for the
#' affected signal only.
#'
#' @param y values at the positions in `t`, in that order
#' @param t observed grid positions
#' @param N grid length for the chromosome
#' @param terms optional precomputed gls_prepare(t, N), reused when complete
#' @param min_observed refuse to fit fewer than this many points
#' @return list(y, terms, n_dropped), or NULL when too little survives
gls_observed <- function(y, t, N, terms = NULL, min_observed = 8L) {
  y <- as.numeric(y)
  ok <- is.finite(y)
  if (all(ok)) {
    return(list(y = y, terms = terms %||% gls_prepare(t, N), n_dropped = 0L))
  }
  if (sum(ok) < min_observed) return(NULL)
  list(y = y[ok], terms = gls_prepare(t[ok], N), n_dropped = sum(!ok))
}

#' Spectral window: what the sampling pattern alone produces.
#'
#' A peak whose frequency coincides with a strong window feature can be an alias
#' of structure at another frequency rather than structure at its own.
spectral_window <- function(t, N) {
  terms <- gls_window_terms(t, N)
  data.frame(k = terms$k, freq = terms$k / N, period = N / terms$k,
             window_power = terms$window_power,
             coverage = terms$n / N,
             stringsAsFactors = FALSE)
}

#' Flag peaks sitting on a strong window feature.
#'
#' `window_rank` is the rank of that frequency in the window spectrum: rank 1
#' means the sampling pattern itself is strongest exactly there, which is the
#' worst case for interpreting the peak as biological.
annotate_window <- function(peaks, terms) {
  wp <- terms$window_power
  peaks$window_power <- wp[match(peaks$k, terms$k)]
  peaks$window_rank <- rank(-wp)[match(peaks$k, terms$k)]
  peaks$window_pct <- 100 * peaks$window_rank / length(wp)
  peaks
}
