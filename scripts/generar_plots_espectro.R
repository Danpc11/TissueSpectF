#!/usr/bin/env Rscript
# generar_plots_espectro.R -- para UNA cohorte, corriendo DIRECTO en el
# servidor sobre el arbol de resultados real (nada de copias locales):
#
#  1. Espectro completo por condicion (todas las k, todos los chr,
#     concatenados por cromosoma), con los picos ya clasificados
#     (condition_invariants_<condicion>.tsv, 60/80/90%, el resultado por
#     default de stage_consensus) marcados con punto rojo. Eje X = k
#     (ciclos por gen) concatenado por cromosoma; eje Y = potencia
#     normalizada en escala log.
#
#  2. Para cada pico marcado: reconstruccion de TODOS los genes del
#     cromosoma en esa onda (eq:loading), usando la fase REAL calculada
#     aqui mismo desde spectra_samples_<condicion>.tsv (Arg(mean(exp(i*
#     phase))), exactamente lo que phase_locking() hace en consensus.R) --
#     no una fase generica de la cohorte, la de la condicion exacta donde
#     el pico se detecto. Plot coloreado (rojo=nodo |carga|<0.2,
#     amarillo=cresta carga>=0.8, azul=valle carga<=-0.8, etiquetados
#     SOLO esos genes) + TSV completo (todos los genes, sin excepcion) +
#     TSV de 3 columnas (nodo/cresta/valle, alineados por orden de
#     aparicion en el cromosoma -- un gen esta en una sola zona a la vez,
#     nunca en las tres, asi que no comparten fila por posicion fisica).
#
#  3. Resta de espectros entre condiciones consecutivas (healthy->F0->F1->
#     ...->F4) y el extremo (healthy vs F4): el log de una diferencia
#     negativa no existe, asi que en vez de log(A-B) se superponen las DOS
#     condiciones en el mismo plot (misma escala log del punto 1) -- la
#     resta se lee como la separacion visual entre curvas. Picos
#     marcados: la union de los picos de ambas condiciones. Reconstruccion
#     igual que el punto 2, para esos picos, con la fase de CADA condicion
#     por separado (dos reconstrucciones por pico si aparece en ambas).
#
# Corre con SOLO base R (sin paquetes) para no depender de instalaciones.
#
# Uso (desde la raiz de TissueSpectF en el servidor):
#   Rscript generar_plots_espectro.R --cohort GSE130970 \
#     --results-dir /home/storres/tejidos_fft/liver_fft/GEO/results_v50.old \
#     --interim-dir /home/storres/tejidos_fft/liver_fft/GEO/interim_v50.old \
#     --out-dir /home/storres/tejidos_fft/liver_fft/GEO/plots_GSE130970

if (any(commandArgs(TRUE) %in% c("-h", "--help"))) {
  self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  if (!is.na(self) && file.exists(self)) {
    hdr <- readLines(self, warn = FALSE); hdr <- hdr[!grepl("^#!", hdr)]
    stop_at <- which(!grepl("^#", hdr) & nzchar(hdr))[1]
    if (is.na(stop_at)) stop_at <- length(hdr) + 1L
    hdr <- hdr[seq_len(stop_at - 1L)]
    cat(paste(sub("^#[ ]?", "", hdr[grepl("^#", hdr)]), collapse = "\n"), "\n")
  }
  quit(save = "no", status = 0)
}

args <- commandArgs(trailingOnly = TRUE)
opt <- list()
i <- 1
while (i <= length(args)) { opt[[args[i]]] <- args[i + 1]; i <- i + 2 }
need <- c("--cohort", "--results-dir", "--interim-dir", "--out-dir")
miss <- setdiff(need, names(opt))
if (length(miss)) stop("faltan argumentos: ", paste(miss, collapse = ", "))

cohort <- opt[["--cohort"]]
results_dir <- opt[["--results-dir"]]
interim_dir <- opt[["--interim-dir"]]
spectra_dir <- file.path(results_dir, cohort, "spectra")
invariants_dir <- file.path(results_dir, cohort, "consensus")
genes_path <- file.path(interim_dir, cohort, "genes.tsv")

dir.create(opt[["--out-dir"]], showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(opt[["--out-dir"]], "espectros"), showWarnings = FALSE)
dir.create(file.path(opt[["--out-dir"]], "restas"), showWarnings = FALSE)

chr_order <- c(as.character(1:22), "X", "Y")

# Servidores headless (sin X11) -- y algunos Mac con X11 roto -- suelen
# "abrir" png(type="cairo") sin lanzar ningun error y aun asi nunca escribir
# el archivo (falla silenciosa: solo un warning). tryCatch() por si solo NO
# lo detecta, porque no hay error que atrapar. La unica prueba confiable es
# verificar, DESPUES de graficar y cerrar el device, que el archivo existe
# y pesa mas de 0 bytes -- por eso plot_fn() recibe la funcion que hace TODO
# el dibujo, para poder reintentarla de cero con cada tipo de device. Si
# ningun tipo de png funciona, cae a pdf() (sin dependencias del sistema
# operativo) como red de seguridad, cambiando la extension del archivo.
graficar_robusto <- function(path_png, width, height, res, plot_fn) {
  tipos <- list(
    function() grDevices::png(path_png, width = width, height = height, res = res, type = "cairo"),
    function() grDevices::png(path_png, width = width, height = height, res = res, type = "cairo-png"),
    function() grDevices::png(path_png, width = width, height = height, res = res, type = "quartz"),
    function() grDevices::png(path_png, width = width, height = height, res = res)
  )
  for (abrir in tipos) {
    unlink(path_png)
    abierto <- tryCatch({ abrir(); TRUE }, error = function(e) FALSE)
    if (!abierto) next
    dibujado <- tryCatch({ plot_fn(); TRUE }, error = function(e) {
      message("  error dibujando ", path_png, ": ", conditionMessage(e)); FALSE
    })
    grDevices::dev.off()
    if (dibujado && file.exists(path_png) && file.size(path_png) > 0) return(invisible(path_png))
  }
  # Red de seguridad final: pdf(), que no depende de cairo/X11/quartz.
  path_pdf <- sub("\\.png$", ".pdf", path_png)
  grDevices::pdf(path_pdf, width = width / res, height = height / res)
  dibujado <- tryCatch({ plot_fn(); TRUE }, error = function(e) {
    message("  error dibujando (pdf) ", path_pdf, ": ", conditionMessage(e)); FALSE
  })
  grDevices::dev.off()
  if (dibujado && file.exists(path_pdf) && file.size(path_pdf) > 0) {
    message("  (png no disponible en este sistema; escrito como pdf: ", path_pdf, ")")
    return(invisible(path_pdf))
  }
  warning("No se pudo generar ni png ni pdf para: ", path_png)
  invisible(NULL)
}

if (!file.exists(genes_path)) stop("no encuentro ", genes_path)
genes_map <- read.delim(genes_path, stringsAsFactors = FALSE)
genes_map$chr <- as.character(genes_map$chr)
needed_genes_cols <- c("gene_id", "gene_name", "chr", "grid_index")
if (!all(needed_genes_cols %in% colnames(genes_map))) {
  stop("genes.tsv no tiene las columnas esperadas (", paste(needed_genes_cols, collapse = ","), ")")
}

# --- condiciones disponibles, en orden biologico -----------------------------
inv_files <- list.files(invariants_dir, pattern = "^condition_invariants_.*\\.tsv$")
if (!length(inv_files)) stop("no encuentro condition_invariants_*.tsv en ", invariants_dir,
                              " -- corre primero la etapa consensus para esta cohorte")
condiciones_presentes <- sub("^condition_invariants_", "", sub("\\.tsv$", "", inv_files))
orden_conocido <- c("Normal_histology", "Control_disease_cohort", "Control_external_study",
                     "F0", "F1", "F2", "F3", "F4")
condiciones <- orden_conocido[orden_conocido %in% condiciones_presentes]
extra <- setdiff(condiciones_presentes, orden_conocido)
if (length(extra)) condiciones <- c(condiciones, extra)
cat("Condiciones de", cohort, "en orden:", paste(condiciones, collapse = " -> "), "\n")

# --- cargar espectro completo (mediana power_normalised por chr,N,k) y ------
# --- fase real (Arg(mean(exp(i*phase)))) por chr,N,k, de UNA condicion ------
cargar_espectro_y_fase <- function(cond) {
  path <- file.path(spectra_dir, sprintf("spectra_samples_%s.tsv", cond))
  if (!file.exists(path)) { cat("  (no existe:", path, ")\n"); return(NULL) }
  d <- read.delim(path, stringsAsFactors = FALSE)
  need <- c("chr", "N", "k", "power_normalised", "phase")
  if (!all(need %in% colnames(d))) {
    stop(path, " no tiene las columnas esperadas (", paste(need, collapse = ","), ")")
  }
  d$chr <- as.character(d$chr)
  espectro <- aggregate(power_normalised ~ chr + N + k, data = d, FUN = median, na.rm = TRUE)
  fase_agg <- do.call(rbind, lapply(
    split(d[, c("chr", "N", "k", "phase")], list(d$chr, d$N, d$k), drop = TRUE),
    function(sub) {
      phi <- sub$phase[is.finite(sub$phase)]
      if (length(phi) < 2) return(NULL)
      data.frame(chr = sub$chr[1], N = sub$N[1], k = sub$k[1],
                 mean_phase = Arg(mean(exp(1i * phi))), stringsAsFactors = FALSE)
    }))
  list(espectro = espectro, fase = fase_agg)
}

cargar_picos <- function(cond) {
  path <- file.path(invariants_dir, sprintf("condition_invariants_%s.tsv", cond))
  if (!file.exists(path)) return(NULL)
  d <- read.delim(path, stringsAsFactors = FALSE)
  if (!nrow(d)) return(NULL)
  d$chr <- as.character(d$chr)
  d
}

# --- eje X concatenado por cromosoma -----------------------------------------
armar_offsets <- function(espectro) {
  chrs_presentes <- chr_order[chr_order %in% unique(espectro$chr)]
  offset <- 0; offset_map <- c()
  for (ch in chrs_presentes) {
    offset_map[ch] <- offset
    offset <- offset + max(espectro$k[espectro$chr == ch], na.rm = TRUE) + 15
  }
  list(offset_map = offset_map, chrs_presentes = chrs_presentes)
}

plot_espectro_uno <- function(espectro, picos, titulo, archivo) {
  arm <- armar_offsets(espectro)
  espectro <- espectro[order(match(espectro$chr, chr_order), espectro$k), ]
  espectro$x_pos <- espectro$k + arm$offset_map[espectro$chr]

  graficar_robusto(archivo, width = 2200, height = 700, res = 130, plot_fn = function() {
    plot(espectro$x_pos, pmax(espectro$power_normalised, 1e-6), type = "l", col = "grey40",
         log = "y", xlab = "Cromosoma (k, ciclos por gen, concatenado)",
         ylab = "Potencia normalizada (escala log)", main = titulo, xaxt = "n", lwd = 0.6)

    if (!is.null(picos) && nrow(picos)) {
      key_e <- paste(espectro$chr, espectro$N, espectro$k)
      key_p <- paste(picos$chr, picos$N, picos$k)
      idx <- match(key_p, key_e)
      ok <- !is.na(idx)
      points(espectro$x_pos[idx[ok]], pmax(espectro$power_normalised[idx[ok]], 1e-6),
             col = "red", pch = 19, cex = 0.6)
    }
    mids <- vapply(arm$chrs_presentes, function(ch) mean(range(espectro$x_pos[espectro$chr == ch])), numeric(1))
    axis(1, at = mids, labels = arm$chrs_presentes, las = 2, cex.axis = 0.65)
  })
}

# --- reconstruccion de genes de un pico ---------------------------------------
reconstruir_pico <- function(chr, N, k, mean_phase) {
  g <- genes_map[genes_map$chr == chr, c("gene_id", "gene_name", "grid_index")]
  g <- g[order(g$grid_index), ]
  g$carga <- round(cos(2 * pi * k * (g$grid_index - 1) / N + mean_phase), 4)
  g$zona <- ifelse(abs(g$carga) < 0.2, "nodo",
             ifelse(g$carga >= 0.8, "cresta",
             ifelse(g$carga <= -0.8, "valle", "intermedio")))
  g
}

plot_reconstruccion_uno <- function(recon, titulo, archivo) {
  graficar_robusto(archivo, width = 2200, height = 750, res = 130, plot_fn = function() {
    plot(recon$grid_index, recon$carga, type = "l", col = "grey65", lwd = 0.7,
         xlab = "Posicion en el cromosoma (grid_index)", ylab = "Carga u(t)",
         main = titulo, ylim = c(-1.15, 1.15))
    abline(h = 0, lty = 2, col = "grey85")
    abline(h = c(0.8, -0.8, 0.2, -0.2), lty = 3, col = "grey85")
    # Sin nombres de gen en el plot -- con cientos/miles de genes por
    # cromosoma se amontonan y quedan ilegibles. El nombre de cada gen
    # sigue completo en los TSV (el _completo.tsv trae carga + zona; el
    # _nodo_cresta_valle.tsv los agrupa por zona), asi que no se pierde,
    # solo se separa donde se puede leer.
    col_map <- c(nodo = "red", cresta = "#C9A227", valle = "blue")
    for (z in names(col_map)) {
      sub <- recon[recon$zona == z, ]
      if (nrow(sub)) {
        points(sub$grid_index, sub$carga, col = col_map[[z]], pch = 19, cex = 0.5)
      }
    }
    legend("bottomright", legend = names(col_map), col = unname(col_map), pch = 19, bty = "n", cex = 0.8)
  })
}

# El TSV completo YA trae carga (numerica) y zona (nodo/cresta/valle/
# intermedio) lado a lado -- "intermedio" solo se entiende bien junto al
# numero real, por eso van juntas. El TSV de 3 columnas ahora tambien lleva
# la carga de cada gen al lado de su nombre, no solo el nombre.
escribir_tsvs_reconstruccion <- function(recon, prefijo) {
  cols_completo <- c("gene_id", "gene_name", "grid_index", "carga", "zona")
  write.table(recon[, cols_completo], paste0(prefijo, "_completo.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  armar_columna <- function(zona_nombre) {
    sub <- recon[recon$zona == zona_nombre, c("gene_name", "carga")]
    list(gene = sub$gene_name, carga = sub$carga)
  }
  nodo <- armar_columna("nodo"); cresta <- armar_columna("cresta"); valle <- armar_columna("valle")
  n <- max(length(nodo$gene), length(cresta$gene), length(valle$gene), 1)
  pad <- function(x, n) { length(x) <- n; x }
  tabla <- data.frame(orden = seq_len(n),
                       nodo = pad(nodo$gene, n), nodo_carga = pad(nodo$carga, n),
                       cresta = pad(cresta$gene, n), cresta_carga = pad(cresta$carga, n),
                       valle = pad(valle$gene, n), valle_carga = pad(valle$carga, n))
  write.table(tabla, paste0(prefijo, "_nodo_cresta_valle.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

# procesa un conjunto de picos, usando la tabla de fase de UNA condicion
# especifica (la condicion en la que fueron detectados, o para restas, cada
# condicion del par por separado). out_subdir: carpeta donde van el PNG y
# los 2 TSV de cada pico (se crea si no existe).
procesar_picos <- function(picos, fase_de_cond, subdir_prefix, out_subdir) {
  if (is.null(picos) || !nrow(picos)) return(invisible(NULL))
  dir.create(out_subdir, showWarnings = FALSE, recursive = TRUE)
  picos_u <- unique(picos[, c("chr", "N", "k")])
  for (i in seq_len(nrow(picos_u))) {
    chr <- picos_u$chr[i]; N <- picos_u$N[i]; k <- picos_u$k[i]
    fr <- fase_de_cond[fase_de_cond$chr == chr & fase_de_cond$N == N & fase_de_cond$k == k, ]
    if (!nrow(fr)) { cat("  (sin fase para chr", chr, "k", k, ", se omite reconstruccion)\n"); next }
    recon <- reconstruir_pico(chr, N, k, fr$mean_phase[1])
    tag <- sprintf("chr%s_N%s_k%s", chr, N, k)
    # subdir_prefix va en el NOMBRE del archivo, no solo en el titulo del
    # plot: dos llamadas a procesar_picos() pueden compartir out_subdir (las
    # reconstrucciones de fase-A y fase-B de una resta van a la misma
    # carpeta) y sin esto se sobrescribirian entre si.
    archivo_base <- paste0(subdir_prefix, "_", tag)
    plot_reconstruccion_uno(recon, sprintf("%s -- %s (periodo=%.1f genes)", subdir_prefix, tag, N / k),
                             file.path(out_subdir, paste0(archivo_base, ".png")))
    escribir_tsvs_reconstruccion(recon, file.path(out_subdir, archivo_base))
  }
}

# ============================== 1. espectros por condicion ===================
espectros <- list(); fases <- list(); picos_por_cond <- list()
for (cond in condiciones) {
  cat("Cargando espectro + fase:", cond, "\n")
  ef <- cargar_espectro_y_fase(cond)
  espectros[[cond]] <- if (!is.null(ef)) ef$espectro else NULL
  fases[[cond]] <- if (!is.null(ef)) ef$fase else NULL
  picos_por_cond[[cond]] <- cargar_picos(cond)
  if (!is.null(espectros[[cond]])) {
    plot_espectro_uno(espectros[[cond]], picos_por_cond[[cond]],
                       sprintf("%s -- %s", cohort, cond),
                       file.path(opt[["--out-dir"]], "espectros", sprintf("espectro_%s.png", cond)))
    procesar_picos(picos_por_cond[[cond]], fases[[cond]], cond,
                   out_subdir = file.path(opt[["--out-dir"]], "espectros", sprintf("reconstrucciones_%s", cond)))
  }
}

# ============================== 2. restas: consecutivas + extremo ===========
pares <- list()
for (i in seq_len(length(condiciones) - 1)) pares[[length(pares) + 1]] <- c(condiciones[i], condiciones[i + 1])
if (length(condiciones) >= 2) pares[[length(pares) + 1]] <- c(condiciones[1], condiciones[length(condiciones)])

for (par in pares) {
  a <- par[1]; b <- par[2]
  if (is.null(espectros[[a]]) || is.null(espectros[[b]])) next
  cat("Resta:", a, "vs", b, "\n")
  tag <- paste0(a, "_vs_", b)

  # offsets desde la UNION de cromosomas de ambas condiciones -- si a le
  # falta un cromosoma que b si tiene (o viceversa), el offset de ese
  # cromosoma debe existir igual para no perder esos puntos ni desalinear
  # el eje entre las dos curvas.
  espectro_union_chr <- rbind(espectros[[a]][, c("chr", "N", "k")], espectros[[b]][, c("chr", "N", "k")])
  arm <- armar_offsets(espectro_union_chr)
  # ambas curvas ORDENADAS por posicion en el eje X antes de graficar --
  # lines() conecta los puntos en el orden de las filas, no por valor de X,
  # asi que sin esto la curva sale en zigzag ("espagueti").
  ea <- espectros[[a]][order(match(espectros[[a]]$chr, chr_order), espectros[[a]]$k), ]
  eb <- espectros[[b]][order(match(espectros[[b]]$chr, chr_order), espectros[[b]]$k), ]
  ea$x_pos <- ea$k + arm$offset_map[ea$chr]
  eb$x_pos <- eb$k + arm$offset_map[eb$chr]

  picos_a <- picos_por_cond[[a]]; picos_b <- picos_por_cond[[b]]
  picos_union <- unique(rbind(picos_a, picos_b)[, c("chr", "N", "k")])

  graficar_robusto(file.path(opt[["--out-dir"]], "restas", sprintf("resta_%s.png", tag)),
      width = 2200, height = 700, res = 130, plot_fn = function() {
    plot(ea$x_pos, pmax(ea$power_normalised, 1e-6), type = "l", col = "steelblue", lwd = 0.6,
         log = "y", xlab = "Cromosoma (k, ciclos por gen, concatenado)",
         ylab = "Potencia normalizada (escala log)",
         main = sprintf("%s -- %s vs %s (superpuestos)", cohort, a, b), xaxt = "n")
    lines(eb$x_pos, pmax(eb$power_normalised, 1e-6), col = "firebrick", lwd = 0.6)
    if (nrow(picos_union)) {
      key_e <- paste(ea$chr, ea$N, ea$k); key_p <- paste(picos_union$chr, picos_union$N, picos_union$k)
      idx <- match(key_p, key_e); ok <- !is.na(idx)
      points(ea$x_pos[idx[ok]], pmax(ea$power_normalised[idx[ok]], 1e-6), col = "black", pch = 19, cex = 0.6)
    }
    mids <- vapply(arm$chrs_presentes, function(ch) mean(range(ea$x_pos[ea$chr == ch])), numeric(1))
    axis(1, at = mids, labels = arm$chrs_presentes, las = 2, cex.axis = 0.65)
    legend("topright", legend = c(a, b, "pico en A o B"), col = c("steelblue", "firebrick", "black"),
           lwd = c(1, 1, NA), pch = c(NA, NA, 19), bty = "n", cex = 0.8)
  })

  # reconstruccion con la fase de CADA condicion por separado (el mismo pico
  # puede tener una reconstruccion ligeramente distinta en A y en B)
  out_sub <- file.path(opt[["--out-dir"]], "restas", sprintf("reconstrucciones_%s", tag))
  procesar_picos(picos_union, fases[[a]], paste0(tag, "__fase_", a), out_subdir = out_sub)
  procesar_picos(picos_union, fases[[b]], paste0(tag, "__fase_", b), out_subdir = out_sub)
}

cat("\nListo. Salidas en:", opt[["--out-dir"]], "\n")
