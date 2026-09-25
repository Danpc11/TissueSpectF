#!/usr/bin/env Rscript
# condition_spectrum_diff.R -- spectrum subtraction between two conditions,
# from their OWN invariants (condition_invariants_<condition>.tsv, the
# default output of stage_consensus: prevalence/median_power_normalised/plv
# per (chr, N, k), computed within each condition's own samples, plus a
# bootstrap CI and a null p-value, both descriptive -- see consensus.R's
# condition_invariant_bootstrap() and null_component_pvalues()).
#
# For each (chr, N, k) that reaches at least the lowest prevalence tier in
# EITHER condition, reports:
#   power_diff       = median_power_normalised(A) - median_power_normalised(B)
#   prevalence_diff  = prevalence(A) - prevalence(B)
#   plv_diff         = plv(A) - plv(B)
# Plain subtraction, not a ratio, so no epsilon is needed anywhere here.
#
# A component reaching the tier in only ONE condition has NA on the other
# side and `presence = "only_in_A"` / `"only_in_B"` -- it is NOT treated as
# zero power there, because condition_invariants_<cond>.tsv only contains
# components that cleared a prevalence floor in THAT condition; a component
# below the floor may still have real, just-not-invariant power, which this
# script does not know and will not invent.
#
# RELIABILITY OF THE DIFFERENCE. (chr, N, k) identity is exact -- the same
# frequency means the same thing in both files, no fuzzy matching needed.
# What is NOT automatic is whether the difference is bigger than resampling
# noise: `power_diff` alone is a point estimate. When both sides carry a
# bootstrap CI (median_power_ci_lower/upper, from condition_invariant_
# bootstrap()), this script adds `power_ci_overlap`: FALSE when condition A's
# and B's power intervals do not touch -- evidence the difference survives
# resampling of each condition's own patients, not just the specific ones
# observed. TRUE means the intervals overlap: the point difference could
# plausibly be resampling noise. NA when a CI is missing on one side (e.g.
# the file predates condition_invariant_bootstrap(), or that condition had
# too few samples to bootstrap) or the component is only_in_A/only_in_B (no
# second interval to compare against). This is NOT a significance test and
# uses no PLV threshold of its own -- plv_A/plv_B are reported alongside so
# a weakly phase-locked component can be judged separately, not silently
# folded into one pass/fail flag nobody asked for.
#
# Usage:
#   Rscript condition_spectrum_diff.R \
#     --a results/.../consensus/condition_invariants_F0.tsv --a-label F0 \
#     --b results/.../consensus/condition_invariants_F4.tsv --b-label F4 \
#     --out condition_spectrum_diff_F0_vs_F4.tsv

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
opt <- list(`--a-label` = "A", `--b-label` = "B", `--out` = "condition_spectrum_diff.tsv")
i <- 1
while (i <= length(args)) { opt[[args[i]]] <- args[i + 1]; i <- i + 2 }
need <- c("--a", "--b")
miss <- setdiff(need, names(opt))
if (length(miss)) stop("faltan argumentos: ", paste(miss, collapse = ", "))

# Optional: only added to the output, and to the reliability check, when
# BOTH files have them. Absent gracefully rather than required, since an
# older condition_invariants_<cond>.tsv (predating condition_invariant_
# bootstrap()) or a condition too small to bootstrap will not carry them.
ci_cols <- c("median_power_ci_lower", "median_power_ci_upper",
             "plv_ci_lower", "plv_ci_upper", "p_null_fwer")

read_condition_invariants <- function(path) {
  d <- read.delim(path, stringsAsFactors = FALSE)
  needed <- c("chr", "N", "k", "period", "prevalence",
              "median_power_normalised", "plv", "condition_invariant_class")
  missing <- setdiff(needed, colnames(d))
  if (length(missing)) {
    stop(path, " no tiene las columnas esperadas (", paste(missing, collapse = ","),
         "); columnas presentes: ", paste(colnames(d), collapse = ","))
  }
  d$chr <- as.character(d$chr)
  d
}

a <- read_condition_invariants(opt[["--a"]])
b <- read_condition_invariants(opt[["--b"]])
lab_a <- opt[["--a-label"]]; lab_b <- opt[["--b-label"]]

has_ci <- all(ci_cols %in% colnames(a)) && all(ci_cols %in% colnames(b))
if (!has_ci) {
  cat("Aviso: no estan las columnas de intervalo de confianza (", paste(ci_cols, collapse=","),
      ") en ambos archivos -- power_ci_overlap saldra vacio. Corre condition_invariant_bootstrap()",
      "(ya integrado por default en stage_consensus) para tenerlas.\n", sep = "")
}

base_cols <- c("chr", "N", "k", "period", "prevalence", "median_power_normalised", "plv")
keep_cols <- if (has_ci) c(base_cols, ci_cols) else base_cols
a <- a[, keep_cols]; b <- b[, keep_cols]
colnames(a)[5:length(keep_cols)] <- paste0(colnames(a)[5:length(keep_cols)], "_", lab_a)
colnames(b)[5:length(keep_cols)] <- paste0(colnames(b)[5:length(keep_cols)], "_", lab_b)

out <- merge(a, b, by = c("chr", "N", "k", "period"), all = TRUE)

pa <- out[[paste0("median_power_normalised_", lab_a)]]
pb <- out[[paste0("median_power_normalised_", lab_b)]]
prev_a <- out[[paste0("prevalence_", lab_a)]]
prev_b <- out[[paste0("prevalence_", lab_b)]]
plv_a <- out[[paste0("plv_", lab_a)]]
plv_b <- out[[paste0("plv_", lab_b)]]

out$power_diff <- pa - pb
out$prevalence_diff <- prev_a - prev_b
out$plv_diff <- plv_a - plv_b

out$presence <- ifelse(is.na(pa), paste0("only_in_", lab_b),
                 ifelse(is.na(pb), paste0("only_in_", lab_a), "both"))

if (has_ci) {
  lo_a <- out[[paste0("median_power_ci_lower_", lab_a)]]
  hi_a <- out[[paste0("median_power_ci_upper_", lab_a)]]
  lo_b <- out[[paste0("median_power_ci_lower_", lab_b)]]
  hi_b <- out[[paste0("median_power_ci_upper_", lab_b)]]
  # No overlap <=> one interval's upper bound sits below the other's lower
  # bound. FALSE (no overlap) is the informative case: the power difference
  # holds up under resampling of each condition's own patients separately,
  # not just a coincidence of the specific ones observed.
  both_present <- out$presence == "both" & is.finite(lo_a) & is.finite(hi_a) &
                   is.finite(lo_b) & is.finite(hi_b)
  out$power_ci_overlap <- NA
  out$power_ci_overlap[both_present] <-
    !(hi_a[both_present] < lo_b[both_present] | hi_b[both_present] < lo_a[both_present])
} else {
  out$power_ci_overlap <- NA
}

out <- out[order(-abs(ifelse(is.na(out$power_diff), 0, out$power_diff))), ]

write.table(out, opt[["--out"]], sep = "\t", quote = FALSE, row.names = FALSE, na = "")
n_no_overlap <- sum(!out$power_ci_overlap, na.rm = TRUE)
cat("Filas:", nrow(out),
    " | en ambas condiciones:", sum(out$presence == "both"),
    " | solo en", lab_a, ":", sum(out$presence == paste0("only_in_", lab_a)),
    " | solo en", lab_b, ":", sum(out$presence == paste0("only_in_", lab_b)), "\n")
if (has_ci) {
  cat("De los presentes en ambas, con intervalos de potencia que NO se solapan (diferencia mas solida): ",
      n_no_overlap, "\n", sep = "")
}
cat("Escrito:", opt[["--out"]], "\n")
