#!/usr/bin/env Rscript
# gtex_subset.R -- 100 samples per tissue from GTEx, one per donor, from the
# full ~17,000-sample matrix without loading it whole.
#
#   Rscript scripts/gtex_subset.R --geo-dir data [--per-tissue 100] [--seed 42]
#
# Writes into --geo-dir:
#   GTEx_8tissue_100each_reads.tsv.gz   genes x selected samples
#   GTEx_8tissue_100each_pheno.tsv      sample_id, tissue, donor
#
# WHY A SUBSET
# ------------
# The full gene_reads GCT is several GB and ~17,000 samples. The question this
# is for -- can the spectrum tell eight maximally different tissues apart --
# does not need all of it, and answering it on 800 samples first means a
# failure costs an afternoon instead of a week.
#
# ONE SAMPLE PER DONOR PER TISSUE
# -------------------------------
# GTEx takes many tissues from each of ~980 donors. Two samples from one donor
# in the same tissue are not two observations of that tissue; and across the
# train/test split, a donor present on both sides lets a model score well by
# learning donor identity. Sampling one per donor removes the within-tissue
# case here, and `donor` is written out so the splitter can handle the
# across-tissue case.
#
# STREAMWISE
# ----------
# The GCT is read in chunks and only the selected columns are kept, so peak
# memory is one chunk rather than the whole matrix.

suppressWarnings({
  source("R/utils_io.R")
  source("R/config.R")
})

args <- commandArgs(trailingOnly = TRUE)
flag <- function(name, default = NULL) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (length(hit)) return(sub(paste0("^", name, "="), "", hit[1]))
  i <- match(name, args)
  if (!is.na(i) && length(args) > i) args[i + 1] else default
}

geo_dir     <- flag("--geo-dir", Sys.getenv("TSF_GEO_DIR", ""))
per_tissue  <- as.integer(flag("--per-tissue", "100"))
seed        <- as.integer(flag("--seed", "42"))
gct         <- flag("--gct",
                    file.path(geo_dir,
                              "GTEx_Analysis_v10_RNASeQCv2.4.2_gene_reads.gct.gz"))
attrs       <- flag("--attributes",
                    file.path(geo_dir,
                              "GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"))
chunk_rows  <- as.integer(flag("--chunk-rows", "5000"))

if (!nzchar(geo_dir)) tsf_abort("Pass --geo-dir <dir>.")
for (f in c(gct, attrs)) if (!file.exists(f)) tsf_abort("No such file: ", f)

# SMTSD is the detailed site ("Brain - Cortex", "Colon - Transverse"); SMTS is
# the broad one ("Brain", "Colon"). The broad column is used, because eight
# maximally different tissues is the question and Brain subregions are not
# eight different tissues.
#
# Intestine is not an SMTS value: GTEx splits the gut into Colon and Small
# Intestine. Both are mapped to Intestine here, which is a decision to state
# rather than a lookup -- a colon sample and a terminal-ileum sample are not
# interchangeable, and pooling them widens that class relative to the others.
TISSUE_MAP <- c(
  "Lung"            = "Lung",
  "Kidney"          = "Kidney",
  "Brain"           = "Brain",
  "Blood"           = "Whole_Blood",
  "Liver"           = "Liver",
  "Muscle"          = "Muscle",
  "Pancreas"        = "Pancreas",
  "Colon"           = "Intestine",
  "Small Intestine" = "Intestine"
)

tsf_log("Reading sample attributes: ", basename(attrs))
sa <- read.delim(attrs, sep = "\t", header = TRUE, quote = "",
                 stringsAsFactors = FALSE, check.names = FALSE)
need <- c("SAMPID", "SMTS")
miss <- setdiff(need, names(sa))
if (length(miss)) {
  tsf_abort(basename(attrs), " lacks: ", paste(miss, collapse = ", "),
            ". Columns present: ", paste(head(names(sa), 12), collapse = ", "))
}

sa$tissue <- unname(TISSUE_MAP[sa$SMTS])
sa <- sa[!is.na(sa$tissue), , drop = FALSE]

# Donor is the first two hyphen-separated fields of SAMPID: GTEX-1117F-0226-...
sa$donor <- sub("^([^-]+-[^-]+).*$", "\\1", sa$SAMPID)

# Whole blood only, when SMTSD is available: SMTS "Blood" also covers cell
# lines, which are not a tissue.
if ("SMTSD" %in% names(sa)) {
  drop <- sa$tissue == "Whole_Blood" & !grepl("whole blood", sa$SMTSD,
                                              ignore.case = TRUE)
  if (any(drop)) {
    tsf_log("  dropping ", sum(drop), " non-whole-blood sample(s) from Blood")
    sa <- sa[!drop, , drop = FALSE]
  }
}

tsf_log("Candidates per tissue:")
for (t in sort(unique(sa$tissue))) {
  n <- sum(sa$tissue == t)
  tsf_log(sprintf("  %-12s %5d sample(s), %5d donor(s)", t, n,
                  length(unique(sa$donor[sa$tissue == t]))))
}

set.seed(seed)
picked <- do.call(rbind, lapply(split(sa, sa$tissue), function(d) {
  # One per donor, then up to per_tissue of those.
  d <- d[!duplicated(d$donor), , drop = FALSE]
  if (nrow(d) > per_tissue) d <- d[sample.int(nrow(d), per_tissue), , drop = FALSE]
  d
}))
picked <- picked[order(picked$tissue, picked$SAMPID), , drop = FALSE]

short <- table(picked$tissue)[table(picked$tissue) < per_tissue]
if (length(short)) {
  tsf_warn("Below --per-tissue after one-per-donor: ",
           paste(sprintf("%s=%d", names(short), as.integer(short)),
                 collapse = ", "),
           ". Reported rather than topped up with a second sample from a donor ",
           "already present, which would not be a new observation of that tissue.")
}
tsf_log("Selected ", nrow(picked), " sample(s) over ",
        length(unique(picked$tissue)), " tissue(s)")

pheno_path <- file.path(geo_dir, "GTEx_8tissue_100each_pheno.tsv")
write_tsv_tsf(data.frame(sample_id = picked$SAMPID,
                         tissue    = picked$tissue,
                         donor     = picked$donor,
                         stringsAsFactors = FALSE), pheno_path)
tsf_log("Wrote ", pheno_path)

# --- the counts, in chunks ----------------------------------------------------
# GCT: line 1 is the version, line 2 is "nrows ncols", line 3 is the header.
con <- gzfile(gct, "rt")
on.exit(close(con), add = TRUE)
invisible(readLines(con, n = 2L))
header <- strsplit(readLines(con, n = 1L), "\t", fixed = TRUE)[[1]]

keep_cols <- match(picked$SAMPID, header)
absent <- picked$SAMPID[is.na(keep_cols)]
if (length(absent)) {
  tsf_abort(length(absent), " selected sample(s) are not columns of the GCT: ",
            paste(head(absent, 5), collapse = ", "),
            ". The attributes file and the expression release must be the same ",
            "GTEx version.")
}
id_col <- match("Name", header)
if (is.na(id_col)) tsf_abort("The GCT has no 'Name' column.")

out_path <- file.path(geo_dir, "GTEx_8tissue_100each_reads.tsv.gz")
out <- gzfile(out_path, "wt")
writeLines(paste(c("Name", picked$SAMPID), collapse = "\t"), out)

n_genes <- 0L
repeat {
  lines <- readLines(con, n = chunk_rows)
  if (!length(lines)) break
  f <- strsplit(lines, "\t", fixed = TRUE)
  writeLines(vapply(f, function(x) paste(c(x[id_col], x[keep_cols]),
                                         collapse = "\t"), character(1)), out)
  n_genes <- n_genes + length(lines)
  if (n_genes %% (chunk_rows * 4L) == 0L) tsf_log("  ", n_genes, " genes...")
}
close(out)
tsf_log("Wrote ", out_path, " (", n_genes, " genes x ", nrow(picked), " samples)")
tsf_log("Next: ./tsf ingest GTEx --geo-dir ", geo_dir, " ...")
