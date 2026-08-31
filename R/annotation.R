# annotation.R -- read a gene annotation into the one shape the pipeline uses.
#
# Every importer produces the same columns, so nothing downstream knows or cares
# where the annotation came from:
#
#   gene_id gene_name entrez_id chr start end strand gene_length gene_type
#
# BIOTYPE NAMES DIFFER BETWEEN SOURCES, AND THE DIFFERENCE IS SILENT
#
#   GENCODE / Ensembl GTF   gene_type "protein_coding"    underscore
#   NCBI annotation table   GeneType  "protein-coding"    hyphen
#
# `gene_universe` is a regular expression matched against that field, so a
# pattern written for one source quietly matches nothing in the other -- and
# "nothing" means an empty grid or, worse, a grid built from whatever else
# happened to match. This is exactly how an earlier run ended up with a grid of
# pure ncRNA: "^PROTEIN_CODING$" matched no NCBI value while "NCRNA" matched
# theirs. Each importer therefore reports the biotypes it found, and
# `default_gene_universe()` gives the right pattern per format.

#' The protein-coding pattern for a given annotation format.
default_gene_universe <- function(format = c("ncbi", "gtf", "gff3")) {
  switch(match.arg(format),
         ncbi = "^protein-coding$",
         "^protein_coding$")
}

#' Parse the attribute column of a GTF/GFF3 line set.
#'
#' GTF is `key "value";` and GFF3 is `key=value;`. Only the keys asked for are
#' extracted, because the attribute column is long and pulling all of it for
#' 60,000 genes is most of the runtime.
parse_gtf_attributes <- function(attr, keys, gff3 = FALSE) {
  out <- lapply(keys, function(k) {
    pat <- if (gff3) paste0("(^|;)\\s*", k, "=([^;]*)") else
      paste0("(^|;)\\s*", k, ' "([^"]*)"')
    m <- regmatches(attr, regexpr(pat, attr, perl = TRUE))
    v <- rep(NA_character_, length(attr))
    hit <- vapply(regexpr(pat, attr, perl = TRUE), function(x) x > 0, logical(1))
    if (any(hit)) {
      v[hit] <- sub(pat, "\\2", m, perl = TRUE)
    }
    trimws(v)
  })
  stats::setNames(as.data.frame(out, stringsAsFactors = FALSE), keys)
}

#' Read a GENCODE/Ensembl GTF (or GFF3) into the common annotation shape.
#'
#' @param length_mode "exonic" sums the union of exon lengths, which is what TPM
#'   needs; "span" uses end - start + 1, which is not a transcript length and
#'   will inflate long genes. The choice is recorded in the provenance so a
#'   later TPM cannot be misread.
read_gtf_annotation <- function(path, chrom_levels, strip_version = TRUE,
                                length_mode = c("exonic", "span"),
                                gff3 = grepl("\\.gff3?(\\.gz)?$", path)) {
  length_mode <- match.arg(length_mode)
  if (!file.exists(path)) tsf_abort("No annotation at ", path)
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)

  tsf_log("Reading ", basename(path), " (", if (gff3) "GFF3" else "GTF", ", ",
          length_mode, " length)")
  raw <- utils::read.delim(con, header = FALSE, comment.char = "#", quote = "",
                           stringsAsFactors = FALSE,
                           col.names = c("seqname", "source", "feature", "start",
                                         "end", "score", "strand", "frame", "attr"))
  keep_feature <- if (length_mode == "exonic") c("gene", "exon") else "gene"
  raw <- raw[raw$feature %in% keep_feature, , drop = FALSE]
  if (!nrow(raw)) tsf_abort("No gene features in ", basename(path))

  keys <- c("gene_id", "gene_name", "gene_type")
  attrs <- parse_gtf_attributes(raw$attr, keys, gff3)
  # Ensembl's own GTF calls it gene_biotype; GENCODE calls it gene_type.
  if (all(is.na(attrs$gene_type))) {
    attrs$gene_type <- parse_gtf_attributes(raw$attr, "gene_biotype", gff3)[[1]]
  }
  raw <- cbind(raw[, c("seqname", "feature", "start", "end", "strand")], attrs)
  raw$chr <- normalise_chrom_names(raw$seqname)
  raw$gene_id_raw <- raw$gene_id
  if (strip_version) raw$gene_id <- sub("\\..*$", "", raw$gene_id)

  # PAR_Y copies duplicate a gene id on the Y chromosome. Keeping both puts one
  # gene at two grid positions; dropping them keeps the X copy, which is what
  # expression is quantified against in almost every pipeline.
  #
  # They are removed from EVERY feature, not only from the gene rows. Stripping
  # the version turns "ENSG...5_PAR_Y" into the same id as its X copy, so a
  # PAR_Y exon left behind is summed into the X gene's exonic length -- silently,
  # and only in the length, which nothing downstream would flag.
  par_y <- grepl("_PAR_Y$", raw$gene_id_raw)
  if (any(par_y)) {
    tsf_log("Dropping ", sum(par_y & raw$feature == "gene"),
            " PAR_Y duplicate gene(s) and their features")
    raw <- raw[!par_y, , drop = FALSE]
  }

  genes <- raw[raw$feature == "gene", , drop = FALSE]

  gene_length <- if (length_mode == "exonic") {
    ex <- raw[raw$feature == "exon", , drop = FALSE]
    if (!nrow(ex)) {
      tsf_warn("No exon features found; falling back to the genomic span, ",
               "which is not a transcript length")
      genes$end - genes$start + 1L
    } else {
      # Union of exons per gene, so overlapping transcripts are not counted twice.
      by_gene <- split(seq_len(nrow(ex)), ex$gene_id)
      union_len <- vapply(by_gene, function(i) {
        s <- ex$start[i]; e <- ex$end[i]
        o <- order(s); s <- s[o]; e <- e[o]
        total <- 0L; cur_s <- s[1]; cur_e <- e[1]
        for (j in seq_along(s)[-1]) {
          if (s[j] <= cur_e + 1L) { cur_e <- max(cur_e, e[j]) }
          else { total <- total + (cur_e - cur_s + 1L); cur_s <- s[j]; cur_e <- e[j] }
        }
        as.integer(total + (cur_e - cur_s + 1L))
      }, integer(1))
      unname(union_len[genes$gene_id])
    }
  } else {
    as.integer(genes$end - genes$start + 1L)
  }

  out <- data.frame(
    gene_id = genes$gene_id, gene_name = genes$gene_name,
    entrez_id = NA_character_,
    chr = genes$chr, start = as.numeric(genes$start),
    end = as.numeric(genes$end), strand = genes$strand,
    gene_length = as.numeric(gene_length), gene_type = genes$gene_type,
    stringsAsFactors = FALSE)

  observed <- sort(unique(out$gene_type[!is.na(out$gene_type)]))
  tsf_log("Biotypes present (", length(observed), "): ",
          paste(utils::head(observed, 8), collapse = ", "),
          if (length(observed) > 8) ", ..." else "")

  off_grid <- !out$chr %in% chrom_levels
  if (any(off_grid)) {
    tsf_log("Dropping ", sum(off_grid), " gene(s) on scaffolds or chromosomes ",
            "outside the declared set")
    out <- out[!off_grid, , drop = FALSE]
  }
  dup <- duplicated(out$gene_id)
  if (any(dup)) {
    tsf_warn(sum(dup), " duplicated gene id(s) after version stripping; ",
             "keeping the first occurrence of each")
    out <- out[!dup, , drop = FALSE]
  }
  out <- out[!is.na(out$start), , drop = FALSE]
  tsf_log("Annotation: ", nrow(out), " genes with chromosome and start")
  attr(out, "length_mode") <- length_mode
  out
}

#' Add identifiers the annotation itself does not carry.
#'
#' A GENCODE GTF has Ensembl ids and symbols but no Entrez ids, and three of the
#' liver cohorts publish counts keyed on Entrez. Rather than guessing, the
#' project declares a mapping file; the NCBI annotation table already in the
#' data directory is one, since it carries GeneID and EnsemblGeneID side by side.
#'
#' Genes with no mapping keep NA. They are then simply unmatched by an
#' Entrez-keyed count table, which is the honest outcome: a gene the mapping
#' does not cover was not measured under that identifier.
attach_id_map <- function(ann, project) {
  spec <- project$id_map
  if (is.null(spec)) return(ann)
  path <- file.path(project$geo_dir, spec$file)
  if (!file.exists(path)) {
    tsf_abort("id_map file not found: ", path)
  }
  map <- read_tsv_tsf(path)
  from <- spec$ensembl_column %||% "EnsemblGeneID"
  to <- spec$entrez_column %||% "GeneID"
  for (cl in c(from, to)) {
    if (!cl %in% colnames(map)) {
      tsf_abort("id_map file has no column '", cl, "'. Columns: ",
                paste(utils::head(colnames(map), 8), collapse = ", "))
    }
  }
  key <- sub("\\..*$", "", trimws(as.character(map[[from]])))
  val <- trimws(as.character(map[[to]]))
  ok <- nzchar(key) & !is.na(key) & nzchar(val) & !is.na(val)
  lookup <- stats::setNames(val[ok], key[ok])
  lookup <- lookup[!duplicated(names(lookup))]

  ann$entrez_id <- unname(lookup[ann$gene_id])
  n <- sum(!is.na(ann$entrez_id))
  tsf_log("Identifier map: ", n, "/", nrow(ann), " genes gained an Entrez id (",
          round(100 * n / nrow(ann)), "%) from ", basename(path))
  if (n == 0L) {
    tsf_abort("The identifier map produced no matches. Check that '", from,
              "' holds Ensembl gene ids in ", basename(path), ".")
  }
  ann
}

# Parsing a GENCODE GTF takes about a minute and a half, and ingest reads the
# annotation once per dataset. Cached for the session, keyed on the file and its
# modification time so an edited annotation is never served from the cache.
.annotation_cache <- new.env(parent = emptyenv())

#' Read whatever annotation the project declares.
#'
#' Dispatches on `annotation_format`. Adding a source means adding a reader
#' here, not touching anything that consumes annotations.
read_annotation <- function(project) {
  path <- file.path(project$geo_dir, project$annotation_file)
  format <- project$annotation_format %||% "ncbi"

  cache_key <- paste(path, file.mtime(path), format,
                     project$gene_universe %||% "all",
                     project$gene_length_mode %||% "exonic",
                     paste(project$chrom_levels, collapse = ","), sep = "|")
  if (!is.null(.annotation_cache[[cache_key]])) {
    tsf_log("Annotation: reusing the parse from earlier in this run")
    return(.annotation_cache[[cache_key]])
  }

  ann <- switch(format,
    ncbi = read_gene_annotation(path, project$chrom_levels),
    gtf  = ,
    gff3 = read_gtf_annotation(path, project$chrom_levels,
                               strip_version = project$strip_gene_version %||% TRUE,
                               length_mode = project$gene_length_mode %||% "exonic",
                               gff3 = identical(format, "gff3")),
    tsf_abort("annotation_format must be ncbi, gtf or gff3; got '", format, "'"))

  ann <- attach_id_map(ann, project)

  # A universe pattern written for the other source matches nothing and would
  # produce an empty or nonsensical grid. Say so here, where the fix is obvious.
  universe <- project$gene_universe
  if (!is.null(universe) && "gene_type" %in% colnames(ann)) {
    if (!any(grepl(universe, ann$gene_type, ignore.case = TRUE))) {
      tsf_abort("gene_universe '", universe, "' matches no biotype in this ",
                "annotation. ", format, " uses values like '",
                paste(utils::head(sort(unique(ann$gene_type)), 3), collapse = "', '"),
                "'. For this format the protein-coding pattern is '",
                default_gene_universe(format), "'.")
    }
  }

  .annotation_cache[[cache_key]] <- ann
  ann
}
