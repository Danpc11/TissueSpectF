# Project-level paths and constants.
#
# Defaults point at the INMEGEN scratch layout; override with environment
# variables (TSF_GEO_DIR, TSF_INTERIM_DIR, TSF_RESULTS_DIR) to run elsewhere.
list(
  # Raw GEO downloads: count tables, series matrices, NCBI annotation.
  geo_dir     = Sys.getenv("TSF_GEO_DIR",
                           "/scratch/home/dperez/GPIB/gene_notes/TissueSpectF/data"),

  # Common format written by the ingest stage. Kept outside the git tree.
  interim_dir = Sys.getenv("TSF_INTERIM_DIR",
                           "/scratch/home/dperez/GPIB/gene_notes/TissueSpectF/interim"),

  # Downstream spectral results.
  results_dir = Sys.getenv("TSF_RESULTS_DIR",
                           "/scratch/home/dperez/GPIB/gene_notes/TissueSpectF/results"),

  annotation_file = "Human.GRCh38.p13.annot.tsv.gz",

  # Which annotated genes define the spectral axis. Regular expression matched
  # against the annotation's gene_type. The grid is what makes N comparable
  # across datasets, so changing this changes what "period in genes" means --
  # run both and report both rather than picking one silently.
  # Anchor the pattern: the NCBI annotation writes biotypes as "protein-coding"
  # and "ncRNA" (hyphen, mixed case), so an unanchored "PROTEIN_CODING|NCRNA"
  # silently matches ncRNA only and builds a grid with no coding genes at all.
  #   "^protein-coding$"            coding only  (default)
  #   "^(protein-coding|ncRNA)$"    coding + non-coding
  #   NULL                          every annotated gene
  gene_universe = Sys.getenv("TSF_GENE_UNIVERSE", "^protein-coding$"),

  # Genes are ordered by start position within a chromosome; a chromosome with
  # fewer than this many genes cannot support a meaningful spectrum.
  min_genes_per_chr = 8L,

  # Expression filter applied before the FFT stage.
  min_tpm      = 1,
  min_fraction = 0.2,

  chrom_levels = c(as.character(1:22), "X", "Y"),

  # maxT permutation settings (kept here so both datasets provably share them).
  maxt = list(
    # TSF_MAXT_B lets a smoke run use a small B; production runs leave it unset.
    B            = as.integer(Sys.getenv("TSF_MAXT_B", "1000")),
    seed         = 42L,
    block_sizes  = c(10L, 20L, 50L),
    alpha        = 0.05,
    stable_frac  = as.numeric(Sys.getenv("TSF_STABLE_FRAC", "0.9"))
  )
)
