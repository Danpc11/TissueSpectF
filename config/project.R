# Project-level paths and constants. Override via environment or a local copy.
list(
  geo_dir     = Sys.getenv("LFFT_GEO_DIR",
                           "/scratch/home/storres/tejidos_fft/liver_fft/GEO"),
  interim_dir = Sys.getenv("LFFT_INTERIM_DIR", "data/interim"),
  results_dir = Sys.getenv("LFFT_RESULTS_DIR", "results"),

  annotation_file = "Human.GRCh38.p13.annot.tsv.gz",

  # Genes are ordered by start position within a chromosome; a chromosome with
  # fewer than this many genes cannot support a meaningful spectrum.
  min_genes_per_chr = 8L,

  # Expression filter applied before the FFT stage.
  min_tpm      = 1,
  min_fraction = 0.2,

  chrom_levels = c(as.character(1:22), "X", "Y"),

  # maxT permutation settings (kept here so both datasets provably share them).
  maxt = list(
    B            = 1000L,
    seed         = 42L,
    block_sizes  = c(10L, 20L, 50L),
    alpha        = 0.05,
    stable_frac  = 0.9
  )
)
