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

  # CLEAN decomposition: greedy deflation with an extended-BIC stopping rule.
  # No number of components is chosen -- EBIC decides per chromosome. gamma
  # scales the cost of searching the frequency grid; 1 is strict and is what
  # keeps pure noise from yielding components. Lower it only as a declared
  # sensitivity analysis.
  clean = list(
    ebic_gamma     = as.numeric(Sys.getenv("TSF_EBIC_GAMMA", "1")),
    penalty_factor = 1,
    max_components = 20L,
    per_sample     = FALSE   # TRUE gives one fingerprint per sample (slower)
  ),

  # What decides which peaks go downstream.
  #   "condition"   condition-level permutation test, q <= 0.05  (default)
  #   "consistency" >= stable_frac of samples individually significant
  # Both are always computed and written; this only picks which one drives
  # is_stable. See the note at the top of R/condition_test.R.
  stability_criterion = Sys.getenv("TSF_STABILITY_CRITERION", "condition"),

  # maxT permutation settings (kept here so both datasets provably share them).
  maxt = list(
    # TSF_MAXT_B lets a smoke run use a small B; production runs leave it unset.
    B            = as.integer(Sys.getenv("TSF_MAXT_B", "1000")),
    seed         = 42L,
    block_sizes  = c(10L, 20L, 50L),
    # Which null decides significance. "full" permutes every observed value and
    # so destroys local autocorrelation as well as long-range structure; "all"
    # additionally requires the peak to survive the block schemes, i.e. to be
    # more than local correlation. Use "all" for any claim about periodicity.
    primary_scheme = Sys.getenv("TSF_PRIMARY_SCHEME", "full"),
    alpha        = 0.05,
    stable_frac  = as.numeric(Sys.getenv("TSF_STABLE_FRAC", "0.9")),
    # Permutations for the condition-level test. It runs once per condition
    # instead of once per sample, so a larger B is affordable and gives a
    # smaller attainable p-value (the floor is 1/(B+1)).
    condition_B  = as.integer(Sys.getenv("TSF_CONDITION_B", "2000"))
  )
)
