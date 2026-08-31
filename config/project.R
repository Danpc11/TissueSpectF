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

  # --- annotation ------------------------------------------------------------
  # Two formats are supported and they name the biotype field differently:
  #
  #   ncbi   GeneType   "protein-coding"   (hyphen)
  #   gtf    gene_type  "protein_coding"   (underscore)
  #
  # gene_universe is a regular expression matched against that field, so a
  # pattern written for one source matches nothing in the other. Ingest aborts
  # with the right pattern rather than building an empty or accidental grid.
  #
  # GENCODE:
  #   annotation_file   = "gencode.v50.basic.annotation.gtf.gz"
  #   annotation_format = "gtf"
  #   gene_universe     = "^protein_coding$"
  #   annotation_release = "GENCODE v50 basic"
  annotation_file   = "Human.GRCh38.p13.annot.tsv.gz",
  annotation_format = Sys.getenv("TSF_ANNOTATION_FORMAT", "ncbi"),

  # Exonic union length, which is what TPM needs. "span" (end - start) is not a
  # transcript length and inflates long genes; it exists only for annotations
  # with no exon features.
  gene_length_mode  = "exonic",
  strip_gene_version = TRUE,

  # A GENCODE GTF carries Ensembl ids and symbols but no Entrez ids, and three
  # of these cohorts publish counts keyed on Entrez. The NCBI annotation table
  # already in the data directory carries both side by side, so it doubles as
  # the mapping. Set to NULL when the annotation already has what the counts use.
  id_map = list(file = "Human.GRCh38.p13.annot.tsv.gz",
                ensembl_column = "EnsemblGeneID",
                entrez_column  = "GeneID"),

  # Provenance of the annotation. Recorded in every output and checked before
  # two datasets are allowed into the same reference: a feature named chrX_k7
  # means a different thing under a different build or gene universe.
  species            = "Homo sapiens",
  genome_build       = "GRCh38.p13",
  annotation_release = "NCBI",

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

  # Consensus spectrum: how strong, how common and how phase-aligned each
  # frequency is across the samples of a condition. See R/consensus.R for why
  # this is not the spectrum of the mean profile.
  consensus = list(
    n_boot         = 500L,
    quantile_cut   = 0.95,   # "stands out" cut when no maxT is available
    min_prevalence = 0.5,
    plv_q          = 0.05,   # BH-adjusted Rayleigh p for phase alignment
    # Draws of n random samples, ignoring condition, whose best consensus score
    # forms the null a component must beat to be called confirmed.
    # Cost is one consensus spectrum per draw, cached per sample size, so
    # conditions of equal size share a null.
    n_null         = 50L,
    null_q         = 0.05,   # family-wise p against the permuted null
    # Column of samples.tsv identifying non-independent samples (subject, batch,
    # tumour-normal pair). When set, the null draws whole blocks. NULL treats
    # samples as independent, which is only right when they are.
    permutation_block = NULL,
    max_components = 50L
  ),

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

  # Fingerprints for the matcher. These are NOT the significant peaks: matching
  # does not require any component to be individually significant, and filtering
  # by q first would discard what makes classes separable.
  fingerprint = list(
    k_max      = 64L,          # frequency indices kept per chromosome
    features   = "amplitude",  # or "amplitude_phase" to keep the phase
    n_features = 500L,         # features the centroid model selects, on train only
    target     = "condition",  # or "tissue" for a cross-tissue reference
    # Coverage calibration. Loss is simulated on GENES of the grid and the
    # fingerprint is recomputed, because that is what a real query loses;
    # masking spectral features instead would measure an easier, wrong quantity.
    # Cost is one GLS fingerprint per (sample, level, mode, mask), hence the cap.
    n_masks              = 10L,
    max_queries_per_mask = 25L,
    # Which per-band threshold is applied.
    #   "pooled"       the quantile over all masks of a band. By construction it
    #                  rejects about `quantile_correct` of true members.
    #   "conservative" the 90th percentile of the per-mask thresholds, so an
    #                  unlucky pattern of missing regions is not judged against
    #                  a lucky mask.
    # Both are computed and both rejection rates are reported. The default is
    # "pooled" because with gene-level masking the conservative threshold was
    # measured to reject 46-79% of true members, which trades far too much
    # sensitivity for its extra safety. Switch only with that number in view.
    threshold_policy     = "pooled"
  ),

  # What decides which peaks go downstream.
  #   "condition"       family-wise (maxT), very strict: a component has to
  #                     dominate its chromosome. Expect single digits.
  #   "condition_fdr"   pointwise p with BH across frequencies. Answers "which
  #                     frequencies carry structure" rather than "which is the
  #                     strongest", and selects far more.
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
