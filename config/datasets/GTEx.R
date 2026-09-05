# GTEx v10 bulk RNA-seq, eight tissues, 100 samples each.
#
# WHAT THIS IS FOR
# ----------------
# Negatives. Every result the fibrosis reference produces today is bounded on
# one side only: the coverage calibration says how often a true member is
# wrongly rejected, and says nothing about how often a sample from OUTSIDE the
# domain is wrongly accepted, because no out-of-domain sample was in the
# validation. A query from kidney currently receives a fibrosis stage with a
# similarity score and no way to notice it does not belong.
#
# It is also the sanity floor. Telling liver from kidney is a far easier task
# than telling F1 from F2. If the spectrum cannot do the first, the second was
# never possible and the problem is the representation, not the biology. If it
# does the first at 90% and the second at 31%, that pair of numbers characterises
# the resolution of the method, which is a result rather than a disappointment.
#
# THE DESIGN THAT MAKES IT A VALIDATION
# -------------------------------------
# GTEx is ONE cohort, so leave-one-cohort-out does not apply within it. The
# validation is across sources instead: train tissue discrimination on GTEx,
# then ask whether a GTEx-trained model calls the 503 GEO liver biopsies
# "Liver". That crosses study, platform and sample type in one step, which is a
# stronger claim than any internal split of GTEx would support.
#
# DONOR LEAKAGE IS THE TRAP
# -------------------------
# GTEx samples ~980 donors and takes many tissues from each, so the same donor
# appears in most tissues. A model split by sample can learn donor identity
# instead of tissue and score well doing it. `donor_column` below makes the
# donor available to the splitter; scripts/gtex_subset.R additionally takes at
# most ONE sample per donor per tissue, so within a tissue no donor repeats.
#
# WHAT GTEX IS NOT
# ----------------
# Postmortem tissue with substantial ischaemic time. GTEx liver is not a biopsy.
# As a negative for the other seven tissues that is irrelevant; as a control
# class for a biopsy-based fibrosis scale it is a confounder, and it must NOT be
# merged into Controles. Use it to test transfer, not to enlarge the healthy
# class.
#
# INPUTS
# ------
# Both files are summary-level and public; no dbGaP application is needed.
#
#   GTEx_Analysis_v10_RNASeQCv2.4.2_gene_reads.gct.gz
#   GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt
#
# The GCT format carries two header lines before the table, hence skip = 2, and
# a `Description` column beside `Name`. Gene ids carry a version suffix
# (ENSG00000223972.5); read_counts() already strips it.
#
# Run scripts/gtex_subset.R first: the full matrix is ~17,000 samples and this
# config expects the 800-sample subset it writes.
list(
  id          = "GTEx",
  title       = "GTEx v10 bulk RNA-seq; eight tissues, donor-disjoint within tissue",
  vocabulary  = "tissue_atlas",

  # Written by scripts/gtex_subset.R, not downloaded.
  counts_file   = "GTEx_8tissue_100each_reads.tsv.gz",
  series_matrix = "GTEx_8tissue_100each_pheno.tsv",

  # The phenotype table the subset script writes: one row per sample, with the
  # tissue and the donor.
  phenotype_format = "tsv",
  sample_id_column = "sample_id",

  counts_spec = list(
    sep       = "\t",
    id_column = "Name"
  ),
  count_id_type = "ENSEMBL",

  # Counts, so the pipeline applies its own normalisation. The TPM release
  # would need expression_unit = "tpm" and would skip that step; keep the two
  # sources on the same footing as the GEO cohorts, which are counts.
  expression_unit = "counts",

  # One level per tissue. Unordered by construction -- liver is not between
  # lung and kidney -- so `compare` produces no transition tables.
  condition_levels = c("Lung", "Kidney", "Brain", "Whole_Blood",
                       "Liver", "Muscle", "Pancreas", "Intestine"),

  condition_rules = list(
    list(id = "tissue", type = "column_match", column = "tissue",
         values = c("Lung", "Kidney", "Brain", "Whole_Blood",
                    "Liver", "Muscle", "Pancreas", "Intestine"),
         assign = NULL)      # assign = NULL: the matched value IS the level
  ),

  # No control cohort: with tissue_atlas every level is its own healthy state
  # and there is no baseline to compare against.
  has_control_cohort = FALSE,

  # Available to the splitter so a donor never lands in both train and test.
  donor_column = "donor",

  declared_n = 800L,
  notes = paste(
    "Negatives and the tissue-level sanity floor, not extra healthy livers.",
    "Postmortem: do not merge GTEx Liver into Controles."
  )
)
