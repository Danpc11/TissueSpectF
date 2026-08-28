# TissueSpectF

Chromosome-ordered Fourier spectral analysis of tissue level transcriptomes.

This repository replaces two ~4,200-line per-dataset scripts with a layered
pipeline. Everything that was dataset-specific now lives in a config file;
everything else is shared code.

## Layout

```
config/
  project.R              paths, filters, maxT settings (shared by all datasets)
  datasets/<GSE>.R       one file per dataset: file names + label rules
R/
  utils_io.R             logging, TSV I/O, phenotype value cleaning, manifests
  config.R               config loading and validation
  labels.R               controlled vocabulary + label harmonisation engine
  ingest.R               GEO -> common format
  spectrum.R             FFT per chromosome + Fisher g-test
  maxt.R                 maxT permutation test
  stability.R            stable peaks per condition
  peaks_genes.R          gene-level reconstruction of a peak
  compare.R              signature, transitions, cross-dataset crossing
  paths.R                where each stage reads and writes
scripts/
  00_check_inputs.R      verify the GEO files are where the configs expect them
  01_ingest_dataset.R    GEO -> common format
  02_spectra.R           FFT, per sample and per condition
  03_maxt.R              permutation test (the slow stage)
  04_stability.R         stable peaks + condition peak tables
  05_peaks_genes.R       gene reconstruction per peak
  06_compare_datasets.R  signature, transitions, replication
  99_validate_against_legacy.R   peak-by-peak diff against an old run
tests/
  test_labels.R          label engine
  test_spectrum.R        FFT, maxT, stability, Wilson
```

## Common format

Every dataset lands in `<interim_dir>/<GSE>/` with the same schema, so nothing
downstream needs to know which GEO series it came from:

| file | contents |
|---|---|
| `samples.tsv` | `sample_id, dataset_id, condition, fibrosis_stage, cohort, label_rule, keep` |
| `genes.tsv` | `gene_id, gene_name, chr, start, gene_order` |
| `counts.tsv` | raw counts, `gene_id` + one column per sample |
| `expression.tsv` | `asinh(TPM)` (or `asinh(CPM)` when gene lengths are unavailable — recorded in the manifest) |
| `label_audit.tsv` | every input sample, resolved or not, with the rule that fired |
| `manifest.tsv` | code version, config fingerprint, timestamp, per-condition counts |

`gene_order` is the FFT input axis, computed once at ingest instead of being
re-derived with `order(start)` at each downstream call site.

## Condition vocabulary

`Control, F0, F1, F2, F3, F4`

- **Control** is a *cohort* statement: liver from a subject outside the disease
  cohort. It is never inferred from a fibrosis stage of 0 or from a description
  of normal histology.
- **F0–F4** are biopsy fibrosis stages *within* the disease cohort.

This is the fix for the main inconsistency in the original scripts. GSE135251
has a real non-NAFLD control group, so `Control` and `F0` (NAFLD without
fibrosis) are different populations. GSE162694 is a single biopsy cohort: its
`N` / "normal liver histology" samples are `F0`, and the config declares
`has_control_cohort = FALSE`. The engine refuses to emit `Control` for a dataset
that declares it has none.

Cross-dataset work runs over `comparable_conditions()`, which reports what each
dataset lacks instead of silently producing an empty intersection.

## Adding a dataset

Write `config/datasets/<GSE>.R` returning a list with `id`, `counts_file`,
`series_matrix`, `has_control_cohort` and an ordered `condition_rules` list.
Three rule types are available:

- `column_match` — exact match on a phenotype column (used for cohort membership)
- `title_token` — stage token parsed from the sample title
- `fibrosis_stage` — numeric stage from a phenotype column

The first rule that resolves a sample wins; samples no rule resolves come back
with `condition = NA` and are reported, never guessed.

## Running

The configs default to the INMEGEN scratch layout:

```
geo_dir     /scratch/home/dperez/GPIB/gene_notes/TissueSpectF/data
interim_dir /scratch/home/dperez/GPIB/gene_notes/TissueSpectF/interim
results_dir /scratch/home/dperez/GPIB/gene_notes/TissueSpectF/results
```

From the repository root:

```bash
make test        # label engine + spectral core, no Bioconductor needed
make check       # are the GEO files where the configs expect them?
make all         # ingest -> spectra -> maxT -> stability -> peaks -> compare
```

Stage by stage, with the same arguments available everywhere
(`<GSE>`, `--cond=F2`, `--branch=median`):

```bash
Rscript scripts/01_ingest_dataset.R
Rscript scripts/02_spectra.R
Rscript scripts/03_maxt.R               # slow: B x samples x chromosomes
Rscript scripts/04_stability.R
Rscript scripts/05_peaks_genes.R
Rscript scripts/06_compare_datasets.R
```

## Re-running one stage

Stages communicate only through files, never through session state, so any
stage can be re-run on its own. After changing `stable_frac` or `alpha` in
`config/project.R`:

```bash
Rscript scripts/04_stability.R
Rscript scripts/05_peaks_genes.R
```

maxT is not recomputed. `03_maxt.R` is the one stage that reuses existing output
by default, because it is the only one that costs hours; pass `--force` to
recompute it. Every other stage always recomputes, since silent reuse is how two
different runs end up interleaved in one output tree.

## Validating the port

```bash
Rscript scripts/99_validate_against_legacy.R GSE135251 /path/to/old/GSE135251
```

Compares the stable-peak key set and per-peak power against a run of the
original script, condition by condition, and exits non-zero on any difference.
The only differences that should be accepted are samples whose condition changed
because of the label fix, and each one should be traceable in `label_audit.tsv`.

`00_check_inputs.R` verifies each expected file, prints the closest names it
actually found when one is missing, and checks that the output directories are
writable. GEO count tables often carry a suffix (`_raw_counts_GRCh38`,
`_norm_counts_TPM`), so if a name differs, edit `counts_file` / `series_matrix`
in `config/datasets/<GSE>.R` rather than renaming the download.

To run somewhere else, override the paths:

```bash
export TSF_GEO_DIR=/some/other/data
export TSF_INTERIM_DIR=/some/other/interim
export TSF_RESULTS_DIR=/some/other/results
```

`GEOquery` is used when installed; a minimal series-matrix parser is used
otherwise, which keeps ingest testable without Bioconductor. On a cluster with
environment modules, load R first (`module load R/4.3` or equivalent).

## Sanity check after ingest

```bash
cut -f3 $TSF_INTERIM_DIR/GSE135251/samples.tsv | sort | uniq -c
cut -f3 $TSF_INTERIM_DIR/GSE162694/samples.tsv | sort | uniq -c
```

GSE162694 must contain no `Control` rows; GSE135251 must show `Control` and `F0`
as separate, non-identical counts. If both hold, the label fix is in effect.

## What changed from the original scripts

1. `save_tsv(overwrite = FALSE)` silently kept files from earlier runs. Writes
   now overwrite by default and each output tree carries a manifest.
2. The per-transition tables named `*_cambian_magnitud.tsv` held peaks *stable in
   both* conditions, with no amplitude change computed anywhere.
   `transition_table()` now reports the amplitude delta, its log2 ratio, and a
   Welch test on per-sample peak power, and the file is named for what it holds.
3. Cross-dataset results are reported side by side, never pooled. A combined
   percentage over summed counts would be dominated by the larger cohort and
   would hide heterogeneity; `replicated` (same direction, FDR <= 0.05 in both)
   is the replication statement instead.
4. Gene order along a chromosome is computed once at ingest, not re-derived at
   each call site.
