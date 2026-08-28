# TissueSpectF

Chromosome-ordered Fourier spectral analysis of tissue level transcriptomes.

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
scripts/
  01_ingest_dataset.R    CLI entry point
tests/
  test_labels.R          dependency-free tests (Rscript tests/test_labels.R)
data/interim/<GSE>/      the common format (git-ignored)
results/                 downstream outputs (git-ignored)
```

## Common format

Every dataset lands in `data/interim/<GSE>/` with the same schema, so nothing
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

```bash
Rscript tests/test_labels.R                       # no Bioconductor needed
Rscript scripts/01_ingest_dataset.R               # all datasets
Rscript scripts/01_ingest_dataset.R GSE135251     # one dataset
```

Paths can be overridden with `LFFT_GEO_DIR`, `LFFT_INTERIM_DIR`,
`LFFT_RESULTS_DIR`.

`GEOquery` is used when installed; a minimal series-matrix parser is used
otherwise, which keeps ingest testable without Bioconductor.

## Still to port

The spectral stages (`02_spectra`, `03_maxt`, `04_stability`, `05_peaks_genes`,
`06_compare_datasets`) are unchanged in intent from the original scripts and
should be lifted into `R/spectrum.R`, `R/maxt.R`, `R/stability.R`,
`R/peaks_genes.R`, `R/compare.R` one at a time, each reading only the common
format. Two things to change while porting:

1. `save_tsv(overwrite = FALSE)` silently kept files from earlier runs. Here
   writes overwrite by default and each output directory carries a manifest.
2. The per-transition tables named `*_cambian_magnitud.tsv` contained peaks
   *stable in both* conditions, with no amplitude change computed. Either rename
   them or add the amplitude delta and its interval.
