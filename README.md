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
  stages.R               each stage as a callable function
tsf                      the command line entry point
scripts/
  tsf.R                  CLI implementation
  00_check_inputs.R      verify the GEO files are where the configs expect them
  selfcheck.R            end-to-end run on synthetic data with a known peak
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

```bash
./tsf check       # are the GEO files where the configs expect them?
./tsf run         # ingest -> spectra -> maxt -> stability -> peaks -> compare
./tsf status      # what exists on disk
```

`run` prints a per-stage timing table and stops at the first failure, so a long
run never leaves you guessing which stage produced which output.

Any stage can be run on its own, over a subset:

```bash
./tsf spectra                        # one stage, all datasets
./tsf maxt GSE135251 --cond=F3       # one dataset, one condition
./tsf run --from=stability           # resume after changing a threshold
./tsf run --to=maxt --log=run.log    # stop early, tee the output to a file
./tsf run --dry-run                  # print the plan, do nothing
```

Options: `--from=`, `--to=`, `--cond=F2,F3`, `--branch=median`, `--force`,
`--dry-run`, `--log=<file>`. Stages in order: `ingest`, `spectra`, `maxt`,
`stability`, `peaks`, `compare`.

## Checking that it works

```bash
make test         # 38 unit tests: labels, FFT, maxT, stability, Wilson
./tsf selfcheck   # the full pipeline on synthetic data with a known answer
```

`selfcheck` builds a GEO-shaped dataset in a temporary directory with one
sinusoid injected on chromosome 1, whose amplitude grows with fibrosis stage,
runs every stage, and asserts that the labels resolve as declared, that the
recovered peak is the injected `(chr, N, k)`, that the phase matches, that
amplitude increases across every transition, and that the increase replicates
across both datasets. It exercises the same code a real run does, so a
regression anywhere in the chain shows up. It says nothing about whether a
signal found in real data is biologically meaningful -- only that the machinery
recovers a signal known to be there.

## Re-running one stage

Stages communicate only through files, never through session state, so any
stage can be re-run on its own. After changing `stable_frac` or `alpha` in
`config/project.R`:

```bash
./tsf run --from=stability
```

maxT is not recomputed. It is the one stage that reuses existing output by
default, because it is the only one that costs hours; pass `--force` to
recompute it. Every other stage always recomputes, since silent reuse is how two
different runs end up interleaved in one output tree.

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
