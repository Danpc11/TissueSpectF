# TissueSpectF

![R](https://img.shields.io/badge/R-4.1%2B-276DC3?logo=r&logoColor=white)
[![Tests](https://github.com/Danpc11/TissueSpectF/actions/workflows/tests.yml/badge.svg)](https://github.com/Danpc11/TissueSpectF/actions/workflows/tests.yml)
[![Full pipeline](https://img.shields.io/badge/Colab-Full%20pipeline-4285F4?logo=googlecolab&logoColor=white)](https://colab.research.google.com/github/Danpc11/TissueSpectF/blob/main/TissueSpectF_colab.ipynb)
![Datasets](https://img.shields.io/badge/GEO-GSE135251%20%7C%20GSE162694-1f6feb)
![Dependencies](https://img.shields.io/badge/dependencies-base%20R-success)

Chromosome-ordered Fourier spectral analysis of tissue level transcriptomes.

Gene expression is read as a signal along the chromosome, ordered by genomic
position, and decomposed into periodic components. The question is whether
transcription carries positional structure at scales of tens to hundreds of
genes, and whether that structure changes with disease stage.

The model, the estimators and the limits of each claim are written up in
[THEORY.md](THEORY.md).

Three things make the result trustworthy rather than merely computable:

- the spectral axis is the **annotation grid**, not the genes that survived a
  filter, so `N` means the same thing in every dataset
- unmeasured genes stay **unobserved**, never zero-filled, and the spectrum is
  fitted by least squares over the observed positions only
- every peak is reported next to its **spectral window**: what the pattern of
  missing genes alone can produce

This repository replaces two ~4,200-line per-dataset scripts with a layered
pipeline. Everything that was dataset-specific now lives in a config file;
everything else is shared code.

## Layout

```
tsf                        the command line entry point
README.md
THEORY.md                  the model, the estimators, what each test licenses
TissueSpectF_colab.ipynb   the whole pipeline on a free Colab VM
Makefile
.github/workflows/
  tests.yml                unit tests + self-check on every push
config/
  project.R                paths, filters, and every tunable parameter
  datasets/<GSE>.R         one per dataset: file names, tissue, label rules
  vocabularies/<name>.R    condition vocabularies (levels, order, baseline)
R/
  utils_io.R               logging, TSV I/O, phenotype cleaning, manifests
  config.R                 config and vocabulary loading, validation
  labels.R                 label harmonisation engine
  fetch.R                  download the GEO inputs
  ingest.R                 GEO -> common format
  paths.R                  where each stage reads and writes
  grid.R                   reference gene grid + least-squares (Lomb-Scargle) spectrum
  spectrum.R               spectra per sample and per condition
  maxt.R                   per-sample permutation test
  condition_test.R         condition-level significance (permutation + Stouffer)
  clean.R                  CLEAN deflation with an extended-BIC stopping rule
  stability.R              which peaks go downstream, and by which criterion
  peaks_genes.R            gene-level reconstruction of a component
  compare.R                signature, transitions, cross-cohort replication
  fingerprint.R            a sample's spectrum as a comparable feature vector
  reference.R              reference library, out-of-cohort validation, matching
  stages.R                 each stage as a callable function
scripts/
  tsf.R                    CLI implementation
  00_check_inputs.R        verify the GEO inputs are where the configs expect
  match_query.R            the match command
  selfcheck.R              end-to-end run on synthetic data with a known peak
app/
  app.R                    local desktop app (the only part that needs shiny)
tests/
  test_labels.R            label engine
  test_spectrum.R          grid, GLS, maxT, CLEAN, condition test, Wilson
```

Outputs go where `config/project.R` says (`interim_dir`, `results_dir`),
outside the repository by default.

## Common format

Every dataset lands in `<interim_dir>/<GSE>/` with the same schema, so nothing
downstream needs to know which GEO series it came from:

| file | contents |
|---|---|
| `samples.tsv` | `sample_id, dataset_id, condition, fibrosis_stage, cohort, label_rule, keep` |
| `genes.tsv` | `gene_id, gene_name, chr, start, grid_index, grid_N` |
| `grid_coverage.tsv` | observed genes vs grid length, per chromosome |
| `counts.tsv` | raw counts, `gene_id` + one column per sample |
| `expression.tsv` | `asinh(TPM)` (or `asinh(CPM)` when gene lengths are unavailable — recorded in the manifest) |
| `label_audit.tsv` | every input sample, resolved or not, with the rule that fired |
| `manifest.tsv` | code version, config fingerprint, timestamp, per-condition counts |

## The spectral axis

The axis is the **reference grid**: every annotated gene of the allowed biotypes
on a chromosome, ordered by start position, indexed `1..N`. `N` is a property of
the annotation, so it is identical across datasets and conditions and the peak
key `(chr, N, k)` means the same thing everywhere. Set the universe with
`gene_universe` in `config/project.R` (default `PROTEIN_CODING|NCRNA`).

Genes without an expression measurement keep their grid slot and are simply
unobserved. They are never zero-filled or imputed: absence of a measurement is
not evidence of zero expression, and a zero at a fixed position injects a
deterministic pattern whose own spectral structure can manufacture a peak.

Spectra are therefore computed by least squares (generalised Lomb-Scargle with a
floating mean) over the observed positions only. This is still Fourier analysis:
on a complete grid the fit reproduces the DFT coefficients exactly, which
`tests/test_spectrum.R` checks against `fft()` and against a brute-force
least-squares fit on a gapped grid. Internally the sums are evaluated with two
FFTs, so the cost stays O(N log N).

What gaps cost: the sinusoids are no longer orthogonal over the observed
positions, so frequency bins are not independent, Parseval no longer holds
exactly, and the Nyquist limit is not uniquely defined. The first is why the
permutation null is mandatory rather than optional. The third is why every peak
carries `window_power` and `window_rank`, and why `./tsf window` exists.

## Spectral window

```bash
./tsf window
```

Computes, per chromosome, the spectrum of the presence indicator alone -- what
the pattern of missing genes can produce with no expression signal at all -- and
places every stable peak in it. A peak in the top 1% of the window sits exactly
where the gaps are most periodic and should be treated as a sampling artefact
until shown otherwise. Run this before interpreting any peak near the Nyquist
limit.

The permutation null permutes values **among the observed positions**, holding
the positions fixed, so the missingness pattern is identical in the data and in
every permutation and cannot by itself produce significance. That is checked
directly in the test suite.

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
./tsf run         # the whole pipeline, in order
./tsf status      # what exists on disk
```

`run` prints a per-stage timing table and stops at the first failure, so a long
run never leaves you guessing which stage produced which output.

Any stage can be run on its own, over a subset:

```bash
./tsf spectra                             # one stage, all datasets
./tsf maxt GSE135251 --cond F3            # one dataset, one condition
./tsf run --from stability                # resume after changing a threshold
./tsf run --to maxt --log run.log         # stop early, tee output to a file
./tsf run --dry-run                       # print the plan, do nothing
```

Paths and parameters are flags, so nothing has to be exported:

```bash
./tsf run --to spectra --results-dir results_pc
./tsf run --gene-universe '^(protein-coding|ncRNA)$' \
          --results-dir results_pc_nc --interim-dir interim_pc_nc
./tsf stability --stable-frac 0.7 --criterion consistency
./tsf reference --target tissue --k-max 96
```

Both `--key value` and `--key=value` work, names are case-insensitive, and `-`
and `_` are interchangeable — `--results-dir`, `--results_dir` and
`--RESULTS_DIR` are the same flag. Precedence is **command line > environment
(`TSF_*`) > `config/project.R`**, and every override is echoed in the log.

`./tsf --help` lists the flags. Stages in order: `ingest`, `spectra`, `maxt`,
`condition`, `clean`, `stability`, `peaks`, `compare`.

## What decides that a peak exists

Two criteria are computed for every peak; `stability_criterion` in
`config/project.R` picks which one drives `is_stable`, and both are written
either way.

**condition** (default) -- a permutation test on the condition's summary signal.
Values are permuted among the observed grid positions, the spectrum recomputed,
and the maximum power over frequencies gives the null. This is the primary test:
it evaluates the signal that is actually analysed downstream, and its power
grows with the number of samples.

**consistency** -- at least `stable_frac` of the samples individually
significant. This is a reproducibility requirement, not a test with aggregated
power: each sample has to reach significance alone under family-wise control, so
a real but moderate periodicity present in every sample stays invisible however
many samples there are. Kept as a descriptor reported next to the result, not as
the thing that decides.

Where `maxt` output exists, a Stouffer combination of the per-sample p-values is
added as a second opinion. It is conservative (its inputs are already family-wise
adjusted) and assumes independence between samples, which replicates only
approximately satisfy.

Multiplicity is corrected over **chromosomes, not frequencies**: maxT already
controls the family-wise error rate across frequencies within a chromosome.
Correcting over frequencies again would be fatal rather than merely conservative,
because a permutation p-value cannot fall below `1/(B+1)` and the smallest
attainable q would exceed 1. This puts a floor on B: `condition_B >= 20 * n_chr`,
about 460 for 23 chromosomes. The stage warns when B is too small to conclude.

## Running on Colab

`TissueSpectF_colab.ipynb` runs the whole thing on a free VM, starting
from `./tsf fetch` so no data needs copying. Two cores make the per-sample
`maxt` stage impractical (roughly 8 hours), so the notebook runs
`--from=condition`, which is about 1/n the cost and answers the primary
question. The stability stage falls back to the condition test on its own when
no per-sample maxT exists.

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

## The app

```bash
Rscript -e 'install.packages("shiny")'   # once
./tsf run --to spectra
./tsf reference
./tsf app
```

Opens in a browser and runs entirely on the machine that started it: no upload
leaves the computer and no server is contacted. Drop in a counts TSV (gene id
column, then one column per sample, Ensembl or Entrez ids) and it returns the
ranked classes.

`shiny` is needed only for `app.R`. The CLI, including `./tsf match`, works on
base R alone.

### What the app will not do

A matcher always returns a best match. Whether it means anything is a separate
question, and the interface answers it before showing the result:

- the banner reports **out-of-cohort accuracy** against the always-guess-the-
  commonest-class baseline, and says plainly when the reference does not beat
  guessing
- a reference built from one cohort is labelled **uncalibrated** — internal
  cross-validation is not evidence here, because a classifier with thousands of
  features reaches a high internal accuracy by learning batch and sequencing
  depth, neither of which transfers between studies
- a top-two margin under 0.02 is reported as **not separable**
- if a randomly shuffled copy of the query scores as well, it says there is no
  usable spectral shape in the file

- a query below the **calibrated rejection threshold** is reported `UNKNOWN`
  rather than assigned to the least distant centroid

`./tsf reference` writes `out_of_cohort_predictions.tsv` and
`confusion_matrix.tsv` so the validation can be inspected rather than trusted.

The reference is **self-contained**: it carries its own grid, gene identifiers,
frequency ceiling and feature representation, so a `reference.rds` is portable
to a machine that has never seen the cohorts it was built from, and a query can
never be scored on a grid other than the one the reference was built with.

A query is fingerprinted on **the genes it actually contains**. Genes the file
lacks are absent from the observed positions, never zero-filled — the same rule
that governs unmeasured genes everywhere else here.

The rejection threshold is calibrated on the similarity that correct held-out
matches reach, per class. It bounds how often a true member is wrongly rejected.
It does **not** bound how often an out-of-domain sample is wrongly accepted:
no out-of-domain sample was in the validation. Open-set specificity needs
negatives — for a tissue reference, other tissues.

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
