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
RUNBOOK.md                 the sequence for a full run, and what to check
requirements-ml.txt        dependencies of the learned layer (optional)
TissueSpectF_colab.ipynb   the whole pipeline on a free Colab VM
Makefile
.github/workflows/
  tests.yml                unit tests + self-check on every push
  release.yml              builds and publishes the app bundle on a version tag
config/
  project.R                paths, filters, and every tunable parameter
  datasets/<GSE>.R         one per dataset: file names, tissue, label rules
  vocabularies/<name>.R    condition vocabularies (levels, states, roles, order)
  autoencoder.yaml         configuration of the learned layer
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
  consensus.R              characteristic spectrum of a condition (power, prevalence, PLV)
  clean.R                  CLEAN deflation with an extended-BIC stopping rule
  stability.R              which peaks go downstream, and by which criterion
  peaks_genes.R            gene-level reconstruction of a component
  compare.R                signature, transitions, cross-cohort replication
  fingerprint.R            a sample's spectrum as a comparable feature vector
  reference.R              reference library, out-of-cohort validation, matching
  bundle.R                 package the app + reference into a portable folder
  stages.R                 each stage as a callable function
ml/                        TissueSpect-AE, the learned layer (Python)
  utils.py                 seeds, device, manifest and compatibility checks
  dataset.py               spectra -> masked tensor, normalisation
  splits.py                leave-one-cohort-out, leak checks, balanced sampling
  prototypes.py            cohort-balanced class prototypes
  evaluate.py              metrics, silhouette, cohort predictability
  baselines.py             the bar the model has to clear
scripts/
  tsf.R                    CLI implementation
  00_check_inputs.R        verify the GEO inputs are where the configs expect
  inspect_series_matrix.R  read a GEO series before writing rules for it
  prepare_ae_data.R        export spectra for the learned layer
  run_baselines.py         leave-one-cohort-out baselines
  match_query.R            the match command
  selfcheck.R              end-to-end run on synthetic data with a known peak
app/
  app.R                    local desktop app (the only part that needs shiny)
tests/
  test_labels.R            label engine, configs, filters, bundle
  test_spectrum.R          grid, GLS, maxT, CLEAN, condition test, Wilson
  ml/                      the learned layer (pytest, numpy/sklearn only)
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

## Classes

`class_id = tissue::state::condition`

Three levels because two are not enough. `state` separates healthy tissue from
diseased tissue, `condition` names the class within it. With two levels, healthy
liver from a non-disease cohort and healthy liver from a biopsy series would
have to share one label or invent unrelated ones; with three, both sit under
`liver::healthy` and stay distinguishable, and TCGA's adjacent-normal will fit
later without renaming anything.

### Three healthy classes, not one

| class | what it is | cohorts |
|---|---|---|
| `Normal_histology` | a biopsy with normal histology and no NAFLD activity (NAS = 0), taken **within** a biopsy series | GSE162694 (31), GSE130970 (6) |
| `Control_disease_cohort` | a subject **outside** the disease cohort of a NAFLD study | GSE135251 (10) |
| `Control_external_study` | the control group of a **different** disease study | GSE142530 (11) |

None absorbs the others. Merging them into one `Control` would assume exactly
what is worth testing — that a healthy liver looks the same whichever study
recruited it — and would do it by choosing a label. If their spectra agree, a
`Healthy_consensus` class can be created and defended afterwards. Equivalence is
demonstrated, not asserted.

`Control_disease_cohort` is **not a strictly healthy extreme**: two of its ten
carry incidental fibrosis of stage 1 and 2, which is what makes it a cohort
statement rather than a histological one. Audit those two out before using the
class as a healthy reference.

### The disease classes

`NAFLD_fibrosis_F0` through `F4`, biopsy fibrosis stages within the disease
cohort. `Normal_histology` is not `F0`: in GSE162694 the phenotype table lists
`normal liver histology` and `0` as separate values, and NAS settles it — all 31
normal-histology samples score 0 while the 35 at stage 0 score 1 to 5.

### Ordering and roles

A vocabulary declares three things beyond its levels:

- **`progression`** — the ordered subset. Transitions and ordinal trends run
  over it alone: `F0 → F1` is a step, `Control_disease_cohort →
  Normal_histology` is not.
- **`states`** — healthy or disease, which becomes the middle term of `class_id`.
- **`cohort_roles`** — `control`, `disease` or `within_disease_normal`. Deriving
  the role by comparing against the baseline alone marked `Control_external_study`
  as disease, because it is not the baseline: a conclusion drawn from a naming
  convention rather than from the biology.

Cross-dataset work runs over `comparable_conditions()`, which reports what each
dataset lacks instead of silently producing an empty intersection.

## The cohorts

| cohort | Control_disease_cohort | Control_external_study | Normal_histology | F0 | F1 | F2 | F3 | F4 |
|---|---|---|---|---|---|---|---|---|
| GSE135251 | 10 | — | — | 38 | 47 | 53 | 54 | 14 |
| GSE162694 | — | — | 31 | 35 | 30 | 27 | 8 | 12 |
| GSE130970 | — | — | 6 | 19 | 28 | 9 | 14 | 2 |
| GSE276114 | — | — | — | — | — | — | 13 | 42 |
| GSE142530 | — | 11 | — | — | — | — | — | — |

GSE276114 spans three etiologies and contributes MASLD only, and only F3/F4: its
`F0-2` bin spans three classes and resolves to none of them. GSE142530 is an
alcohol study and contributes only its controls. Neither restriction is a
filter applied for convenience — pooling etiologies would put one class label on
a mixture of diseases.

## Adding a dataset

Read the series before writing rules for it:

```bash
Rscript scripts/inspect_series_matrix.R data/<GSE>_series_matrix.txt.gz "field" "other"
```

It prints every phenotype field with its counts and cross-tabulates two of them.
That last view is what settles whether a class is what its name suggests — it is
how GSE162694's 31 "normal liver histology" samples turned out to have NAS = 0
while its 35 stage-0 samples had NAS 1–5, two groups the pipeline had been
merging. **Write the config from that output, not from the paper's description
of it.**

Then `config/datasets/<GSE>.R` returns a list with `id`, `tissue`, `vocabulary`,
`counts_file`, `series_matrix`, `has_control_cohort` and an ordered
`condition_rules` list. Four rule types:

- `column_match` — exact match on a phenotype column
- `title_token` — stage token parsed from the sample title
- `fibrosis_stage` — numeric stage from a phenotype column
- `compound` — every sub-condition must hold, for a class defined by a
  combination of scores rather than one field

The first rule that resolves a sample wins; samples no rule resolves come back
with `condition = NA` and are reported, never guessed.

### Restricting what a cohort contributes

- `sample_filter` — keep only samples matching a field, or exclude those that
  do. A sample whose filter field is missing is **dropped**, never admitted: in
  a mixed cohort an absent label is not evidence that it is the one wanted. A
  missing filter column **aborts** rather than admitting everything.
- `keep_conditions` — the classes this cohort contributes at all.
- `expected_n_samples` — ingest stops if the count changes.

### Reading count matrices that are not plain TSVs

GEO publishes count matrices in whatever shape the submitters chose, and those
are properties of a file rather than of the analysis, so they go in
`counts_spec`: `sep`, `id_column`, `symbol_column`, `skip`, `sample_map_row`,
`map_transform`, `exclude_columns`. `count_id_type` is `ENTREZID`, `ENSEMBL` or
`SYMBOL`.

`sample_map_row` names a second header row carrying the real sample names. It is
extracted as a map, used to rename the columns, and removed **before** anything
is coerced to numeric — left in place it turns every count column into character
and the matrix silently becomes text. `map_transform` reconciles labels that
differ from the phenotype table (`Control_Lille 389` against `Control_389`), and
every substitution is recorded next to the raw value in `count_column_map.tsv`
with whether it matched. A regular expression that rewrites sample names and
leaves no table behind is not a mapping, it is a guess.

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
`condition`, `consensus`, `clean`, `stability`, `peaks`, `compare`.

## The characteristic spectrum of a condition

```bash
./tsf consensus
```

The spectrum of the mean profile is not a summary of the samples' spectra. The
transform is linear, so averaging the profiles averages the complex
coefficients: a component present in every sample at the same frequency but with
scattered phases cancels and disappears, however reproducible it is.

`consensus` therefore works from the per-sample spectra and reports each
frequency on three axes — median normalised power (**how strong**), prevalence
(**how common**, judged within each sample and chromosome so that chromosomes do
not compete), and the phase-locking value (**how aligned**) — plus heterogeneity
of power and phase, bootstrap intervals over samples, and the number of valid
samples. The consensus score is their product, so a component has to satisfy all
three.

It is a ranking statistic, not a test: it has no null and no error rate. Every
PLV comes with a Rayleigh p-value, because under random phases the expected PLV
is about `sqrt(pi)/(2*sqrt(n))`, not zero — with 8 samples a PLV of 0.3 means
nothing. That p-value cannot fall below `exp(-n)`, so a condition with fewer
than roughly `log(n_frequencies/q)` samples cannot establish alignment at all;
the stage says so and ranks by prevalence and score instead of returning an
empty signature.

`signature_<condition>.tsv` is the exported characteristic signature. Each
component carries `signature_class`: **confirmed** when prevalence holds, the
adjusted Rayleigh p holds, and the component beats its own label-permuted null
(BH-adjusted across components); **exploratory** otherwise. The permuted
comparison uses the rank-based score on both sides — the maxT-based score exists
only where the per-sample test was run, so it is reported as confirmatory
evidence rather than used as the test statistic — including every component from a condition with too few samples to
test phase alignment at all. An exploratory component is a lead, not a
signature, and the stage counts only confirmed ones in its summary.

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

No cohort should show a class it does not contribute: GSE162694 has no
`Control_*`, GSE135251 has no `Normal_histology`, and `Control_disease_cohort`
and `F0` are separate with different counts. Compare against the cohort table
above; a count that differs means a label went somewhere unexpected, and nothing
downstream is worth reading until it is reconciled.

## TissueSpect-AE, the learned layer

A complementary layer, not a replacement. The statistical results — permutation
p-values, consensus spectra, replication across cohorts — remain the evidence.

> **A peak reconstructed by TissueSpect-AE is not statistical or causal evidence
> on its own. Components only the model finds are reported as AI candidates.**

What exists today is the data layer and the evaluation harness. Neither needs
`torch`; the R pipeline and its tests do not need Python at all.

```bash
pip install -r requirements-ml.txt          # optional
./tsf ae-prepare --interim-dir interim --results-dir results_pc
python3 scripts/run_baselines.py --data results_pc/autoencoder/data \
                                 --out  results_pc/autoencoder/baselines
python3 -m pytest tests/ml/ -q
```

`ae-prepare` exports the per-sample spectra as plain TSVs plus a manifest, and
refuses when cohorts disagree on annotation, gene universe, grid digest or
expression unit — a model trained across two grids would be learning the grids.
Every frequency is exported, not a peak list: a fingerprint does not need any
component to be individually significant, and filtering first would remove
exactly what makes classes separable.

Missing frequencies are carried as a mask with four levels, distinguishing what
is above a chromosome's Nyquist limit (structurally absent, identical for every
sample) from what this sample could not estimate (varies by sample, and so
correlates with cohort). Arrays are rectangular, so unobserved positions hold
zero — but only after the mask exists, and normalisation re-zeroes them, because
subtracting a centre would otherwise turn every hole into a constant the network
could read as a missingness flag.

### The baselines come first, deliberately

`run_baselines.py` runs a cohort-balanced nearest centroid, elastic net and
random forest through leave-one-cohort-out, and prints a per-class report saying
which classes were actually evaluated and which exist in one cohort only and so
contribute nothing to any average. A single macro-F1 presents every class as if
it had been tested equally; with these cohorts it never is.

The column that matters is `lift` over the training fold's majority class, not
accuracy. This is the bar the autoencoder has to clear, fixed before the model
exists so it cannot be set to fit the answer. If nothing lifts over the baseline,
the spectra do not distinguish stages across cohorts and no architecture fixes
that — a network with a domain adversary trained on a few hundred samples would
find something, and that something would be the cohorts.

## The app

Locally, from the repository:

```bash
Rscript -e 'install.packages("shiny")'   # once
./tsf run --to spectra
./tsf reference
./tsf app
```

### Giving it to someone else

The usual route is a Release: tag a version and GitHub Actions builds the
bundle, so the recipient downloads a zip and nobody has to run anything.

```bash
git tag v0.2.0 && git push --tags
```

`.github/workflows/release.yml` fetches the GEO inputs, ingests, builds the
reference, runs the out-of-cohort validation, packages the bundle, and attaches
it to the Release with its SHA-256 and the provenance file as the release notes.
It **refuses to publish** a reference that does not beat its majority-class
baseline: shipping one would mean distributing confident-looking output with no
information in it. The validation numbers, including the per-coverage-band
table, land on the run's summary page.

The reference needs only `ingest` — fingerprints come from the expression
matrices — so the slow spectral stages are not on the release path. What does
cost time there is the coverage calibration, which recomputes a fingerprint per
(sample, coverage level, loss mode, mask); `workflow_dispatch` exposes
`n_masks` for a cheaper trial build.

Binary references do not belong in git, which is why the bundle is a release
asset rather than a committed file.

Locally, the same thing:

```bash
./tsf bundle --out TissueSpectF-app
```

Produces a folder (and a zip) that depends on nothing else — not this
repository, not `config/project.R`, not the interim or results directories, not
the GEO downloads. It contains the reference, the four R modules a query needs,
a launcher for macOS/Linux and one for Windows, a README and a provenance file
recording exactly what the reference was built from and how well it validated.

The recipient unzips it and runs `run.sh` (or double-clicks `run.bat`). R is the
only prerequisite; `shiny` installs itself on first run. Nothing is uploaded.

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
- thresholds are calibrated **per coverage band** (90–100%, 75–90%, 50–75%,
  below 50%), because the similarity distribution shifts as coverage falls; a
  band with no calibrated threshold gives `UNCALIBRATED_COVERAGE`, and below 50%
  nothing is classified

`./tsf reference` writes `out_of_cohort_predictions.tsv` and
`confusion_matrix.tsv` so the validation can be inspected rather than trusted.

The reference is **self-contained**: it carries the canonical annotation grid
(not the genes one cohort happened to observe), gene identifiers, species,
genome build, annotation release, frequency ceiling and feature representation.
A `reference.rds` is portable to a machine that has never seen the cohorts it
was built from. Datasets are refused entry into one reference unless species,
build, release, gene universe and grid digest all match — a feature named
`chr7_k12` means a different thing under a different build.

A query is fingerprinted on **the genes it actually contains**. Genes the file
lacks are absent from the observed positions, never zero-filled — the same rule
that governs unmeasured genes everywhere else here. The fingerprint is
normalised *after* intersecting with the reference feature space, so query and
centroids are standardised over the same features, and similarity is computed
over the shared features only: an unobserved frequency contributes nothing
rather than contributing the mean.

The reference records the expression unit it was built on, and a TPM query
against a CPM reference is refused: length normalisation changes the relative
height of every gene, so it changes the spectrum.

Declare what the values are with `--input-unit` (`counts`, `cpm`, `tpm`,
`logged`). Duplicate identifiers are summed for counts and refused for anything
already normalised; negative values are refused. A query covering less than 20%
of the grid, or less than 50% of the features the model uses, is reported
`LOW_COVERAGE` and **not scored** — in the CLI and in the app alike.

Coverage is calibrated rather than assumed, and measured against the canonical
grid rather than against whatever the validation cohort happened to observe: a
cohort covering 70% of the grid, masked to 80% of its own genes, is a 56% query,
not an 80% one. Validation masks **genes** five ways — scattered, one retained
block, several missing blocks, whole chromosomes, and dropout of the least
expressed — recomputes the fingerprint from the survivors, and reports accuracy and threshold per band
with the spread between masks. Masking spectral features instead would measure
an easier quantity: half a chromosome's genes can go missing and the GLS still
estimates nearly every frequency, so feature coverage stays near 100% while gene
coverage is 50%.

Two thresholds are computed per band, pooled and conservative, and both
rejection rates are reported; `fingerprint$threshold_policy` declares which is
applied. The default is pooled: with gene-level masking the conservative
threshold rejected 46–79% of true members.

The rejection threshold is calibrated on the similarity that correct held-out
matches reach, per class and per coverage band. It bounds how often a true member is wrongly rejected.
It does **not** bound how often an out-of-domain sample is wrongly accepted:
no out-of-domain sample was in the validation. Open-set specificity needs
negatives — for a tissue reference, other tissues.
