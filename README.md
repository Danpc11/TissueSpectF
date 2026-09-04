# TissueSpectF

![R](https://img.shields.io/badge/R-4.1%2B-276DC3?logo=r&logoColor=white)
![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)
[![Tests](https://github.com/Danpc11/TissueSpectF/actions/workflows/tests.yml/badge.svg)](https://github.com/Danpc11/TissueSpectF/actions/workflows/tests.yml)
[![Full pipeline](https://img.shields.io/badge/Colab-Full%20pipeline-4285F4?logo=googlecolab&logoColor=white)](https://colab.research.google.com/github/Danpc11/TissueSpectF/blob/main/TissueSpectF_colab.ipynb)
![Cohorts](https://img.shields.io/badge/GEO-5%20liver%20cohorts-1f6feb)
![Classes](https://img.shields.io/badge/classes-3%20healthy%20%7C%20F0--F4-1f6feb)
![Dependencies](https://img.shields.io/badge/pipeline-base%20R%20only-success)
![ML](https://img.shields.io/badge/TissueSpect--AE-optional-lightgrey)

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

The pipeline is layered: everything dataset-specific lives in a config file and
everything else is shared code, so adding a cohort means writing one config
rather than copying an analysis. Stages communicate only through files on disk,
never through session state, so any stage can be re-run on its own.

A fresh clone runs without configuration. Every default path in
`config/project.R` is relative to the repository, and all of them are
gitignored:

```bash
git clone https://github.com/Danpc11/TissueSpectF && cd TissueSpectF
make test         # base R only, no Bioconductor, no install step
./tsf selfcheck   # synthetic data, known answer, needs no paths of its own

# Every command that touches your data states where it reads and writes.
# config/project.R names no paths, so there is no default to inherit.
export TSF_GEO_DIR=data
export TSF_INTERIM_DIR=interim
export TSF_RESULTS_DIR=run_2026_09_02

./tsf fetch       # the five GEO cohorts into $TSF_GEO_DIR
./tsf check       # confirm every input is where the configs expect
./tsf run
```

The three paths can equally be given per command as `--geo-dir`,
`--interim-dir` and `--results-dir`, which take precedence. There is
deliberately no default: a run's output location should be visible in the
command that produced it, not decided by a file nobody read.

Point the outputs elsewhere with `TSF_ROOT`, with `TSF_GEO_DIR` /
`TSF_INTERIM_DIR` / `TSF_RESULTS_DIR`, or with `--geo-dir` / `--interim-dir` /
`--results-dir`. Nothing has to be edited to run on a cluster:

```bash
TSF_ROOT=/scratch/$USER/TissueSpectF ./tsf run
```

## Layout

```
tsf                        the command line entry point
README.md
THEORY.md                  the model, the estimators, what each test licenses
LICENSE                    CC BY-NC 4.0
requirements-ml.txt        dependencies of the learned layer (optional)
requirements-sonify.txt    dependencies of the sonification (optional)
TissueSpectF_colab.ipynb   the whole pipeline on a free Colab VM
Makefile
.github/workflows/
  tests.yml                every test suite + self-check on every push
  release.yml              builds and publishes the app bundle on a version tag
config/
  project.R                paths, filters, and every tunable parameter
  datasets/<GSE>.R         one per dataset: file names, tissue, label rules
  vocabularies/<name>.R    condition vocabularies (levels, states, roles, order)
  autoencoder.yaml         the learned layer's parameters -- EARLY, see below
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
  clean_results.R          empty results_dir, with a guard on what it may delete
  sonify_tissuespectf.py   condition spectra -> MIDI
app/
  app.R                    local desktop app (the only part that needs shiny)
tests/
  test_labels.R            label engine, configs, filters, bundle
  test_spectrum.R          grid, GLS, maxT, CLEAN, condition test, Wilson
  ml/                      the learned layer (pytest, numpy/sklearn only)
```

## Where things go

`config/project.R` names three directories and defaults all three to the
repository, so a clone runs with no configuration:

| | default | holds |
|---|---|---|
| `geo_dir` | `<repo>/data` | raw GEO downloads and the annotation |
| `interim_dir` | `<repo>/interim` | the common format, one folder per cohort |
| `results_dir` | `<repo>/results` | every spectral output |

All three are gitignored: they hold downloads and derived output, never inputs
that need versioning. No default names a particular machine — a default pointing
at one person's scratch directory makes a fresh clone fail with paths the user
has never seen, and makes a config fingerprint record a location instead of a
choice. `tests/test_labels.R` asserts this rather than leaving it to review.

Four ways to point them elsewhere, highest precedence first:

```bash
./tsf run --results-dir results_other            # one flag, one invocation
./tsf run --config config/other_tree.R           # a named set of settings
TSF_RESULTS_DIR=/scratch/$USER/results ./tsf run # per-shell
TSF_ROOT=/scratch/$USER/TissueSpectF ./tsf run   # all three at once
```

Paths on the command line may be relative: `./tsf` changes into the repository
first, so `--results-dir results_other` is enough and no absolute path is
needed.

`--config` is the one to reach for when a run tree differs in more than its
output directory — a gene universe, a period floor, a permutation count.
Copy `config/project.R`, or inherit from it and override only what changes:

```r
# config/other_tree.R
base <- local({
  .tsf_root <- Sys.getenv("TSF_ROOT", unset = getwd())
  source(file.path(.tsf_root, "config", "project.R"), local = TRUE)$value
})

modifyList(base, list(
  results_dir = "results_other"
))
```

Commit it. The file is the record of how a result was produced; a flag typed
into a terminal six weeks ago is not. Every run logs the config it read and the
results tree it resolved, so a stage never writes somewhere unnoticed.

`TSF_ROOT` is the one to reach for on a cluster where code and output belong on
different filesystems. Every override is echoed in the log, so a run never
leaves you guessing which tree it wrote to.

To empty the results tree:

```bash
make clean-dry    # list what would go
make clean        # asks you to type the directory name first
make clean-force  # no prompt; for scripts that mean it
```

`results_dir` points at the **active** run tree, so `make clean` deletes real
work — a consensus stage is tens of minutes and is not reproducible from what
survives. It reports the entry count and size, then requires the directory name
typed back. An unanswered or closed stdin counts as a refusal, so a CI job
cannot delete a tree by accident.

The guard is structural, not a length heuristic: `scripts/clean_results.R`
refuses the filesystem root, the home directory and the repository itself, and
verifies afterwards that the removal actually happened rather than assuming it.

### Every environment variable

| variable | effect | default |
|---|---|---|
| `TSF_ROOT` | base for all three directories below | the working directory |
| `TSF_GEO_DIR` | raw GEO downloads | `$TSF_ROOT/data` |
| `TSF_INTERIM_DIR` | the common format | `$TSF_ROOT/interim` |
| `TSF_RESULTS_DIR` | spectral outputs | `$TSF_ROOT/results` |
| `TSF_LIBRARY_DIR` | condition library, for the peak-gene and sonification scripts | `$TSF_RESULTS_DIR/library_domains` |
| `TSF_MAXT_B` | permutations for the per-sample maxT test | as `config/project.R` says |
| `TSF_CONDITION_B` | permutations for the condition-level test | as `config/project.R` says |
| `TSF_APP_MAX_UPLOAD_MB` | upload cap in the app | 512 |

`TSF_MAXT_B` and `TSF_CONDITION_B` exist so a smoke run finishes: they lower the
permutation count, which raises the floor on the smallest reportable p-value to
`1/(B+1)`. Use them for a self-check or a Colab demo, never for a result — the
Colab notebook and CI both set them for exactly that reason and say so.

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

### Three healthy groups, one class

| raw label | what it is | cohorts |
|---|---|---|
| `Normal_histology` | a biopsy with normal histology and no NAFLD activity (NAS = 0), taken **within** a biopsy series | GSE162694 (31), GSE130970 (6) |
| `Control_disease_cohort` | a subject **outside** the disease cohort of a NAFLD study | GSE135251 (8, see below) |
| `Control_external_study` | the control group of a **different** disease study | GSE142530 (11) |

All three map to one condition, **`Controles`**, 56 samples over four cohorts.
The raw labels stay distinct, so ingest, the audit trail and `cohort_roles`
still record which cohort a control came from; only the class is shared.

An earlier version kept them apart, on the argument that a healthy liver might
not look the same whichever study recruited it. Two things decided against it,
neither of which depends on any accuracy figure. The source publications place
all three in the same group, and `states` already marked all three `healthy` —
the split lived in the label. And the split was not validable: two of the three
existed in a single cohort each, so leave-one-cohort-out had no second cohort to
learn them from. On a real run GSE142530 was skipped as a hold-out entirely,
ten samples were dropped from validation, and both classes sat among the
centroids competing for every prediction while winning zero.

**What it costs, stated:** the class is heterogeneous in clinical context.
`Normal_histology` is a NAFLD-cohort patient with no fibrosis yet, not a donor
liver; `Control_external_study` is the control arm of an alcohol study under a
different protocol.

**F0 is not folded in.** F0 is a histological stage in a NAFLD patient — a
measurement — while these three are the absence of the disease context. Merging
them would assume what the `Controles` vs `F0` contrast exists to test.

Two of GSE135251's ten controls carried incidental fibrosis of stage 1 and 2
(`GSM3998224`, `GSM3998341`). They are excluded by `exclude_samples` in that
dataset's config, which is why the count above is 8 and the class is 56 rather
than 58. `exclude_samples` aborts if a named accession is not in the series — a
typo would otherwise exclude nothing and leave the run looking corrected — and
warns instead of aborting when none of them match, which is the legitimate case
for synthetic or subset data.
Excluding by accession rather than by a rule on `fibrosis stage` matters:
`sample_filter` evaluates one column across all samples, so a filter on stages 1
and 2 would also delete the 47 F1 and 53 F2 of the NAFLD cohort. This departs
from the publication, which groups all ten as controls — the control class here
is histologically clean, not the paper's recruitment group.

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
./tsf run --to spectra --results-dir results_proteincoding
./tsf run --gene-universe '^(protein-coding|ncRNA)$' \
          --results-dir results_pc_lnc --interim-dir interim_pc_lnc
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

The spectrum of the mean profile is not a summary of the per-sample spectra. The
transform is linear, so the spectrum of the mean is the *vector* mean of the
complex coefficients: a component present in every sample at scattered phases
cancels and disappears. `./tsf consensus` therefore works from the per-sample
spectra and reports three things the mean cannot separate — how strong (median
normalised power), how common (prevalence), how aligned (phase-locking value) —
plus a permutation null built by drawing samples at random from the whole
dataset, ignoring condition.

```bash
./tsf consensus --n-null 999 --n-boot 200
```

A component is **confirmed** when prevalence holds, phase alignment is unlikely
under uniform phases, and it beats the label-permuted null family-wise.
Everything else is **exploratory**. Clearing zero is not evidence: the score is a
product of non-negative quantities, so any signal at all clears it.

Two warnings from this stage are the pipeline saying a claim is *not reachable*,
which is different from absent — a permutation p cannot go below 1/(B+1), so a
family with more members than draws can allow has a smallest attainable q above
any threshold. The logs print that number.

## The condition library

`scripts/build_final_condition_spectra.R` turns the per-cohort consensus spectra
into one table per class, aggregating by the **median across cohorts** rather
than by pooling samples: with 42 of 70 F4 samples from one series, a mean would
place the F4 spectrum where that series is.

```bash
Rscript scripts/build_final_condition_spectra.R \
  --results-dir results --out-dir results/library_all

Rscript scripts/build_final_condition_spectra.R \
  --results-dir results --out-dir results/library_domains \
  --min-period auto --period-margin 2 --min-period-biological 10
```

The signature is not a top-N list. Each cohort's family-wise permutation
p-values combine by Stouffer and are BH-adjusted across frequencies; membership
is `q_meta_null <= 0.05`. Combining across cohorts is what makes that cut
reachable — with four cohorts it works at 99 draws, with three at 199, with two
it needs many more, and **with one it is impossible at any number of draws**.
Those classes are reported `single_cohort` and provisional.

`--min-period` removes short periods from the analysis *and from the testing
family*. That is legitimate because the period is a property of the grid alone,
so the filter can be fixed before a spectrum is seen; filtering by enrichment or
prevalence would not be. The effective floor is the larger of a technical one
(`2/coverage + margin`, the Nyquist limit corrected for the gaps actually
present) and a biological one (the scale the mechanisms of interest could
produce). Report results with and without it: a component that appears only with
the floor is not a finding of the floor.

chrY and the mitochondrion are excluded by default. chrY's expression tracks the
sex composition of a group, which differs between conditions and cohorts, so a
component there reports who was recruited.

Three figures per class: `power_spectrum_` (every frequency, selected ones
marked), `signature_peaks_` (selected only, one panel per chromosome, on a period
axis so chromosomes are comparable), and `genome_spectrum_` (all chromosomes on
one axis; that axis is a running frequency index, not a genomic coordinate).

## What decides that a peak exists

Two criteria are computed for every peak; `stability_criterion` in
`config/project.R` picks which one drives `is_stable`, and both are written
either way.

**condition** (default) -- a permutation test on the condition's summary signal.
Values are permuted among the observed grid positions, the spectrum recomputed,
and the *maximum* power over frequencies gives the null. Family-wise: it asks
whether a frequency beats the strongest frequency of a permuted spectrum, which
on 500-800 frequencies is close to asking whether it dominates its chromosome.
Expect single digits, and believe what passes.

**condition_fdr** -- the same permutations, read pointwise: does this frequency
beat its own null? BH across frequencies then controls the false discovery rate
instead of the probability of any error. This is the question a signature is
about, and it selects far more, a fraction of them false by construction.

The pointwise p needs a *pooled* null, because against its own frequency it
inherits the `1/(B+1)` floor and BH over ~10,000 frequencies would need
B = 200,000. Each frequency's null is standardised by its own mean and standard
deviation and the standardised values are pooled across the chromosome. The
assumption is that the standardised null is exchangeable across frequencies;
what standardising does *not* remove is the sampling window, which is why
`./tsf window` is not optional.

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

For the family-wise criterion, multiplicity is corrected over **chromosomes, not
frequencies**: maxT already controls the family-wise error rate across
frequencies within a chromosome.

### The permutation floor

A permutation p-value cannot go below `1/(B+1)`. Any procedure that multiplies it
by the size of a family therefore has a **smallest attainable q**, and if that
exceeds the threshold, nothing can pass however strong the signal — a result of
"zero components" would then say nothing about the data. This has bitten in five
different places in this pipeline, so every one of them now computes the
attainable floor and prints it:

```
Pointwise null: with 20 draws over 297 frequencies the smallest reachable
BH q is 14 > 0.05, so q_null cannot confirm anything.

calibrated cut: 999 draws, 3 cohort(s), 5013 frequencies ->
smallest reachable BH q = 0.0004 (usable)
```

**"Not reachable" and "not present" are different findings.** Read the floor
before reading a zero.

## Running on Colab

`TissueSpectF_colab.ipynb` is where to try a change online: it runs the tests
and the self-check with no data at all, then fetches the five cohorts, checks the
labels against the cohort table, and runs the spectral stages restricted with
`--chromosomes` so they finish. It ends with the baselines, which are worth
running there even when nothing else is.

Two cores make the per-sample `maxt` stage impractical, so the notebook skips
it; the rank-based statistics carry the analysis without it. **Nothing produced
on a partial genome is a result** — it is a rehearsal that proves the code runs
and the labels are right. Two cores make the per-sample
`maxt` stage impractical (roughly 8 hours), so the notebook runs
`--from=condition`, which is about 1/n the cost and answers the primary
question. The stability stage falls back to the condition test on its own when
no per-sample maxT exists.

## Checking that it works

```bash
make test         # every suite: labels, config, FFT, maxT, CLEAN, stability, Wilson
make test-r       # the R suites alone (base R, no install step)
make test-ml      # the learned layer's tests (skipped if Python is absent)
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

**Early phase.** This is the one part of the repository that is not finished, so
what exists is stated plainly rather than described as though it were done:

| | status |
|---|---|
| `ml/dataset.py`, `ml/splits.py` | done — export to masked tensors, leave-one-cohort-out, leak checks |
| `ml/prototypes.py`, `ml/evaluate.py`, `ml/baselines.py` | done — cohort-balanced prototypes, metrics, the baselines |
| `tests/ml/` | done — 46 tests, run in CI |
| the model itself | **not written** |
| `config/autoencoder.yaml` | **read by nothing yet** — it records the intended parameters, not a configuration in use |

So the order is deliberate: the data layer and the bar the model has to clear
come first, and the model comes last. Nothing here is on the path of any result
the statistical pipeline produces.

A complementary layer, not a replacement. The statistical results — permutation
p-values, consensus spectra, replication across cohorts — remain the evidence.

> **A peak reconstructed by TissueSpect-AE is not statistical or causal evidence
> on its own. Components only the model finds are reported as AI candidates.**

Neither the data layer nor the harness needs `torch`; the R pipeline and its
tests do not need Python at all.

```bash
pip install -r requirements-ml.txt          # optional
./tsf ae-prepare --interim-dir interim --results-dir results
python3 scripts/run_baselines.py --data results/autoencoder/data \
                                 --out  results/autoencoder/baselines
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

## Sonification

```bash
pip install -r requirements-sonify.txt
make sonify
```

One MIDI composition per class. Core invariant peaks — the components present in
every condition — become a shared accompaniment; each condition's consensus
peaks become a melody over it. So what you hear changing between `F0` and `F4`
is the part of the spectrum that changes, against a bed that does not.

The mapping:

| spectral quantity | musical quantity |
|---|---|
| frequency `k/N` | pitch, logarithmically — higher spectral frequency, higher pitch |
| phase | a sub-beat timing displacement |
| peak strength (`invariant_score`, `meta_score`, `final_power`) | velocity |
| period | note duration, inversely — long genomic periods, longer notes |
| evidence class (`robust` / `candidate`) | articulation; candidates are softer |

Pitches are snapped to D Dorian. That is the one frankly aesthetic decision in
the chain: it makes the result listenable, and it does discard information.
Frequency still sets pitch height; snapping only quantises it.

Every mapping table is written next to the MIDI, so the transformation stays
inspectable — a note can be traced back to the `(chr, N, k)` it came from. The
run also writes a manifest and a `README_sonification.txt`.

```bash
python3 scripts/sonify_tissuespectf.py \
  --library-dir results/library_domains \
  --conditions Normal_histology,F0,F1,F2,F3,F4 \
  --bpm 82 --bars 16
```

`--library-dir` defaults to `$TSF_LIBRARY_DIR`, then
`$TSF_RESULTS_DIR/library_domains`, then `<repo>/results/library_domains`. It
needs the condition library, so run `build_final_condition_spectra.R` first.

**This is not evidence.** A sonification is a presentation of a result, not a
test of one. Nothing audible here supports a claim that the statistical
pipeline does not already support on its own, and a component that sounds
striking is not thereby more real. The evidence is the permutation p-values,
the consensus spectra and the cross-cohort replication.

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
`n_masks`, which is passed through to `./tsf reference --n-masks`, for a cheaper
trial build.

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

### Which spectrum the fingerprint uses

`--features` selects the representation. All of them go through the same
leave-one-cohort-out validation, the same training-only feature selection and
the same thresholds, so comparing them isolates the representation rather than
confounding it with the harness.

| `--features` | features | in one line |
|---|---|---|
| `amplitude` | ~1,500 | log-amplitude per `(chromosome, k)`, the default |
| `amplitude_phase` | ~3,000 | the above, plus amplitude-weighted phase |
| `period_bins` | ~960 | per `(chromosome, period)` on a common log grid |
| `period_bins_genomic` | 40 | one curve over period, averaged across chromosomes |
| `band_ratios` | 780 | ratios between period bands; invariant to global scale |
| `expression_baseline` | ~14,000 | the **control**: raw expression, no transform |

`(chromosome, k)` mixes incomparable scales. `k` is cycles per chromosome, so
with N = 2066 on chr1 the index `k = 64` is a period of 32 genes and with
N = 759 on chr21 it is 12: `chr1_k64` and `chr21_k64` are columns the model
treats as parallel while they describe different things. Indexing by period
fixes that, and only then can chromosomes be averaged into one curve.

`k_max` does not apply to the period representations — it caps cycles per
chromosome, which would empty every short-period bin on the long chromosomes.

Run the sweep before trusting any single number:

```bash
for f in amplitude period_bins period_bins_genomic band_ratios expression_baseline; do
  ./tsf reference GSE130970 GSE135251 GSE142530 GSE162694 GSE276114 \
    --config config/project.R --geo-dir data --interim-dir interim \
    --results-dir results_feat_$f --cores 24 --seed 42 --features $f
done
```

`expression_baseline` is the one that makes the others interpretable. If raw
expression classifies these cohorts better than the spectrum, the transform is
discarding information and the write-up has to say so; if it classifies worse,
the spectrum is a real compression. Keep `--n-features` identical across the
sweep — changing it for one invalidates the comparison.

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

## Licence

**CC BY-NC 4.0.** See [LICENSE](LICENSE). All of it — the pipeline, the app,
the learned layer and the sonification. One licence, no exceptions.

Use, modify and share it for any non-commercial purpose with attribution;
academic research, teaching, peer review and reproducing a published result
all qualify and need no permission. Commercial use needs written permission.

Two consequences worth knowing rather than discovering: CC is not a software
licence and carries no patent or warranty provisions, and NC is not
OSI-approved, so a repository that has to be deposited under an open licence
cannot be deposited as-is.
