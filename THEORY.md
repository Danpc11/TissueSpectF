# TissueSpectF — model and theory

What the pipeline assumes, what it estimates, what each test licenses you to
claim, and where the whole construction can fail. Implementation lives in `R/`;
every section below names the file that implements it.

---

## 1. The object of study

Transcription is not spatially random along a chromosome. Co-regulated genes
cluster, chromatin is organised into domains, and replication timing and lamina
association vary along the arm. If any of that leaves a trace in steady-state
mRNA abundance, expression read *in genomic order* should carry positional
structure, and some of that structure should be quasi-periodic.

The question the pipeline asks is therefore:

> Ordered along the chromosome, does gene expression contain periodic components
> at scales of tens to hundreds of genes; are those components characteristic of
> a tissue and of a disease state; and do they replicate across independent
> cohorts?

This is a claim about *organisation*, not about individual genes. A component
with a period of 200 genes says nothing about any single gene in the crest.

---

## 2. The axis

### 2.1 Rank, not base pairs

Position is the **index of a gene on a reference grid**, not its coordinate in
base pairs. Gene density varies by more than an order of magnitude along a
chromosome, so a periodicity in base pairs and a periodicity in genes are
different hypotheses. We test the second: "every ~200 genes", not "every ~5 Mb".

The choice is defensible on mechanism — co-regulation acts on genes, and the
number of genes between two co-regulated blocks is the quantity a domain model
predicts — but it is a choice, and a result in gene units does not transfer to
base-pair units.

### 2.2 The reference grid

Let $\mathcal{G}_c$ be the set of annotated genes on chromosome $c$ whose
biotype is in the declared universe (`gene_universe`, default
`^protein-coding$`), ordered by start coordinate and indexed

$$t = 1, \dots, N_c .$$

$N_c$ is a property of the **annotation**, not of any experiment. Two datasets
processed with the same annotation therefore live in the same index space, and
a peak key $(c, N_c, k)$ means the same thing in both. This is what makes
cross-dataset comparison possible at all.

Implemented in `build_reference_grid()` (`R/grid.R`).

> **Why this matters, concretely.** An earlier version used the rank of a gene
> among those that passed the expression filter. That axis is dataset-specific
> ($N$ = 14,495 vs 16,269 in our two cohorts) and, worse, *non-uniform with
> respect to biology*: dropping a gene makes its neighbours adjacent, so the
> axis compresses wherever coverage is poor. A component that is smooth in grid
> space folds into high frequency in compacted space. On real data this produced
> a chromosome-17 peak at a period of 2.18 genes, present in every condition
> including healthy controls, which vanished entirely once the axis was fixed to
> the annotation. It was an artefact of the axis, not of the biology.

### 2.3 Missing genes

A gene with no usable measurement — filtered for low expression, absent from the
count matrix, lacking a valid length — keeps its grid slot and is **unobserved**.
It is never zero-filled and never imputed. Absence of a measurement is not
evidence of zero expression, and a zero placed at a fixed position injects a
deterministic pattern whose own spectral structure can manufacture a peak.

Write $T_c \subseteq \{1,\dots,N_c\}$ for the observed positions,
$n_c = |T_c|$, and $\rho_c = n_c/N_c$ for coverage. In our liver cohorts
$\rho$ ranges from 0.16 to 0.78 with medians of 0.66 and 0.73.

---

## 3. The signal

For sample $s$ and gene at grid position $t \in T_c$,

$$
y_s(t)=\mathrm{asinh}\!\left(\mathrm{TPM}_s(t)\right).
$$

`asinh` rather than `log(x+1)`: it is defined at zero, behaves linearly near
zero and logarithmically in the tail, and needs no pseudocount whose value would
have to be justified. TPM when transcript lengths are available, CPM otherwise —
the manifest records which, so nothing downstream assumes length normalisation
it did not get.

A condition $g$ is summarised by the pointwise mean or median over its samples:

$$\bar y_g(t) = \frac{1}{|S_g|}\sum_{s \in S_g} y_s(t).$$

Both branches (`average`, `median`) are carried through the whole pipeline. The
condition spectrum is the spectrum *of the summary signal*, never the average of
per-sample spectra — averaging spectra discards phase and would conflate
components that are aligned across samples with components that are not.

---

## 4. Spectral estimation on an incomplete grid

### 4.1 Model

For frequency index $k$, with $\omega_k = 2\pi k/N$ and positions taken as
$t_0 = t-1$ (the DFT's zero-based convention):

$$y(t) = \mu + a_k \cos(\omega_k t_0) + b_k \sin(\omega_k t_0) + \varepsilon(t),
\qquad t \in T,$$

fitted by ordinary least squares over the observed positions. This is the
generalised Lomb–Scargle periodogram with a floating mean (Zechmeister &
Kürster 2009). Reported quantities:

$$
A_k=\sqrt{a_k^2+b_k^2}, \qquad
\phi_k=\mathrm{atan2}\!\left(-b_k,a_k\right), \qquad
P_k=\frac{1}{2}A_k^2.
$$

so that the fitted component is $A_k\cos(\omega_k t_0 + \phi_k)$, and

$$f_k = k/N \ \text{(cycles per gene)}, \qquad \Pi_k = N/k \ \text{(genes per cycle)}.$$

### 4.2 This is still Fourier analysis

On a complete grid ($T = \{1,\dots,N\}$) the sinusoids are orthogonal, the least
squares solution is the DFT coefficient, and $A_k, \phi_k, P_k$ coincide with the
one-sided FFT values term by term. `tests/test_spectrum.R` checks this against
`fft()` to $10^{-9}$, and checks the gapped case against an explicit `qr.solve`
of the design matrix. Lomb–Scargle is not an alternative to Fourier analysis; it
is the least-squares spectral estimator, of which the DFT is the complete-grid
special case.

### 4.3 Computing it in $O(N\log N)$

The positions are integers on a regular grid, which the general Lomb–Scargle
case does not enjoy, and that buys an exact fast algorithm. Let $z$ be the
signal zero-filled at missing positions and $w$ the 0/1 presence indicator.
Then, with $Z = \mathrm{FFT}(z)$ and $W = \mathrm{FFT}(w)$:

$$\textstyle\sum_T y\cos = \Re Z_k, \qquad \sum_T y\sin = -\Im Z_k,$$
$$\textstyle\sum_T \cos = \Re W_k, \qquad \sum_T \sin = -\Im W_k,$$
$$\textstyle\sum_T \cos^2 = \tfrac{1}{2}(n + \Re W_{2k}), \quad
\sum_T \sin^2 = \tfrac{1}{2}(n - \Re W_{2k}), \quad
\sum_T \cos\sin = -\tfrac{1}{2}\Im W_{2k},$$

using $\cos^2\theta = (1+\cos 2\theta)/2$ etc., with the $2k$ index taken
modulo $N$. The floating-mean (hatted) sums follow, and

$$\begin{pmatrix} a_k \\ b_k \end{pmatrix} =
\frac{1}{D_k}\begin{pmatrix} \hat{SS}_k & -\hat{CS}_k \\ -\hat{CS}_k & \hat{CC}_k \end{pmatrix}
\begin{pmatrix} \hat{YC}_k \\ \hat{YS}_k \end{pmatrix},
\qquad D_k = \hat{CC}_k\hat{SS}_k - \hat{CS}_k^2 .$$

**The zero-filling here is a computational device, not imputation.** The missing
positions enter the normalisation through $w$ — that is precisely what naive
zero-filling omits, and why naive zero-filling is biased. The result is
identical to fitting only the observed points, at the cost of two FFTs instead
of $O(N \cdot K)$ inner products. Without this, the permutation tests below
would be roughly a thousand times slower and the analysis would not be
computable on the hardware available.

Implemented in `gls_window_terms()` and `gls_spectrum()` (`R/grid.R`).

### 4.4 What gaps cost

Three properties of the complete-grid transform are lost, and each has a
consequence the pipeline has to handle explicitly.

**Orthogonality.** The sinusoids are not orthogonal over $T$, so frequency bins
are not independent and power leaks between them. *Consequence:* per-bin
p-values would be badly calibrated; the family-wise permutation null of §5 is
mandatory, not optional.

**Parseval.** $\sum_k P_k \neq \mathrm{Var}(y)$ exactly, so normalised
power is not a clean fraction of variance. *Consequence:* `power_normalised` is
reported as the classical GLS statistic in $[0,1]$, not as a variance share.

**Nyquist.** With irregular sampling the Nyquist limit is not uniquely defined:
structure above the naive limit is partly recoverable, and aliasing appears
where it could not before. *Consequence:* the spectral window, next.

### 4.5 The spectral window

$$W(f_k) = \frac{|\sum_{t\in T} e^{-i\omega_k t_0}|^2}{n^2}$$

is the spectrum of the sampling pattern alone — what the gaps produce with no
expression signal whatsoever. A complete grid gives $W \equiv 0$ off zero; a
gapped grid does not. If an observed peak sits at a frequency where $W$ is
large, that peak may be an alias of structure at another frequency rather than
structure at its own.

Every peak is therefore reported with `window_power` and `window_rank`, and
`./tsf window` places all selected peaks inside the window spectrum. **A peak in
the top 1% of the window is treated as a sampling artefact until shown
otherwise.** Because coverage differs between cohorts, the window differs too —
so a component that replicates across datasets replicates *under different
sampling patterns*, which is stronger evidence than replication under the same
one.

---

## 5. Inference: is a component real?

### 5.1 Null hypothesis and exchangeability

$H_0$: the values are exchangeable across the observed grid positions — the
expression profile carries no positional structure. Under $H_0$, permuting the
values among $T$ leaves the distribution unchanged.

Two properties of this null matter:

- **Positions are held fixed.** The sampling pattern is identical in the data
  and in every permutation, so the pattern of missing genes cannot by itself
  produce significance. `tests/test_spectrum.R` verifies this directly: white
  noise on a grid with 50% gaps yields significance at the nominal rate.
- **It is distribution-free.** No assumption about the marginal distribution of
  `asinh` TPM, which is skewed and heavy-tailed.

### 5.2 Family-wise control by maxT

For each permutation $b$, record $M_b = \max_k P_k^{(b)}$, and

$$
p_k=\frac{1+\lvert\{b:M_b\ge P_k\}\rvert}{B+1}.
$$

Comparing every observed power against the null *maximum* controls the
family-wise error rate across all frequencies within a chromosome
(Westfall & Young). The $+1$ in numerator and denominator keeps the p-value
valid and bounds it below by $1/(B+1)$ — a floor that has consequences in §5.5.

### 5.3 Two permutation schemes

**full** permutes all observed values: it destroys structure at every scale,
including local autocorrelation, and is therefore the most permissive null.
A peak significant under `full` may still be nothing more than local
correlation.

**block** permutes intervals of the **reference grid** — blocks defined by
$\lceil t / \text{block size} \rceil$, not by runs of consecutive entries in the
observed list. With gaps the two differ badly: at 20% coverage, ten consecutive
observed genes span about fifty grid positions, so list-index blocks would treat
distant genes as neighbours and preserve no locality at all. Blocks hold a
variable number of observed genes; the permutation shuffles block order and
writes the concatenated values back onto the same observed positions.

*Known compromise:* when block sizes differ, within-block spacing is not
preserved exactly after relocation. No permutation preserves both the position
set and every within-block distance; holding positions fixed is the property
that keeps the null valid, so that is the one kept.

A scheme needs at least 10 blocks. With 4 blocks there are only 24 orderings, so
a sizeable fraction of null draws leave the signal nearly intact, no p-value can
get small, and the scheme would veto every peak — not because the peak is local
structure but because the test cannot see. Schemes below the floor are skipped
and recorded in `block_schemes_skipped`.

`maxt$primary_scheme` declares which null decides: `full`, or `all` (the maximum
across schemes, requiring the peak to survive the block nulls too). **For any
claim about periodicity, `all` is the defensible choice**; `full` remains the
default only so that changing it is a deliberate act.

### 5.4 Condition-level test — the primary inference

The per-sample criterion originally used — "significant in $\ge 90\%$ of
samples" — is a *consistency requirement*, not a test with aggregated power.
Each sample must reach significance alone under family-wise control, so evidence
is never pooled: a periodicity that is real but moderate, present in every
sample and individually short of significance in each, stays invisible however
many samples the condition has. Adding samples does not increase power, which is
backwards.

The primary test is therefore applied to the summary signal $\bar y_g$ — the
signal that is actually analysed downstream — with the same maxT logic, once per
condition instead of once per sample. Cost falls by a factor of $|S_g|$, which
is also what makes a complete run feasible on two cores.

Where per-sample results exist, Stouffer's combination

$$Z = \frac{1}{\sqrt{m}}\sum_{i=1}^{m}\Phi^{-1}(1-p_i)$$

is reported as a second opinion. It aggregates evidence and gains power with
$m$, but it is conservative (its inputs are already family-wise adjusted) and
assumes independence between samples, which biological replicates satisfy only
approximately. It is not the primary.

The consistency fraction survives as a **reproducibility descriptor** reported
next to the test — never as the thing that decides what exists.

### 5.5 Multiplicity: the family is chromosomes

maxT already controls the family-wise error rate across frequencies *within* a
chromosome. The remaining multiplicity is the number of chromosomes tested, so
that is what is corrected:

$$q = \min(1,\ p \cdot n_{\mathrm{chr}}).$$

Applying BH across all ~9,500 frequencies on top of maxT would not be merely
conservative, it would be **fatal**: with $B = 2000$ the smallest attainable
p-value is $5\times10^{-4}$, so the smallest attainable $q$ would be
$5\times10^{-4} \times 9500 = 4.7$, capped at 1. No peak could reach any
threshold whatever the data. (This was a real bug, caught because the test
returned zero peaks on synthetic data with a known injected signal.)

This puts a floor on $B$: $q \le 0.05$ requires
$B \ge 20\,n_{\mathrm{chr}} - 1$, about 460 for 23 chromosomes. The stage warns
when $B$ is too small to be able to conclude.

---

### 5.6 The characteristic spectrum of a condition

The spectrum of the mean profile is not the same object as a summary of the
per-sample spectra. The transform is linear, so the spectrum of the mean is the
**vector** mean of the complex coefficients: a component present in every sample
at the same frequency but with scattered phases cancels in that sum and vanishes
from the condition spectrum, however reproducible it is. One extreme sample can
equally carry a peak no other sample has.

A condition is therefore summarised on three axes the mean profile cannot
separate (`R/consensus.R`):

| axis | statistic | reads as |
|---|---|---|
| how strong | median normalised power across samples | robust to one outlier sample |
| how common | prevalence: fraction of samples where the frequency stands out within its own sample and chromosome | independent of depth and of chromosome |
| how aligned | phase-locking value | 1 when every sample puts the crest in the same place |

$$
\mathrm{PLV}=\left\lvert \frac{1}{n}\sum_{s=1}^{n} e^{i\phi_s} \right\rvert,
\qquad
S = \tilde P_{\mathrm{norm}} \cdot \pi \cdot \mathrm{PLV}.
$$

$S$ is a product of median normalised power, prevalence $\pi$ and $\mathrm{PLV}$,
so a component must satisfy all three: any factor near zero sends the score to
zero. Normalised power rather than raw, because raw power differs by orders of
magnitude between chromosomes and the score would otherwise rank chromosomes
instead of components. Prevalence is judged within a sample **and** a
chromosome, for the same reason.

Two cautions. Under uniformly random phases

$$
\mathbb{E}\left[\mathrm{PLV}\right] \approx \frac{\sqrt{\pi}}{2\sqrt{n}},
$$

which is not zero — with 8 samples a $\mathrm{PLV}$ of $0.3$ is unremarkable —
so a Rayleigh p-value accompanies every value. And that p-value cannot fall
below $e^{-n}$, so after BH over $n_f$ frequencies the smallest reachable $q$ is
$e^{-n} n_f$: below roughly $\log(n_f/q)$ samples, about 9 for 5,000
frequencies, perfect alignment could not pass. The stage says so and falls back
to prevalence and score rather than returning an empty signature.

Prevalence has two definitions and both are reported, never interchanged. The
**rank** definition — above this sample's own 95th percentile for this
chromosome — needs nothing but the spectra. The **maxT** definition — significant
in that sample — is stronger evidence but exists only where the per-sample
permutation test was run, which a permuted null cannot assume. Comparing an
observed score built on one against a null built on the other is not a
permutation test, whatever the direction of the bias, so
`consensus_score_rank` is used on both sides of the permuted comparison and
`consensus_score_maxt` is reported beside it as confirmatory evidence.

The null draws $n$ samples at random from the whole dataset, ignoring condition,
and scores them the same way. Two summaries of it are kept:

- **per component**, an empirical p-value at each $(\mathrm{chr}, k)$ against
  its own null, BH-adjusted across components. This is what `confirmed` rests
  on, and it is what lets a signature be *localised*.
- **global maximum**, the 95th percentile of the best score over all frequencies
  of each draw. This controls error over the whole signature and is very strict;
  in practice it confirms almost nothing, which is the right answer to "is there
  any component at all" and useless for "which ones". Kept as
  `beats_global_null`.

Clearing zero is not evidence: the score is a product of non-negative
quantities, so any signal at all clears it. Clearing the null means the
component belongs to *this condition* rather than to any group of samples of
that size — the confound that matters, since tissue-wide structure has high
prevalence and high phase locking too. Everything else is **exploratory**.

The score is a **ranking statistic, not a test**. Bootstrap intervals over
samples say how stable it is; the Rayleigh p-value covers phase alignment alone.
Neither makes it a significance claim — §5.4 is for that. Stouffer stays a
secondary analysis; Fisher's method is not used, because its sensitivity to a
single small p-value is the wrong behaviour when the question is whether a
component is *shared*.

## 6. Decomposition: which components are there?

### 6.1 Why thresholding the periodogram is wrong

Because the sinusoids are not orthogonal, a strong component leaks into
neighbouring frequencies and into window aliases. Those sidelobes are local
maxima of the periodogram without being components of the signal. Taking the
$M$ highest peaks returns the real component *together with its own artefacts*
and counts them as independent structure. On a 50%-coverage grid with one
injected component, the top five raw peaks contain four false maxima.

### 6.2 CLEAN / orthogonal matching pursuit

Iteratively: take the strongest frequency in the current residual, refit **all**
selected components jointly by least squares, subtract the fit, search the new
residual. Each component is chosen against what remains rather than against what
the previous one contaminated. Amplitudes and phases reported are those of the
final joint fit.

### 6.3 Stopping: extended BIC

The number of components is not chosen. Extraction continues while a component
lowers the criterion:

$$\mathrm{EBIC} = n\log\frac{\mathrm{RSS}}{n} + p\log n + 2\gamma\log\binom{K}{m},
\qquad p = 2m+1 .$$

The third term is essential and is the whole reason this is not a heuristic in
disguise. Plain BIC assumes a model fixed in advance, whereas each component is
the maximum over $K \approx N/2$ candidate frequencies; that search is unpaid
for. Measured on pure noise, $\gamma = 0$ selects at least one component in
about 65% of draws; $\gamma = 1$ selects one in none. Chen & Chen's extended BIC
charges the selection (Chen & Chen 2008).

$\gamma$ remains a declared parameter, but it is of a different kind from a
chosen $M$: it has a justified range, and its false-positive rate under the null
is measurable and is measured (`TSF_EBIC_GAMMA` exists for exactly that
sensitivity analysis).

**EBIC selection is not a significance test and carries no error rate.** "The
data prefer to keep this component" and "this periodicity is unlikely under a
null" are different claims. The two procedures disagree by design: CLEAN reports
components the permutation test calls non-significant. That is correct — it is
what a fingerprint needs and what an inferential claim must not use.

---

## 7. From component to genes

A component $(c, N, k, A, \phi)$ is projected back onto the chromosome as

$$r(t) = A\cos(\omega_k t_0 + \phi), \qquad t \in T_c,$$

evaluated at the observed positions. Genes where $r > 0$ sit on a crest of that
periodic component, genes where $r < 0$ in a trough.

This is a projection, not a per-gene statement. A gene on a crest is not
"periodically expressed"; it occupies a position where the fitted component is
positive. Gene-set analysis of crest membership is a legitimate downstream
question — it is what would connect a component to biology — but it inherits no
error control from the spectral test and needs its own.

---

## 8. Replication across cohorts

Datasets are compared **only within a tissue and a vocabulary**, and results are
reported side by side, never pooled. Pooling proportions over summed counts
would be dominated by the larger cohort and would hide heterogeneity; the
question is replication, not a bigger $n$.

For an ordered vocabulary, each transition $g \to g'$ reports the amplitude
delta, its $\log_2$ ratio, and a Welch test on per-sample peak power. A peak is
`replicated` when the direction of change agrees across all datasets and
FDR $\le 0.05$ in each.

Two independent cohorts with different coverage — hence different spectral
windows — agreeing on both the frequency and the direction of change is the
strongest evidence this design can produce.

---

## 9. Generalisation beyond liver fibrosis

Nothing in `R/` knows what a fibrosis stage is. The condition vocabulary is
data: `config/vocabularies/` declares levels, whether they are ordered, and
which is the baseline.

- `liver_fibrosis` — `Control, F0..F4`, ordered, baseline `Control`
- `case_control` — two groups, any tissue, unordered
- `tissue_atlas` — one level per tissue, unordered, no baseline

Unordered vocabularies produce no transition tables: "liver → lung" is not a
transition, and the pipeline knows it. Adding a tissue is a vocabulary file and
a dataset file; no code changes.

One invariant carries the generalisation: **`Control` is a cohort statement, not
a histology statement.** It is never inferred from a stage of 0 or from a
description of normal tissue. A dataset with no non-disease cohort declares
`has_control_cohort = FALSE`, and the label engine refuses to emit the baseline
for it. Without that rule, `Control` in one cohort (healthy donors) and
`Control` in another (biopsies with normal histology from the disease cohort)
would merge into one label meaning two different things.

---

## 10. Fingerprinting: a different question

"Is this component real?" and "does this spectrum identify the tissue?" are
distinct, and conflating them costs power in both directions.

A fingerprint does **not** require any individual component to be significant. A
set of frequencies can discriminate classes while each one sits in the noise.
Filtering by $q \le 0.05$ before building a fingerprint discards precisely the
information that would make it work. The per-sample CLEAN output
(`clean$per_sample = TRUE`) is the sparse representation for this purpose:
components with amplitude *and phase*, no threshold.

Phase deserves emphasis. The classical audio-fingerprinting approach discards
it. Here it should not be discarded: amplitude alone is invariant to translation
along the chromosome, telling you a 200-gene structure exists but not where its
crest falls; phase anchors the structure to genomic coordinates. For tissue
identity, phase is plausibly the more discriminative half.

Coverage is calibrated, not assumed. A query rarely covers the whole grid, and
the similarity distribution shifts with coverage, so a threshold calibrated on
complete fingerprints rejects members of a 60%-covered query that should be
accepted. Validation therefore re-scores every held-out sample under simulated loss and
calibrates accuracy and rejection threshold per band: 90–100%, 75–90%, 50–75%,
and below 50%, where nothing is classified.

What is masked are **genes on the grid**, and the fingerprint is recomputed from
the survivors. Masking spectral features instead would measure a different and
easier quantity: dropping half the genes of a chromosome still leaves the GLS
able to estimate almost every frequency, so feature coverage stays near 100%
while gene coverage is 50% — every estimate merely becomes noisier and the
spectral window changes. Bands are keyed on gene coverage **relative to the canonical grid**, which is
what a query can report about itself before anything is computed. The
denominator matters: a held-out cohort covering 70% of the grid, masked to keep
80% of its own genes, has 56% coverage of the reference, not 80%. Calibrating
against the dataset's own denominator would place every simulated query in an
easier band than the real query it stands for. Each record keeps
`baseline_dataset_coverage`, `mask_retention` and
`absolute_reference_coverage` so the three are never confused again.

Loss is simulated five ways — scattered genes, one retained contiguous block,
several missing disjoint blocks, whole chromosomes, and dropout of the least
expressed genes — repeated over many independent masks (`fingerprint$n_masks`),
since *which* regions are missing matters as much as how many. The dropout mode
is the most realistic: coverage loss is not independent of expression, and the
genes that drop out cluster in the tissue-specific families that carry much of
the between-condition signal.

Thresholds are calibrated per band **and per predicted class** where a class has
enough correct held-out matches to support one, since classes differ in how
tight their centroids are; the band-level threshold is the fallback. Two
band-level thresholds are computed and both rejection rates are reported: the
pooled quantile over all masks, and a conservative one at the 90th percentile of
the per-mask thresholds. Which is applied is declared in
`fingerprint$threshold_policy`. The default is pooled, because with gene-level
masking the conservative threshold was measured to reject 46–79% of true
members — safety bought at a price that is not worth paying, and a number that
only became visible once the masking was done on genes rather than on
frequencies. A query whose band has no calibrated threshold is reported
`UNCALIBRATED_COVERAGE` rather than judged against the wrong number.

The reference is self-contained — it carries the canonical annotation grid, gene
identifiers, species, genome build, annotation release, frequency ceiling,
feature representation and expression unit — so a query is never scored against
a grid or a scale other than the one the reference was built on. Similarity is
computed over the shared features only: a frequency the query never observed
contributes nothing, rather than contributing the mean of a z-scored vector.

**The decisive experiment**, implemented in `R/reference.R` as leave-one-cohort-
out validation: train on one cohort's fingerprints and predict the other
cohort's labels, with no tuning against the second. Internal cross-validation cannot substitute for this — with
thousands of features and hundreds of samples, a classifier will achieve a high
internal AUC by learning batch, library protocol and sequencing depth. Batch
does not transfer between studies; a real spectral signature does.

The fingerprint is built by `R/fingerprint.R` (log-amplitude, optionally
amplitude-weighted phase, per `(chromosome, k)` up to `k_max`), normalised per
sample so that library depth cannot drive a match, and classified by nearest
centroid on features selected from training data only. The model is deliberately
inflexible: with thousands of features, anything richer would fit cohort
idiosyncrasy and the out-of-cohort number would stop meaning anything.

`reference_status()` turns the validation into one sentence, and both the CLI
and the app print it above every result. A reference that does not beat the
majority-class baseline is reported as carrying no information, rather than
silently returning a confident label.

Anticipated outcome worth naming in advance: classifying *tissue* is likely to
work far better than classifying *stage*, because the identity of a tissue is a
much larger difference than five grades of fibrosis within one. That result
would still be the useful one.

---

## 11. Assumptions, in the order they could break

1. **Gene order is the right axis.** If the relevant scale is physical distance,
   a gene-rank analysis will find nothing, or find something with no mechanistic
   reading.
2. **The gene universe is a choice.** Coding-only and coding-plus-noncoding
   grids give different $N$, hence different periods for the same structure. Run
   both; a component surviving both is robust to the definition, one that does
   not tells you the result depends on which genes you count.
3. **Coverage is not random.** Filtered-out genes are the lowly expressed ones,
   which cluster (tissue-specific families, gene deserts). The permutation null
   is immune to this by construction — positions are fixed — but the *aliasing*
   it induces is not, which is what the spectral window exists to expose.
4. **Samples within a condition are exchangeable.** Sex, age, batch and ancestry
   are not modelled. A component tracking any of these would present as a
   condition-level component if the covariate is unbalanced across stages.
5. **Stage labels are ordinal, not metric.** F1 → F2 and F3 → F4 are not
   comparable increments; the transition tables treat them as steps, not
   distances.
6. **Steady-state mRNA, one tissue, bulk.** Cell-type composition changes with
   fibrosis. A component appearing with stage may reflect a shift in the
   cellular mixture rather than a change in spatial organisation within any
   cell type. Bulk data cannot separate these; deconvolution or single-cell
   data would be required.

Point 6 is, in our judgement, the most serious threat to a biological reading
of any positive result.

---

## 12. What each output licenses

| Output | Supports | Does **not** support |
|---|---|---|
| `p_condition`, `q_condition` | this frequency shows structure unlikely under positional exchangeability | that the structure is biological rather than compositional |
| `window_rank` | whether sampling alone could produce the frequency | anything about biology |
| CLEAN components | the components the data prefer to retain | a per-component error rate |
| consensus score | how strong, common and phase-aligned a frequency is | a significance claim; it has no null |
| `signature_class` | `confirmed` when prevalence, bootstrap bound and adjusted Rayleigh all hold; `exploratory` otherwise | that an exploratory component is a signature — it is a lead |
| PLV + Rayleigh q | whether the crest falls in the same place across samples | that the component is biological |
| band accuracy | how well a query of that coverage can be classified | that an out-of-domain sample will be rejected |
| `pct_samples_significant` | reproducibility across samples | evidence strength; it does not aggregate |
| `replicated` | agreement in frequency and direction across cohorts | a mechanism |
| crest gene lists | genes positioned on a crest | that those genes are periodically regulated |

---

## References

- Lomb, N. R. (1976). Least-squares frequency analysis of unequally spaced data. *Astrophysics and Space Science* 39, 447–462.
- Scargle, J. D. (1982). Studies in astronomical time series analysis II. *The Astrophysical Journal* 263, 835–853.
- Zechmeister, M. & Kürster, M. (2009). The generalised Lomb–Scargle periodogram. *Astronomy & Astrophysics* 496, 577–584.
- Press, W. H. & Rybicki, G. B. (1989). Fast algorithm for spectral analysis of unevenly sampled data. *The Astrophysical Journal* 338, 277–280.
- Högbom, J. A. (1974). Aperture synthesis with a non-regular distribution of interferometer baselines. *Astronomy & Astrophysics Supplement* 15, 417–426.
- Westfall, P. H. & Young, S. S. (1993). *Resampling-Based Multiple Testing*. Wiley.
- Chen, J. & Chen, Z. (2008). Extended Bayesian information criteria for model selection with large model spaces. *Biometrika* 95, 759–771.
- Fisher, R. A. (1929). Tests of significance in harmonic analysis. *Proceedings of the Royal Society A* 125, 54–59.
- Hurst, L. D., Pál, C. & Lercher, M. J. (2004). The evolutionary dynamics of eukaryotic gene order. *Nature Reviews Genetics* 5, 299–310.
