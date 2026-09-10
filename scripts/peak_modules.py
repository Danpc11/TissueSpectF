#!/usr/bin/env python3
"""peak_modules.py -- the gene set behind a spectral peak, and whether it is
co-expressed.

    python3 scripts/peak_modules.py \
        --expr GSE162694_chr20_expr.tsv.gz \
        --grid GSE162694_chr20_grid.tsv \
        --labels GSE162694_labels.tsv \
        --chr 20 --period-mb 2.92 --phase-mb 2.57 \
        [--self-test]

THE QUESTION
------------
A spectral peak is not a location: a period of 2.92 Mb means power at the
frequency of a 2.92 Mb repeat, spread over the chromosome. What it does have is
a PHASE, which says where the crests fall. The genes near those crests are the
ones contributing to the component.

So: does that set co-express? If the genes of a peak are correlated with each
other more than chance, the peak is capturing a real module that is organised
positionally -- which is the thing a per-gene differential test cannot see,
because it treats every gene as independent.

THE NULL, AND WHY IT IS THE HARD PART
-------------------------------------
Six earlier attempts in this analysis failed the same way: the null did not
destroy the structure the observed statistic had by construction, so it
compared a value with itself or with something geometrically guaranteed to be
smaller.

The specific trap here: a crest set is a union of narrow windows, so its genes
sit in tight groups. Neighbouring genes on a chromosome are co-expressed for
reasons that have nothing to do with any period -- shared enhancers, operon-like
domains, co-duplication. A null drawing genes uniformly from the chromosome
would produce dispersed sets with lower correlation, and the crest set would
look significant purely because it is clumped.

`shift` is therefore the default: the SAME crest geometry, slid along the
chromosome by a random offset. Spacing, clump sizes and set size are all
preserved exactly, and the only thing that changes is where the crests land.
That isolates the phase, which is the claim.

`--self-test` checks the null can tell the two cases apart before any real
data is read: a planted module at the crests must come out significant, and a
clumped-but-uncorrelated set must not.
"""

import argparse
import sys

import numpy as np
import pandas as pd


def crest_positions(period, phase0, lo, hi):
    """Crest centres of cos(2*pi*(x - phase0)/period) over [lo, hi]."""
    first = phase0 + period * np.ceil((lo - phase0) / period)
    return np.arange(first, hi + period, period)


def crest_mask(pos, crests, half_width):
    """Which genes fall within half_width of any crest."""
    if not len(crests):
        return np.zeros(len(pos), bool)
    d = np.abs(pos[:, None] - crests[None, :]).min(axis=1)
    return d <= half_width


def residualise_by_condition(X, cond):
    """Quitar la media de cada gen DENTRO de cada condición.

    Sin esto, una diferencia de medias entre condiciones crea correlación
    aparente: si cien genes suben en F4, entre todas las muestras juntas
    covarían --las de F4 altas, las demás bajas-- aunque dentro de F4 no
    covaríen en absoluto. Eso es exactamente lo que un pico de expresión
    diferencial produciría, y llamarlo módulo de coexpresión sería confundir
    dos cosas distintas.

    Residualizar deja sólo la covarianza INTRA-condición, que es la que define
    un módulo. El costo: un módulo que existiera únicamente como diferencia
    entre condiciones se vuelve invisible, y eso es correcto -- ese caso es
    expresión diferencial y ya lo mide `differential`.
    """
    R = X.astype(float).copy()
    for c in np.unique(cond):
        sel = cond == c
        if sel.sum() < 2:
            R[:, sel] = np.nan
            continue
        R[:, sel] -= np.nanmean(R[:, sel], axis=1, keepdims=True)
    return R


def mean_abs_corr(X):
    """Mean |correlation| over gene pairs. Absolute value because a module can
    contain anti-correlated members -- a repressor and its target belong to the
    same module -- and a signed mean would cancel them out."""
    if X.shape[0] < 2:
        return np.nan
    C = np.corrcoef(X)
    iu = np.triu_indices(C.shape[0], 1)
    v = C[iu]
    v = v[np.isfinite(v)]
    return float(np.mean(np.abs(v))) if len(v) else np.nan


def peak_module(pos, X, period, phase0, half_width, n_null=2000, seed=42):
    """Observed co-expression of the crest set, against a shifted-crest null."""
    lo, hi = float(pos.min()), float(pos.max())
    obs_mask = crest_mask(pos, crest_positions(period, phase0, lo, hi), half_width)
    n_in = int(obs_mask.sum())
    if n_in < 3:
        return {"n": n_in, "obs": np.nan, "null_median": np.nan, "p": np.nan}
    obs = mean_abs_corr(X[obs_mask])

    rng = np.random.default_rng(seed)
    null = []
    for _ in range(n_null):
        # Slide the crest comb by a random offset. Same period, same width,
        # same comb -- only the phase differs.
        m = crest_mask(pos, crest_positions(period, phase0 + rng.uniform(0, period),
                                            lo, hi), half_width)
        if m.sum() >= 3:
            null.append(mean_abs_corr(X[m]))
    null = np.array([v for v in null if np.isfinite(v)])
    if not len(null):
        return {"n": n_in, "obs": obs, "null_median": np.nan, "p": np.nan}
    p = (1 + int((null >= obs).sum())) / (len(null) + 1)
    return {"n": n_in, "obs": obs, "null_median": float(np.median(null)),
            "null_q95": float(np.quantile(null, 0.95)), "p": p,
            "n_null": len(null)}


def self_test():
    """The null has to separate a planted module from mere clumping.

    Without this the whole analysis is unfalsifiable: a null that cannot fail
    reports a discovery either way.
    """
    rng = np.random.default_rng(0)
    period, phase0, hw = 2.92, 2.57, 0.5
    pos = np.sort(rng.uniform(0, 64, 350))
    n_s = 143
    crests = crest_positions(period, phase0, pos.min(), pos.max())
    at_crest = crest_mask(pos, crests, hw)

    ok = True

    # CASE 1: a module planted AT the crests must be found.
    X = rng.normal(size=(len(pos), n_s))
    shared = rng.normal(size=n_s)
    X[at_crest] += 1.2 * shared
    r = peak_module(pos, X, period, phase0, hw, n_null=400, seed=1)
    print(f"  planted at crests   n={r['n']:3d} obs={r['obs']:.3f} "
          f"null={r['null_median']:.3f} p={r['p']:.4f}", end="")
    if r["p"] <= 0.05:
        print("   PASS")
    else:
        print("   FAIL: cannot detect a real module")
        ok = False

    # CASE 2: clumped but uncorrelated must NOT be found. The genes sit in
    # tight groups, so a uniform-draw null would call this significant.
    X = rng.normal(size=(len(pos), n_s))
    r = peak_module(pos, X, period, phase0, hw, n_null=400, seed=2)
    print(f"  clumped, no module  n={r['n']:3d} obs={r['obs']:.3f} "
          f"null={r['null_median']:.3f} p={r['p']:.4f}", end="")
    if r["p"] > 0.05:
        print("   PASS")
    else:
        print("   FAIL: calls clumping a module")
        ok = False

    # CASE 3: a module planted OFF phase must not be credited to this phase.
    X = rng.normal(size=(len(pos), n_s))
    off = crest_mask(pos, crest_positions(period, phase0 + period / 2,
                                          pos.min(), pos.max()), hw)
    X[off] += 1.2 * rng.normal(size=n_s)
    r = peak_module(pos, X, period, phase0, hw, n_null=400, seed=3)
    print(f"  module in antiphase n={r['n']:3d} obs={r['obs']:.3f} "
          f"null={r['null_median']:.3f} p={r['p']:.4f}", end="")
    if r["p"] > 0.05:
        print("   PASS")
    else:
        print("   FAIL: wrong phase credited")
        ok = False

    print("\nself-test", "PASSED" if ok else "FAILED")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--expr")
    ap.add_argument("--grid")
    ap.add_argument("--labels")
    # NOMBRES EXPLÍCITOS, no la primera y la última columna.
    #
    # samples.tsv tiene 17 columnas y termina en `keep`, no en `condition`:
    # sample_id, dataset_id, tissue, vocabulary, state, class_id, condition,
    # fibrosis_stage, ..., keep, filtered_out. Tomar `columns[-1]` residualiza
    # por TRUE/FALSE y el análisis por condición sale mal sin avisar.
    ap.add_argument("--sample-column", default="sample_id")
    ap.add_argument("--condition-column", default="condition")
    ap.add_argument("--chr", default="20")
    ap.add_argument("--period-mb", type=float, default=2.92)
    ap.add_argument("--phase-mb", type=float, default=2.57)
    ap.add_argument("--half-width-mb", type=float, default=0.5)
    ap.add_argument("--n-null", type=int, default=2000)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--out")
    a = ap.parse_args()

    if a.self_test:
        sys.exit(0 if self_test() else 1)
    if not (a.expr and a.grid):
        ap.error("--expr and --grid are required unless --self-test")

    print("Validating the null before touching the data:")
    if not self_test():
        sys.exit("null does not separate the cases; nothing below is meaningful")
    print()

    grid = pd.read_csv(a.grid, sep="\t", dtype={"chr": str})
    grid = grid[grid.chr == str(a.chr)]
    X = pd.read_csv(a.expr, sep="\t", index_col=0)
    g = grid[grid.gene_id.isin(X.index)].sort_values("start")
    X = X.loc[g.gene_id].to_numpy(dtype=float)
    pos = (g.start.to_numpy(dtype=float)) / 1e6
    print(f"chr{a.chr}: {len(g)} expressed genes over "
          f"{pos.min():.1f}-{pos.max():.1f} Mb, {X.shape[1]} samples")

    rows = []
    # RESIDUALIZADO por condición cuando las etiquetas están. Sin residualizar,
    # el resultado sobre todas las muestras juntas mezcla coexpresión con
    # diferencias de medias entre condiciones.
    Xr = X
    if a.labels:
        lab0 = pd.read_csv(a.labels, sep="\t")
        for c in (a.sample_column, a.condition_column):
            if c not in lab0.columns:
                sys.exit(f"--labels no tiene la columna '{c}'. Presentes: "
                         f"{', '.join(lab0.columns)}. Usá --sample-column y "
                         f"--condition-column; adivinar la primera y la última "
                         f"residualiza por la columna equivocada sin avisar.")
        cols0 = pd.read_csv(a.expr, sep="\t", index_col=0, nrows=0).columns
        m0 = dict(zip(lab0[a.sample_column], lab0[a.condition_column]))
        cnd = np.array([m0.get(c, "NA") for c in cols0])
        if len(set(cnd) - {"NA"}) >= 2:
            Xr = residualise_by_condition(X, cnd)
            keep = np.isfinite(Xr).all(axis=0)
            Xr = Xr[:, keep]
            print(f"residualizado por condición: {len(set(cnd) - {'NA'})} "
                  f"condiciones, {Xr.shape[1]} muestras usables")
        else:
            print("una sola condición: no se residualiza")

    r = peak_module(pos, Xr, a.period_mb, a.phase_mb, a.half_width_mb,
                    n_null=a.n_null)
    r.update(period_mb=a.period_mb, phase_mb=a.phase_mb,
             scope="all_samples_residualised" if Xr is not X else "all_samples")
    rows.append(r)
    print(f"\nperiod {a.period_mb} Mb, phase {a.phase_mb} Mb:")
    print(f"  crest genes      {r['n']}")
    print(f"  mean |r|         {r['obs']:.4f}")
    print(f"  null median      {r['null_median']:.4f}   q95 {r['null_q95']:.4f}")
    print(f"  p                {r['p']:.4f}")

    # Same peak within each condition: a module present only in disease is a
    # different claim from one present throughout.
    if a.labels:
        lab = pd.read_csv(a.labels, sep="\t")
        cols = pd.read_csv(a.expr, sep="\t", index_col=0, nrows=0).columns
        m = dict(zip(lab[a.sample_column], lab[a.condition_column]))
        cond = np.array([m.get(c, "NA") for c in cols])
        print("\nby condition:")
        for c in sorted(set(cond) - {"NA"}):
            sel = cond == c
            if sel.sum() < 8:
                print(f"  {c:22} n={sel.sum():3d} samples -- too few, skipped")
                continue
            rc = peak_module(pos, X[:, sel], a.period_mb, a.phase_mb,
                             a.half_width_mb, n_null=max(400, a.n_null // 4))
            rc.update(period_mb=a.period_mb, phase_mb=a.phase_mb, scope=c)
            rows.append(rc)
            print(f"  {c:22} n={sel.sum():3d}  |r|={rc['obs']:.4f}  "
                  f"null={rc['null_median']:.4f}  p={rc['p']:.4f}")

    if a.out:
        pd.DataFrame(rows).to_csv(a.out, sep="\t", index=False)
        print(f"\nwrote {a.out}")

    n_tests = len(rows)
    if n_tests > 1:
        print(f"\nMULTIPLICIDAD: {n_tests} pruebas en esta corrida (global mas "
              f"una por condicion). Y a lo largo del proyecto se prueban varios "
              f"picos, cromosomas y anchos de cresta: los p de arriba son SIN "
              f"corregir, y hay que corregirlos por el total de combinaciones "
              f"evaluadas, no por las de una sola invocacion.")

    print("\nWhat a small p does and does not mean: the crest set is more "
          "co-expressed than the same comb at another phase. It does not "
          "establish that the period is the cause -- only that this phase is "
          "special among phases of this period.")


if __name__ == "__main__":
    main()
