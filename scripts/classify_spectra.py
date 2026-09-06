#!/usr/bin/env python3
"""classify_spectra.py -- can the spectrum say which condition a sample is?

    python scripts/classify_spectra.py \
        --results-dir results_ncbi \
        --datasets GSE130970,GSE135251,GSE162694,GSE276114 \
        [--n-bins 200] [--models all]

WHAT IT COMPARES
----------------
The same leave-one-cohort-out folds, the same tensor, four models:

  centroid    nearest shrunken centroid. What the R pipeline does today, as the
              reference point rather than as a contender.
  catch22     22 canonical time-series features per chromosome, then a forest.
              Compact enough for this sample size and every feature has a name,
              which matters when the question is *which periods*.
  minirocket  fixed convolutional kernels + ridge. Among the top performers in
              recent time-series classification surveys. Read as the reachable
              accuracy CEILING, not as the final model: it does not say which
              periods carry the signal.
  expression  the CONTROL: the gene expression itself, same folds, no spectral
              transform. Without it none of the three numbers above can be read
              as good or bad.

WHY THE TENSOR IS (sample, chromosome, period)
----------------------------------------------
A spectrum indexed by (chromosome, k) is not one series: k is cycles per
chromosome, so k = 64 is a period of 32 genes on chr1 and 12 on chr21, and the
columns are not comparable. Re-indexed by period on a common log grid, each
chromosome becomes a channel of the SAME series, and the whole multivariate
time-series toolbox applies without adaptation.

WHY MORE BINS HERE THAN IN THE DIFFERENTIAL TEST
------------------------------------------------
The differential test collapses to ~40 bands because multiplicity is its
binding constraint and m is the only lever. Classification corrects for no
multiplicity -- it feeds a model rather than testing m hypotheses -- so
resolution is information, and 40 points is too short for dilated kernels to
have any range of scales to work with. Hence --n-bins 200 by default.
"""

import argparse
import glob
import importlib.util
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ml.splits import leave_one_cohort_out, assert_no_cohort_leak  # noqa: E402


def period_grid(n_bins, lo=10.0, hi=500.0):
    """Log-spaced, because period resolution is multiplicative: the gap between
    20 and 21 genes is not the gap between 400 and 401."""
    return np.exp(np.linspace(np.log(lo), np.log(hi), n_bins + 1))


def load_tensor(results_dir, datasets, n_bins):
    """(samples, chromosomes, bins), plus labels and cohort ids."""
    frames = []
    for ds in datasets:
        pat = os.path.join(results_dir, ds, "spectra", "spectra_samples_*.tsv")
        for f in sorted(glob.glob(pat)):
            cond = os.path.basename(f)[len("spectra_samples_"):-len(".tsv")]
            d = pd.read_csv(f, sep="\t")
            need = {"chr", "k", "sample", "power_normalised", "period"}
            if not need.issubset(d.columns):
                print(f"  skip {f}: missing {sorted(need - set(d.columns))}")
                continue
            d["condition"] = cond
            d["dataset"] = ds
            frames.append(d)
    if not frames:
        raise SystemExit("No per-sample spectra found. Run the spectra stage.")
    sp = pd.concat(frames, ignore_index=True)

    br = period_grid(n_bins)
    sp = sp[(sp["period"] >= br[0]) & (sp["period"] <= br[-1])]
    sp["bin"] = np.digitize(sp["period"], br) - 1
    sp["bin"] = sp["bin"].clip(0, n_bins - 1)

    agg = (sp.groupby(["sample", "chr", "bin"])["power_normalised"]
             .mean().reset_index())

    samples = sorted(sp["sample"].unique())
    chroms = sorted(sp["chr"].astype(str).unique())
    s_ix = {s: i for i, s in enumerate(samples)}
    c_ix = {c: i for i, c in enumerate(chroms)}

    X = np.full((len(samples), len(chroms), n_bins), np.nan)
    X[agg["sample"].map(s_ix).values,
      agg["chr"].astype(str).map(c_ix).values,
      agg["bin"].values] = agg["power_normalised"].values

    # A bin no frequency of that chromosome falls into is genuinely unmeasured.
    # Filling with zero would assert the spectrum has no power there; the column
    # mean keeps the sample comparable without inventing a measurement, and the
    # models below cannot take NaN.
    col_mean = np.nanmean(X, axis=0)
    inds = np.where(np.isnan(X))
    X[inds] = np.take(col_mean, inds[1] * n_bins + inds[2])
    X = np.nan_to_num(X, nan=0.0)

    meta = sp[["sample", "condition", "dataset"]].drop_duplicates("sample")
    meta = meta.set_index("sample").loc[samples].reset_index()
    return X, meta["condition"].values, meta["dataset"].values, chroms


def run_folds(X, y, cohorts, fit_predict, name):
    """One accuracy per fold, plus the training-fold majority baseline.

    The baseline is the majority of the TRAINING fold, never of the held-out
    truth: taking the majority of what is being predicted would let the
    baseline see the test labels, which is the one thing the held-out design
    exists to prevent.
    """
    folds = leave_one_cohort_out(cohorts, y)
    if not folds:
        print(f"  {name}: no usable fold")
        return None
    acc, base, n = [], [], 0
    for fold in folds:
        assert_no_cohort_leak(fold, cohorts)
        tr, te = fold.train_idx, fold.test_idx
        if len(np.unique(y[tr])) < 2 or not len(te):
            continue
        try:
            pred = fit_predict(X[tr], y[tr], X[te])
        except Exception as exc:                      # noqa: BLE001
            print(f"  {name}: fold {fold.held_out} failed: {exc}")
            continue
        acc.append(float(np.mean(pred == y[te])))
        vals, cnt = np.unique(y[tr], return_counts=True)
        base.append(float(np.mean(y[te] == vals[np.argmax(cnt)])))
        n += len(te)
    if not acc:
        print(f"  {name}: every fold failed")
        return None
    print(f"  {name:<12} {np.mean(acc):6.1%}  vs {np.mean(base):5.1%} baseline"
          f"   (n={n}, {len(acc)} folds, sd {np.std(acc):.3f})")
    return {"model": name, "accuracy": float(np.mean(acc)),
            "baseline": float(np.mean(base)), "n": n, "folds": len(acc),
            "sd": float(np.std(acc))}


# --- models -------------------------------------------------------------------

def m_centroid(Xtr, ytr, Xte):
    from sklearn.neighbors import NearestCentroid
    from sklearn.preprocessing import StandardScaler
    a, b = Xtr.reshape(len(Xtr), -1), Xte.reshape(len(Xte), -1)
    sc = StandardScaler().fit(a)
    return NearestCentroid().fit(sc.transform(a), ytr).predict(sc.transform(b))


def m_catch22(Xtr, ytr, Xte):
    import pycatch22
    from sklearn.ensemble import RandomForestClassifier

    def feats(X):
        out = np.empty((X.shape[0], X.shape[1] * 22))
        for i in range(X.shape[0]):
            row = []
            for c in range(X.shape[1]):
                row.extend(pycatch22.catch22_all(
                    list(map(float, X[i, c])))["values"])
            out[i] = row
        return np.nan_to_num(out, nan=0.0, posinf=0.0, neginf=0.0)

    return RandomForestClassifier(
        n_estimators=300, random_state=42, class_weight="balanced"
    ).fit(feats(Xtr), ytr).predict(feats(Xte))


def m_minirocket(Xtr, ytr, Xte):
    from sktime.transformations.panel.rocket import MiniRocketMultivariate
    from sklearn.linear_model import RidgeClassifierCV
    mr = MiniRocketMultivariate(random_state=42).fit(Xtr)
    # Ridge handles ~10k features against a few hundred samples; that is what
    # the regularisation is for, and RidgeClassifierCV picks alpha itself.
    return RidgeClassifierCV(alphas=np.logspace(-3, 3, 10)).fit(
        mr.transform(Xtr), ytr).predict(mr.transform(Xte))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--datasets", required=True)
    ap.add_argument("--n-bins", type=int, default=200)
    ap.add_argument("--models", default="all")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    X, y, cohorts, chroms = load_tensor(
        a.results_dir, a.datasets.split(","), a.n_bins)
    print(f"tensor: {X.shape[0]} samples x {X.shape[1]} chromosomes "
          f"x {X.shape[2]} period bins")
    print(f"classes: {dict(zip(*np.unique(y, return_counts=True)))}")
    print(f"cohorts: {sorted(set(cohorts))}\n")

    want = a.models.split(",") if a.models != "all" else \
        ["centroid", "catch22", "minirocket"]
    rows = []
    # A missing package is not a fold failure. Reported once, by name, with the
    # install line -- the previous version printed the same ImportError once per
    # fold and buried the reason.
    needs = {"catch22": ["pycatch22"], "minirocket": ["sktime"]}
    for name, fn in [("centroid", m_centroid), ("catch22", m_catch22),
                     ("minirocket", m_minirocket)]:
        if name not in want:
            continue
        missing = [m for m in needs.get(name, [])
                   if importlib.util.find_spec(m) is None]
        if missing:
            print(f"  {name:<12} skipped: pip install {' '.join(missing)}")
            continue
        r = run_folds(X, y, cohorts, fn, name)
        if r:
            rows.append(r)

    if rows and a.out:
        pd.DataFrame(rows).to_csv(a.out, sep="\t", index=False)
        print(f"\nwrote {a.out}")

    print("\nRun the expression control before reading any of these: an "
          "accuracy has no meaning without it.\n"
          "  ./tsf reference ... --features expression_baseline")


if __name__ == "__main__":
    main()
