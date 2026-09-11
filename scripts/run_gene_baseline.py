#!/usr/bin/env python3
"""run_gene_baseline.py -- the trivial baseline the spectrum has to beat.

Leave-one-cohort-out over the ingested datasets, on GENES, not spectra: the
top-K most variable genes of the training fold, median-imputed and
standardised on the training fold, multinomial logistic regression. Feature
selection and imputation live inside the sklearn pipeline, so in the
within-cohort cross-validation they are refit per fold and never see the
validation samples. Same folds, same target, same
metric as `./tsf reference` (accuracy on the held-out cohort over the classes
it shares with training, against the training majority class).

If the spectral fingerprint does not beat this, it is transforming information,
not adding it. The defensible claim for the spectrum is a different one: that
its accuracy falls LESS when the cohort changes. This script reports, per fold,
the within-cohort cross-validated accuracy next to the out-of-cohort one, so
that drop can be compared with the spectral drop.

Usage:
  python3 scripts/run_gene_baseline.py --interim-dir $TSF_INTERIM_DIR \
      --datasets GSE135251,GSE130970,R3_LIVER [--target class_id|condition] \
      [--top-genes 1000] [--out results/gene_baseline.tsv]

Reads whatever expression.tsv is on disk: if apply_reference_profile.R ran, the
baseline is on deviations too, so the comparison stays like for like.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


class TopVariance(BaseEstimator, TransformerMixin):
    """Keep the k most variable columns, chosen on the data given to fit().

    Inside the pipeline, so that in cross-validation the selection sees the
    training fold only. Selecting on the whole cohort first leaks the held-out
    samples' variance into the choice of features.
    """

    def __init__(self, k=1000):
        self.k = k

    def fit(self, X, y=None):
        var = np.nanvar(X, axis=0)
        var = np.where(np.isfinite(var), var, -np.inf)
        self.idx_ = np.argsort(-var)[: min(self.k, X.shape[1])]
        return self

    def transform(self, X):
        return X[:, self.idx_]


def make_model(top_genes, seed):
    # Imputation is explicit and per fold: the training-fold MEDIAN of each
    # gene. Zero is not a neutral fill on the deviation scale -- it means
    # "equal to the GTEx median", a measurement the sample never made.
    return make_pipeline(
        TopVariance(k=top_genes),
        SimpleImputer(strategy="median"),
        StandardScaler(),
        LogisticRegression(max_iter=3000, C=1.0, random_state=seed),
    )


def load(interim: Path, ds: str, target: str):
    d = interim / ds
    expr = pd.read_csv(d / "expression.tsv", sep="\t", index_col=0)
    expr.index = expr.index.astype(str).str.replace(r"\..*$", "", regex=True)
    expr = expr[~expr.index.duplicated()]
    samples = pd.read_csv(d / "samples.tsv", sep="\t")
    samples = samples[samples["sample_id"].isin(expr.columns)]
    y = samples[target].astype(str).values
    x = expr[samples["sample_id"]].T  # samples x genes
    return x, y, samples["sample_id"].values


def fit_predict(xtr, ytr, xte, top_genes, seed):
    model = make_model(top_genes, seed)
    model.fit(xtr.values, ytr)
    return model.predict(xte.values), model


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--interim-dir", required=True)
    ap.add_argument("--datasets", required=True)
    ap.add_argument("--target", default="class_id")
    ap.add_argument("--top-genes", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    interim = Path(a.interim_dir)
    ids = [s.strip() for s in a.datasets.split(",") if s.strip()]

    data = {}
    for ds in ids:
        try:
            data[ds] = load(interim, ds, a.target)
        except FileNotFoundError as e:
            print(f"[skip] {ds}: {e}", file=sys.stderr)
    if len(data) < 2:
        sys.exit("need at least two ingested datasets")

    rows = []
    for held in data:
        xte, yte, _ = data[held]
        train = [d for d in data if d != held]
        genes = set.intersection(*[set(data[d][0].columns) for d in train]) & set(xte.columns)
        genes = sorted(genes)
        xtr = pd.concat([data[d][0][genes] for d in train])
        ytr = np.concatenate([data[d][1] for d in train])
        shared = sorted(set(ytr) & set(yte))
        if len(shared) < 2:
            rows.append(dict(held_out=held, n_test=len(yte), n_shared_classes=len(shared),
                             note="fewer than 2 shared classes; skipped"))
            continue
        m_tr = np.isin(ytr, shared); m_te = np.isin(yte, shared)
        pred, _ = fit_predict(xtr[m_tr], ytr[m_tr], xte[genes][m_te], a.top_genes, a.seed)
        acc = float((pred == yte[m_te]).mean())
        maj = pd.Series(ytr[m_tr]).mode()[0]
        base = float((yte[m_te] == maj).mean())

        # within-cohort CV on the held-out cohort alone, same feature count,
        # to measure how much accuracy the cohort change costs
        within = np.nan
        counts = pd.Series(yte[m_te]).value_counts()
        if (counts >= 3).all() and len(counts) >= 2:
            k = int(min(5, counts.min()))
            skf = StratifiedKFold(n_splits=k, shuffle=True, random_state=a.seed)
            xw = xte[genes][m_te]
            # the whole pipeline -- selection, imputation, scaling -- is refit
            # inside each fold by cross_val_predict
            pw = cross_val_predict(make_model(a.top_genes, a.seed), xw.values, yte[m_te], cv=skf)
            within = float((pw == yte[m_te]).mean())

        rows.append(dict(held_out=held, n_train=int(m_tr.sum()), n_test=int(m_te.sum()),
                         n_shared_classes=len(shared), n_genes=len(genes),
                         top_genes=a.top_genes, out_of_cohort_acc=round(acc, 4),
                         majority_baseline=round(base, 4),
                         within_cohort_cv_acc=(round(within, 4) if np.isfinite(within) else np.nan),
                         cohort_drop=(round(within - acc, 4) if np.isfinite(within) else np.nan),
                         note=""))
        print(f"{held:>14}  out-of-cohort {acc:.3f}  majority {base:.3f}  "
              f"within-cohort CV {within:.3f}  classes {shared}")

    tab = pd.DataFrame(rows)
    out = a.out or str(interim.parent / "results" / "gene_baseline.tsv")
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    tab.to_csv(out, sep="\t", index=False)
    ok = tab["out_of_cohort_acc"].dropna()
    if len(ok):
        print(f"\nmean out-of-cohort accuracy {ok.mean():.3f} over {len(ok)} fold(s) -> {out}")
        print("Compare with results/reference/validation (spectral LOCO): the number to "
              "report is the cohort_drop of each, not the accuracy alone.")


if __name__ == "__main__":
    main()
