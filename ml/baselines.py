"""Baselines the autoencoder has to beat before it means anything.

The specification is explicit that TissueSpect-AE should not be adopted unless
it improves on or complements these under leave-one-cohort-out, and that is the
right rule. These are built first, deliberately: a bar set after seeing the
model's score is not a bar.

All of them consume the same flattened, masked spectrum the model consumes, and
all of them are fitted inside a fold — feature selection included. Selecting
features on the full data and then cross-validating is the most common way to
report an inflated number while every other step looks correct.
"""
from __future__ import annotations

import numpy as np

from .evaluate import evaluate_fold
from .prototypes import assign, build_prototypes, prototype_probabilities
from .splits import (assert_no_cohort_leak, class_evaluation_report,
                     leave_one_cohort_out)


def flatten(x: np.ndarray, mask: np.ndarray, channels=(0, 1, 2, 3)) -> np.ndarray:
    """(n, chr, k, ch) -> (n, chr*k*len(channels)), padding already zeroed."""
    sub = x[..., list(channels)]
    return sub.reshape(len(sub), -1)


class NearestCentroidBaseline:
    """Cohort-balanced prototypes on the raw spectrum. No learning at all.

    This is the one to watch. If a linear centroid on the spectra matches the
    autoencoder out of cohort, the encoder is not earning its complexity.
    """

    def __init__(self, metric: str = "cosine", temperature: float = 0.1):
        self.metric, self.temperature = metric, temperature

    def fit(self, xtr, ytr, dstr):
        self.prototypes_, self.classes_, self.info_ = build_prototypes(xtr, ytr, dstr)
        return self

    def predict(self, xte):
        return assign(xte, self.prototypes_, self.classes_, self.metric)[0]

    def predict_proba(self, xte):
        return prototype_probabilities(xte, self.prototypes_, self.metric, self.temperature)


class SklearnBaseline:
    """Thin wrapper so every baseline has the same three methods."""

    def __init__(self, estimator, name: str):
        self.est, self.name = estimator, name

    def fit(self, xtr, ytr, dstr=None):
        self.est.fit(xtr, ytr)
        self.classes_ = list(self.est.classes_)
        return self

    def predict(self, xte):
        return self.est.predict(xte)

    def predict_proba(self, xte):
        return self.est.predict_proba(xte)


def default_baselines(seed: int = 42) -> dict:
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    out = {
        "nearest_centroid": NearestCentroidBaseline(),
        # Elastic net: saga is the solver that supports the l1_ratio mix. Scaling
        # is inside the pipeline so it is fitted on the training fold only.
        "elastic_net": SklearnBaseline(
            make_pipeline(StandardScaler(),
                          # l1_ratio alone selects the elastic-net mix; passing
                          # penalty= as well is deprecated from sklearn 1.8.
                          LogisticRegression(solver="saga", l1_ratio=0.5, C=0.1,
                                             max_iter=3000, random_state=seed)),
            "elastic_net"),
        "random_forest": SklearnBaseline(
            RandomForestClassifier(n_estimators=400, min_samples_leaf=2,
                                   random_state=seed, n_jobs=-1),
            "random_forest"),
    }
    try:  # optional dependency, absent by default
        from xgboost import XGBClassifier  # noqa: F401
        out["xgboost"] = SklearnBaseline(
            XGBClassifier(n_estimators=300, max_depth=3, learning_rate=0.05,
                          subsample=0.8, random_state=seed), "xgboost")
    except ImportError:
        pass
    return out


def run_baselines(x: np.ndarray, mask: np.ndarray, class_ids, dataset_ids,
                  order=None, models: dict | None = None, seed: int = 42) -> dict:
    """Every baseline through every leave-one-cohort-out fold.

    Returns per-fold rows and a per-model summary. The summary averages folds
    unweighted: a fold is a cohort, and the question is how the model behaves on
    a new cohort, not on a new sample.
    """
    models = models or default_baselines(seed)
    flat = flatten(x, mask)
    class_ids, dataset_ids = np.asarray(class_ids), np.asarray(dataset_ids)
    folds = leave_one_cohort_out(dataset_ids, class_ids)

    rows = []
    for fold in folds:
        assert_no_cohort_leak(fold, dataset_ids)
        xtr, xte = flat[fold.train_idx], flat[fold.test_idx]
        ytr, yte = class_ids[fold.train_idx], class_ids[fold.test_idx]
        dtr = dataset_ids[fold.train_idx]

        for name, model in models.items():
            try:
                model.fit(xtr, ytr, dtr) if isinstance(model, NearestCentroidBaseline) \
                    else model.fit(xtr, ytr)
                pred = model.predict(xte)
                proba = model.predict_proba(xte)
                classes = list(model.classes_)
            except Exception as exc:  # a baseline failing must not stop the rest
                rows.append({"model": name, "held_out": fold.held_out,
                             "error": str(exc)[:200]})
                continue
            metrics = evaluate_fold(yte, pred, proba, classes, ytr, order=order)
            rows.append({"model": name, "held_out": fold.held_out, **metrics})

    summary = {}
    for name in models:
        vals = [r for r in rows if r["model"] == name and "error" not in r]
        if not vals:
            continue
        summary[name] = {
            k: float(np.nanmean([v[k] for v in vals]))
            for k in vals[0] if isinstance(vals[0][k], (int, float))
        }
        summary[name]["n_folds"] = len(vals)
    return {"folds": rows, "summary": summary,
            "class_report": class_evaluation_report(folds, dataset_ids, class_ids)}


def format_summary(result: dict) -> str:
    """A table you can paste into a decision, not a dict you have to squint at.

    The per-class report is printed with it, always. Every average below is over
    the classes that were evaluated out of cohort, and which those are is not
    obvious from the number.
    """
    keys = ["balanced_accuracy", "accuracy", "majority_baseline",
            "lift_over_baseline", "macro_f1", "stage_mae"]
    header = f"{'model':<18}" + "".join(f"{k[:14]:>16}" for k in keys)
    lines = [header, "-" * len(header)]
    for name, s in sorted(result["summary"].items(),
                          key=lambda kv: -kv[1].get("balanced_accuracy", 0)):
        lines.append(f"{name:<18}" + "".join(
            f"{s.get(k, float('nan')):>16.3f}" for k in keys))
    if "class_report" in result:
        from .splits import format_class_report
        lines += ["", format_class_report(result["class_report"])]
    return "\n".join(lines)
