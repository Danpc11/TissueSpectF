"""Metrics for cohort-held-out evaluation.

Everything here is computed on a held-out cohort. Nothing in this file should
ever be handed a random split: the numbers would be higher and would mean
nothing, and that is the specific failure this project is built to avoid.

Two families are reported side by side because they answer different questions.
Classification metrics ask whether the right class was chosen. Representation
metrics ask what the space is organised by — a model can classify passably while
its latent space clusters by cohort, and that model will not survive contact
with a sixth series.
"""
from __future__ import annotations

import numpy as np


def _confusion(y_true, y_pred, classes) -> np.ndarray:
    idx = {c: i for i, c in enumerate(classes)}
    m = np.zeros((len(classes), len(classes)), dtype=int)
    for t, p in zip(y_true, y_pred):
        if t in idx and p in idx:
            m[idx[t], idx[p]] += 1
    return m


def balanced_accuracy(y_true, y_pred, classes=None) -> float:
    """Mean per-class recall. Plain accuracy would reward predicting F4."""
    classes = classes or sorted(set(map(str, y_true)))
    cm = _confusion(y_true, y_pred, classes)
    present = cm.sum(1) > 0
    if not present.any():
        return float("nan")
    return float(np.mean(cm.diagonal()[present] / cm.sum(1)[present]))


def macro_f1(y_true, y_pred, classes=None) -> float:
    classes = classes or sorted(set(map(str, y_true)))
    cm = _confusion(y_true, y_pred, classes)
    f1s = []
    for i in range(len(classes)):
        tp, fp, fn = cm[i, i], cm[:, i].sum() - cm[i, i], cm[i, :].sum() - cm[i, i]
        if tp + fn == 0:
            continue
        prec = tp / (tp + fp) if tp + fp else 0.0
        rec = tp / (tp + fn)
        f1s.append(0.0 if prec + rec == 0 else 2 * prec * rec / (prec + rec))
    return float(np.mean(f1s)) if f1s else float("nan")


def majority_baseline(y_train, y_test) -> float:
    """Always predict the training fold's commonest class.

    The baseline has to be computed on the TRAINING fold. Taking the majority of
    the held-out truth lets the baseline see the test labels, which is exactly
    what the held-out design exists to prevent.
    """
    vals, counts = np.unique(np.asarray(y_train), return_counts=True)
    guess = vals[np.argmax(counts)]
    return float(np.mean(np.asarray(y_test) == guess))


def top_k_accuracy(y_true, proba: np.ndarray, classes, k: int = 2) -> float:
    order = np.argsort(-proba, axis=1)[:, :k]
    hits = [y in np.asarray(classes)[row] for y, row in zip(y_true, order)]
    return float(np.mean(hits)) if hits else float("nan")


# --- ordinal ----------------------------------------------------------------
def _ordinal_position(labels, order) -> np.ndarray:
    pos = {c: i for i, c in enumerate(order)}
    return np.array([pos.get(str(l), np.nan) for l in labels], dtype=float)


def stage_mae(y_true, y_pred, order) -> float:
    """Mean absolute error in ordinal steps, over samples that have a rank.

    A step is a step, not a distance: F1 to F2 and F3 to F4 are not comparable
    increments, and this number should be read as "how many stages off", never
    as a clinical magnitude.
    """
    t, p = _ordinal_position(y_true, order), _ordinal_position(y_pred, order)
    ok = np.isfinite(t) & np.isfinite(p)
    return float(np.mean(np.abs(t[ok] - p[ok]))) if ok.any() else float("nan")


def quadratic_weighted_kappa(y_true, y_pred, order) -> float:
    t, p = _ordinal_position(y_true, order), _ordinal_position(y_pred, order)
    ok = np.isfinite(t) & np.isfinite(p)
    t, p = t[ok].astype(int), p[ok].astype(int)
    n = len(order)
    if len(t) == 0:
        return float("nan")
    o = np.zeros((n, n))
    for a, b in zip(t, p):
        o[a, b] += 1
    w = np.array([[((i - j) ** 2) / ((n - 1) ** 2 or 1) for j in range(n)] for i in range(n)])
    ht, hp = np.bincount(t, minlength=n), np.bincount(p, minlength=n)
    e = np.outer(ht, hp) / max(len(t), 1)
    denom = (w * e).sum()
    return float(1 - (w * o).sum() / denom) if denom else float("nan")


# --- calibration ------------------------------------------------------------
def expected_calibration_error(y_true, proba: np.ndarray, classes, n_bins: int = 10) -> float:
    """How far the stated confidence is from the observed hit rate."""
    classes = np.asarray(classes)
    conf = proba.max(axis=1)
    pred = classes[proba.argmax(axis=1)]
    correct = (pred == np.asarray(y_true)).astype(float)
    edges = np.linspace(0, 1, n_bins + 1)
    ece = 0.0
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (conf > lo) & (conf <= hi)
        if m.any():
            ece += m.mean() * abs(correct[m].mean() - conf[m].mean())
    return float(ece)


def brier_score(y_true, proba: np.ndarray, classes) -> float:
    classes = np.asarray(classes)
    onehot = (classes[None, :] == np.asarray(y_true)[:, None]).astype(float)
    return float(np.mean(np.sum((proba - onehot) ** 2, axis=1)))


# --- representation ---------------------------------------------------------
def silhouette(z: np.ndarray, labels) -> float:
    """Silhouette by whatever labels are given. Run it twice.

    By condition, high is good. By cohort, high is bad — it means the space is
    organised by the batch. Reporting only the first would hide the failure the
    adversarial head exists to prevent.
    """
    labels = np.asarray(labels)
    if len(set(labels.tolist())) < 2 or len(z) < 3:
        return float("nan")
    from scipy.spatial.distance import cdist

    d = cdist(z, z)
    scores = []
    for i in range(len(z)):
        same = (labels == labels[i]) & (np.arange(len(z)) != i)
        if not same.any():
            continue
        a = d[i, same].mean()
        b = min(d[i, labels == other].mean()
                for other in set(labels.tolist()) if other != labels[i])
        scores.append((b - a) / max(a, b) if max(a, b) > 0 else 0.0)
    return float(np.mean(scores)) if scores else float("nan")


def cohort_predictability(z: np.ndarray, dataset_ids, n_splits: int = 3) -> dict:
    """Can a simple classifier read the cohort off the representation?

    The number to watch. If cohort is as predictable from the latent space as
    the condition is, the model has learned the batch, whatever its accuracy
    says. Compared against the majority-cohort rate, because with an imbalanced
    set of cohorts a high raw accuracy can be trivial.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold, cross_val_predict

    y = np.asarray(dataset_ids)
    if len(set(y.tolist())) < 2:
        return {"accuracy": float("nan"), "majority": float("nan"), "lift": float("nan")}
    counts = np.bincount(np.unique(y, return_inverse=True)[1])
    n_splits = int(min(n_splits, counts.min()))
    if n_splits < 2:
        return {"accuracy": float("nan"), "majority": float(counts.max() / len(y)),
                "lift": float("nan")}
    clf = LogisticRegression(max_iter=2000)
    pred = cross_val_predict(clf, z, y, cv=StratifiedKFold(n_splits, shuffle=True,
                                                           random_state=0))
    acc = float(np.mean(pred == y))
    maj = float(counts.max() / len(y))
    return {"accuracy": acc, "majority": maj, "lift": acc - maj}


def evaluate_fold(y_true, y_pred, proba, classes, y_train, order=None,
                  z=None, dataset_ids=None) -> dict:
    """Everything a fold reports, in one dictionary."""
    out = {
        "n_test": int(len(y_true)),
        "balanced_accuracy": balanced_accuracy(y_true, y_pred, classes),
        "macro_f1": macro_f1(y_true, y_pred, classes),
        "accuracy": float(np.mean(np.asarray(y_true) == np.asarray(y_pred))),
        "majority_baseline": majority_baseline(y_train, y_true),
    }
    if proba is not None:
        out.update({
            "top2_accuracy": top_k_accuracy(y_true, proba, classes, 2),
            "expected_calibration_error": expected_calibration_error(y_true, proba, classes),
            "brier": brier_score(y_true, proba, classes),
        })
    if order:
        out["stage_mae"] = stage_mae(y_true, y_pred, order)
        out["quadratic_weighted_kappa"] = quadratic_weighted_kappa(y_true, y_pred, order)
    if z is not None:
        out["silhouette_condition"] = silhouette(z, y_true)
        if dataset_ids is not None:
            out["silhouette_cohort"] = silhouette(z, dataset_ids)
            out.update({f"cohort_{k}": v for k, v in
                        cohort_predictability(z, dataset_ids).items()})
    out["lift_over_baseline"] = out["accuracy"] - out["majority_baseline"]
    return out
