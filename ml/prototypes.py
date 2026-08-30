"""Class prototypes, balanced across cohorts.

A prototype is the centre of a class in whatever space it is computed in —
latent for the model, raw spectrum for the baseline. The arithmetic is the same,
and so is the trap:

    mu_c = mean(z_i | y_i = c)

is dominated by whichever cohort contributed most of class c. Forty-two of the
seventy F4 samples come from one series, so a plain mean puts the F4 prototype
where that series is, and a query then matches on cohort. The fix is two-stage:
centre each cohort first, then combine those centres robustly, so every cohort
contributes one vote regardless of size.

This module is deliberately numpy-only. It is used by the baselines, by the
statistical prototype, and later by the model, and none of those should have to
agree on a tensor library to agree on what a prototype is.
"""
from __future__ import annotations

import numpy as np


def cohort_balanced_prototype(z: np.ndarray, dataset_ids, combine: str = "median",
                              min_per_cohort: int = 2) -> tuple[np.ndarray, int]:
    """Centre per cohort, then combine the centres.

    Returns the prototype and how many cohorts contributed. A cohort with fewer
    than `min_per_cohort` samples of this class is skipped: a centre computed
    from one sample is that sample, and giving it a full vote would be worse
    than the imbalance being corrected.
    """
    dataset_ids = np.asarray(dataset_ids)
    centres = []
    for ds in sorted(set(dataset_ids.tolist())):
        rows = z[dataset_ids == ds]
        if len(rows) < min_per_cohort:
            continue
        centres.append(np.median(rows, axis=0) if combine == "median" else rows.mean(0))
    if not centres:
        # No cohort qualifies: fall back to the plain mean and say so through
        # the count, so a caller can flag the prototype as unbalanced.
        return z.mean(axis=0), 0
    centres = np.stack(centres)
    proto = np.median(centres, axis=0) if combine == "median" else centres.mean(0)
    return proto, len(centres)


def build_prototypes(z: np.ndarray, class_ids, dataset_ids, combine: str = "median",
                     min_per_cohort: int = 2) -> tuple[np.ndarray, list[str], dict]:
    """One prototype per class. Returns (prototypes, class order, provenance)."""
    class_ids = np.asarray(class_ids)
    classes = sorted(set(class_ids.tolist()))
    protos, info = [], {}
    for c in classes:
        m = class_ids == c
        p, n_cohorts = cohort_balanced_prototype(z[m], np.asarray(dataset_ids)[m],
                                                 combine, min_per_cohort)
        protos.append(p)
        info[c] = {"n_samples": int(m.sum()), "n_cohorts_used": n_cohorts,
                   "balanced": n_cohorts > 0}
    return np.stack(protos), classes, info


def _l2_normalise(a: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    n = np.linalg.norm(a, axis=-1, keepdims=True)
    return a / np.maximum(n, eps)


def distances(z: np.ndarray, prototypes: np.ndarray, metric: str = "cosine") -> np.ndarray:
    """Distance from every sample to every prototype. Smaller is closer."""
    if metric == "cosine":
        return 1.0 - _l2_normalise(z) @ _l2_normalise(prototypes).T
    if metric == "euclidean":
        return np.linalg.norm(z[:, None, :] - prototypes[None, :, :], axis=-1)
    raise ValueError(f"metric must be cosine or euclidean; got {metric!r}")


def prototype_probabilities(z: np.ndarray, prototypes: np.ndarray,
                            metric: str = "cosine", temperature: float = 0.1) -> np.ndarray:
    """Softmax over negative distances. Temperature is a calibration knob.

    These are not calibrated probabilities out of the box — the temperature sets
    how sharp they are, and it has to be chosen on training folds like any other
    parameter. Reported alongside expected calibration error for that reason.
    """
    if temperature <= 0:
        raise ValueError("temperature must be positive")
    logits = -distances(z, prototypes, metric) / temperature
    logits -= logits.max(axis=1, keepdims=True)
    e = np.exp(logits)
    return e / e.sum(axis=1, keepdims=True)


def prototype_separation(prototypes: np.ndarray, metric: str = "cosine") -> dict:
    """How far apart the classes are, and which two are closest.

    A model can reach a good loss with prototypes almost on top of each other:
    the classifier head compensates and the latent space carries no structure. If
    the minimum separation is near zero the prototypes are decorative.
    """
    d = distances(prototypes, prototypes, metric)
    np.fill_diagonal(d, np.inf)
    i, j = np.unravel_index(np.argmin(d), d.shape)
    finite = d[np.isfinite(d)]
    return {"min": float(d[i, j]), "mean": float(finite.mean()),
            "closest_pair": (int(i), int(j))}


def assign(z: np.ndarray, prototypes: np.ndarray, classes: list[str],
           metric: str = "cosine") -> tuple[np.ndarray, np.ndarray]:
    """Nearest-prototype assignment plus the margin to the runner-up."""
    d = distances(z, prototypes, metric)
    order = np.argsort(d, axis=1)
    best = order[:, 0]
    margin = (d[np.arange(len(z)), order[:, 1]] - d[np.arange(len(z)), best]
              if d.shape[1] > 1 else np.full(len(z), np.nan))
    return np.asarray(classes)[best], margin
