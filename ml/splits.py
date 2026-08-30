"""Cohort-aware splitting and balanced sampling.

WHY THIS FILE EXISTS BEFORE THE MODEL
-------------------------------------
Leave-one-cohort-out is not one evaluation option among several here; it is the
only one that answers the question. With thousands of features and a few hundred
samples, a random split gives a high accuracy by learning batch, library
protocol and sequencing depth — none of which transfer to a new cohort, all of
which are perfectly available inside a cohort. A model that looks excellent
under random splitting and useless under LOCO has not failed; it has been
measured correctly for the first time.

So the splitter comes first, and it refuses to produce a split that leaks.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class Fold:
    """One leave-one-cohort-out fold."""

    held_out: str
    train_idx: np.ndarray
    test_idx: np.ndarray
    shared_classes: tuple[str, ...]
    dropped_classes: tuple[str, ...] = field(default_factory=tuple)
    dropped_test_samples: int = 0

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return (f"Fold(held_out={self.held_out!r}, n_train={len(self.train_idx)}, "
                f"n_test={len(self.test_idx)}, classes={len(self.shared_classes)})")


def leave_one_cohort_out(dataset_ids, class_ids, min_classes: int = 2) -> list[Fold]:
    """One fold per cohort, restricted to classes present on both sides.

    A class the training folds never saw cannot be predicted, and a class absent
    from the test cohort cannot be scored. Both are dropped from the fold and
    counted, rather than being silently carried as impossible errors that drag
    the reported accuracy down for a reason that has nothing to do with the
    model.
    """
    dataset_ids = np.asarray(dataset_ids)
    class_ids = np.asarray(class_ids)
    if dataset_ids.shape != class_ids.shape:
        raise ValueError("dataset_ids and class_ids must have the same length")

    folds: list[Fold] = []
    for held in sorted(set(dataset_ids.tolist())):
        test_mask = dataset_ids == held
        train_mask = ~test_mask
        if not train_mask.any():
            continue
        shared = sorted(set(class_ids[train_mask]) & set(class_ids[test_mask]))
        if len(shared) < min_classes:
            continue
        dropped = sorted(set(class_ids[test_mask]) - set(shared))
        keep_tr = train_mask & np.isin(class_ids, shared)
        keep_te = test_mask & np.isin(class_ids, shared)
        folds.append(Fold(
            held_out=str(held),
            train_idx=np.where(keep_tr)[0],
            test_idx=np.where(keep_te)[0],
            shared_classes=tuple(shared),
            dropped_classes=tuple(dropped),
            dropped_test_samples=int((test_mask & ~keep_te).sum()),
        ))
    if not folds:
        raise ValueError(
            "No usable leave-one-cohort-out fold: every cohort shares fewer than "
            f"{min_classes} classes with the rest. Cross-cohort evaluation is not "
            "possible with these data."
        )
    return folds


def assert_no_cohort_leak(fold: Fold, dataset_ids) -> None:
    """A sample of the held-out cohort must never appear in training.

    Cheap to check, catastrophic to get wrong, and easy to reintroduce by
    reindexing a matrix somewhere upstream. Asserted rather than assumed.
    """
    dataset_ids = np.asarray(dataset_ids)
    if set(dataset_ids[fold.train_idx]) & {fold.held_out}:
        raise AssertionError(
            f"Cohort leak: samples from {fold.held_out} appear in the training "
            "set of its own fold."
        )
    if set(dataset_ids[fold.test_idx]) != {fold.held_out}:
        raise AssertionError("The test set contains samples from another cohort")
    if set(fold.train_idx) & set(fold.test_idx):
        raise AssertionError("Train and test indices overlap")


def class_cohort_counts(dataset_ids, class_ids) -> dict[tuple[str, str], int]:
    keys, counts = np.unique(
        np.stack([np.asarray(class_ids), np.asarray(dataset_ids)], axis=1),
        axis=0, return_counts=True)
    return {(str(c), str(d)): int(n) for (c, d), n in zip(keys, counts)}


def balanced_sample_weights(dataset_ids, class_ids, idx=None,
                            max_per_cell: int | None = None) -> np.ndarray:
    """Weights that balance condition x cohort cells, not classes alone.

    Balancing by class alone still lets one cohort supply most of a class: 42 of
    the 70 F4 samples come from a single series, so a class-balanced sampler
    would show the model that cohort's F4 most of the time and it would learn
    the cohort. Weighting each (class, cohort) cell equally is what makes the
    class the thing in common.

    `max_per_cell` additionally caps a cell's total influence, for the case
    where one cell is so large that even equal weighting gives it every draw
    within the class.
    """
    dataset_ids = np.asarray(dataset_ids)
    class_ids = np.asarray(class_ids)
    idx = np.arange(len(class_ids)) if idx is None else np.asarray(idx)

    cells: dict[tuple[str, str], list[int]] = {}
    for i in idx:
        cells.setdefault((str(class_ids[i]), str(dataset_ids[i])), []).append(int(i))

    weights = np.zeros(len(class_ids), dtype=np.float64)
    n_classes = len({c for c, _ in cells})
    for (cls, _ds), members in cells.items():
        cells_in_class = sum(1 for c, _ in cells if c == cls)
        effective = len(members) if max_per_cell is None else min(len(members), max_per_cell)
        # Each class gets equal total mass; within a class, each cohort gets an
        # equal share of it; within a cohort, samples share that equally.
        per_sample = 1.0 / (n_classes * cells_in_class * max(effective, 1))
        weights[members] = per_sample * (effective / len(members))
    total = weights.sum()
    return weights / total if total > 0 else weights


def describe_folds(folds: list[Fold], dataset_ids, class_ids) -> str:
    """A human-readable summary, printed before training rather than after."""
    lines = []
    for f in folds:
        n_by_class = {c: int((np.asarray(class_ids)[f.test_idx] == c).sum())
                      for c in f.shared_classes}
        lines.append(
            f"hold out {f.held_out}: train {len(f.train_idx)}, test {len(f.test_idx)}"
            f" over {len(f.shared_classes)} class(es) "
            + ", ".join(f"{c.split('::')[-1]}={n}" for c, n in n_by_class.items()))
        if f.dropped_classes:
            lines.append(
                f"    {f.dropped_test_samples} test sample(s) dropped: "
                f"{', '.join(c.split('::')[-1] for c in f.dropped_classes)} "
                "absent from training")
    return "\n".join(lines)
