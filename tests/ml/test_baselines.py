"""The bar the autoencoder has to clear, and the metrics that define it."""
import numpy as np

from ml.baselines import NearestCentroidBaseline, flatten, format_summary, run_baselines
from ml.evaluate import (balanced_accuracy, brier_score, cohort_predictability,
                         expected_calibration_error, macro_f1, majority_baseline,
                         quadratic_weighted_kappa, silhouette, stage_mae)

ORDER = ["F0", "F1", "F2", "F3"]


def synthetic(cohort_shift=2.0, seed=0, per_cell=8):
    """Three cohorts, four ordered classes, one real peak, one cohort offset."""
    rng = np.random.default_rng(seed)
    n_chr, k_max = 3, 16
    X, ds, cl = [], [], []
    for ci, coh in enumerate("ABC"):
        for si, cls in enumerate(ORDER):
            for _ in range(per_cell):
                x = rng.normal(0, 0.3, (n_chr, k_max, 6))
                x[..., 5] = 1.0
                x[0, 5, 0] += 1.5 * si            # the class signal
                x[..., 0] += cohort_shift * ci    # the batch offset
                X.append(x); ds.append(coh); cl.append(cls)
    X = np.array(X)
    return X, np.ones(X.shape[:3], bool), np.array(cl), np.array(ds)


def test_metrics_agree_with_hand_computation():
    y = ["a", "a", "b", "b"]
    p = ["a", "b", "b", "b"]
    assert balanced_accuracy(y, p, ["a", "b"]) == 0.75
    assert 0 < macro_f1(y, p, ["a", "b"]) < 1


def test_majority_baseline_uses_the_training_fold():
    # Taking the majority of the held-out truth would let the baseline see the
    # test labels, which is what the held-out design exists to prevent.
    y_train = ["F0"] * 9 + ["F4"]
    y_test = ["F4"] * 10
    assert majority_baseline(y_train, y_test) == 0.0


def test_stage_error_and_kappa_are_ordinal():
    assert stage_mae(["F0"], ["F3"], ORDER) == 3.0
    assert quadratic_weighted_kappa(ORDER, ORDER, ORDER) == 1.0


def test_calibration_penalises_overconfidence():
    y = ["a"] * 10
    sure_and_wrong = np.tile([0.01, 0.99], (10, 1))
    sure_and_right = np.tile([0.99, 0.01], (10, 1))
    assert (expected_calibration_error(y, sure_and_wrong, ["a", "b"])
            > expected_calibration_error(y, sure_and_right, ["a", "b"]))
    assert brier_score(y, sure_and_right, ["a", "b"]) < brier_score(y, sure_and_wrong, ["a", "b"])


def test_silhouette_detects_what_a_space_is_organised_by():
    rng = np.random.default_rng(3)
    # Two tight clusters that correspond to cohort, not to condition.
    z = np.vstack([rng.normal(0, 0.1, (20, 2)), rng.normal(5, 0.1, (20, 2))])
    cohort = np.array(["A"] * 20 + ["B"] * 20)
    condition = np.array(["x", "y"] * 20)
    assert silhouette(z, cohort) > 0.8
    assert silhouette(z, condition) < 0.2


def test_cohort_predictability_flags_a_batch_dominated_space():
    rng = np.random.default_rng(4)
    z = np.vstack([rng.normal(0, 0.1, (30, 3)), rng.normal(5, 0.1, (30, 3))])
    ds = np.array(["A"] * 30 + ["B"] * 30)
    r = cohort_predictability(z, ds)
    assert r["accuracy"] > 0.95 and r["lift"] > 0.4


def test_baselines_run_over_every_cohort_fold():
    X, mask, cl, ds = synthetic()
    res = run_baselines(X, mask, cl, ds, order=ORDER)
    assert {r["held_out"] for r in res["folds"]} == {"A", "B", "C"}
    assert "nearest_centroid" in res["summary"]
    assert res["summary"]["nearest_centroid"]["n_folds"] == 3


def test_cohort_balanced_centroid_survives_a_batch_offset():
    # The offset is additive and identical within a cohort, so a cosine centroid
    # built one-vote-per-cohort should still find the class signal. This is the
    # result that sets the bar: if a linear centroid does this, an encoder has
    # to do better to be worth its complexity.
    X, mask, cl, ds = synthetic(cohort_shift=2.0)
    res = run_baselines(X, mask, cl, ds, order=ORDER)
    s = res["summary"]["nearest_centroid"]
    assert s["balanced_accuracy"] > s["majority_baseline"]


def test_a_failing_baseline_does_not_stop_the_others():
    class Broken:
        classes_ = []
        def fit(self, *a, **k): raise RuntimeError("boom")
        def predict(self, x): ...
        def predict_proba(self, x): ...

    X, mask, cl, ds = synthetic(per_cell=4)
    res = run_baselines(X, mask, cl, ds, order=ORDER,
                        models={"broken": Broken(),
                                "nearest_centroid": NearestCentroidBaseline()})
    assert any("error" in r for r in res["folds"])
    assert "nearest_centroid" in res["summary"]


def test_flatten_keeps_biological_channels_only():
    X, mask, cl, ds = synthetic(per_cell=2)
    flat = flatten(X, mask)
    assert flat.shape == (len(X), 3 * 16 * 4)   # coverage and mask excluded


def test_summary_table_renders():
    X, mask, cl, ds = synthetic(per_cell=4)
    text = format_summary(run_baselines(X, mask, cl, ds, order=ORDER))
    assert "nearest_centroid" in text and "balanced_accur" in text
