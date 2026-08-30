"""Cohort-aware splitting. The one evaluation that answers the question."""
import numpy as np
import pytest

from ml.splits import (assert_no_cohort_leak, balanced_sample_weights,
                       class_cohort_counts, class_evaluation_report,
                       describe_folds, format_class_report, leave_one_cohort_out)


def toy():
    ds = np.array(["A"] * 6 + ["B"] * 6 + ["C"] * 4)
    cl = np.array(["F0", "F0", "F1", "F1", "F2", "F2"] * 2 + ["F0", "F0", "F1", "F1"])
    return ds, cl


def test_one_fold_per_cohort():
    ds, cl = toy()
    folds = leave_one_cohort_out(ds, cl)
    assert [f.held_out for f in folds] == ["A", "B", "C"]


def test_no_sample_of_the_held_out_cohort_is_in_training():
    ds, cl = toy()
    for f in leave_one_cohort_out(ds, cl):
        assert_no_cohort_leak(f, ds)


def test_a_leak_is_detected():
    ds, cl = toy()
    f = leave_one_cohort_out(ds, cl)[0]
    f.train_idx = np.append(f.train_idx, f.test_idx[0])
    with pytest.raises(AssertionError, match="Cohort leak"):
        assert_no_cohort_leak(f, ds)


def test_classes_absent_from_training_are_dropped_and_counted():
    # C has an F9 that nobody else has: it cannot be predicted, so scoring it as
    # an error would blame the model for the design of the cohorts.
    ds = np.array(["A"] * 4 + ["C"] * 4)
    cl = np.array(["F0", "F0", "F1", "F1", "F0", "F1", "F9", "F9"])
    fold = [f for f in leave_one_cohort_out(ds, cl) if f.held_out == "C"][0]
    assert fold.dropped_classes == ("F9",)
    assert fold.dropped_test_samples == 2
    assert set(cl[fold.test_idx]) == {"F0", "F1"}


def test_a_cohort_sharing_too_little_is_skipped():
    # C contributes only a class nobody else has, so holding C out leaves no
    # scorable class. It yields no fold rather than a fold of zero samples.
    ds = np.array(["A"] * 6 + ["B"] * 4 + ["C"] * 2)
    cl = np.array(["F0", "F0", "F1", "F1", "F2", "F2",
                   "F0", "F0", "F1", "F1", "F9", "F9"])
    folds = leave_one_cohort_out(ds, cl)
    assert [f.held_out for f in folds] == ["A", "B"]


def test_no_usable_fold_raises_rather_than_returning_nothing():
    ds = np.array(["A", "A", "B", "B"])
    cl = np.array(["F0", "F0", "F9", "F9"])
    with pytest.raises(ValueError, match="not possible"):
        leave_one_cohort_out(ds, cl)


def test_balanced_weights_equalise_class_mass():
    ds = np.array(["A"] * 10 + ["B"] * 2)
    cl = np.array(["F4"] * 10 + ["F0"] * 2)
    w = balanced_sample_weights(ds, cl)
    assert np.isclose(w[cl == "F4"].sum(), w[cl == "F0"].sum())
    assert np.isclose(w.sum(), 1.0)


def test_balanced_weights_equalise_cohorts_within_a_class():
    # 42 of 70 F4 samples from one series is the real case: class-only balancing
    # would still show the model that series most of the time.
    ds = np.array(["A"] * 42 + ["B"] * 28)
    cl = np.array(["F4"] * 70)
    w = balanced_sample_weights(ds, cl)
    assert np.isclose(w[ds == "A"].sum(), w[ds == "B"].sum())


def test_class_cohort_counts():
    ds, cl = toy()
    counts = class_cohort_counts(ds, cl)
    assert counts[("F0", "A")] == 2 and counts[("F1", "C")] == 2


def test_describe_folds_mentions_dropped_classes():
    ds = np.array(["A"] * 4 + ["C"] * 4)
    cl = np.array(["F0", "F0", "F1", "F1", "F0", "F1", "F9", "F9"])
    text = describe_folds(leave_one_cohort_out(ds, cl), ds, cl)
    assert "dropped" in text and "F9" in text


def test_class_report_names_classes_never_evaluated():
    # A class living in one cohort is dropped from every fold and contributes
    # nothing to any average. Reporting a single macro-F1 without saying so
    # presents it as if it had been tested like the rest.
    ds = np.array(["A"] * 6 + ["B"] * 4)
    cl = np.array(["F0", "F0", "F1", "F1", "F2", "F2", "F0", "F0", "F1", "F1"])
    folds = leave_one_cohort_out(ds, cl)
    rows = {r["class_id"]: r for r in class_evaluation_report(folds, ds, cl)}
    assert rows["F0"]["evaluated_out_of_cohort"] and rows["F0"]["n_folds_evaluated"] == 2
    assert not rows["F2"]["evaluated_out_of_cohort"]
    assert rows["F2"]["n_samples_tested"] == 0
    text = format_class_report(list(rows.values()))
    assert "NEVER EVALUATED OUT OF COHORT" in text


def test_class_report_counts_cohorts_per_class():
    ds = np.array(["A"] * 4 + ["B"] * 4 + ["C"] * 2)
    cl = np.array(["F0", "F0", "F1", "F1"] * 2 + ["F0", "F1"])
    rows = {r["class_id"]: r for r in
            class_evaluation_report(leave_one_cohort_out(ds, cl), ds, cl)}
    assert rows["F0"]["n_cohorts"] == 3 and rows["F0"]["cohorts"] == "A,B,C"
