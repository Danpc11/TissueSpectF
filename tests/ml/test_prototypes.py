"""Prototypes must describe a class, not the cohort that supplied most of it."""
import numpy as np

from ml.prototypes import (assign, build_prototypes, cohort_balanced_prototype,
                           distances, prototype_probabilities, prototype_separation)


def test_a_large_cohort_does_not_capture_the_prototype():
    # Cohort A contributes 40 samples around +1 and B contributes 4 around -1.
    # A plain mean lands near +0.8, i.e. on A. One vote per cohort lands near 0.
    rng = np.random.default_rng(0)
    z = np.vstack([rng.normal(1.0, 0.05, (40, 4)), rng.normal(-1.0, 0.05, (4, 4))])
    ds = np.array(["A"] * 40 + ["B"] * 4)
    plain = z.mean(axis=0)
    balanced, n = cohort_balanced_prototype(z, ds)
    assert n == 2
    assert abs(plain.mean()) > 0.7
    assert abs(balanced.mean()) < 0.2


def test_a_cohort_with_too_few_samples_gets_no_vote():
    rng = np.random.default_rng(1)
    z = np.vstack([rng.normal(1.0, 0.05, (10, 3)), rng.normal(-5.0, 0.05, (1, 3))])
    ds = np.array(["A"] * 10 + ["B"])
    proto, n = cohort_balanced_prototype(z, ds, min_per_cohort=2)
    assert n == 1 and abs(proto.mean() - 1.0) < 0.2


def test_provenance_records_whether_a_prototype_is_balanced():
    rng = np.random.default_rng(2)
    z = rng.normal(size=(12, 3))
    cl = np.array(["a"] * 6 + ["b"] * 6)
    ds = np.array(["A"] * 3 + ["B"] * 3 + ["A"] * 6)   # class b comes from A only
    _, classes, info = build_prototypes(z, cl, ds)
    assert info["a"]["n_cohorts_used"] == 2 and info["a"]["balanced"]
    assert info["b"]["n_cohorts_used"] == 1


def test_distances_and_assignment_agree():
    protos = np.array([[1.0, 0.0], [0.0, 1.0]])
    z = np.array([[0.9, 0.1], [0.1, 0.9]])
    labels, margin = assign(z, protos, ["x", "y"])
    assert list(labels) == ["x", "y"]
    assert (margin > 0).all()


def test_cosine_ignores_magnitude_and_euclidean_does_not():
    protos = np.array([[1.0, 0.0], [0.0, 1.0]])
    z = np.array([[10.0, 0.0]])
    assert distances(z, protos, "cosine")[0, 0] < 1e-9
    assert distances(z, protos, "euclidean")[0, 0] > 8


def test_probabilities_sum_to_one_and_sharpen_with_temperature():
    protos = np.array([[1.0, 0.0], [0.0, 1.0]])
    z = np.array([[0.8, 0.2]])
    warm = prototype_probabilities(z, protos, temperature=1.0)
    cold = prototype_probabilities(z, protos, temperature=0.05)
    assert np.isclose(warm.sum(), 1.0) and np.isclose(cold.sum(), 1.0)
    assert cold.max() > warm.max()


def test_separation_flags_collapsed_prototypes():
    apart = prototype_separation(np.array([[1.0, 0.0], [0.0, 1.0]]))
    collapsed = prototype_separation(np.array([[1.0, 0.0], [1.0, 1e-9]]))
    assert apart["min"] > 0.5 and collapsed["min"] < 1e-9
