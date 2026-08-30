"""Absence must never become signal.

These tests exist because the failure they guard against is invisible: a model
fed structured zeros trains normally, scores well internally, and has learned
the pattern of missingness — which correlates with cohort — instead of biology.
"""
import numpy as np
import pandas as pd

from ml.dataset import (MASK_ABOVE_NYQUIST, MASK_CHROMOSOME_ABSENT,
                        MASK_NOT_ESTIMABLE, apply_normalisation,
                        assert_no_padding_leak, build_tensor, normalisation_stats)
from tests.ml.test_dataset import write_export


def test_missing_values_never_become_zero_expression(tmp_path):
    # Two exports identical except that one sample lost a chromosome. The values
    # that remain must be untouched: masking removes positions, it does not
    # rescale or shift what is left.
    (tmp_path / "a").mkdir(); (tmp_path / "b").mkdir()
    full = build_tensor(write_export(tmp_path / "a", seed=1), k_max=32)
    part = build_tensor(write_export(tmp_path / "b", drop={"S1": ["2"]}, seed=1), k_max=32)
    si = list(full.sample_ids).index("S1")
    ci = full.manifest.chromosome_order.index("1")
    assert np.allclose(full.x[si, ci], part.x[si, ci])


def test_padding_survives_normalisation(tmp_path):
    t = build_tensor(write_export(tmp_path, drop={"S1": ["2"]}), k_max=32)
    stats = normalisation_stats(t, np.arange(t.n_samples))
    x = apply_normalisation(t, stats)
    assert_no_padding_leak(x, t.mask)


def test_a_leak_is_detected(tmp_path):
    t = build_tensor(write_export(tmp_path, drop={"S1": ["2"]}), k_max=32)
    x = t.x.copy()
    si = list(t.sample_ids).index("S1")
    ci = t.manifest.chromosome_order.index("2")
    x[si, ci, 0, 0] = 0.001            # one value written where there is none
    try:
        assert_no_padding_leak(x, t.mask)
    except AssertionError:
        return
    raise AssertionError("a padding leak went undetected")


def test_mask_reasons_are_mutually_exclusive(tmp_path):
    t = build_tensor(write_export(tmp_path, drop={"S1": ["2"]}), k_max=32)
    reasons = set(np.unique(t.mask_reason).tolist())
    assert reasons <= {2, MASK_ABOVE_NYQUIST, MASK_CHROMOSOME_ABSENT, MASK_NOT_ESTIMABLE}
    # observed_mask is true exactly where the reason says observed
    assert (t.mask == (t.mask_reason == 2)).all()


def test_structural_absence_is_identical_across_samples(tmp_path):
    # k above a chromosome's Nyquist limit does not exist for anyone. If it
    # varied by sample it could carry cohort identity; it must not.
    t = build_tensor(write_export(tmp_path, drop={"S1": ["2"]}), k_max=32)
    struct = t.mask_reason == MASK_ABOVE_NYQUIST
    assert (struct == struct[0]).all()
