"""Tensor assembly, masks and normalisation. No torch needed."""
import numpy as np
import pandas as pd
import pytest

from ml.dataset import (CHANNELS, MASK_ABOVE_NYQUIST, MASK_CHROMOSOME_ABSENT,
                        MASK_NOT_ESTIMABLE, MASK_OBSERVED, apply_normalisation,
                        assert_no_padding_leak, build_tensor, normalisation_stats)


def write_export(tmp, chroms=(("1", 200), ("2", 60)), samples=4, k_max=32,
                 drop=None, seed=0):
    """A minimal but valid ae-prepare export."""
    rng = np.random.default_rng(seed)
    drop = drop or {}
    rows, srows = [], []
    for s in range(samples):
        sid = f"S{s}"
        ds = "A" if s < samples // 2 else "B"
        cls = "liver::disease::NAFLD_fibrosis_F2" if s % 2 else "liver::healthy::Control"
        srows.append(dict(sample_id=sid, dataset_id=ds, tissue="liver",
                          state="disease", condition="F2", class_id=cls, keep=True))
        for chrom, N in chroms:
            if chrom in drop.get(sid, []):
                continue
            for k in range(1, min(k_max, (N - 1) // 2) + 1):
                rows.append(dict(sample_id=sid, dataset_id=ds, tissue="liver",
                                 state="disease", condition="F2", class_id=cls,
                                 chr=chrom, N=N, k=k, freq=k / N, period=N / k,
                                 amplitude=abs(rng.normal(1, 0.2)),
                                 power=abs(rng.normal(1, 0.2)),
                                 phase=rng.uniform(-np.pi, np.pi), coverage=0.8))
    pd.DataFrame(rows).to_csv(tmp / "spectra.tsv", sep="\t", index=False)
    pd.DataFrame(srows).to_csv(tmp / "samples.tsv", sep="\t", index=False)
    pd.DataFrame([{"chr": c, "N": n} for c, n in chroms]).to_csv(
        tmp / "chromosomes.tsv", sep="\t", index=False)
    pd.DataFrame({"key": ["k_max", "grid_digest", "expression_unit", "seed",
                          "genome_build", "gene_universe"],
                  "value": [k_max, "sha256:abc", "asinh(CPM)", 42,
                            "GRCh38", "^protein-coding$"]}).to_csv(
        tmp / "manifest.tsv", sep="\t", index=False)
    return tmp


def test_shape_and_channel_order(tmp_path):
    t = build_tensor(write_export(tmp_path), k_max=32)
    assert t.shape == (4, 2, 32, len(CHANNELS))
    assert t.manifest.channels == CHANNELS
    assert t.manifest.chromosome_order == ("1", "2")


def test_padding_is_exactly_zero(tmp_path):
    t = build_tensor(write_export(tmp_path), k_max=32)
    assert_no_padding_leak(t.x, t.mask)


def test_above_nyquist_is_structural_not_missing(tmp_path):
    # chr2 has N=60, so k > 29 does not exist. That is not a failed measurement
    # and must not be confused with one: it is identical for every sample and
    # can never carry cohort information.
    t = build_tensor(write_export(tmp_path), k_max=32)
    ci = t.manifest.chromosome_order.index("2")
    assert (t.mask_reason[:, ci, 29:] == MASK_ABOVE_NYQUIST).all()
    assert not t.mask[:, ci, 29:].any()
    ci1 = t.manifest.chromosome_order.index("1")
    assert (t.mask_reason[:, ci1, :32] == MASK_OBSERVED).all()


def test_absent_chromosome_is_distinguished(tmp_path):
    t = build_tensor(write_export(tmp_path, drop={"S1": ["2"]}), k_max=32)
    si = list(t.sample_ids).index("S1")
    ci = t.manifest.chromosome_order.index("2")
    assert (t.mask_reason[si, ci, :29] == MASK_CHROMOSOME_ABSENT).all()
    # still structural above Nyquist, not "chromosome absent"
    assert (t.mask_reason[si, ci, 29:] == MASK_ABOVE_NYQUIST).all()
    assert not t.mask[si, ci].any()


def test_coverage_excludes_structural_absence(tmp_path):
    # A sample observing everything it could must report coverage 1.0, even
    # though a third of the array is above chr2's Nyquist limit.
    t = build_tensor(write_export(tmp_path), k_max=32)
    assert np.allclose(t.sample_coverage(), 1.0)


def test_coverage_falls_when_a_chromosome_is_missing(tmp_path):
    t = build_tensor(write_export(tmp_path, drop={"S1": ["2"]}), k_max=32)
    cov = dict(zip(t.sample_ids, t.sample_coverage()))
    assert cov["S1"] < cov["S0"]


def test_phase_enters_as_cos_and_sin(tmp_path):
    t = build_tensor(write_export(tmp_path), k_max=32)
    c, s = t.channel("cos_phase"), t.channel("sin_phase")
    m = t.mask
    assert np.allclose(c[m] ** 2 + s[m] ** 2, 1.0, atol=1e-6)


def test_normalisation_uses_training_samples_only(tmp_path):
    t = build_tensor(write_export(tmp_path), k_max=32)
    train = np.array([0, 1])
    stats = normalisation_stats(t, train)
    # Corrupting a held-out sample must not move the statistics at all.
    t2 = build_tensor(write_export(tmp_path), k_max=32)
    t2.x[2:] *= 1000.0
    assert stats == normalisation_stats(t2, train)


def test_normalisation_re_zeroes_padding(tmp_path):
    # Standardising before re-zeroing would turn every unobserved slot into
    # -center, a constant the network could read as a missingness flag on
    # channels where the mask is meant to be the only route.
    t = build_tensor(write_export(tmp_path, drop={"S1": ["2"]}), k_max=32)
    stats = normalisation_stats(t, np.arange(t.n_samples))
    x = apply_normalisation(t, stats)
    assert np.abs(x[~t.mask]).max() == 0.0


def test_inconsistent_export_is_refused(tmp_path):
    d = write_export(tmp_path)
    sp = pd.read_csv(d / "spectra.tsv", sep="\t", dtype={"chr": str})
    sp.loc[0, "chr"] = "99"
    sp.to_csv(d / "spectra.tsv", sep="\t", index=False)
    with pytest.raises(ValueError, match="absent from chromosomes.tsv"):
        build_tensor(d, k_max=32)


def test_missing_columns_are_refused(tmp_path):
    d = write_export(tmp_path)
    sp = pd.read_csv(d / "spectra.tsv", sep="\t").drop(columns=["phase"])
    sp.to_csv(d / "spectra.tsv", sep="\t", index=False)
    with pytest.raises(ValueError, match="missing columns"):
        build_tensor(d, k_max=32)
