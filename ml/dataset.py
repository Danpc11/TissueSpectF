"""Turn the exported spectra into a masked tensor.

THE CENTRAL RULE OF THIS FILE
-----------------------------
A frequency that was not observed is not zero. It is unobserved, and every loss
must ignore it. Arrays are rectangular, so unobserved positions do get the
number 0.0 written into them — but only after the mask exists, and only as
computational padding. Anything downstream that reads a value without consulting
`observed_mask` is reading padding as if it were biology.

This is the same rule the R pipeline enforces on the genomic grid, where absent
genes stay out of the observed positions rather than being filled with zeros. It
is repeated here because the failure mode is the same and the consequence is
worse: a network handed structured zeros will happily learn the pattern of
missingness, which correlates with cohort, and report it as signal.

MASK LEVELS
-----------
`observed_mask` is binary because losses need it that way, but the reason a
position is missing is kept separately in `mask_reason` for diagnostics:

    2  observed
    1  frequency above this chromosome's Nyquist limit (structurally absent:
       k > (N-1)/2 simply does not exist, it was never measurable)
    0  chromosome absent for this sample (dropped at ingest for coverage)
   -1  frequency not estimable for this sample though the chromosome is present

Distinguishing 1 from -1 matters: the first is a property of the grid and is
identical for every sample, so a model could never learn cohort identity from
it; the second varies by sample and is exactly the kind of thing an adversary
would exploit.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .utils import Manifest, digest_array, software_versions

CHANNELS = ("log_power", "amplitude", "cos_phase", "sin_phase", "coverage", "observed_mask")

MASK_OBSERVED = 2
MASK_ABOVE_NYQUIST = 1
MASK_CHROMOSOME_ABSENT = 0
MASK_NOT_ESTIMABLE = -1


@dataclass
class SpectraTensor:
    """sample x chromosome x frequency x channel, plus labels and mask reasons."""

    x: np.ndarray            # (n_samples, n_chr, k_max, n_channels) float32
    mask: np.ndarray         # (n_samples, n_chr, k_max) bool -- observed
    mask_reason: np.ndarray  # (n_samples, n_chr, k_max) int8
    sample_ids: np.ndarray
    dataset_ids: np.ndarray
    class_ids: np.ndarray
    coverage: np.ndarray     # (n_samples, n_chr) float32
    manifest: Manifest

    @property
    def n_samples(self) -> int:
        return self.x.shape[0]

    @property
    def shape(self) -> tuple[int, ...]:
        return self.x.shape

    def sample_coverage(self) -> np.ndarray:
        """Fraction of the grid a sample actually observed, per sample.

        Structurally absent frequencies (above Nyquist) are excluded from the
        denominator: they are not missing, they do not exist.
        """
        usable = self.mask_reason != MASK_ABOVE_NYQUIST
        seen = (self.mask & usable).reshape(self.n_samples, -1).sum(1)
        total = usable.reshape(self.n_samples, -1).sum(1)
        return (seen / np.maximum(total, 1)).astype(np.float32)

    def channel(self, name: str) -> np.ndarray:
        return self.x[..., self.manifest.channels.index(name)]

    def digest(self) -> str:
        return digest_array(self.x)


def _read_manifest_tsv(path: Path) -> dict[str, str]:
    df = pd.read_csv(path, sep="\t", dtype=str)
    return dict(zip(df["key"], df["value"]))


def build_tensor(data_dir: str | Path, k_max: int | None = None) -> SpectraTensor:
    """Read the R export and assemble the tensor.

    Raises rather than guesses whenever the export is internally inconsistent:
    a silently mis-assembled tensor produces a model that trains fine and means
    nothing.
    """
    data_dir = Path(data_dir)
    spectra = pd.read_csv(data_dir / "spectra.tsv", sep="\t")
    samples = pd.read_csv(data_dir / "samples.tsv", sep="\t")
    chroms = pd.read_csv(data_dir / "chromosomes.tsv", sep="\t", dtype={"chr": str})
    meta = _read_manifest_tsv(data_dir / "manifest.tsv")

    k_max = int(k_max or meta["k_max"])
    chrom_order = tuple(str(c) for c in chroms["chr"])
    grid_n = {str(r.chr): int(r.N) for r in chroms.itertuples()}

    required = {"sample_id", "dataset_id", "class_id", "chr", "N", "k",
                "amplitude", "power", "phase"}
    missing = required - set(spectra.columns)
    if missing:
        raise ValueError(f"spectra.tsv is missing columns: {sorted(missing)}")

    spectra = spectra.copy()
    spectra["chr"] = spectra["chr"].astype(str)
    spectra = spectra[spectra["k"].between(1, k_max)]
    if spectra.empty:
        raise ValueError("No rows with 1 <= k <= k_max; nothing to build")

    keep = samples[samples.get("keep", True).astype(bool)] if "keep" in samples else samples
    sample_ids = np.array(sorted(set(spectra["sample_id"]) & set(keep["sample_id"])))
    if sample_ids.size == 0:
        raise ValueError("No sample appears in both spectra.tsv and samples.tsv")

    s_index = {s: i for i, s in enumerate(sample_ids)}
    c_index = {c: i for i, c in enumerate(chrom_order)}
    n_s, n_c, n_ch = len(sample_ids), len(chrom_order), len(CHANNELS)

    x = np.zeros((n_s, n_c, k_max, n_ch), dtype=np.float32)
    mask = np.zeros((n_s, n_c, k_max), dtype=bool)
    reason = np.full((n_s, n_c, k_max), MASK_NOT_ESTIMABLE, dtype=np.int8)
    coverage = np.zeros((n_s, n_c), dtype=np.float32)

    # Structural absence first, identical for every sample: k above (N-1)/2 is
    # not a measurement that failed, it is a frequency the grid cannot carry.
    for c, ci in c_index.items():
        nyq = (grid_n[c] - 1) // 2
        if nyq < k_max:
            reason[:, ci, nyq:] = MASK_ABOVE_NYQUIST

    rows = spectra[spectra["sample_id"].isin(s_index)]
    si = rows["sample_id"].map(s_index).to_numpy()
    ci = rows["chr"].map(c_index)
    valid = ci.notna().to_numpy()
    if not valid.all():
        dropped = sorted(set(rows["chr"][~valid]))
        raise ValueError(
            f"spectra.tsv contains chromosomes absent from chromosomes.tsv: {dropped}. "
            "The export is inconsistent; re-run ae-prepare."
        )
    ci = ci.to_numpy().astype(int)
    ki = rows["k"].to_numpy().astype(int) - 1

    power = np.asarray(rows["power"], dtype=np.float64)
    amp = np.asarray(rows["amplitude"], dtype=np.float64)
    phase = np.asarray(rows["phase"], dtype=np.float64)
    cov = (np.asarray(rows["coverage"], dtype=np.float64)
           if "coverage" in rows else np.ones(len(rows)))

    finite = np.isfinite(power) & np.isfinite(amp) & np.isfinite(phase)

    x[si[finite], ci[finite], ki[finite], 0] = np.log1p(np.maximum(power[finite], 0))
    x[si[finite], ci[finite], ki[finite], 1] = amp[finite]
    x[si[finite], ci[finite], ki[finite], 2] = np.cos(phase[finite])
    x[si[finite], ci[finite], ki[finite], 3] = np.sin(phase[finite])
    x[si[finite], ci[finite], ki[finite], 4] = cov[finite]
    x[si[finite], ci[finite], ki[finite], 5] = 1.0
    mask[si[finite], ci[finite], ki[finite]] = True
    reason[si[finite], ci[finite], ki[finite]] = MASK_OBSERVED
    coverage[si[finite], ci[finite]] = cov[finite]

    # A chromosome with no observed frequency for a sample was dropped at
    # ingest for that sample; say so rather than leaving it as "not estimable".
    empty_chr = ~mask.any(axis=2)
    for s, c in zip(*np.where(empty_chr)):
        keep_struct = reason[s, c] == MASK_ABOVE_NYQUIST
        reason[s, c][~keep_struct] = MASK_CHROMOSOME_ABSENT

    labels = keep.set_index("sample_id").reindex(sample_ids)
    manifest = Manifest(
        annotation=meta.get("annotation", "?"),
        species=meta.get("species", "?"),
        genome_build=meta.get("genome_build", "?"),
        annotation_release=meta.get("annotation_release", "?"),
        gene_universe=meta.get("gene_universe", "?"),
        grid_digest=meta.get("grid_digest", "?"),
        expression_unit=meta.get("expression_unit", "?"),
        k_max=k_max,
        channels=CHANNELS,
        chromosome_order=chrom_order,
        grid_n=grid_n,
        dataset_ids=tuple(sorted(set(labels["dataset_id"].astype(str)))),
        class_ids=tuple(sorted(set(labels["class_id"].astype(str)))),
        seed=int(meta.get("seed", 42)),
        software=software_versions(),
    )

    return SpectraTensor(
        x=x, mask=mask, mask_reason=reason, sample_ids=sample_ids,
        dataset_ids=labels["dataset_id"].astype(str).to_numpy(),
        class_ids=labels["class_id"].astype(str).to_numpy(),
        coverage=coverage, manifest=manifest,
    )


def assert_no_padding_leak(x: np.ndarray, mask: np.ndarray) -> None:
    """Every unobserved position must hold exactly zero in every channel.

    Guards the invariant this module exists to protect. If a non-zero value ever
    reaches a masked position, some code path wrote data where there is none,
    and the model would be free to read it.
    """
    bad = np.abs(x[~mask]).max(initial=0.0)
    if bad > 0:
        raise AssertionError(
            f"Unobserved positions carry non-zero values (max {bad:g}). "
            "Padding is being written as if it were data."
        )


def normalisation_stats(tensor: SpectraTensor, train_idx: np.ndarray) -> dict:
    """Per-channel robust location and scale, from the TRAINING samples only.

    Median and IQR rather than mean and sd: spectra are heavy-tailed and a
    single dominant component would otherwise set the scale for everything.

    Computing these on anything but the training fold is the most common way to
    leak a held-out cohort into a model while every other precaution looks
    correct — the test cohort's own scale would be baked into the inputs.
    """
    stats = {}
    m = tensor.mask[train_idx]
    for i, name in enumerate(tensor.manifest.channels):
        if name in ("observed_mask",):
            stats[name] = {"center": 0.0, "scale": 1.0}
            continue
        vals = tensor.x[train_idx, ..., i][m]
        vals = vals[np.isfinite(vals)]
        if vals.size == 0:
            stats[name] = {"center": 0.0, "scale": 1.0}
            continue
        centre = float(np.median(vals))
        q1, q3 = np.percentile(vals, [25, 75])
        scale = float(q3 - q1) or float(np.std(vals)) or 1.0
        stats[name] = {"center": centre, "scale": scale}
    return stats


def apply_normalisation(tensor: SpectraTensor, stats: dict) -> np.ndarray:
    """Standardise, then re-zero the padding.

    Order matters: subtracting a centre would otherwise turn every unobserved
    zero into -centre, which is a constant the network can detect and use as a
    missingness flag on channels where the mask was meant to be the only route.
    """
    x = tensor.x.copy()
    for i, name in enumerate(tensor.manifest.channels):
        if name in ("observed_mask", "cos_phase", "sin_phase"):
            continue  # a mask is a mask; cos/sin are already bounded
        s = stats[name]
        x[..., i] = (x[..., i] - s["center"]) / (s["scale"] or 1.0)
    x[~tensor.mask] = 0.0
    return x
