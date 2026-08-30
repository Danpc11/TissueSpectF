"""Shared helpers: seeding, config loading, manifests."""
from __future__ import annotations

import hashlib
import json
import os
import platform
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import yaml


def load_config(path: str | Path) -> dict:
    with open(path) as fh:
        return yaml.safe_load(fh)


def set_seed(seed: int) -> None:
    """Seed every generator that can affect a run, and record that we did.

    Reproducibility here means: same config plus same seed gives the same
    numbers within tolerance. It is not a nicety — a peak that moves between
    runs of the same configuration is not a finding.
    """
    random.seed(seed)
    np.random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)
    try:  # torch is optional for the data layer
        import torch

        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        torch.use_deterministic_algorithms(True, warn_only=True)
    except ImportError:
        pass


def resolve_device(requested: str = "auto") -> str:
    """CPU always works; CUDA is used when present and never required."""
    if requested not in ("auto", "cpu", "cuda"):
        raise ValueError(f"device must be auto, cpu or cuda; got {requested!r}")
    if requested == "cpu":
        return "cpu"
    try:
        import torch

        if torch.cuda.is_available():
            return "cuda"
    except ImportError:
        pass
    if requested == "cuda":
        raise RuntimeError("CUDA was requested but is not available")
    return "cpu"


def software_versions() -> dict[str, str]:
    versions = {"python": platform.python_version(), "numpy": np.__version__}
    for name in ("torch", "pandas", "scipy", "sklearn"):
        try:
            mod = __import__(name)
            versions[name] = getattr(mod, "__version__", "?")
        except ImportError:
            versions[name] = "absent"
    return versions


def digest_array(arr: np.ndarray) -> str:
    """Stable fingerprint of an array, for manifests and reproducibility tests."""
    return "sha256:" + hashlib.sha256(np.ascontiguousarray(arr).tobytes()).hexdigest()[:32]


@dataclass(frozen=True)
class Manifest:
    """What a tensor was built from. Travels with the model bundle.

    Compatibility is checked against this before a query is scored: a model
    trained on one annotation, gene universe, grid or expression unit cannot
    score a sample prepared under another. The pipeline refuses rather than
    producing a number nobody can interpret.
    """

    annotation: str
    species: str
    genome_build: str
    annotation_release: str
    gene_universe: str
    grid_digest: str
    expression_unit: str
    k_max: int
    channels: tuple[str, ...]
    chromosome_order: tuple[str, ...]
    grid_n: dict[str, int]
    dataset_ids: tuple[str, ...]
    class_ids: tuple[str, ...]
    seed: int
    software: dict[str, str]

    CRITICAL = (
        "species",
        "genome_build",
        "annotation_release",
        "gene_universe",
        "grid_digest",
        "expression_unit",
        "k_max",
    )

    def incompatibilities(self, other: "Manifest") -> list[str]:
        problems = [
            f"{f}: {getattr(self, f)!r} vs {getattr(other, f)!r}"
            for f in self.CRITICAL
            if getattr(self, f) != getattr(other, f)
        ]
        if self.channels != other.channels:
            problems.append(f"channels: {self.channels} vs {other.channels}")
        if self.chromosome_order != other.chromosome_order:
            problems.append("chromosome_order differs")
        if self.grid_n != other.grid_n:
            problems.append("grid_N differs for at least one chromosome")
        return problems

    def assert_compatible(self, other: "Manifest") -> None:
        problems = self.incompatibilities(other)
        if problems:
            raise ValueError(
                "Incompatible inputs; a model cannot be applied across these:\n  "
                + "\n  ".join(problems)
            )

    def to_json(self, path: str | Path) -> None:
        payload = {k: (list(v) if isinstance(v, tuple) else v) for k, v in self.__dict__.items()}
        Path(path).write_text(json.dumps(payload, indent=2, sort_keys=True))

    @classmethod
    def from_json(cls, path: str | Path) -> "Manifest":
        payload: dict[str, Any] = json.loads(Path(path).read_text())
        payload["channels"] = tuple(payload["channels"])
        payload["chromosome_order"] = tuple(payload["chromosome_order"])
        payload["dataset_ids"] = tuple(payload["dataset_ids"])
        payload["class_ids"] = tuple(payload["class_ids"])
        return cls(**payload)
