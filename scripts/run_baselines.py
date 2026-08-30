#!/usr/bin/env python3
"""Run the baselines over the exported spectra, leave-one-cohort-out.

    python3 scripts/run_baselines.py --data results_pc/autoencoder/data \
        --out results_pc/autoencoder/baselines

This is usable the moment `./tsf ae-prepare` has run, before any model exists,
and it is what the autoencoder will have to beat. Needs numpy, pandas and
scikit-learn; torch is not involved.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ml.baselines import format_summary, run_baselines            # noqa: E402
from ml.dataset import apply_normalisation, build_tensor, normalisation_stats  # noqa: E402
from ml.splits import describe_folds, leave_one_cohort_out        # noqa: E402
from ml.utils import set_seed, software_versions                  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True, help="output of ./tsf ae-prepare")
    ap.add_argument("--out", default=None)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--ordinal", default=None,
                    help="comma-separated ordered class ids for stage metrics")
    args = ap.parse_args()
    set_seed(args.seed)

    tensor = build_tensor(args.data)
    print(f"tensor {tensor.shape}  classes {len(tensor.manifest.class_ids)}  "
          f"cohorts {len(tensor.manifest.dataset_ids)}")

    folds = leave_one_cohort_out(tensor.dataset_ids, tensor.class_ids)
    print(describe_folds(folds, tensor.dataset_ids, tensor.class_ids))

    # Normalisation is fitted per fold inside run_baselines' models where it
    # matters; here the tensor is only re-zeroed so padding cannot leak.
    stats = normalisation_stats(tensor, np.arange(tensor.n_samples))
    x = apply_normalisation(tensor, stats)

    order = args.ordinal.split(",") if args.ordinal else [
        c for c in sorted(tensor.manifest.class_ids) if "::disease::" in c]

    result = run_baselines(x, tensor.mask, tensor.class_ids, tensor.dataset_ids,
                           order=order, seed=args.seed)
    print()
    print(format_summary(result))

    if args.out:
        out = Path(args.out)
        out.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(result["folds"]).to_csv(out / "baseline_folds.tsv",
                                             sep="\t", index=False)
        pd.DataFrame(result["summary"]).T.to_csv(out / "baseline_summary.tsv", sep="\t")
        (out / "provenance.json").write_text(json.dumps(
            {"seed": args.seed, "data": str(args.data),
             "software": software_versions(),
             "manifest": {k: str(v) for k, v in tensor.manifest.__dict__.items()
                          if k not in ("grid_n", "software")}}, indent=2))
        print(f"\nwritten to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
