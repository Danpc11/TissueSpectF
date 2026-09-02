#!/usr/bin/env python3
"""Report whether the ml/ test suite can be collected, and what is missing.

Exits 0 when every module is importable, 1 otherwise. On failure it names the
missing modules and the pip packages that provide them, which are not always the
same string -- `import yaml` comes from `pyyaml`, and guessing wrong sends the
reader to a package that does not exist.

Why this exists: the Makefile guard used to check `numpy, sklearn, pytest`. It
demanded sklearn, which ml/ imports lazily inside functions and does not need
to be collected, and it did not check yaml or pandas, which are imported at
module level and do. So the guard passed and pytest then died with
`ModuleNotFoundError: No module named 'yaml'`. A guard that passes and then lets
the thing it guards fail is worse than no guard at all.
"""

import importlib.util
import sys

# import name -> pip name, where they differ.
PIP_NAME = {"yaml": "pyyaml", "sklearn": "scikit-learn"}

modules = sys.argv[1:]
if not modules:
    print("usage: check_ml_deps.py <module> [<module> ...]", file=sys.stderr)
    sys.exit(2)

missing = [m for m in modules if importlib.util.find_spec(m) is None]

if not missing:
    sys.exit(0)

packages = " ".join(PIP_NAME.get(m, m) for m in missing)
print(f"skipping tests/ml: cannot import {', '.join(missing)}")
print(f"  pip install {packages}")
print("  (or: pip install -r requirements-ml.txt)")
print("  The R pipeline and its tests do not depend on any of this.")
sys.exit(1)
