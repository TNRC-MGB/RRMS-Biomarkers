#!/usr/bin/env bash
# ==============================================================================
# Record the Python environment that produced MS.score.gz.
# ==============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

python - <<'PY'
import sys
from importlib.metadata import version, PackageNotFoundError

pkgs = ["scdrs", "scanpy", "anndata", "numpy", "pandas", "scipy",
        "scikit-learn", "statsmodels", "h5py", "matplotlib"]
print(f"{'python':<15s} {sys.version.split()[0]}")
for p in pkgs:
    try:
        print(f"{p:<15s} {version(p)}")
    except PackageNotFoundError:
        print(f"{p:<15s} (not installed)")
PY

pip freeze > "$HERE/requirements.txt"
echo
echo "full environment written -> $HERE/requirements.txt"
