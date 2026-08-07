#!/usr/bin/env python3
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham / Harvard Medical School
#
# Updated: 7-25-2026
#
#
# Usage:
#   python scdrs/05_make_covariates.py <input.h5ad> <output.cov>
#
#   python scdrs/05_make_covariates.py scdrs/work/ms_pbmc.h5ad   scdrs/work/ms_pbmc.cov
#   python scdrs/05_make_covariates.py scdrs/work/ms_bcells.h5ad scdrs/work/ms_bcells.cov
# ==============================================================================

import argparse
import sys


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Build a scDRS covariate file from an AnnData object."
    )
    ap.add_argument("h5ad", help="input .h5ad (raw counts in .X)")
    ap.add_argument("out", help="output covariate file (tab separated)")
    ap.add_argument(
        "--mt-prefix",
        default="MT-",
        help="gene-name prefix identifying mitochondrial genes (default: MT-)",
    )
    args = ap.parse_args()

    # imported here so that --help works without the scientific stack installed
    import numpy as np
    import pandas as pd
    import scanpy as sc

    print(f"reading {args.h5ad} ...", flush=True)
    adata = sc.read_h5ad(args.h5ad)
    print(f"  {adata.n_obs} cells x {adata.n_vars} genes")

    adata.var["mt"] = adata.var_names.str.startswith(args.mt_prefix)
    n_mt = int(adata.var["mt"].sum())
    print(f"  {n_mt} mitochondrial genes matched '{args.mt_prefix}*'")
    if n_mt == 0:
        print(
            "  WARNING: no mitochondrial genes matched; pct_counts_mt will be 0 "
            "for every cell.",
            file=sys.stderr,
        )

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, inplace=True
    )

    cov = pd.DataFrame(index=adata.obs_names)
    cov["log1p_total_counts"] = np.log1p(adata.obs["total_counts"].values)
    cov["n_genes_by_counts"] = adata.obs["n_genes_by_counts"].values
    cov["pct_counts_mt"] = adata.obs["pct_counts_mt"].values

    cov.to_csv(args.out, sep="\t")
    print(f"wrote {args.out}  ({cov.shape[0]} cells x {cov.shape[1]} covariates)")
    print(cov.describe().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
