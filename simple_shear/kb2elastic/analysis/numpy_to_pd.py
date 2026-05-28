#!/usr/bin/env python3
# saves .npy files from compute_elastic_moduli.py to one dataset; workaround for ovito package incompatibility 

import os
import argparse
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description="Build pandas DataFrame (cases mode only) from minimal .npy arrays")
    ap.add_argument("--base-dir", required=True, help="Folder containing the .npy files")
    ap.add_argument("--out-name", default="GE_results.parquet", help="Output parquet filename")
    args = ap.parse_args()

    base = os.path.abspath(args.base_dir)

    # required arrays
    kb = np.load(os.path.join(base, "kb.npy"))
    eps = np.load(os.path.join(base, "eps.npy"))
    Gm = np.load(os.path.join(base, "G_mean.npy"))
    Em = np.load(os.path.join(base, "E_mean.npy"))

    n = kb.shape[0]
    for arr, nm in [(eps, "eps"), (Gm, "G_mean"), (Em, "E_mean")]:
        if arr.shape[0] != n:
            raise ValueError(f"{nm} length {arr.shape[0]} != kb length {n}")

    df = pd.DataFrame({
        "kb": kb.astype(float),
        "eps": eps.astype(float),
        "G_mean": Gm.astype(float),
        "E_mean": Em.astype(float),
    }).sort_values("kb").reset_index(drop=True)

    out_path = os.path.join(base, args.out_name)
    df.to_parquet(out_path, index=False)

    print(f"[OK] wrote {out_path}")
    print(f"[OK] rows={len(df)} cols={list(df.columns)}")


if __name__ == "__main__":
    main()
