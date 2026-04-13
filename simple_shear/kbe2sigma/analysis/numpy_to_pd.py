#!/usr/bin/env python3
# saves .npy files from compute_max_strengths.py to one dataset

import os
import argparse
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser(
        description="Build pandas DataFrame from minimal .npy arrays"
    )
    ap.add_argument("--base-dir", required=True, help="Folder containing the .npy files")
    ap.add_argument(
        "--out-name",
        default="strength_results.parquet",
        help="Output parquet filename"
    )
    args = ap.parse_args()

    base = os.path.abspath(args.base_dir)

    # required arrays
    eps = np.load(os.path.join(base, "eps.npy"))
    critical_strain = np.load(os.path.join(base, "critical_strain.npy"))
    strength = np.load(os.path.join(base, "strength.npy"))
    max_shear_stress = np.load(os.path.join(base, "max_shear_stress.npy"))

    n = eps.shape[0]
    for arr, nm in [
        (critical_strain, "critical_strain"),
        (strength, "strength"),
        (max_shear_stress, "max_shear_stress"),
    ]:
        if arr.shape[0] != n:
            raise ValueError(f"{nm} length {arr.shape[0]} != eps length {n}")

    df = pd.DataFrame({
        "eps": eps.astype(float),
        "critical_strain": critical_strain.astype(float),
        "strength": strength.astype(float),
        "max_shear_stress": max_shear_stress.astype(float),
    }).sort_values("eps").reset_index(drop=True)

    out_path = os.path.join(base, args.out_name)
    df.to_parquet(out_path, index=False)

    print(f"[OK] wrote {out_path}")
    print(f"[OK] rows={len(df)} cols={list(df.columns)}")


if __name__ == "__main__":
    main()