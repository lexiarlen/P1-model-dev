#!/usr/bin/env python3
import argparse
from pathlib import Path
import os
import numpy as np
import pandas as pd


def load_theta_y_sig(state_dir: Path):
    theta_path = state_dir / "theta.tsv"
    y_path = state_dir / "Y.npy"
    sig_path = state_dir / "SIG.npy"

    if not theta_path.exists():
        raise FileNotFoundError(f"Missing {theta_path}")
    if not y_path.exists():
        raise FileNotFoundError(f"Missing {y_path}")
    if not sig_path.exists():
        raise FileNotFoundError(f"Missing {sig_path}")

    THETA = np.array(np.loadtxt(theta_path), dtype=float).T   # shape (n_params, J)
    Y = np.load(y_path).astype(float)                         # shape (2, J)
    SIG = np.load(sig_path).astype(float)                     # shape (2, J)

    if Y.shape[0] != 2:
        raise ValueError(f"Expected Y shape (2, J), got {Y.shape}")
    if SIG.shape != Y.shape:
        raise ValueError(f"Expected SIG shape {Y.shape}, got {SIG.shape}")
    if THETA.shape[1] != Y.shape[1]:
        raise ValueError(
            f"Mismatch: THETA has {THETA.shape[1]} samples but Y has {Y.shape[1]}"
        )

    return THETA, Y, SIG


def classify_labels(y1, y2):
    """
    y1, y2 are 1D arrays (same length), may contain NaNs.

    0 = (NaN, NaN)  no fracture
    1 = (x,   NaN)  stable arch
    2 = (NaN, x)    no arch
    3 = (x,   x)    unstable arch
    """
    m1 = np.isnan(y1)
    m2 = np.isnan(y2)

    label = np.full(y1.shape, 3, dtype=np.int8)
    label[m1 & m2] = 0
    label[(~m1) & m2] = 1
    label[m1 & (~m2)] = 2
    return label


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--state_dir",
        type=str,
        required=True,
        help="Directory containing theta.tsv, Y.npy, SIG.npy"
    )
    ap.add_argument("--kb_idx", type=int, default=0)
    ap.add_argument("--eps_idx", type=int, default=1)
    args = ap.parse_args()

    state_dir = Path(args.state_dir)
    out_path = state_dir / "regimes.parquet"

    THETA, Y, SIG = load_theta_y_sig(state_dir)

    kb = THETA[args.kb_idx, :].astype(float).ravel()
    eps = THETA[args.eps_idx, :].astype(float).ravel()

    y1_s = Y[0, :].astype(float).ravel()
    y2_s = Y[1, :].astype(float).ravel()

    sig1 = SIG[0, :].astype(float).ravel()
    sig2 = SIG[1, :].astype(float).ravel()

    member_id = np.arange(1, len(kb) + 1, dtype=int)

    # Keep only valid parameters for plotting / analysis
    valid = np.isfinite(kb) & np.isfinite(eps) & (kb > 0) & (eps > 0)

    kb = kb[valid]
    eps = eps[valid]
    y1_s = y1_s[valid]
    y2_s = y2_s[valid]
    sig1 = sig1[valid]
    sig2 = sig2[valid]
    member_id = member_id[valid]

    y1_h = y1_s / 3600.0
    y2_h = y2_s / 3600.0
    dt_h = np.abs(y2_h - y1_h)

    label = classify_labels(y1_s, y2_s)

    # only meaningful when both peaks exist
    dt_h[label != 3] = np.nan

    df = pd.DataFrame({
        "member_id": member_id,
        "kb": kb,
        "eps": eps,
        "label": label,
        "y1_s": y1_s,
        "y2_s": y2_s,
        "y1_h": y1_h,
        "y2_h": y2_h,
        "dt_h": dt_h,
        "sig1": sig1,
        "sig2": sig2,
    })

    df.to_parquet(out_path, index=False)

    counts = df["label"].value_counts().sort_index()
    print(f"Saved parquet: {out_path}")
    print(f"Rows: {len(df)}")
    print("Label counts:")
    for k, v in counts.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()