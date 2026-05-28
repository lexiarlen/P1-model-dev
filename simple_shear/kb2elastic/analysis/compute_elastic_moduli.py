#!/usr/bin/env python3
# extracts per atom stress from lammps dump and computes elastic properties by 
# taking the slope through the origin of the shear stress v shear strain plot 
# up to a few dumps before the ice initially fails

import os, sys, glob, re
import numpy as np
import argparse
from ovito.io import import_file

def moving_average(x, w):
    x = np.asarray(x, float)
    if w <= 1:
        return x
    k = np.ones(w, dtype=float) / w
    return np.convolve(x, k, mode="same")


def first_sustained_decrease(y, tol_frac=0.005, min_grow_len=10, lookahead=1, buffer = 10):
    # get failure from stress data; subtract buffer safely to make sure 
    # getting slope from purely elastic regime
    y = np.asarray(y, float)
    if len(y) < min_grow_len + lookahead + 2:
        return len(y) - 1 - buffer

    finite = np.isfinite(y)
    if not np.any(finite):
        return len(y) - 1 - buffer 

    ymax = np.nanmax(np.abs(y[finite]))
    tol = tol_frac * (ymax if ymax > 0 else 1.0)

    for i in range(min_grow_len, len(y) - lookahead - 1):
        if np.all(y[i + 1 : i + 1 + lookahead] < y[i] - tol):
            if i - buffer >= 0:
                return i - buffer
            else: 
                return i

    return len(y) - 1 - buffer

def slope_through_origin(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if not np.any(m):
        return np.nan
    x = x[m]
    y = y[m]
    den = np.dot(x, x)
    return np.dot(x, y) / den if den > 0 else np.nan


def parse_kb_eps(path):
    # file parsing
    m = re.search(r"kb_([0-9.eE+-]+)_eps_([0-9.eE+-]+)", path)
    if not m:
        return None, None
    return float(m.group(1)), float(m.group(2))

def tau_xy_from_dump(dump_file, area):
    # get shear stress from dumped data
    pipeline = import_file(dump_file)
    totals = []
    nframes = pipeline.source.num_frames

    for fr in range(nframes):
        data = pipeline.compute(fr)
        sxy = data.particles["c_peratom_stress[4]"]
        totals.append(np.sum(sxy))

    totals = np.asarray(totals, float)

    return totals / area


def main():
    ap = argparse.ArgumentParser(description="Compute G, E and save minimal .npy arrays")
    ap.add_argument("--base-dir", required=True, help="Folder containing kb_*_eps_* subfolders")
    ap.add_argument("--nu", type=float, default=0.3333)
    ap.add_argument("--shear-rate", type=float, default=1e-7)
    ap.add_argument("--L0", type=float, default=100e3)
    ap.add_argument("--run-time", type=float, default=3600.0)
    ap.add_argument("--num-dumps", type=int, default=500)
    ap.add_argument("--smooth-window", type=int, default=7)
    ap.add_argument("--tol-frac", type=float, default=0.005)
    ap.add_argument("--min-grow-len", type=int, default=10)
    ap.add_argument("--lookahead", type=int, default=3)
    args = ap.parse_args()

    base = os.path.abspath(args.base_dir)
    area = args.L0 ** 2
    t = np.linspace(0.0, args.run_time, args.num_dumps + 1)
    gamma = args.shear_rate * t

    rows = []
    case_dirs = sorted(
        d for d in glob.glob(os.path.join(base, "kb_*_eps_*"))
        if os.path.isdir(d)
    )
    if not case_dirs:
        print(f"No kb_*_eps_* dirs in {base}")
        sys.exit(0)

    for cdir in case_dirs:
        kb, eps = parse_kb_eps(cdir)
        if kb is None:
            print(f"Skipping {cdir}")
            continue

        dumps = sorted(glob.glob(os.path.join(cdir, "N*.lammps")))
        if not dumps:
            print(f"No dumps in {cdir}")
            continue

        G_list = []
        E_list = []

        for dump in dumps:
            try:
                tau = tau_xy_from_dump(dump, area)
            except Exception as e:
                print(f"Erron with {dump}: {e}")
                continue

            n = min(len(tau), len(gamma))
            if n < 3:
                continue

            tau = np.asarray(tau[:n], float)
            g = np.asarray(gamma[:n], float)

            tau_s = moving_average(tau, args.smooth_window)

            idx = first_sustained_decrease(
                tau_s,
                tol_frac=args.tol_frac,
                min_grow_len=args.min_grow_len,
                lookahead=args.lookahead,
            )
            if idx < 2:
                idx = max(2, len(tau_s) // 10)

            x = g[: idx + 1]
            y = tau[: idx + 1]
            G = slope_through_origin(x, y)
            E = 2.0 * G * (1.0 + args.nu)

            G_list.append(G)
            E_list.append(E)

        if not G_list:
            continue

        rows.append(dict(
            kb=kb,
            eps=eps,
            G_mean=float(np.nanmean(G_list)),
            E_mean=float(np.nanmean(E_list)),
        ))

    if not rows:
        print("No results to save.")
        sys.exit(0)

    kb_arr = np.array([r["kb"] for r in rows], float)
    eps_arr = np.array([r["eps"] for r in rows], float)
    Gm_arr = np.array([r["G_mean"] for r in rows], float)
    Em_arr = np.array([r["E_mean"] for r in rows], float)

    np.save(os.path.join(base, "kb.npy"), kb_arr)
    np.save(os.path.join(base, "eps.npy"), eps_arr)
    np.save(os.path.join(base, "G_mean.npy"), Gm_arr)
    np.save(os.path.join(base, "E_mean.npy"), Em_arr)

if __name__ == "__main__":
    main()
