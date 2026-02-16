# for each kb, save a (3 x 10) numpy array with kb, G, E as the rows
#!/usr/bin/env python3
"""
Compute shear (G) and elastic (E) moduli from LAMMPS shear runs,
then save results as NumPy arrays (.npy) so we can build a NetCDF later
from a different (clean) Python environment.

Outputs under --base-dir:
  kb.npy, eps.npy, product_kb_eps.npy,
  G_mean.npy, G_std.npy, E_mean.npy, E_std.npy,
  peak_gamma_mean.npy, peak_tau_mean.npy, n_files.npy
  metadata.txt  (simple human-readable metadata)

"""

import os, sys, glob, re
import numpy as np
import argparse
from ovito.io import import_file

# ---------- helpers ----------
def moving_average(x, w):
    if w <= 1:
        return x
    k = np.ones(w, dtype=float) / w
    return np.convolve(x, k, mode='same')

def first_sustained_decrease(y, tol_frac=0.005, min_grow_len=10, lookahead=3):
    if len(y) < min_grow_len + lookahead + 2:
        return len(y) - 1
    ymax = np.nanmax(np.abs(y)) if np.any(np.isfinite(y)) else 0.0
    tol = tol_frac * (ymax if ymax > 0 else 1.0)
    for i in range(min_grow_len, len(y) - lookahead - 1):
        if np.all(y[i+1:i+1+lookahead] < y[i] - tol):
            return i
    return len(y) - 1

def slope_through_origin(x, y):
    den = np.dot(x, x)
    return np.dot(x, y) / den if den > 0 else np.nan

def parse_kb_eps(path):
    m = re.search(r"kb_([0-9.eE+-]+)_eps_([0-9.eE+-]+)", path)
    if not m: return None, None
    return float(m.group(1)), float(m.group(2))

def tau_xy_from_dump(dump_file, area, phi):
    pipeline = import_file(dump_file)
    totals = []
    nframes = pipeline.source.num_frames
    for fr in range(nframes):
        data = pipeline.compute(fr)
        if 'c_peratom_stress[4]' not in data.particles.keys():
            raise RuntimeError(f"c_peratom_stress[4] missing in {dump_file}")
        sxy = data.particles['c_peratom_stress[4]']
        totals.append(np.sum(sxy))
    totals = np.asarray(totals, float)
    return phi * totals / area

# ---------- main ----------
def main():
    ap = argparse.ArgumentParser(description="Compute G, E and save as .npy arrays")
    ap.add_argument("--base-dir", required=True, help="Folder containing kb_*_eps_* subfolders")
    ap.add_argument("--nu", type=float, default=0.3)
    ap.add_argument("--phi", type=float, default=0.71)
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
    area = args.L0**2
    t = np.linspace(0.0, args.run_time, args.num_dumps + 1)
    gamma = args.shear_rate * t

    rows = []
    case_dirs = sorted(d for d in glob.glob(os.path.join(base, "kb_*_eps_*")) if os.path.isdir(d))
    if not case_dirs:
        print(f"[WARN] No kb_*_eps_* dirs in {base}")
        sys.exit(0)

    for cdir in case_dirs:
        kb, eps = parse_kb_eps(cdir)
        if kb is None:
            print(f"[WARN] skip {cdir}")
            continue
        dumps = sorted(glob.glob(os.path.join(cdir, "N*.lammps")))
        if not dumps:
            print(f"[WARN] no dumps in {cdir}")
            continue

        G_list, E_list, pg_list, pt_list = [], [], [], []
        for dump in dumps:
            try:
                tau = tau_xy_from_dump(dump, area, args.phi)
            except Exception as e:
                print(f"[ERROR] {dump}: {e}")
                continue

            n = min(len(tau), len(gamma))
            tau = tau[:n]
            g = gamma[:n]

            tau_s = moving_average(tau, args.smooth_window)
            idx = first_sustained_decrease(
                tau_s, tol_frac=args.tol_frac,
                min_grow_len=args.min_grow_len,
                lookahead=args.lookahead,
            )
            if idx < 2:
                idx = max(2, len(tau_s)//10)

            x = g[:idx+1]
            y = tau[:idx+1]
            G = slope_through_origin(x, y)
            E = 2.0 * G * (1.0 + args.nu)

            G_list.append(G)
            E_list.append(E)
            pg_list.append(g[idx])
            pt_list.append(tau[idx])

        if not G_list:
            continue

        rows.append(dict(
            kb=kb, eps=eps, product=kb*eps,
            G_mean=np.nanmean(G_list), G_std=np.nanstd(G_list),
            E_mean=np.nanmean(E_list), E_std=np.nanstd(E_list),
            peak_gamma_mean=np.nanmean(pg_list), peak_tau_mean=np.nanmean(pt_list),
            n_files=len(G_list)
        ))

    if not rows:
        print("[WARN] No results to save.")
        sys.exit(0)

    # pack into arrays and save
    kb_arr   = np.array([r["kb"] for r in rows], float)
    eps_arr  = np.array([r["eps"] for r in rows], float)
    prod_arr = np.array([r["product"] for r in rows], float)
    Gm_arr   = np.array([r["G_mean"] for r in rows], float)
    Gs_arr   = np.array([r["G_std"] for r in rows], float)
    Em_arr   = np.array([r["E_mean"] for r in rows], float)
    Es_arr   = np.array([r["E_std"] for r in rows], float)
    pg_arr   = np.array([r["peak_gamma_mean"] for r in rows], float)
    pt_arr   = np.array([r["peak_tau_mean"] for r in rows], float)
    nfiles   = np.array([r["n_files"] for r in rows], int)

    # ensure deterministic case order (already sorted by folder name)
    np.save(os.path.join(base, "kb.npy"), kb_arr)
    np.save(os.path.join(base, "eps.npy"), eps_arr)
    np.save(os.path.join(base, "product_kb_eps.npy"), prod_arr)
    np.save(os.path.join(base, "G_mean.npy"), Gm_arr)
    np.save(os.path.join(base, "G_std.npy"), Gs_arr)
    np.save(os.path.join(base, "E_mean.npy"), Em_arr)
    np.save(os.path.join(base, "E_std.npy"), Es_arr)
    np.save(os.path.join(base, "peak_gamma_mean.npy"), pg_arr)
    np.save(os.path.join(base, "peak_tau_mean.npy"), pt_arr)
    np.save(os.path.join(base, "n_files.npy"), nfiles)

    # tiny metadata file for human inspection
    with open(os.path.join(base, "metadata.txt"), "w") as f:
        f.write(
            "Saved arrays:\n"
            "  kb.npy, eps.npy, product_kb_eps.npy,\n"
            "  G_mean.npy, G_std.npy, E_mean.npy, E_std.npy,\n"
            "  peak_gamma_mean.npy, peak_tau_mean.npy, n_files.npy\n\n"
            f"nu={args.nu}, phi={args.phi}, shear_rate={args.shear_rate}, "
            f"L0={args.L0}, run_time={args.run_time}, num_dumps={args.num_dumps}\n"
            f"smooth_window={args.smooth_window}, tol_frac={args.tol_frac}, "
            f"min_grow_len={args.min_grow_len}, lookahead={args.lookahead}\n"
        )

    print(f"[OK] .npy arrays written under {base}")

if __name__ == "__main__":
    main()
