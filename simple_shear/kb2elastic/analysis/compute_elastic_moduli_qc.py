#!/usr/bin/env python3
"""
Compute shear (G) and elastic (E) moduli from LAMMPS shear runs,
saving ONLY the necessary arrays for downstream pandas/parquet:

Outputs under --base-dir:
  kb.npy, eps.npy, G_mean.npy, E_mean.npy
  metadata.txt (lightweight run info)

QC:
  Saves a stress-vs-strain plot for EACH dump (ensemble member) with best-fit line tau = G*gamma:
    --base-dir/kb_<kb>_eps_<eps>/qc_plots/stress_strain_<dump>.png

Assumes directory structure:
  --base-dir/kb_<kb>_eps_<eps>/N*.lammps
"""

import os, sys, glob, re
import numpy as np
import argparse
from ovito.io import import_file

# QC plotting
import matplotlib
matplotlib.use("Agg")  # safe on HPC
import matplotlib.pyplot as plt


# ---------- helpers ----------
def moving_average(x, w):
    x = np.asarray(x, float)
    if w <= 1:
        return x
    k = np.ones(w, dtype=float) / w
    return np.convolve(x, k, mode="same")


def first_sustained_decrease(y, tol_frac=0.005, min_grow_len=10, lookahead=1):
    y = np.asarray(y, float)
    if len(y) < min_grow_len + lookahead + 2:
        return len(y) - 1

    finite = np.isfinite(y)
    if not np.any(finite):
        return len(y) - 1

    ymax = np.nanmax(np.abs(y[finite]))
    tol = tol_frac * (ymax if ymax > 0 else 1.0)

    for i in range(min_grow_len, len(y) - lookahead - 1):
        if np.all(y[i + 1 : i + 1 + lookahead] < y[i] - tol):
            return i

    return len(y) - 1


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
    m = re.search(r"kb_([0-9.eE+-]+)_eps_([0-9.eE+-]+)", path)
    if not m:
        return None, None
    return float(m.group(1)), float(m.group(2))



def tau_xy_from_dump(dump_file, area, phi):
    pipeline = import_file(dump_file)
    totals = []
    nframes = pipeline.source.num_frames

    for fr in range(nframes):
        data = pipeline.compute(fr)
        sxy = data.particles["c_peratom_stress[4]"]
        totals.append(np.sum(sxy))

    totals = np.asarray(totals, float)
    return phi * totals / area


def save_qc_plot(out_png, gamma, tau, idx, G, title=None):
    """
    Save a QC plot:
      - full curve tau(gamma)
      - highlighted fit segment (0..idx)
      - fit line tau = G*gamma on fit segment
    """
    gamma = np.asarray(gamma, float)
    tau = np.asarray(tau, float)

    # convert for display only
    tau_kpa = tau / 1e3

    # fit window
    g_fit = gamma[: idx + 1]
    tau_fit = tau[: idx + 1]
    tau_fitline_kpa = (G * g_fit) / 1e3

    fig = plt.figure(figsize=(6, 4))

    # full curve (thin)
    plt.plot(gamma, tau_kpa, linewidth=1.0)

    # fit window (thicker)
    plt.plot(g_fit, tau_fit / 1e3, linewidth=2.0)

    # fit line
    plt.plot(g_fit, tau_fitline_kpa, linestyle="--", linewidth=2.0)

    plt.xlabel("shear strain, $\\gamma$ [-]")
    plt.ylabel("shear stress, $\\tau$ [kPa]")

    if title:
        plt.title(title, fontsize=10)

    # annotate G in plot (Pa) and also kPa per strain
    plt.text(
        0.02, 0.98,
        f"G = {G:.3g} Pa\n(G/1e3 = {G/1e3:.3g} kPa)",
        transform=plt.gca().transAxes,
        va="top",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
    )

    plt.grid(linestyle="--", linewidth=0.5, color="0.7", alpha=0.7)

    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    plt.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ---------- main ----------
def main():
    ap = argparse.ArgumentParser(description="Compute G, E and save minimal .npy arrays + QC plots")
    ap.add_argument("--base-dir", required=True, help="Folder containing kb_*_eps_* subfolders")
    ap.add_argument("--nu", type=float, default=0.3)
    ap.add_argument("--phi", type=float, default=0.71)
    ap.add_argument("--shear-rate", type=float, default=1e-7)
    ap.add_argument("--L0", type=float, default=100e3)
    ap.add_argument("--run-time", type=float, default=3600.0)
    ap.add_argument("--num-dumps", type=int, default=500)
    ap.add_argument("--smooth-window", type=int, default=2)
    ap.add_argument("--tol-frac", type=float, default=0.005)
    ap.add_argument("--min-grow-len", type=int, default=10)
    ap.add_argument("--lookahead", type=int, default=3)
    ap.add_argument("--qc-plots", action="store_true", help="Save per-dump stress-strain QC plots")
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

        G_list = []
        E_list = []

        qc_dir = os.path.join(cdir, "qc_plots")

        for dump in dumps:
            try:
                tau = tau_xy_from_dump(dump, area, args.phi)
            except Exception as e:
                print(f"[ERROR] {dump}: {e}")
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

            if args.qc_plots:
                dump_base = os.path.splitext(os.path.basename(dump))[0]
                out_png = os.path.join(qc_dir, f"stress_strain_{dump_base}.png")
                title = f"kb={kb:.3g}, eps={eps:.3g} | {dump_base}"
                try:
                    save_qc_plot(out_png, g, tau, idx, G, title=title)
                except Exception as e:
                    print(f"[WARN] QC plot failed for {dump}: {e}")

        if not G_list:
            continue

        rows.append(dict(
            kb=kb,
            eps=eps,
            G_mean=float(np.nanmean(G_list)),
            E_mean=float(np.nanmean(E_list)),
        ))

        if args.qc_plots:
            print(f"[OK] QC plots written under {qc_dir}")

    if not rows:
        print("[WARN] No results to save.")
        sys.exit(0)

    kb_arr = np.array([r["kb"] for r in rows], float)
    eps_arr = np.array([r["eps"] for r in rows], float)
    Gm_arr = np.array([r["G_mean"] for r in rows], float)
    Em_arr = np.array([r["E_mean"] for r in rows], float)

    np.save(os.path.join(base, "kb.npy"), kb_arr)
    np.save(os.path.join(base, "eps.npy"), eps_arr)
    np.save(os.path.join(base, "G_mean.npy"), Gm_arr)
    np.save(os.path.join(base, "E_mean.npy"), Em_arr)

    # lightweight metadata for reproducibility
    meta_path = os.path.join(base, "metadata.txt")
    with open(meta_path, "w") as f:
        f.write("compute_elastic_moduli (minimal + qc)\n")
        f.write(f"base_dir: {base}\n")
        f.write(f"nu: {args.nu}\n")
        f.write(f"phi: {args.phi}\n")
        f.write(f"shear_rate: {args.shear_rate}\n")
        f.write(f"L0: {args.L0}\n")
        f.write(f"run_time: {args.run_time}\n")
        f.write(f"num_dumps: {args.num_dumps}\n")
        f.write(f"smooth_window: {args.smooth_window}\n")
        f.write(f"tol_frac: {args.tol_frac}\n")
        f.write(f"min_grow_len: {args.min_grow_len}\n")
        f.write(f"lookahead: {args.lookahead}\n")
        f.write(f"qc_plots: {bool(args.qc_plots)}\n")
        f.write(f"n_cases_saved: {len(rows)}\n")

    print(f"[OK] wrote kb.npy, eps.npy, G_mean.npy, E_mean.npy under {base}")
    print(f"[OK] wrote {meta_path}")


if __name__ == "__main__":
    main()
