#!/usr/bin/env python3
# get maximum strengths acheived for each shear stress simulation

import os
import sys
import glob
import re
import argparse
import numpy as np
from ovito.io import import_file


def parse_eps(path):
    base = os.path.basename(os.path.normpath(path))
    m = re.search(r"eps_([0-9.eE+-]+)", base)
    if not m:
        return None
    return float(m.group(1))


def parse_seed_from_dump(dump_file):
    base = os.path.basename(dump_file)
    m = re.search(r"_s([0-9]+)_eps", base)
    if not m:
        return None
    return int(m.group(1))


def tau_xy_from_dump(dump_file, area):
    pipeline = import_file(dump_file)
    totals = []
    nframes = pipeline.source.num_frames

    for fr in range(nframes):
        data = pipeline.compute(fr)
        sxy = data.particles["c_peratom_stress[4]"]
        totals.append(np.sum(sxy))

    totals = np.asarray(totals, float)
    return totals / area


def first_peak_before_drop_idx(stress_array, strain_array):
    """
    Return (idx, stress_at_peak, strain_at_peak).

    Peak = stress at the first negative diff ("first drop").
    If there is no drop, return the global max stress and its strain.
    """
    s = np.asarray(stress_array, float)
    e = np.asarray(strain_array, float)

    if s.size == 0:
        return np.nan, np.nan, np.nan
    if s.size == 1:
        return 0, s[0], e[0]

    diffs = np.diff(s)
    drop_idxs = np.where(diffs < 0)[0]

    if drop_idxs.size == 0:
        idx = int(np.argmax(s))
        return idx, s[idx], e[idx]

    idx = int(drop_idxs[0])
    return idx, s[idx], e[idx]


def make_qc_plot(gamma, tau, idx, outpath):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(gamma, tau, label="tau", linewidth=1.5)
    ax.axvline(gamma[idx], linestyle="--", linewidth=1.2, label="critical strain")
    ax.plot(gamma[idx], tau[idx], "o", markersize=5, label="first peak before drop")

    ax.set_xlabel("shear strain")
    ax.set_ylabel("shear stress")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(
        description="Compute critical strain, strength, and max shear stress from eps_* cases for all seeds"
    )
    ap.add_argument("--base-dir", required=True, help="Folder containing eps_* subfolders")
    ap.add_argument("--shear-rate", type=float, default=1e-7)
    ap.add_argument("--L0", type=float, default=100e3)
    ap.add_argument("--run-time", type=float, default=3600.0)
    ap.add_argument("--num-dumps", type=int, default=500)
    ap.add_argument("--qc-plots", action="store_true", help="Save tau-strain QC plot for each simulation")
    args = ap.parse_args()

    base = os.path.abspath(args.base_dir)
    area = args.L0 ** 2
    t = np.linspace(0.0, args.run_time, args.num_dumps + 1)
    gamma = args.shear_rate * t

    rows = []

    case_dirs = sorted(
        d for d in glob.glob(os.path.join(base, "eps_*"))
        if os.path.isdir(d)
    )

    if not case_dirs:
        print(f"No eps_* dirs found in {base}")
        sys.exit(0)

    qc_dir = None
    if args.qc_plots:
        qc_dir = os.path.join(base, "qc_plots")
        os.makedirs(qc_dir, exist_ok=True)

    for cdir in case_dirs:
        eps = parse_eps(cdir)
        if eps is None:
            print(f"Skipping {cdir}")
            continue

        dumps = sorted(glob.glob(os.path.join(cdir, "N*.lammps")))
        if not dumps:
            print(f"No dumps in {cdir}")
            continue

        for dump in dumps:
            seed = parse_seed_from_dump(dump)
            if seed is None:
                print(f"Could not parse seed from {dump}; skipping")
                continue

            try:
                tau = tau_xy_from_dump(dump, area)
            except Exception as e:
                print(f"Error with {dump}: {e}")
                continue

            n = min(len(tau), len(gamma))
            if n < 3:
                print(f"Too few frames in {dump}")
                continue

            tau = np.asarray(tau[:n], float)
            g = np.asarray(gamma[:n], float)

            idx, tau_peak, critical_strain = first_peak_before_drop_idx(tau, g)

            if not np.isfinite(idx):
                print(f"Could not determine peak for {dump}")
                continue

            idx = int(max(0, min(idx, n - 1)))

            max_shear_stress = tau[idx]
            strength = critical_strain * 1e9

            rows.append(
                dict(
                    eps=eps,
                    seed=seed,
                    critical_strain=critical_strain,
                    strength=strength,
                    max_shear_stress=max_shear_stress,
                )
            )

            if args.qc_plots:
                plot_name = f"eps_{eps:.6e}_s{seed}.png"
                make_qc_plot(
                    g, tau, idx,
                    os.path.join(qc_dir, plot_name)
                )

    if not rows:
        print("No results to save.")
        sys.exit(0)

    rows = sorted(rows, key=lambda r: (r["eps"], r["seed"]))

    eps_arr = np.array([r["eps"] for r in rows], float)
    seed_arr = np.array([r["seed"] for r in rows], int)
    critical_strain_arr = np.array([r["critical_strain"] for r in rows], float)
    strength_arr = np.array([r["strength"] for r in rows], float)
    max_tau_arr = np.array([r["max_shear_stress"] for r in rows], float)

    np.save(os.path.join(base, "eps.npy"), eps_arr)
    np.save(os.path.join(base, "seed.npy"), seed_arr)
    np.save(os.path.join(base, "critical_strain.npy"), critical_strain_arr)
    np.save(os.path.join(base, "strength.npy"), strength_arr)
    np.save(os.path.join(base, "max_shear_stress.npy"), max_tau_arr)

    print("Saved:")
    print(f"  {os.path.join(base, 'eps.npy')}")
    print(f"  {os.path.join(base, 'seed.npy')}")
    print(f"  {os.path.join(base, 'critical_strain.npy')}")
    print(f"  {os.path.join(base, 'strength.npy')}")
    print(f"  {os.path.join(base, 'max_shear_stress.npy')}")
    if args.qc_plots:
        print(f"  QC plots in {qc_dir}")


if __name__ == "__main__":
    main()