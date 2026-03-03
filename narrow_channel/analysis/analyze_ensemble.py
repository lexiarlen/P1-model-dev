#!/usr/bin/env python3
import argparse
import re
import os
import glob
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks
from ovito.io import import_file
import matplotlib.pyplot as plt


# ----------------------------- helpers -----------------------------

def get_ts_from_data(fname, dt, y_cutoff=120e3):
    """
    Load LAMMPS dump and compute dcodt and avg vy in the upper half.
    Returns:
      dcodt: |d(avg_co)/dt_frame|  (units per frame; consistent with your thresholds)
      avg_vy: mean vy in upper half (m/s)
      time: seconds
    Notes:
      - Always appends per-frame values to keep arrays aligned.
      - Missing data is recorded as np.nan.
    """
    pipeline = import_file(str(fname))
    nframes = pipeline.source.num_frames

    timesteps = np.empty(nframes, dtype=float)
    avg_co = np.empty(nframes, dtype=float)
    avg_vy = np.empty(nframes, dtype=float)

    for frame in range(nframes):
        data = pipeline.compute(frame)

        ts = data.attributes.get("Timestep", frame)
        timesteps[frame] = float(ts)

        # default to NaN each frame; fill if we can compute
        avg_co[frame] = np.nan
        avg_vy[frame] = np.nan

        # mean coordination number (active particles only)
        if ('c_nbond' in data.particles.keys()) and ('v_ingroup' in data.particles.keys()):
            co = data.particles['c_nbond']
            ingroup = data.particles['v_ingroup']
            mask = (ingroup == 1)
            filtered = co[mask]
            if len(filtered) > 0:
                avg_co[frame] = float(np.mean(filtered))

        # velocity in upper half (active particles only)
        if ('Velocity.Y' in data.particles.keys()) and ('v_ingroup' in data.particles.keys()) and ('Position.Y' in data.particles.keys()):
            vy = data.particles['Velocity.Y']
            ypos = data.particles['Position.Y']
            ingroup = data.particles['v_ingroup']
            mask = (ingroup == 1) & (ypos >= y_cutoff)
            filtered = vy[mask]
            if len(filtered) > 0:
                avg_vy[frame] = float(np.mean(filtered))

    # Build a NaN-safe avg_co for gradient.
    # If all NaN (pathological), fallback to zeros -> dcodt zeros.
    if np.all(np.isnan(avg_co)):
        avg_co_filled = np.zeros_like(avg_co)
    else:
        # Fill NaNs using linear interpolation over frame index
        x = np.arange(nframes)
        ok = np.isfinite(avg_co)
        avg_co_filled = avg_co.copy()
        avg_co_filled[~ok] = np.interp(x[~ok], x[ok], avg_co[ok])

    dcodt = np.abs(np.gradient(avg_co_filled))
    time = timesteps * dt
    return dcodt, avg_vy, time


def compute_y(dcodt, vy, time, threshold=0.002, second_ratio=0.15, none_val=np.nan, epsilon=0.01):
    """
    Returns array([y1, y2]) in seconds, with np.nan representing "none".
    """
    # NaN-safe peak finding: treat NaNs as 0 so they don't create peaks.
    dcodt_safe = np.nan_to_num(dcodt, nan=0.0, posinf=0.0, neginf=0.0)

    peaks, _ = find_peaks(dcodt_safe)
    heights = dcodt_safe[peaks]

    # no peaks
    if len(peaks) == 0:
        if np.abs(np.nanmean(vy)) > epsilon:
            return np.array([0.0, none_val], dtype=float)  # fail immediately case (kept as you had it)
        return np.array([none_val, none_val], dtype=float)  # no fracturing

    # sort peaks by height (tallest first)
    sort_idx = np.argsort(heights)[::-1]
    peaks_sorted = peaks[sort_idx]
    heights_sorted = heights[sort_idx]

    # tallest peak check
    if heights_sorted[0] < threshold:
        if np.abs(np.nanmean(vy)) > epsilon:
            return np.array([0.0, none_val], dtype=float)  # fail immediately case
        return np.array([none_val, none_val], dtype=float)

    peak_indices = [int(peaks_sorted[0])]

    # Search for a second peak among the next 4 highest
    if len(peaks_sorted) > 1:
        max_check = min(5, len(peaks_sorted))
        p1 = int(peaks_sorted[0])

        for j in range(1, max_check):
            p2 = int(peaks_sorted[j])
            h2 = heights_sorted[j]

            # amplitude criterion
            if h2 < second_ratio * heights_sorted[0]:
                continue

            # valley criterion between p1 and p2
            lo, hi = (p1, p2) if p1 < p2 else (p2, p1)
            if hi - lo <= 1:
                continue

            valley = np.nanmin(dcodt_safe[lo + 1:hi])
            has_valley = valley < threshold
            if not has_valley:
                continue

            # velocity logic
            pre = np.nanmean(vy[:p2]) if p2 > 0 else np.nanmean(vy)
            post = np.nanmean(vy[p2:]) if p2 < len(vy) else np.nanmean(vy)

            if (np.abs(pre) < epsilon) and (np.abs(post) > epsilon):
                peak_indices.append(p2)
                break

    # Build y
    if len(peak_indices) == 2:
        y = time[np.sort(np.array(peak_indices, dtype=int))]
        return y.astype(float)

    # single peak classification
    y = np.array([none_val, none_val], dtype=float)
    idx = int(peak_indices[0])

    # stable arch
    if (idx > 1) and (idx < len(vy) - 2):
        local_mean = np.nanmean(vy[idx - 2:idx + 2])
        if np.abs(local_mean) < epsilon:
            y[0] = float(time[idx])
            y[1] = none_val

    # no arch
    if np.abs(np.nanmean(vy)) > epsilon:
        y[0] = none_val
        y[1] = float(time[idx])

    # stable arch but fail quick
    if (idx <= 1) and (np.abs(np.nanmean(vy[:2])) < epsilon):
        y[0] = float(time[idx])
        y[1] = none_val

    return y


def get_u_vec(time_array, u_max, max_time):
    """NaN-propagating vectorized wind speed profile."""
    return u_max * time_array / max_time


def get_wind_stress_vec(u_array, rho_a, C_a):
    """NaN-propagating vectorized wind stress."""
    return rho_a * C_a * np.abs(u_array) * u_array


# ----------------------------- main -----------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--J", type=int, required=True)
    parser.add_argument("--base_dir", type=str, required=True)

    # optional knobs
    parser.add_argument("--threshold", type=float, default=0.002)
    parser.add_argument("--second_ratio", type=float, default=0.15)
    parser.add_argument("--epsilon", type=float, default=0.01)
    parser.add_argument("--y_cutoff", type=float, default=120e3)

    args = parser.parse_args()

    J = args.J
    base_dir = Path(args.base_dir)

    # hard coding simulation setup
    u_max = 20
    rho_a = 1.3
    C_a = 0.0012
    T_max = 26
    max_time = T_max * 3600.0  # [s]

    radius = 400
    rho_i = 920
    radius_small = radius * 3 / 5
    m_ice = np.pi * radius_small**2 * rho_i

    # Paths
    lammps_dir = base_dir / "lammps_output"
    state_dir = base_dir / "state"
    state_dir.mkdir(parents=True, exist_ok=True)

    # Load THETA and extract kb for each ensemble member
    theta_path = state_dir / "theta.tsv"
    if not theta_path.exists():
        raise FileNotFoundError(f"Missing {theta_path}. Expected theta.tsv under {state_dir}")

    THETA = np.array(np.loadtxt(theta_path), dtype=float).T  # shape (n_params, J)
    if THETA.shape[1] != J:
        raise ValueError(f"theta.tsv has J={THETA.shape[1]} columns but you passed --J {J}")

    kb_values = THETA[0, :]  # k_b row

    # analysis
    pattern = re.compile(r'^dump(\d+)\.lammps$')

    Y = np.full((2, J), np.nan, dtype=float)
    SIG = np.full((2, J), np.nan, dtype=float)

    files = sorted(glob.glob(str(lammps_dir / "*.lammps")))
    print(f"Found {len(files)} dump files under {lammps_dir}")

    for dumpfile in files:
        dumpfile = Path(dumpfile)
        fname = dumpfile.name
        m = pattern.match(fname)
        if not m:
            print(f"Skipped: {dumpfile}")
            continue

        ensemble_member_id = int(m.group(1))
        if not (1 <= ensemble_member_id <= J):
            print(f"Skipping {dumpfile}: member id {ensemble_member_id} out of 1...{J}")
            continue

        print(f"File: {dumpfile} → ID: {ensemble_member_id}")

        kb_local = float(kb_values[ensemble_member_id - 1])
        dt_local = 0.1 * np.sqrt(m_ice / kb_local)

        dcodt, vy, time = get_ts_from_data(dumpfile, dt_local, y_cutoff=args.y_cutoff)
        y = compute_y(
            dcodt, vy, time,
            threshold=args.threshold,
            second_ratio=args.second_ratio,
            none_val=np.nan,
            epsilon=args.epsilon,
        )

        # -------------------- figure generation --------------------
        fig, ax = plt.subplots(figsize=(10, 4))
        line1, = ax.plot(time / 3600.0, dcodt, lw=4, color='gray', label='damage rate')
        ax.set_xlabel("time (hrs)", size='large')
        ax.tick_params(labelsize='large')
        ax.set_ylabel(r"$\frac{d\text{Co#}}{dt}$", size='xx-large', color='gray')
        ax.set_ylim(-0.01, 0.35)
        ax.grid(True)

        ax2 = ax.twinx()
        line2, = ax2.plot(time / 3600.0, vy, lw=4, color='r', label=r'$v_{y}^{up}$')
        ax2.set_ylabel(r"mean $v_y$ (m/s) in upper domain", size='large', color='tab:blue')
        ax2.tick_params(labelsize='large')
        ax2.set_ylim(-0.3, 0.1)

        lines = [line1, line2]

        if np.isfinite(y[0]):
            line3 = ax.axvline(y[0] / 3600.0, color='k', ls='-.',
                               label=rf'$y_1 = {(y[0]/3600.0):.1f}$ h')
            lines.append(line3)

        if np.isfinite(y[1]):
            line4 = ax.axvline(y[1] / 3600.0, color='k', ls='--',
                               label=rf'$y_2 = {(y[1]/3600.0):.1f}$ h')
            lines.append(line4)

        ax.legend(lines, [ln.get_label() for ln in lines], loc='best')
        ax.set_xlim(0, 36)

        fig_path = state_dir / f"diagnostic_{ensemble_member_id:03d}.png"
        plt.savefig(fig_path, bbox_inches='tight', dpi=300)
        plt.close(fig)
        # ------------------------------------------------------------

        u = get_u_vec(y, u_max, max_time)
        sig = get_wind_stress_vec(u, rho_a, C_a)

        SIG[:, ensemble_member_id - 1] = sig
        Y[:, ensemble_member_id - 1] = y

    np.save(state_dir / 'Y.npy', Y)
    np.save(state_dir / 'SIG.npy', SIG)

    print(f"Saved Y.npy and SIG.npy under {state_dir}")
    print(f"Saved diagnostic figures under {state_dir}")


if __name__ == "__main__":
    main()
