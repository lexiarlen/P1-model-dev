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

    Convention:
      - If no qualifying particles exist for a frame, store 0.0.
      - No NaN handling here; NaN is reserved for output labels in compute_y().
    """
    pipeline = import_file(str(fname))
    nframes = pipeline.source.num_frames

    timesteps = np.empty(nframes, dtype=float)
    avg_co = np.zeros(nframes, dtype=float)
    avg_vy = np.zeros(nframes, dtype=float)

    for frame in range(nframes):
        data = pipeline.compute(frame)

        ts = data.attributes.get("Timestep", frame)
        timesteps[frame] = float(ts)

        # mean coordination number (active particles only)
        if ('c_nbond' in data.particles.keys()) and ('v_ingroup' in data.particles.keys()):
            co = data.particles['c_nbond']
            ingroup = data.particles['v_ingroup']
            mask = (ingroup == 1)
            filtered = co[mask]
            if len(filtered) > 0:
                avg_co[frame] = float(np.mean(filtered))

        # velocity in upper half (active particles only)
        if ('Velocity' in data.particles.keys()) and ('v_ingroup' in data.particles.keys()) and ('Position' in data.particles.keys()):
            vy = np.asarray(data.particles['Velocity'])[:, 1]
            ypos = np.asarray(data.particles['Position'])[:, 1]
            ingroup = data.particles['v_ingroup']
            mask = (ingroup == 1) & (ypos >= y_cutoff)
            filtered = vy[mask]
            if len(filtered) > 0:
                avg_vy[frame] = float(np.mean(filtered))

    dcodt = np.abs(np.gradient(avg_co, 36/500)) # fracture per hour - 36 hours, 500 dumps
    time = timesteps * dt
    return dcodt, avg_vy, time


def compute_y(dcodt, vy, time, threshold=0.05, second_ratio=0.1, none_val=np.nan, epsilon=0.01):
    """
    Returns array([y1, y2]) in seconds, with np.nan representing "none".

    Interpretation:
      y1 = stable arch formation time
      y2 = failure time

    Assumes dcodt and vy are finite arrays (no NaN preprocessing needed).
    """
    peaks, _ = find_peaks(dcodt)
    heights = dcodt[peaks]

    mean_vy = np.mean(vy)

    # no peaks
    if len(peaks) == 0:
        if np.abs(mean_vy) > epsilon:
            return np.array([0.0, none_val], dtype=float)   # fail immediately
        return np.array([none_val, none_val], dtype=float)  # no fracturing

    # sort peaks by height (tallest first)
    sort_idx = np.argsort(heights)[::-1]
    peaks_sorted = peaks[sort_idx]
    heights_sorted = heights[sort_idx]

    # tallest peak check
    if heights_sorted[0] < threshold:
        if np.abs(mean_vy) > epsilon:
            return np.array([0.0, none_val], dtype=float)   # fail immediately
        return np.array([none_val, none_val], dtype=float)  # no fracturing

    p1 = int(peaks_sorted[0])

    # search for second peak among next 4 highest
    max_check = min(5, len(peaks_sorted))
    for j in range(1, max_check):
        p2 = int(peaks_sorted[j])
        h2 = heights_sorted[j]

        # amplitude criterion
        if h2 < second_ratio * heights_sorted[0]:
            continue

        # valley criterion between p1 and p2
        lo, hi = sorted((p1, p2))
        if hi - lo <= 1:
            continue

        valley = np.min(dcodt[lo + 1:hi])
        if valley >= threshold:
            continue

        # velocity logic
        pre = np.mean(vy[:p2]) if p2 > 0 else mean_vy
        post = np.mean(vy[p2:]) if p2 < len(vy) else mean_vy

        if (np.abs(pre) < epsilon) and (np.abs(post) > epsilon):
            y = time[np.sort(np.array([p1, p2], dtype=int))]
            return y.astype(float)

    # single peak classification
    y = np.array([none_val, none_val], dtype=float)
    idx = p1

    # no arch
    if np.abs(mean_vy) > epsilon:
        y[1] = float(time[idx])
        return y

    # stable arch
    if (idx > 1) and (idx < len(vy) - 2):
        local_mean = np.mean(vy[idx - 2:idx + 2])
        if np.abs(local_mean) < epsilon:
            y[0] = float(time[idx])
            return y

    # stable arch but fail quick
    if (idx <= 1) and (np.abs(np.mean(vy[:2])) < epsilon):
        y[0] = float(time[idx])

    return y


#### testing

def compute_y2(dcodt, vy, time, threshold=0.02, second_ratio=0.1, none_val=np.nan, epsilon=0.01, local_window=2, min_peak_separation_h=0.75):
    """
    Returns array([y1, y2]) in seconds.

    y1 = stable arch formation time
    y2 = failure time

    none_val returned if arch never forms/collapses

    Same basic logic as compute_y, but safer for:
      - failure peak taller than arch peak
      - failure near the end of the time series
      - velocity logic needing to be based on the later peak
      - rejecting two-peak arch/failure cases if peaks are too close in time
    """

    dcodt = np.asarray(dcodt, dtype=float)
    vy = np.asarray(vy, dtype=float)
    time = np.asarray(time, dtype=float)

    if not (len(dcodt) == len(vy) == len(time)):
        raise ValueError("dcodt, vy, and time must have the same length.")

    n = len(time)
    y = np.array([none_val, none_val], dtype=float)

    if n == 0:
        return y

    peaks, _ = find_peaks(dcodt)
    heights = dcodt[peaks]

    min_peak_separation_s = min_peak_separation_h * 3600.0

    # ------------------------------------------------------------
    # helper functions
    # ------------------------------------------------------------
    def mean_window(x, center, before=local_window, after=local_window):
        lo = max(0, center - before)
        hi = min(len(x), center + after + 1)

        if hi <= lo:
            return np.nan

        return np.mean(x[lo:hi])

    def has_valley(p_a, p_b):
        lo, hi = sorted((int(p_a), int(p_b)))

        if hi - lo <= 1:
            return False

        return np.min(dcodt[lo + 1:hi]) < threshold

    def peaks_far_enough_apart(p_a, p_b):
        lo, hi = sorted((int(p_a), int(p_b)))
        return (time[hi] - time[lo]) >= min_peak_separation_s

    def post_mean(idx):
        # use all data after idx, edge safe
        if idx < n - 1:
            return np.mean(vy[idx:])
        return vy[idx]

    def pre_mean(idx):
        # use all data before idx, edge safe
        if idx > 0:
            return np.mean(vy[:idx])
        return vy[idx]

    # ------------------------------------------------------------
    # No peaks
    # ------------------------------------------------------------
    if len(peaks) == 0:
        if np.abs(np.mean(vy)) > epsilon:
            return np.array([none_val, 0.0], dtype=float)

        return y

    # Sort peaks by height, tallest first
    sort_idx = np.argsort(heights)[::-1]
    peaks_sorted = peaks[sort_idx]
    heights_sorted = heights[sort_idx]

    # Tallest peak too small
    if heights_sorted[0] < threshold:
        if np.abs(np.mean(vy)) > epsilon:
            return np.array([none_val, 0.0], dtype=float)

        return y
    # ------------------------------------------------------------
    # Two-peak logic
    # ------------------------------------------------------------

    # Collapse/failure peak must be the tallest significant peak
    p_fail = int(peaks_sorted[0])
    h_fail = dcodt[p_fail]

    # Arch formation must be an earlier, strictly smaller, significant peak
    arch_candidates = peaks[
        (peaks < p_fail) &
        (dcodt[peaks] >= second_ratio * h_fail) &
        (dcodt[peaks] < h_fail) &
        (dcodt[peaks] >= threshold)
    ]

    # Search earlier arch candidates in time order
    for p_arch in np.sort(arch_candidates):
        p_arch = int(p_arch)
        h_arch = dcodt[p_arch]

        valley_ok = has_valley(p_arch, p_fail)
        distance_ok = peaks_far_enough_apart(p_arch, p_fail)

        if not (valley_ok and distance_ok):
            continue

        pre = pre_mean(p_fail)
        post = post_mean(p_fail)

        if (np.abs(pre) < epsilon) and (np.abs(post) > epsilon):
            return np.array(
                [float(time[p_arch]), float(time[p_fail])],
                dtype=float,
            )

    # ------------------------------------------------------------
    # Single-peak / no-arch logic
    # ------------------------------------------------------------
    p1 = int(peaks_sorted[0])

    # If velocity is generally moving, interpret tallest peak as failure / no arch.
    if np.abs(np.mean(vy)) > epsilon:
        y[1] = float(time[p1])
        return y

    # Otherwise, if velocity is locally still around the peak, stable arch.
    local_v = mean_window(vy, p1)

    if np.abs(local_v) < epsilon:
        y[0] = float(time[p1])
        return y

    return y
    

def compute_y3(
    dcodt,
    vy,
    time,
    threshold=0.1,
    second_ratio=0.1,
    none_val=np.nan,
    epsilon=0.01,
    local_window=2,
    min_valley_duration_h=1/3,
):
    """
    Returns array([y1, y2]) in seconds.

    y1 = arch formation time
    y2 = failure / arch collapse time

    none_val returned if arch never forms/collapses.

    Logic:
      - unstable arch:
          arch formation is an earlier, strictly smaller peak
          collapse is the tallest peak
          damage rate stays below threshold between the two peaks
          for at least min_valley_duration_h
      - stable arch:
          tallest peak is arch formation
          no valid later collapse
      - no arch:
          moving case with no valid arch formation before collapse
    """

    dcodt = np.asarray(dcodt, dtype=float)
    vy = np.asarray(vy, dtype=float)
    time = np.asarray(time, dtype=float)

    if not (len(dcodt) == len(vy) == len(time)):
        raise ValueError("dcodt, vy, and time must have the same length.")

    n = len(time)
    y = np.array([none_val, none_val], dtype=float)

    if n == 0:
        return y

    peaks, _ = find_peaks(dcodt)

    if len(peaks) == 0:
        if np.abs(np.mean(vy)) > epsilon:
            return np.array([none_val, 0.0], dtype=float)

        return y

    heights = dcodt[peaks]

    min_valley_duration_s = min_valley_duration_h * 3600.0

    if n > 1:
        dt_sample = np.nanmedian(np.diff(time))
    else:
        dt_sample = 0.0

    # ------------------------------------------------------------
    # helper functions
    # ------------------------------------------------------------
    def mean_window(x, center, before=local_window, after=local_window):
        lo = max(0, center - before)
        hi = min(len(x), center + after + 1)

        if hi <= lo:
            return np.nan

        return np.mean(x[lo:hi])

    def has_sustained_valley(p_a, p_b):
        """
        Returns True if dcodt is below threshold for at least
        min_valley_duration_s between p_a and p_b.

        This looks for a continuous run of samples below threshold.
        """
        lo, hi = sorted((int(p_a), int(p_b)))

        if hi - lo <= 1:
            return False

        idx = np.arange(lo + 1, hi)
        below = dcodt[idx] < threshold

        if not np.any(below):
            return False

        # Find contiguous True runs in `below`
        padded = np.r_[False, below, False]
        starts = np.flatnonzero((~padded[:-1]) & padded[1:])
        ends = np.flatnonzero(padded[:-1] & (~padded[1:])) - 1

        for start, end in zip(starts, ends):
            i0 = idx[start]
            i1 = idx[end]

            # Add one sample spacing because a one-sample run represents
            # roughly one saved-output interval.
            duration = time[i1] - time[i0] + dt_sample

            if duration >= min_valley_duration_s:
                return True

        return False

    def post_mean(idx):
        # use all data after idx, edge safe
        if idx < n - 1:
            return np.mean(vy[idx:])
        return vy[idx]

    def pre_mean(idx):
        # use all data before idx, edge safe
        if idx > 0:
            return np.mean(vy[:idx])
        return vy[idx]

    # ------------------------------------------------------------
    # Sort peaks by height, tallest first
    # ------------------------------------------------------------
    sort_idx = np.argsort(heights)[::-1]
    peaks_sorted = peaks[sort_idx]
    heights_sorted = heights[sort_idx]

    # Tallest peak too small
    if heights_sorted[0] < threshold:
        if np.abs(np.mean(vy)) > epsilon:
            return np.array([none_val, 0.0], dtype=float)

        return y

    # ------------------------------------------------------------
    # Two-peak logic
    # ------------------------------------------------------------

    # Collapse/failure peak must be the tallest significant peak
    p_fail = int(peaks_sorted[0])
    h_fail = dcodt[p_fail]

    # Arch formation must be an earlier, strictly smaller, significant peak
    arch_candidates = peaks[
        (peaks < p_fail) &
        (dcodt[peaks] >= second_ratio * h_fail) &
        (dcodt[peaks] < h_fail) &
        (dcodt[peaks] >= threshold)
    ]

    # Search earlier arch candidates in time order
    for p_arch in np.sort(arch_candidates):
        p_arch = int(p_arch)

        valley_ok = has_sustained_valley(p_arch, p_fail)

        if not valley_ok:
            continue

        pre = pre_mean(p_fail)
        post = post_mean(p_fail)

        if (np.abs(pre) < epsilon) and (np.abs(post) > epsilon):
            return np.array(
                [float(time[p_arch]), float(time[p_fail])],
                dtype=float,
            )

    # ------------------------------------------------------------
    # Single-peak / no-arch logic
    # ------------------------------------------------------------
    p1 = int(peaks_sorted[0])

    # If velocity is generally moving, interpret tallest peak as failure / no arch.
    if np.abs(np.mean(vy)) > epsilon:
        y[1] = float(time[p1])
        return y

    # Otherwise, if velocity is locally still around the peak, stable arch.
    local_v = mean_window(vy, p1)

    if np.abs(local_v) < epsilon:
        y[0] = float(time[p1])
        return y

    return y
##### testing

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

    # optional
    parser.add_argument("--threshold", type=float, default=0.1)
    parser.add_argument("--second_ratio", type=float, default=0.1)
    parser.add_argument("--epsilon", type=float, default=0.01)
    parser.add_argument("--y_cutoff", type=float, default=120e3)

    args = parser.parse_args()

    J = args.J
    base_dir = Path(args.base_dir)

    # hard coding simulation setup
    u_max = 20
    rho_a = 1.3
    C_a = 0.0012
    T_max = 36
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
        y = compute_y3(
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
        ax.set_ylabel(r"$|\frac{d\text{Co#}}{dt}|$", size='xx-large', color='gray')
        ax.set_ylim(-0.01, 1.8)
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
