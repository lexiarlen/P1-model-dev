import argparse
import re
import os
import glob
from pathlib import Path
import numpy as np
from scipy.signal import find_peaks
from ovito.io import import_file
import matplotlib.pyplot as plt  

### some functions for getting y from lammps dump files ###

def get_ts_from_data(fname, dt, y_cutoff=120e3):
    """
    Load LAMMPS dump and compute dcodt and avg vy in the upper half.
    Returns: dcodt [arbitrary units per frame], avg_vy [m/s], time [s].
    """
    pipeline = import_file(fname)

    nframes = pipeline.source.num_frames
    timesteps = []
    avg_co = []
    avg_vy = []

    for frame in range(nframes):
        data = pipeline.compute(frame)

        # record timestep
        ts = data.attributes.get("Timestep", frame)
        timesteps.append(float(ts))

        # mean coordination number (active particles only)
        if 'c_nbond' in data.particles.keys() and 'v_ingroup' in data.particles.keys():
            co = data.particles['c_nbond']
            ingroup = data.particles['v_ingroup']
            mask = (ingroup == 1)
            filtered_co = co[mask]
            if len(filtered_co) > 0:
                avg_co.append(np.mean(filtered_co))
            else:
                print(f"No particles with v_ingroup == 1 at frame {frame}")
        else:
            print(f"Properties 'c_nbond' or 'v_ingroup' not found at frame {frame}")

        # velocity in upper half (active particles only)
        if 'Velocity' in data.particles.keys() and 'v_ingroup' in data.particles.keys():
            vy = data.particles['Velocity.Y']
            ypos = data.particles['Position.Y']
            v_ingroup = data.particles['v_ingroup']
            mask = (v_ingroup == 1) & (ypos >= y_cutoff)
            filtered_vy = vy[mask]
            if len(filtered_vy) > 0:
                avg_vy.append(np.mean(filtered_vy))
            else:
                print(f"No particles in upper half (y>={y_cutoff}) with v_ingroup == 1 at frame {frame}")
        else:
            print(f"Properties 'Velocity' or 'v_ingroup' not found at frame {frame}")

    dcodt = np.abs(np.gradient(np.array(avg_co)))  # keep as-is w.r.t. your thresholds
    return dcodt, np.array(avg_vy), np.array(timesteps) * dt

def compute_y(dcodt, vy, time, threshold=0.005, second_ratio=0.05, none_val=999, epsilon=0.01):
    peaks, _ = find_peaks(dcodt)
    heights = dcodt[peaks]

    # if there are no peaks: 
    if len(peaks) == 0:
        if (np.abs(np.mean(vy)) > epsilon):
            return np.array([0, none_val]) # fail immediately case
        else:
            return np.array([none_val, none_val]) # no fracturing case

    # sort peaks by height (tallest first)
    sort_idx = np.argsort(heights)[::-1]
    peaks_sorted   = peaks[sort_idx]
    heights_sorted = heights[sort_idx]

    # tallest peak check
    if heights_sorted[0] < threshold:
        if (np.abs(np.mean(vy)) > epsilon):
            print('fail immediate case')
            return np.array([0, none_val]) # fail immediately case
        else:
            return np.array([none_val, none_val]) # no fracturing case

    peak_indices = [peaks_sorted[0]]
    print(time[peaks_sorted] / 3600)

    # ---------------------------------------------
    # Search for a second peak among the next 4 highest
    # ---------------------------------------------
    if len(peaks_sorted) > 1:
        print('in potential second peak loop')

        # we will check up to the 5 highest peaks total (1 primary + next up to 4)
        max_check = min(5, len(peaks_sorted))

        for j in range(1, max_check):
            p1 = peaks_sorted[0]       # tallest peak
            p2 = peaks_sorted[j]       # j-th tallest peak
            h2 = heights_sorted[j]

            # amplitude criterion
            if h2 < second_ratio * heights_sorted[0]:
                continue

            print(f'testing candidate p2 (idx={p2}) at t={time[p2]}')

            # valley criterion between p1 and p2
            lo, hi = (p1, p2) if p1 < p2 else (p2, p1)
            has_valley = (hi - lo > 1) and (np.nanmin(dcodt[lo+1:hi]) < threshold)

            if not has_valley:
                continue

            print('  -> has valley')

            # velocity logic
            if (np.abs(np.mean(vy[:p2])) < epsilon) and (np.abs(np.mean(vy[p2:])) > epsilon):
                peak_indices.append(p2)
                print('  -> accepted as second peak')
                break  # stop once we find a valid second peak

    # ---------------------------------------------
    # Rest of your logic (unchanged)
    # ---------------------------------------------
    if len(peak_indices) == 2:
        y = time[np.sort(np.array(peak_indices))]
    else:
        y = np.ones(2) * none_val
        idx = peak_indices[0]
        if (idx > 1) & (idx < len(vy)-2) & (np.abs(np.mean(vy[idx-2:idx+2])) < epsilon):
            print('stable arch')
            y[0] = time[idx]
            y[1] = none_val
        if (np.abs(np.mean(vy)) > epsilon):
            print('no arch')
            y[0] = none_val
            y[1] = time[idx]
        if (idx <= 1) & (np.abs(np.mean(vy[:2])) < epsilon):
            print('stable arch but fail quick')
            y[0] = time[idx]
            y[1] = none_val

    return y

def get_u(time, u_max, max_time):
    if time == 999:
        return 999
    else:
        return u_max * time / max_time

get_u_vec = np.vectorize(get_u)

def get_wind_stress(u, rho_a, C_a):
    if u == 999:
        return 999
    else:
        return rho_a * C_a * np.abs(u) * u

get_wind_stress_vec = np.vectorize(get_wind_stress)

# -------------------- CLI --------------------
parser = argparse.ArgumentParser()
parser.add_argument("--iter", type=int, default=0)   # was required; made optional to match launcher
parser.add_argument("--J", type=int, required=True)
parser.add_argument("--max_wind_speed", type=float, required=False)
parser.add_argument("--parent_dir_name", type=str, required=False)
parser.add_argument("--rho_a", type=float, required=False)
parser.add_argument("--C_a", type=float, required=False)
parser.add_argument("--T_max", type=float, required=False)
parser.add_argument("--dt", type=float, required=False)  # kept for compatibility, not used now

# NEW: ice radius and density (for mass and dt from k_b)
parser.add_argument("--radius", type=float, required=False, default = 500)
parser.add_argument("--rho_i", type=float, required=False, default = 920)

args = parser.parse_args()

ITER = args.iter
J = args.J
u_max = args.max_wind_speed if args.max_wind_speed is not None else 20.0
parent_dir_name = args.parent_dir_name if args.parent_dir_name is not None else "ensemble"
rho_a = args.rho_a if args.rho_a is not None else 1
C_a = args.C_a if args.C_a is not None else 0.0012
T_max = args.T_max if args.T_max is not None else 36
max_time = T_max*3600 # convert [hrs] to [s]

radius = args.radius
rho_i = args.rho_i

# ice mass (per floe), used for dt computation
radius_small = radius*3/5
m_ice = np.pi * radius_small**2 * rho_i

# Paths consistent with the rest of the pipeline
base = Path(f"/scratch/groups/earlew/arlenlex/P1-model-dev/nares/{parent_dir_name}")
lammps_dir = base / "lammps_output"
state_dir = base / "state"
state_dir.mkdir(parents=True, exist_ok=True)

# ---- NEW: load THETA and extract k_b for each ensemble member ----
theta_path = state_dir / "theta.tsv"
THETA = np.array(np.loadtxt(theta_path), dtype=float).T   # shape (n_params, J); k_b is row 0
kb_values = THETA[0, :]   # k_b for each ensemble member, 0-based index

# -------------------- analysis --------------------
pattern = re.compile(r'^dump(\d+)\.lammps$')

Y = 9999 * np.ones((2, J))
SIG = 9999 * np.ones((2, J))

files = sorted(glob.glob(str(lammps_dir / "*.lammps")))
print(f"Found {len(files)} dump files under {lammps_dir}")

for dumpfile in files:
    fname = os.path.basename(dumpfile)
    m = pattern.match(fname)
    if not m:
        print(f"Skipped: {dumpfile}")
        continue

    ensemble_member_id = int(m.group(1))
    if not (1 <= ensemble_member_id <= J):
        print(f"Skipping {dumpfile}: member id {ensemble_member_id} out of 1..{J}")
        continue

    print(f'File: {dumpfile} → ID: {ensemble_member_id}')

    # k_b for this ensemble member (column index ensemble_member_id-1)
    kb_local = kb_values[ensemble_member_id - 1]

    # dt = 0.1 * sqrt(m / k_b)
    dt_local = 0.1 * np.sqrt(m_ice / kb_local)

    dcodt, vy, time = get_ts_from_data(dumpfile, dt_local)
    y = compute_y(dcodt, vy, time)

    # -------------------- figure generation (per dumpfile) --------------------
    vy_up = vy  # alias so we can use your plotting code as-is

    fig, ax = plt.subplots(figsize=(10,4))
    line1, = plt.plot(time/3600, dcodt, lw = 4, color = 'gray', label = 'damage rate')
    plt.xlabel("time (hrs)", size = 'large')
    plt.xticks(size='large'); plt.yticks(size='large')

    plt.ylabel(r"$\frac{d\text{Co#}}{dt}$", size = 'xx-large', color = 'gray')
    plt.ylim(-0.01,0.35)
    plt.grid()

    ax2 = ax.twinx()
    plt.sca(ax2)
    plt.ylabel(r"mean $v_y$ (m/s) in upper domain", size = 'large', color = 'tab:blue')
    plt.yticks(size='large')
    line2, = plt.plot(time/3600, vy_up, lw = 4, color = 'r', label = r'$v_{y}^{up}$')
    plt.ylim(-0.3, 0.1)

    line3 = plt.axvline(y[0]/3600, color = 'k', ls = '-.', label = rf'$y_1 = {(y[0]/3600):.1f}$ h')
    line4 = plt.axvline(y[1]/3600, color = 'k', ls = '--', label = rf'$y_2 = {(y[1]/3600):.1f}$ h')
    lines = [line1,line2,line3,line4]
    labels = [line.get_label() for line in lines]
    ax.legend(lines, labels, loc='best')
    plt.xlim(0,36)

    # save figure to the same directory where Y.npy and SIG.npy go
    fig_path = state_dir / f"diagnostic_{ensemble_member_id:03d}.png"
    plt.savefig(fig_path, bbox_inches='tight', dpi = 300)
    plt.close(fig)
    # -------------------------------------------------------------------------

    u = get_u_vec(y, u_max, max_time)
    sig = get_wind_stress_vec(u, rho_a, C_a)

    SIG[:, ensemble_member_id - 1] = sig
    Y[:, ensemble_member_id - 1] = y


np.save(state_dir / 'Y.npy', Y)
np.save(state_dir / 'SIG.npy', SIG)

print(f"Saved Y.npy and SIG.npy under {state_dir}")
print(f"Saved diagnostic figures (one per dumpfile) under {state_dir}")
