# obtain domain averaged shear stress from lammps dump files

import xarray as xr
import numpy as np
import sys
import glob
import os

def first_peak_before_drop_idx(stress_array, strain_array):
    """
    Return (idx, stress_at_peak, strain_at_peak).

    Peak = stress at the first negative diff ("first drop").
    If there is no drop, return the global max stress and its strain.
    """
    s = np.asarray(stress_array)
    e = np.asarray(strain_array)

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


def main():
    # expecting 2 arguments: path_to_dumped_atoms output_directory
    if len(sys.argv) != 2:
        print("Usage: compute_macroscopic_stress.py path_to_dumped_atoms_dir")
        sys.exit(1)

    _, datadir = sys.argv

    outputdir = os.path.join(datadir, 'figures')
    os.makedirs(outputdir, exist_ok=True)

    # ----- hard coding file parsing and strain calculation -----
    Ls = np.array(['100', '200', '300', '400', '500', '600'])
    gs = 1/Ls.astype(float)
    seeds = np.array([8302, 3094, 3894, 1092, 4999])
    natoms_list = np.array(['3496', '13987', '31473', '55952', '87426', '125893'])
    run_time = 1800
    time = np.linspace(0, run_time, len(macro))
    L = L_ini * (1+shear_rate*time)
    strain = (L-L_ini)/L_ini

    # ------ iterate ------
    all_stress_data = np.empty([len(strain), len(gs), len(seeds)])
    for j, L in enumerate(Ls):
        directory = os.path.join(datadir, f'L{L}')
        if os.path.isdir(directory):
            for k, seed in enumerate(seeds):
                fname = os.path.join(directory, f'N{natoms_list[j]}Dump_seed{seed}.lammps')
                base = os.path.splitext(os.path.basename(fname))[0][:-7]
                fpath = os.path.join(directory, f'{base}_shear.npy')
                stress = np.load(fpath)
                all_stress_data[:, j, k] = stress[:500]

    # ------- create dataset ------- 
    ds = xr.Dataset(
        data_vars = dict(
                stress=(['strain', 'granularity', 'seed'], all_stress_data)
                ), 
            coords = dict(
                strain = strain,
                granularity = gs,
                seed = seeds
            )
        )

if __name__ == "__main__":
    main()




# ---------- paths ----------
base_dir = "data/size_effect"          # folder containing 'processed' subdir
processed_dir = os.path.join(base_dir, "processed")

# ---------- simulation constants (for strain axis) ----------
shear_rate = 1e-7      # [1/s]
run_time  = 3600.0     # [s]

# ---------- system sizes / mapping ----------
Ls = np.array([50., 100., 200., 300., 400., 500.])          # same L_s as in processing script
natoms_list = np.array([873, 3496, 13987, 31473, 55952, 87426])

granularity = 1.0 / Ls   # 1/L_s
order = np.argsort(granularity)   # increasing 1/L_s (low g -> high g)

# ---------- infer number of timesteps & build strain ----------
example_tag = f"L{int(Ls[0])}"
example_path = os.path.join(processed_dir,
                            f"macros_{example_tag}_runsxsteps.npy")
example_macros = np.load(example_path)
n_steps = example_macros.shape[1]

time = np.linspace(0.0, run_time, n_steps)
strain = shear_rate * time   # gamma = dotgamma * t

# ---------- plotting ----------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
fig.tight_layout(pad=2.0, w_pad=4.0)

colors = ['indianred', 'sandybrown', 'gold', 'yellowgreen', 'skyblue', 'slateblue']

peak_means = []
peak_stderrs = []
gran_arr = []
peak_strain_means = []

# NEW: store peak frame stats
peak_frame_means = []
peak_frame_stderrs = []
peak_frame_mean_curve = []  # peak frame from ensemble-mean curve

for j, idx in enumerate(order):
    Lval = Ls[idx]
    Nval = natoms_list[idx]
    tag = f"L{int(Lval)}"

    # stacked runs array: shape (n_seeds, n_steps)
    macros_path = os.path.join(processed_dir,
                               f"macros_{tag}_runsxsteps.npy")
    macros = np.load(macros_path)
    n_seeds = macros.shape[0]

    # mean and stderr
    mean_macro = macros.mean(axis=0)
    if n_seeds > 1:
        stderr_macro = macros.std(axis=0, ddof=1) / np.sqrt(n_seeds)
    else:
        stderr_macro = np.zeros_like(mean_macro)

    # ---- left subplot: stress–strain with shaded stderr ----
    c = colors[j % len(colors)]
    label = rf"$g = {sci_fmt(granularity[idx], 2)}$"
    ax1.plot(strain, mean_macro/1e3, color=c, lw=2, label=label, zorder=6-j)
    ax1.fill_between(
        strain,
        (mean_macro - stderr_macro)/1e3,
        (mean_macro + stderr_macro)/1e3,
        color=c, alpha=0.2, zorder=6-j
    )

    # ---- per-run peaks (frame + stress + strain) ----
    peak_stresses = []
    peak_strains = []
    peak_frames = []

    for m in macros:
        i_peak, p_stress, p_strain = first_peak_before_drop_idx(m, strain)
        peak_frames.append(i_peak)
        peak_stresses.append(p_stress)
        peak_strains.append(p_strain)

    peak_stresses = np.array(peak_stresses)
    peak_strains = np.array(peak_strains)
    peak_frames = np.array(peak_frames, dtype=float)

    peak_mean = peak_stresses.mean()
    peak_strain_mean = peak_strains.mean()

    if n_seeds > 1:
        peak_stderr = peak_stresses.std(ddof=1) / np.sqrt(n_seeds)
        frame_stderr = peak_frames.std(ddof=1) / np.sqrt(n_seeds)
    else:
        peak_stderr = 0.0
        frame_stderr = 0.0

    # ---- peak frame on ENSEMBLE-MEAN curve (often the most "representative") ----
    i_mean, s_mean_peak, e_mean_peak = first_peak_before_drop_idx(mean_macro, strain)

    # store for later
    peak_means.append(peak_mean)
    peak_stderrs.append(peak_stderr)
    peak_strain_means.append(peak_strain_mean)
    gran_arr.append(granularity[idx])

    peak_frame_means.append(np.nanmean(peak_frames))
    peak_frame_stderrs.append(frame_stderr)
    peak_frame_mean_curve.append(i_mean)

# ---------- print peak frame # by granularity ----------
gran_arr = np.array(gran_arr)
peak_strain_means = np.array(peak_strain_means)
peak_frame_means = np.array(peak_frame_means)
peak_frame_stderrs = np.array(peak_frame_stderrs)
peak_frame_mean_curve = np.array(peak_frame_mean_curve)

print("\nPeak FRAME # by g (same order as gran_arr):")
for g_val, f_run, f_run_se, f_mean in zip(gran_arr, peak_frame_means, peak_frame_stderrs, peak_frame_mean_curve):
    print(f"g={g_val:.6f}  peak_frame(runs)= {f_run:.2f} ± {f_run_se:.2f}   peak_frame(mean_curve)= {int(f_mean)}")

print("\nMean peak strains by g (same order as gran_arr):")
for g_val, eps in zip(gran_arr, peak_strain_means):
    print(g_val, eps)

# ---------- cosmetics ----------
# left
ax1.set_xlabel(r"$\gamma$", size='x-large')
ax1.set_ylabel(r"$\tau$ [kPa]", size='x-large')
ax1.legend(fontsize=8, loc='upper left')
ax1.set_ylim(-2, 25)

plt.sca(ax1)
formatter = ScalarFormatter(useMathText=True)
formatter.set_powerlimits((-2, 2))
ax1.xaxis.set_major_formatter(formatter)
plt.xticks(rotation=45)

# ---- Add shear modulus annotation (triangle) to ax1 ----
x0 = 5.35e-5     # left bottom x of the triangle
x1 = 7.5e-5      # right top x of the triangle
y0 = 7           # bottom y
y1 = 9.3         # top y

ax1.plot([x1, x1], [y0, y1], color='k', lw=1)        # vertical
ax1.plot([x0, x1], [y0, y0], color='k', lw=1)        # top horizontal
ax1.text((x0+x1)/2-1e-6, y0-0.6, r"1", fontsize=8, va='center')
ax1.text(x1+3e-6, (y0 + y1)/2, r"$G$", fontsize=8, va='center')

# right
peak_means = np.array(peak_means)
peak_stderrs = np.array(peak_stderrs)

ax2.errorbar(gran_arr, peak_means/1e3,
             yerr=peak_stderrs/1e3, color='k',
             fmt="o-", capsize=4, lw=1.5)
plt.sca(ax2)
plt.xticks(rotation=45)

formatter = ScalarFormatter(useMathText=True)
formatter.set_powerlimits((-2, 2))
ax2.xaxis.set_major_formatter(formatter)

ax2.set_xlabel(r"$g$", size='x-large')
ax2.axvline(0.005, color='indianred', lw=1, ls='--', zorder=0)
ax2.fill_betweenx([14200/1e3, 16800/1e3], -0.002, 0.005,
                  color='indianred', alpha=0.1, zorder=1)
ax2.set_ylim(14200/1e3, 16800/1e3)
ax2.set_xlim(0.001, 0.021)
plt.xticks(rotation=45)

ax2.grid(which="major", linestyle="-", linewidth=0.5, color="0.85", alpha=0.6)
ax2.set_axisbelow(True)
ax1.grid(which="major", linestyle="-", linewidth=0.5, color="0.85", alpha=0.6)
ax1.set_axisbelow(True)

ax2.set_ylabel(r"$\tau_{\text{max}}$ [kPa]", size='x-large')

plt.savefig(
    '/Users/arlenlex/Library/CloudStorage/OneDrive-Stanford/Stanford_Research/2025 Autumn/first_draft_paper_figs/size_effect.png',
    dpi=300, bbox_inches='tight'
)

plt.show()
