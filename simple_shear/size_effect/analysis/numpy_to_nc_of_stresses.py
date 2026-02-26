# obtain domain averaged shear stress from lammps dump files
import numpy as np
import sys
import os
import pandas as pd

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

    outputdir = os.path.join(datadir, 'output_files')
    os.makedirs(outputdir, exist_ok=True)

    # ----- hard coding file parsing and strain calculation -----
    Ls = np.array([50, 60, 80, 100, 200, 300, 400, 500, 600, 700, 800, 900], dtype=int)
    gs = 1.0 / Ls.astype(float)
    seeds = np.array([8302, 3094, 3894, 1092, 4999, 5673, 1927, 2380, 8095, 8920], dtype=int)
    natoms_list = np.array([873, 1257, 2237, 3496, 13987, 31475, 55954, 87431, 125900, 171366, 223824, 283279], dtype=int)

    # ---------- simulation constants (for strain axis) ----------
    shear_rate = 1e-7      # [1/s]
    run_time  = 1800       # [s]

    # ------ infer number of timesteps from the first file ------
    first_L = Ls[0]
    first_seed = seeds[0]
    first_dir = os.path.join(outputdir, f'L{first_L}')
    first_N = f"{natoms_list[0]:06d}"
    first_base = f"N{first_N}Dump_seed{first_seed}"
    first_fpath = os.path.join(first_dir, f"{first_base}_shear.npy")

    example = np.load(first_fpath)
    n_steps = example.shape[0]

    time = np.linspace(0.0, run_time, n_steps)
    strain = shear_rate * time

    # ------ iterate ------
    all_stress_data = np.empty([len(strain), len(gs), len(seeds)])
    for j, L in enumerate(Ls):
        directory = os.path.join(outputdir, f'L{L}')
        if os.path.isdir(directory):
            for k, seed in enumerate(seeds):
                Nstr = f"{natoms_list[j]:06d}"   # leading zeros: 6 digits
                base = f"N{Nstr}Dump_seed{seed}"
                fpath = os.path.join(directory, f'{base}_shear.npy')
                stress = np.load(fpath)
                all_stress_data[:, j, k] = stress[:len(strain)]
        else:
            print(f'path {directory} does not exist.', flush = True)

    # ------- save as a tidy dataframe -------
    # columns: strain, granularity, seed, stress
    df = (
        pd.DataFrame(
            all_stress_data.reshape(len(strain), len(gs) * len(seeds)),
            index=pd.Index(strain, name="strain"),
            columns=pd.MultiIndex.from_product([gs, seeds], names=["granularity", "seed"]),
        )
        .stack(["granularity", "seed"])
        .rename("stress")
        .reset_index()
    )

    outpath = os.path.join(outputdir, "macroscopic_shear_stress.parquet")
    df.to_parquet(outpath, index=False)
    print(f"Saved: {outpath}", flush=True)

if __name__ == "__main__":
    main()
